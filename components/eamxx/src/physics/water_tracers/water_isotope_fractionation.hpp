
#ifndef WATER_ISOTOPE_FRACTIONATION
#define WATER_ISOTOPE_FRACTIONATION

#include "share/eamxx_types.hpp"
#include "physics_constants.hpp"

#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/logging/ekat_logger.hpp"

#include <string>
#include <vector>
#include <array>

namespace scream { 
namespace WaterIsotopes {

/*

Functions to calculate equilibrium fractionation factors 

*/ 

template <typename Scalar>
Scalar AlphaEqIceVapor(const std::string& iso, const Scalar& tk) {

  using Wiso = WaterIsotopologues<Scalar>;
  // calculate equilibrium alpha for ice<->vapor transitions
  Scalar wiso_alpi = 1.0;
  if (iso != "H2O") {
    // Calculate fractionation factors after Merlivat & Nief, 1967 for HDO
    // and Majoube 1971 for H218O and H217O. Need to modify for H217O and HTO.
    wiso_alpi = exp(Wiso::alpai[IsotopologueToIndex.at(iso)]*pow(tk,-2) +
                Wiso::alpbi[IsotopologueToIndex.at(iso)]/tk +
                Wiso::alpci[IsotopologueToIndex.at(iso)]) ;
  };

  // update for H217O
  if (iso == "H217O") {
    wiso_alpi = pow(wiso_alpi, 0.528);
  } else if (iso == "HTO") { // update for HTO
    wiso_alpi = pow(wiso_alpi, 2.0);
  };
  assert(wiso_alpi >= 1.0);
  return wiso_alpi;
};

template <typename Scalar>
Scalar AlphaEqLiquidVapor(const std::string& isosp, const Scalar& tk) {

  using Wiso = WaterIsotopologues<Scalar>;
  // calculate equilibrium alpha for liquid<->vapor transitions
  Scalar wiso_alpl = 1.0;
  if (isosp != "H2O") {
    if (isosp == "HDO" || isosp == "HTO") { // these have a different structure
      wiso_alpl = exp(Wiso::alpal[IsotopologueToIndex.at(isosp)]*pow(tk,3) + 
                  Wiso::alpbl[IsotopologueToIndex.at(isosp)]*pow(tk,2) +
                  Wiso::alpcl[IsotopologueToIndex.at(isosp)]*tk + 
                  Wiso::alpdl[IsotopologueToIndex.at(isosp)] + 
                  Wiso::alpel[IsotopologueToIndex.at(isosp)]*pow(tk,-3));
    } else {
      wiso_alpl = exp(Wiso::alpal[IsotopologueToIndex.at(isosp)]*pow(tk,-3) + 
                  Wiso::alpbl[IsotopologueToIndex.at(isosp)]*pow(tk,-2) +
                  Wiso::alpcl[IsotopologueToIndex.at(isosp)]/tk + 
                  Wiso::alpdl[IsotopologueToIndex.at(isosp)]);
    };
  };
  // update for H217O
  if (isosp == "H217O") {
    wiso_alpl = pow(wiso_alpl, 0.528);
  } else if (isosp == "HTO") { // update for HTO
    wiso_alpl = pow(wiso_alpl, 2.0);
  };
  assert(wiso_alpl >= 1.0);
  return wiso_alpl;
};

/* 
Functions to calculate fractionations including kinetic effects (e.g., partial evaporation of
hydrometeors, formation of ice crystals under supersaturation, evaporation from ocean surface)
*/

template <typename Scalar>
Scalar AlphaKineticEvap(const std::string& isosp, const Scalar& tk, const Scalar& hum0, const Scalar& alpeq) {
  
  using Wiso = WaterIsotopologues<Scalar>;
  // return fractionation factor modified for kinetic effects during
  // cloud/hydrometeor liquid evaporation into unsaturated air
  // note that this is configured such that R_dest = wiso_akel*R_source

  //declare local variables
  const Scalar h0 = min(1.0, hum0);
  const Scalar difrmj = Wiso::difrm[IsotopologueToIndex.at(isosp)];
  Scalar heff;
  Scalar dondi;
  Scalar wiso_akel;

  heff = h0; // tempoerary, needs to be updated
  dondi = pow((1.0/difrmj),Wiso::dkfac);
  //TODO: Rederive this equation
  wiso_akel = alpeq*heff / (alpeq*dondi*(heff-1.0) + 1.0);
  assert(wiso_akel >= 1.0)
  return wiso_akel;

};

template <typename Scalar>
Scalar AlphaKineticDepo(const std::string& isosp, const Scalar& tk, const Scalar& rh, const Scalar& alpeq) {
  
  using Wiso = WaterIsotopologues<Scalar>;
  // return fractionation factor modified for kinetic effects during
  // deposition to ice. Optional use of vapor supersaturation.

  //declare local variables
  Scalar sat1;
  Scalar difrmj;
  Scalar dondi;
  Scalar wiso_akci;

  if (tk < Wiso::tkini) {
    sat1 = min(max(1.0, rh), 0); //this needs updating
    difrmj = Wiso::difrm[IsotopologueToIndex.at(isosp)];
    dondi = (1.0/difrmj);
    wiso_akci = alpeq*sat1 / (alpeq*dondi*(sat1-1.0) + 1.0);
  } else {
    wiso_akci = alpeq;
  };
  return wiso_akci;
};


/*
* Functions to calculate kinetic fractionation factors
*/

template <typename Scalar>
Scalar AlphaKMol(const std::string& isosp, const Scalar& rbot, const Scalar& zbot, const Scalar& ustar) {

  /* INPUTS:
  isosp: isotopologue species
  rbot: density of lowest layer [kg/m3]
  zbot: height of lowest model layer [m]
  ustar: friction velocity [m/s]
  */
  using Wiso = WaterIsotopes::WaterIsotopologues<Scalar>;
  using Constants = physics::Constants<Scalar>;

  /* compute kinetic modifier for drag coefficient.
  uses Brutsaert equations for the turbulent layer using GCM 
  quantities. This calculates kinetic fractionation using the 
  Merlivat and Jouzel 1979 closure assumption */

  const Scalar difair = 2.36e-5;    // molecular diffusivity of air
  const Scalar muair = 1.789e-5;    // dynamic viscosity of air [Pa*s]

  // local variables
  const Scalar z0 = pow(ustar, 2.0)/(81.1*Constants::gravit); // roughness length via Charnock's equation
  const Scalar vmu = muair / rbot; // kinematic viscocity of air [m2/s]
  const Scalar reno = ustar*z0 / vmu; // Reynolds number
  const Scalar sc = vmu/difair; // Schmidt number (momentum to mass diffusivity)
  
  if (reno < Wiso::renocrit) { // Smooth regime
    const Scalar diffpow = 2.0/3.0;
    const scalar tmr = ((1.0/Constants::Karman)*log(zbot*ustar/(30.0*vmu)))/(13.6 * pow(sc,diffpow));
  } else {
    const Scalar diffpow = 1.0/2.0;
    const Scalar tmr = ((1.0/Constants::Karman)*log(zbot/z0)-5.0)/(7.3*pow(reno,0.25)*pow(sc,diffpow));
  }

  const difn = pow(1.0/Wiso::difrm[isosp];, diffpow);
  const kmol = (difn - 1.0)/(difn + tmr);

  const Scalar alphakn = 1.0 - kmol;

  assert(alphakn <= 1.0);

  return alphakn;

};

void WtypeGetAlpha() {
  /* return fractionation factor given a source water type and 
  a destination water type.
  */
};


}
}

#endif // #ifndef