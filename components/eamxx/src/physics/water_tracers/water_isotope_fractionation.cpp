
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

// Functions to calculate equilibrium fractionation factors 

template <typename Scalar>
Scalar AlphaEqIceVapor(const std::string& iso, const Scalar& tk) {

  using Wiso = WaterIsotopologues<Scalar>;\
  // calculate equilibrium alpha for ice<->vapor transitions
  Scalar wiso_alpi = 1.0;
  if (iso != "H2O") {
    // Calculate fractionation factors after Merlivat & Nief, 1967 for HDO
    // and Majoube 1971 for H218O and H217O. Need to modify for H217O and HTO.
    wiso_alpi = exp(Wiso::alpai[IsotopologueToIndex.at(iso)]*pow(tk,-2) +
                Wiso::alpbi[IsotopologueToIndex.at(iso)]/tk +
                Wiso::alpci[IsotopologueToIndex.at(iso)]) ;
  } 

  // update for H217O
  if (iso == "H217O") {
    wiso_alpi = pow(wiso_alpi, 0.528);
  } else if (iso == "HTO") { // update for HTO
    wiso_alpi = pow(wiso_alpi, 2.0);
  } 
  return wiso_alpi;
};

template <typename Scalar>
Scalar AlphaEqLiquidVapor(const std::string& isosp, const Scalar tk) {

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
    }
    // apply H217O fractionation factor adjustment after Schonemann et al. 2014
    if (isosp == "H217O") {
      wiso_alpl = pow(wiso_alpl, 2.0);
    }
  }
  return wiso_alpl;
};

template <typename Scalar>
void AlphaKineticEvap() {
  // return fractionation factor modified for kinetic effects during
  // liquid evaporation into unsaturated air

};

template <typename Scalar>
void AlphaKineticDepo() {
  // return fractionation factor modified for kinetic effects during
  // deposition to ice. Optional use of vapor supersaturation.
};


/*
* Functions to calculate kinetic fractionation factors
*/

template <typename Scalar>
Scalar AlphaKMol(const std::string& isosp, const Scalar rbot, const Scalar zbot, const Scalar ustar) {

  using Wiso = WaterIsotopes::WaterIsotopologues<Scalar>;
  using Constants = physics::Constants<Scalar>;

  /* compute kinetic modifier for drag coefficient.
  uses Brutsaert equations for the turbulent layer using GCM 
  quantities. This calculates kinetic fractionation using the 
  Merlivat and Jouzel 1979 closure assumption */

  const Scalar difair = 2.36e-5;    // molecular diffusivity of air
  const Scalar muair = 1.7e-5;        // dynamic viscosity of air
  
  // need gravitational constant and karman constant - are these present in share?

  // local variables
  Scalar z0;
  Scalar reno;
  Scalar tmr;
  Scalar sc;
  Scalar vmu; 
  Scalar difn;
  Scalar difrmj;
  Scalar diffpow;
  Scalar kmol = 0.5;
  Scalar alphakn = 1.0;
  
  difrmj = Wiso::difrm[isosp];
  z0 = pow(ustar, 2.0)/(81.1*Constants::gravit); // charnock's equation
  vmu = muair / rbot;
  sc = vmu/difair;
  reno = ustar*z0 / vmu;

  Scalar renocrit = 1.0; // critical reynolds number to distinguish between smooth/rough regimes

  if (reno < renocrit) { // Smooth regime
    diffpow = 2.0/3.0;
    tmr = ((1.0/Constants::Karman)*log(zbot*ustar/(30.0*vmu)))/(13.6 * pow(sc,diffpow));
  } else {
    diffpow = 1.0/2.0;
    tmr = ((1.0/Constants::Karman)*log(zbot/z0)-5.0)/(7.3*pow(reno,0.25)*pow(sc,diffpow));
  }

  difn = pow(1.0/difrmj, diffpow);
  kmol = (difn - 1.0)/(difn + tmr);

  alphakn = 1.0 - kmol;

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