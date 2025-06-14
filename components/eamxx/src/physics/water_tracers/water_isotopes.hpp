#ifndef WATER_ISOTOPES
#define WATER_ISOTOPES

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

// Define water isotopologue species and parameters
template <typename Scalar>
struct WaterIsotopologues {
  // this structure defines the various water isotopologues considered
  static constexpr int isospec  = 6;    // maximum number of water species

   
  static constexpr int isiundef = 0;    // is water species undefined? (Needed?)
  static constexpr int isih2o   = 1;    // is isotopologue = h2o
  static constexpr int isih216o = 2;    // is isotopologue = h216o
  static constexpr int isihdo   = 3;    // is isotopologue = hdo
  static constexpr int isih218o = 4;    // is isotopologue = h218o
  static constexpr int isih217o = 5;    // is isotopologue = h217o
  static constexpr int isihto   = 6;    // is isotopologue = hto
 
  enum class WaterIsotopologue : int { // define a class of water isotopologues
    Undefined = 0,
    H2O = 1,
    H216O = 2,
    HDO = 3,
    H218O = 4,
    H217O = 5,
    HTO = 6
  }; //future development?
   

  static constexpr std::array<Scalar,isospec> fisub = {1.0, 1.0, 2.0, 1.0, 1.0, 2.0}; // not sure what this is?
  static constexpr std::array<Scalar,isospec> mwiso = {18.0, 18.0, 19.0, 20.0, 19.0, 20.0}; // molecular weights of water isotope species
  static constexpr std::array<Scalar,isospec> mwratiso = {1.0, 1.0, 19.0/18.0, 20.0/18.0, 19.0/18.0, 20.0/18.0};   // molecular weight ratios with respect to h216o
  static constexpr std::array<Scalar,isospec> rnat = {1.0, 0.9976, 155.76e-6, 2005.2e-6, 402e-6, 77.88e-6}; // VSMOW ratios
  static constexpr std::array<Scalar,isospec> difrm = {1.0, 1.0, 0.9757, 0.9727, 0.9858, 0.9679}; // Merlivat 1978 + Schoenemann 2014
  static constexpr std::array<Scalar,isospec> rstd = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; // Prescribed isotope ratio wrt H216O (mostly for numerics) 
  static constexpr std::array<Scalar,isospec> boce = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; // Mean ocean surface enrichment relative to VSMOW
// Parameters for kinetic fractionation using Merlivat & Jouzel method
  static constexpr std::array<Scalar,isospec> aksmc = {0.0, 0.0, 0.00528, 0.006, 0.00314, 0.01056}; // not sure! RPF
  static constexpr std::array<Scalar,isospec> akrfa = {0.0, 0.0, 0.2508e-3, 0.285e-3, 0.1495e-3, 0.5016e-3};  // not sure! RPF
  static constexpr std::array<Scalar,isospec> akrfb = {0.0, 0.0, 0.7216e-3, 0.82e-3, 0.430e-3, 1.4432e-3};   // not sure! RPF
// define coefficients for liquid-vapor fractionations
// using polynomial regression. Coefficients are from
// Horita and Wesolowski, 1994. https://doi.org/10.1016/0016-7037(94)90096-5
  static constexpr std::array<Scalar,isospec> alpal = {0.0, 0.0, 1158.8e-12, 0.35041e6, 0.35041e6, 1158.8e-12}; // Liquid/vapor fractionation factor polynomial coefficient
  static constexpr std::array<Scalar,isospec> alpbl = {0.0, 0.0, -1620.1e-9, -1.6664e-3, -1.6664e-3, -1620.1e-9};        // Liquid/vapor fractionation factor polynomial coefficient
  static constexpr std::array<Scalar,isospec> alpcl = {0.0, 0.0, 794.84e-6, 6.7123, 6.7123, 794.84e-6};        // Liquid/vapor fractionation factor polynomial coefficient
  static constexpr std::array<Scalar,isospec> alpdl = {0.0, 0.0, -161.04e-3, -7.683e-3, -7.685e-3, -161.04e-3};        // Liquid/vapor fractionation factor polynomial coefficient
  static constexpr std::array<Scalar,isospec> alpel = {0.0, 0.0, 2.9992e6, 0.0, 0.0, 2.9992e6};     // Liquid/vapor fractionation factor polynomial coefficient
  // define coefficients for ice-vapor fractionation
  // from Merlivat 1978
  static constexpr std::array<Scalar,isospec> alpai = {0.0, 0.0, 16289.0, 0.0, 0.0, 16289.0};       // Ice/vapor fractionation factor polynomial coefficient
  static constexpr std::array<Scalar,isospec> alpbi = {0.0, 0.0, 0.0, 11.839, 11.839, 0.0};      // Ice/vapor fractionation factor polynomial coefficient
  static constexpr std::array<Scalar,isospec> alpci = {0.0, 0.0, -9.45e-2, -28.224e-3, -28.224e-3, -9.45e-2};        // Ice/vapor fractionation factor polynomial coefficient
  // decay rates
  static constexpr Scalar hlhto = 3.888e8;

  // tunable parameters for fractionation scheme
  static constexpr Scalar dkfac = 0.58; // diffusive ivaporartion kinetic power law (Stewart 1975?)
  static constexpr Scalar fkhum = 0.25; // effective humidity factor
  static constexpr Scalar recrit = 1.0; // critical Reynolds number for kmol calculation


};

static const std::vector<std::string> isoname = {"H2O","H216O","HD16O","H218O","H217O","HTO"};


// Functions to calculate equilibrium fractionation factors 

template <typename Scalar>
Scalar AlphaEqIceVapor(const typename WaterIsotopologues<Scalar>::WaterIsotopologue iso, const Scalar& tk) {

  using Wiso = typename WaterIsotopologues<Scalar>::WaterIsotopologue;
  // calculate equilibrium alpha for ice<->vapor transitions
  Scalar wiso_alpi = 1.0;
  if (iso != Wiso::H2O) {
    // Calculate fractionation factors after Merlivat & Nief, 1967 for HDO
    // and Majoube 1971 for H218O and H217O. Need to modify for H217O and HTO.
    wiso_alpi = exp(Wiso::alpai[iso]*pow(tk,-2) +
                Wiso::alpbi[iso]/tk +
                Wiso::alpci[iso]) ;
  } 

  // update for H217O
  if (iso == Wiso::H217O) {
    wiso_alpi = pow(wiso_alpi, 0.528);
  } else if (iso == Wiso::HTO) { // update for HTO
    wiso_alpi = pow(wiso_alpi, 2.0);
  } 
  return wiso_alpi;
};

template <typename Scalar>
Scalar AlphaEqLiquidVapor(const std::string& isosp, const Scalar tk) {

  using Wiso = typename WaterIsotopologues<Scalar>::WaterIsotopologue; 
  // calculate equilibrium alpha for liquid<->vapor transitions
  Scalar wiso_alpl = 1.0;
  if (isosp != Wiso::isih2o) {
    if (isosp == Wiso::isihdo || isosp == Wiso::isihto) { // these have a different structure
      wiso_alpl = exp(Wiso::alpal[isosp]*pow(tk,3) + 
                  Wiso::alpbl[isosp]*pow(tk,2) +
                  Wiso::alpcl[isosp]*tk + 
                  Wiso::alpdl[isosp] + 
                  Wiso::alpel[isosp]*pow(tk,-3));
    } else {
      wiso_alpl = exp(Wiso::alpal[isosp]*pow(tk,-3) + 
                  Wiso::alpbl[isosp]*pow(tk,-2) +
                  Wiso::alpcl[isosp]/tk + 
                  Wiso::alpdl[isosp]);
    }
    // apply H217O fractionation factor adjustment after Schonemann et al. 2014
    if (isosp == Wiso::isih217o) {
      wiso_alpl = pow(wiso_alpl, 2.0);
    }
  }
  return wiso_alpl;
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

template <typename Scalar>
void AlphaKineticEvap() {
  // return fractionation factor modified for kinetic effects during
  // liquid evaporation into unsaturated air

};

void AlphaKineticDepo() {
  // return fractionation factor modified for kinetic effects during
  // deposition to ice. Optional use of vapor supersaturation.
};


void WtypeGetAlpha() {
  /* return fractionation factor given a source water type and 
  a destination water type .
  */
};

/*template <typename Scalar>
Scalar WisoDecay(const Int isosp, const Scalar dtime, const Scalar q) {
  // calculate change in tracer mass due to decay over time step
  Scalar dqdt;
  dqdt = q * pow(0.5, dtime/WaterIsotopologues<Scalar>::hlhto) / dtime;
  return dqdt;
}; */

}
}

#endif // WATER_ISOTOPES
