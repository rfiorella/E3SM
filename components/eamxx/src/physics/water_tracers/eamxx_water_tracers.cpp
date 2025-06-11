#include "water_tracers.hpp"

#include <vector>
#include <string>

namespace scream { 
namespace WaterTracers { 
// =========================================================================================
// WaterTracers::WaterTracers (const ekat::Comm& comm, const ekat::ParameterList& params)
//   : AtmosphereProcess(comm, params)
// {
//   // Nothing to do here
// }

// void WtypeGetItype(name) {
//   // Retrieve WaterType based on type name
// };


void WtypeGetAlpha() {
  // retrieve fractionation factor for process that converts
  // source water type (isrctype) to desination water type (idsttype)

};

// // Set of functions to determine the type of water tracer
// bool trc_is_wtrc(tracer_num) {
//   // returns true if tracer is a water tracer
// };

template <typename Scalar>
Real AlphaEqIceVapor(Int isosp, Real tk) {

  using Wiso = WaterIsotopologues<Scalar>;
  // calculate equilibrium alpha for ice<->vapor transitions
  Real wiso_alpi = 1.0;
  if (isosp != Wiso::isih2o) {
    // Calculate fractionation factors after Merlivat & Nief, 1967 for HDO
    // and Majoube 1971 for H218O and H217O. Need to modify for H217O and HTO.
    wiso_alpi = exp(Wiso::alpai[isosp]/pow(tk,2) +
                Wiso::alpbi[isosp]/tk +
                Wiso::alpci[isosp]) ;
  } 

  // update for H217O
  if (isosp == Wiso::isih217o) {
    wiso_alpi = pow(wiso_alpi, 0.528);
  } else if (isosp == Wiso::isihto) { // update for HTO
    wiso_alpi = pow(wiso_alpi, 2.0);
  } 
  return wiso_alpi;
};

template <typename Scalar>
Real AlphaEqLiquidVapor(Int isosp, Real tk) {

  using Wiso = WaterIsotopologues<Scalar>;
  // calculate equilibrium alpha for liquid<->vapor transitions
  Real wiso_alpl = 1.0;
  if (isosp != Wiso::isih2o) {
    if (isosp == Wiso::isihdo || isosp == Wiso::isihto) { // these have a different structure
      wiso_alpl = exp(Wiso::alpal[isosp]*pow(tk,3) + 
                  Wiso::alpbl[isosp]*pow(tk,2) +
                  Wiso::alpcl[isosp]*tk + 
                  Wiso::alpdl[isosp] + 
                  Wiso::alpel[isosp]/pow(tk,3));
    } else {
      wiso_alpl = exp(Wiso::alpal[isosp]/pow(tk,3) + 
                  Wiso::alpbl[isosp]/pow(tk,2) +
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

template <typename Scalar>
void AlphaKineticEvap() {
  // return fractionation factor modified for kinetic effects during
  // liquid evaporation into unsaturated air

};

void AlphaKineticDepo() {
  // return fractionation factor modified for kinetic effects during
  // deposition to ice. Optional use of vapor supersaturation.
};

void OceanTracerFlux() {
  // compute water tracer exchange from ocean
  // see L. 531-635 of wiso_flxoce in CAM
};

template <typename Scalar>
Real AlphaKMol(Int isosp, Real rbot, Real zbot, Real ustar) {

  using Wiso = WaterIsotopologues<Scalar>;
  // compute kinetic modifier for drag coefficient.
  // uses Brutsaert equations for the turbulent layer using GCM 
  // quantities.

  const Real difair = 2.36e-5;    // molecular diffusivity of air
  const Real muair = 1.7e-5;        // dynamic viscosity of air
  const Real gravit = 9.81;       // gravity constant
  const Real karman = 0.6;        // von Karman constant
  // need gravitational constant and karman constant - are these present in share?

  // local variables
  Real z0;
  Real reno;
  Real tmr;
  Real sc;
  Real vmu; 
  Real difn;
  Real difrmj;
  Real diffpow;
  Real kmol = 1.0;
  Real alphakn = 1.0;
  
  difrmj = Wiso::difrm[isosp];
  z0 = pow(ustar, 2.0)/(81.1*gravit); // charnock's equation
  vmu = muair / rbot;
  sc = vmu/difair;
  reno = ustar*z0 / vmu;

  Real renocrit = 1.0; // critical reynolds number to distinguish between smooth/rough regimes

  if (reno < renocrit) { // Smooth regime
    diffpow = 2.0/3.0;
    tmr = ((1.0/karman)*log(zbot*ustar/(30.0*vmu)))/(13.6 * pow(sc,diffpow));
  } else {
    diffpow = 1.0/2.0;
    tmr = ((1.0/karman)*log(zbot/z0)-5.0)/(7.3*pow(reno,0.25)*pow(sc,diffpow));
  }

  difn = pow(1.0/difrmj, diffpow);
  kmol = (difn - 1.0)/(difn + tmr);

  alphakn = 1.0 - kmol;

  return alphakn;

};

}
}
