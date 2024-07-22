#include "physics/p3/eamxx_p3_process_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"
// Needed for p3_init, the only F90 code still used.
#include "physics/p3/p3_functions.hpp"
#include "physics/share/physics_constants.hpp"
#include "physics/p3/p3_f90.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>
#include <vector>
#include <string>

namespace scream
{

// =========================================================================================
// WaterTracers::WaterTracers (const ekat::Comm& comm, const ekat::ParameterList& params)
//   : AtmosphereProcess(comm, params)
// {
//   // Nothing to do here
// }

struct WaterTypes {
  // this defines the various types water types that will be queried by other functions
  static const int pwtype = 7;      // number of water types

  static const int iwtundef = 0;    // tracer is not a water type
  static const int iwtvap   = 1;    // water type is vapor (qv)
  static const int iwtliq   = 2;    // water type is cloud liquid
  static const int iwtice   = 3;    // water type is cloud ice
  static const int iwtstrain = 4;   // water type is stratiform rain
  static const int iwtstsnow = 5;   // water type is stratiform snow
  static const int iwtcvrain = 6;   // water type is convective rain
  static const int iwtcvsnow = 7;   // water type is convective snow
};

struct WaterIsotopologues {
  // this structure defines the various water isotopologues considered
  static const int isospec  = 6;    // maximum number of water species
  static const int isiundef = 0;    // is water species undefined? (Needed?)
  static const int isih2o   = 1;    // is isotopologue = h2o
  static const int isih216o = 2;    // is isotopologue = h216o
  static const int isihdo   = 3;    // is isotopologue = hdo
  static const int isih218o = 4;    // is isotopologue = h218o
  static const int isih217o = 5;    // is isotopologue = h2o
  static const int isihto   = 6;    // is isotopologue = h2o

  static const std::vector<std::string> isoname; // names of water isotope species
  static const std::vector<double> fisub;        // not sure what this is?
  static const std::vector<double> mwiso;        // molecular weights of water isotope species
  static const std::vector<double> mwratiso;     // molecular weight ratios with respect to h216o
  static const std::vector<double> rnat;         // heavy-to-light isotope ratios for species / SMOW
  static const std::vector<double> difrm;        // atomic diffusivity ratios (e.g., D/H, not HDO/H216O)
  static const std::vector<double> rstd;         // Prescribed isotope ratio wrt H216O (mostly for numerics) 
  static const std::vector<double> boce;         // Mean ocean surface enrichment relative to VSMOW
  static const std::vector<double> aksmc;        // not sure! RPF
  static const std::vector<double> akrfa;        // not sure! RPF
  static const std::vector<double> akrfb;        // not sure! RPF
  static const std::vector<double> alpal;        // Liquid/vapor fractionation factor polynomial coefficient
  static const std::vector<double> alpbl;        // Liquid/vapor fractionation factor polynomial coefficient
  static const std::vector<double> alpcl;        // Liquid/vapor fractionation factor polynomial coefficient
  static const std::vector<double> alpdl;        // Liquid/vapor fractionation factor polynomial coefficient
  static const std::vector<double> alpel;        // Liquid/vapor fractionation factor polynomial coefficient
  static const std::vector<double> alpai;        // Ice/vapor fractionation factor polynomial coefficient
  static const std::vector<double> alpbi;        // Ice/vapor fractionation factor polynomial coefficient
  static const std::vector<double> alpci;        // Ice/vapor fractionation factor polynomial coefficient

};

// now fill the vectors outside of the WaterIsotoplogues struct definition:
const std::vector<std::string> WaterIsotopologues::isoname = {"H2O","H216O","HD16O","H218O","H217O","HTO"};
const std::vector<double> WaterIsotopologues::fisub = {1.0, 1.0, 2.0, 1.0, 1.0, 2.0};
const std::vector<double> WaterIsotopologues::mwiso = {18.0, 18.0, 19.0, 20.0, 19.0, 20.0}; // correct to 3 sig figs.
const std::vector<double> WaterIsotopologues::mwratiso = {1.0, 1.0, 19.0/18.0, 20.0/18.0, 19.0/18.0, 20.0/18.0};
const std::vector<double> WaterIsotopologues::rnat  = {1.0, 0.9976, 155.76e-6, 2005.2e-6, 402e-6, 77.88e-6}; // VSMOW ratios
const std::vector<double> WaterIsotopologues::difrm = {1.0, 1.0, 0.9757, 0.9727, 0.9858, 0.9679}; // Merlivat 1978 + Schoenemann 2014
const std::vector<double> WaterIsotopologues::rstd  = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
const std::vector<double> WaterIsotopologues::boce  = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
// Parameters for kinetic fractionation using Merlivat & Jouzel method
const std::vector<double> WaterIsotopologues::aksmc = {0.0, 0.0, 0.00528, 0.006, 0.00314, 0.01056};
const std::vector<double> WaterIsotopologues::akrfa = {0.0, 0.0, 0.2508e-3, 0.285e-3, 0.1495e-3, 0.5016e-3};
const std::vector<double> WaterIsotopologues::akrfb = {0.0, 0.0, 0.7216e-3, 0.82e-3, 0.430e-3, 1.4432e-3}; 
// define coefficients for liquid-vapor fractionations
// using polynomial regression. Coefficients are from
// Horita and Wesolowski, 1994. https://doi.org/10.1016/0016-7037(94)90096-5
const std::vector<double> WaterIsotopologues::alpal = {0.0, 0.0, 1158.8e-12, 0.35041e6, 0.35041e6, 1158.8e-12}; // not sure where the 17O and hto values come from
const std::vector<double> WaterIsotopologues::alpbl = {0.0, 0.0, -1620.1e-9, -1.6664e-3, -1.6664e-3, -1620.1e-9};
const std::vector<double> WaterIsotopologues::alpcl = {0.0, 0.0, 794.84e-6, 6.7123, 6.7123, 794.84e-6};
const std::vector<double> WaterIsotopologues::alpdl = {0.0, 0.0, -161.04e-3, -7.683e-3, -7.685e-3, -161.04e-3};
const std::vector<double> WaterIsotopologues::alpel = {0.0, 0.0, 2.9992e6, 0.0, 0.0, 2.9992e6};
// define coefficients for ice-vapor fractionation
// from Merlivat 1978
const std::vector<double> WaterIsotopologues::alpai = {0.0, 0.0, 16289.0, 0.0, 0.0, 16289.0};
const std::vector<double> WaterIsotopologues::alpbi = {0.0, 0.0, 0.0, 11.839, 11.839, 0.0};
const std::vector<double> WaterIsotopologues::alpci = {0.0, 0.0, -9.45e-2, -28.224e-3, -28.224e-3, -9.45e-2};



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

void AlphaEqIceVapor() {
  // calculate equilibrium alpha for ice<->vapor transitions

};

void AlphaEqLiquidVapor(Int isosp, Real tk) {
  real wiso_alpl = 1.0
  // calculate equilibrium alpha for liquid<->vapor transitions
  if (isosp != WaterIsotopologues::isih2o) {
    if (isosp == WaterIsotopologues::isihdo || isosp == WaterIsotpologues::isihto) { // these have a different structure
      wiso_alpl = exp(WaterIsotopologues::alpal[isosp]*pow(tk,3) + 
                  WaterIsotopologues::alpbl[isosp]*pow(tk,2) +
                  WaterIsotopologues::alpcl[isosp]*tk + 
                  WaterIsotopologues::alpdl[isosp] + 
                  WaterIsotopologues::alpel[isosp]/pow(tk,3))
    } else {
      wiso_alpl = exp(WaterIsotopologues::alpal[isosp]/pow(tk,3) + 
                  WaterIsotpologues::alpbl[isosp]/pow(tk,2) +
                  WaterIsotoplogues::alpcl[isosp]/tk + 
                  WaterIsotopologues::alpdl[isosp])
    }
    // apply H217O fractionation factor adjustment after Schonemann et al. 2014
    if (isosp == WaterIsotopologues::isih217o) {
      wiso_alpl = pow(wiso_alpl, 2.0)
    }
  }

  return wiso_alpl
};

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

void AlphaKMol() {
  // compute kinetic modifier for drag coefficient.
  // uses Brutsaert equations for the turbulent layer using GCM 
  // quantities.

  const real difair = 2.36e-5;    // molecular diffusivity of air
  const real muair = 1.7e-5;        // dynamic viscosity of air
  // need gravitational constant and karman constant - are these present in share?

};

}
