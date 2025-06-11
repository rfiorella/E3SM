#ifndef WATER_TRACERS
#define WATER_TRACERS

#include "share/eamxx_types.hpp"

#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/logging/ekat_logger.hpp"

#include <string>
#include <vector>
#include <array>

namespace scream {
namespace WaterTracers {

/* 
* Isotopic constants used to determine fractionation factors, etc.
*/

// Define all water types that can be specified
struct WaterTypes {
  // this defines the various types water types that will be queried by other functions
  static constexpr int pwtype = 7;      // number of water types

  static constexpr int iwtundef = 0;    // tracer is not a water type
  static constexpr int iwtvap   = 1;    // water type is vapor (qv)
  static constexpr int iwtliq   = 2;    // water type is cloud liquid
  static constexpr int iwtice   = 3;    // water type is cloud ice
  static constexpr int iwtstrain = 4;   // water type is stratiform rain
  static constexpr int iwtstsnow = 5;   // water type is stratiform snow
  static constexpr int iwtcvrain = 6;   // water type is convective rain
  static constexpr int iwtcvsnow = 7;   // water type is convective snow
};

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

  static const std::vector<std::string> isoname; // names of water isotope species
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

};

static const std::vector<std::string> isoname = {"H2O","H216O","HD16O","H218O","H217O","HTO"};

}
}

#endif // WATER_TRACERS