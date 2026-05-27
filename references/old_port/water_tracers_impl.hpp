#ifndef SCREAM_WATER_TRACERS_IMPL_HPP
#define SCREAM_WATER_TRACERS_IMPL_HPP

//======================================================================================
// water_tracers_impl.hpp
//
// Implementation details for water_tracers.hpp
// Contains template implementations for device-compatible physics functions
//
// Original Fortran: water_tracers.F90 (7,813 lines)
// C++ Port: Rich Fiorella (2026)
//
//======================================================================================

#include "water_tracers.hpp"
#include <cmath>

namespace scream {
namespace water_tracers {

//======================================================================================
// GROUP 5: RATIO CALCULATION (Most heavily used function - ~80+ calls)
//======================================================================================

// Calculate isotope ratio from masses with numerical checks
// Fortran: wtrc_ratio (lines 6966-6995)
// Called extensively throughout the code for ratio calculations
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT wtrc_ratio(int ispec,
                   const ScalarT& qtrc,
                   const ScalarT& qtot)
{
  // Purpose: Compute tracer ratio from masses with numerical checks
  // Author: David Noone <dcn@colorado.edu> - Sat Jul  3 18:52:40 MDT 2004
  //
  // Args:
  //   ispec - isotopic species index
  //   qtrc  - tracer water mass (kg/kg)
  //   qtot  - "base" water mass (kg/kg)
  //
  // Returns:
  //   Isotope ratio (dimensionless)
  //
  // Note: Uses wtrc_qmin threshold to avoid division by very small numbers.
  //       When qtot is below threshold, returns standard ratio for the species.
  
  using ekat::impl::abs;
  
  constexpr ScalarT zero = 0.0;
  
  // Check if denominator is below minimum threshold
  if (abs(qtot) < wtrc_qmin) {
    // Return standard isotopic ratio for this species
    return wtrc_get_rstd(ispec);
  } else {
    // Calculate ratio directly
    return qtrc / qtot;
  }
}

//======================================================================================
// GROUP 4: FRACTIONATION & EQUILIBRIUM
//======================================================================================

// Get fractionation factor for phase transition
// Fortran: wtrc_get_alpha (lines 7042-7140)
// Wrapper function that combines functionality from water_isotopes and water_types
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT wtrc_get_alpha(const ScalarT& qvap,
                       const ScalarT& temp,
                       int ispec,
                       int isrctype,
                       int idsttype,
                       bool do_kinetic,
                       const ScalarT& pmid,
                       bool use_ice_supsat)
{
  // Purpose: Retrieve the fractionation for a process that goes from
  //          the source water type to the destination water type.
  //
  // Author: Chuck Bardeen
  // Modified for RH use by Jesse Nusbaumer
  // Modified for RH-ice (RHi) use by Marina Dutsch
  //
  // Args:
  //   qvap - water vapor mixing ratio (kg/kg)
  //   temp - temperature (K)
  //   ispec - isotope species index  
  //   isrctype - source water type index
  //   idsttype - destination water type index
  //   do_kinetic - use kinetic fractionation?
  //   pmid - pressure (Pa)
  //   use_ice_supsat - use ice supersaturation?
  //
  // Returns:
  //   Fractionation factor (dimensionless)
  
  using namespace water_isotopes;
  using namespace water_types;
  
  constexpr ScalarT zero = 0.0;
  constexpr ScalarT one = 1.0;
  
  // Determine if kinetic fractionation should be used
  bool alpkin = wtrc_alpha_kinetic || do_kinetic;
  
  // Check if isotopes are enabled
  if (!wisotope) {
    // Return fixed fractionation factor
    return wtrc_fixed_alpha[ispec];
  }
  
  // Calculate saturation specific humidity
  // Note: This requires access to saturation vapor pressure functions
  // For now, we'll use a simplified approach - full implementation would
  // call qsat_water or qsat_ice from wv_saturation module
  
  ScalarT qs = zero;  // Placeholder - needs proper saturation calculation
  ScalarT rh = zero;  // Relative humidity
  ScalarT rhi = zero; // Relative humidity w.r.t. ice
  
  // TODO: Implement proper saturation calculations
  // For full implementation, need to call:
  // - qsat_water(temp, pmid, es, qs) for liquid saturation
  // - qsat_ice(temp, pmid, esi, qsi) for ice saturation
  // - wtrc_eff_sat for effective saturation
  
  if (use_ice_supsat) {
    // Calculate effective RH with respect to ice
    // rhi = wtrc_eff_sat(qvap, temp, pmid, true);
    rhi = one;  // Placeholder
  }
  
  // Calculate RH with respect to liquid
  if (qs > zero) {
    rh = qvap / qs;
  } else {
    rh = one;
  }
  
  // Call water_types function to get fractionation factor
  return wtype_get_alpha(ispec, isrctype, idsttype, temp, rh, alpkin, rhi);
}

// Equilibrate liquid and vapor isotopes
// Fortran: wtrc_liqvap_equil (lines 6289-6341)
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
void wtrc_liqvap_equil(const ScalarT& alpha,
                       const ScalarT& feq0,
                       const ScalarT& vaptot,
                       const ScalarT& liqtot,
                       ScalarT& vapiso,
                       ScalarT& liqiso,
                       ScalarT& dliqiso)
{
  // Purpose: Equilibrate isotopic composition between liquid and vapor phases
  //
  // Algorithm based on instantaneous reservoir coupling with optional
  // fractional equilibration (feq0 parameter).
  //
  // Args:
  //   alpha - equilibrium fractionation factor (dimensionless)
  //   feq0 - fraction of equilibration to apply (0-1)
  //   vaptot - total vapor mass (H2O) (kg/kg)
  //   liqtot - total liquid mass (H2O) (kg/kg)
  //   vapiso - isotope vapor mass (kg/kg) [IN/OUT]
  //   liqiso - isotope liquid mass (kg/kg) [IN/OUT]
  //   dliqiso - change in liquid isotope mass (kg/kg) [OUT]
  //
  // Note: Uses implicit solution to avoid instabilities
  
  using ekat::impl::abs;
  
  constexpr ScalarT zero = 0.0;
  constexpr ScalarT one = 1.0;
  constexpr ScalarT tiny = 1.e-20;
  
  // Calculate equilibration factor
  // Based on wtrc_efac function (lines 6721-6748)
  ScalarT feq = zero;
  ScalarT denom = liqtot + alpha * vaptot;
  
  if (abs(denom) > tiny) {
    feq = (liqtot * alpha * vaptot) / denom;
  }
  
  // Apply fractional equilibration
  feq = feq0 * feq;
  
  // Calculate change needed for equilibration
  // Based on wtrc_dqequil function (lines 6344-6387)
  ScalarT vapnew = vaptot;
  ScalarT liqnew = liqtot;
  ScalarT visoold = vapiso;
  ScalarT lisoold = liqiso;
  
  // Total isotope mass (conserved)
  ScalarT isotot = vapiso + liqiso;
  
  // Calculate equilibrated isotope distribution
  ScalarT liqeq = zero;
  if (abs(denom) > tiny) {
    liqeq = (liqnew * isotot + feq * alpha * (liqnew * visoold - vapnew * lisoold)) / denom;
  }
  
  // Limit to physically reasonable values
  liqeq = ekat::impl::max(zero, ekat::impl::min(isotot, liqeq));
  
  // Calculate change in liquid isotope
  dliqiso = liqeq - lisoold;
  
  // Update isotope masses
  liqiso = lisoold + dliqiso;
  vapiso = visoold - dliqiso;
}

// Distill vapor for some increment
// Fortran: wtrc_vap_distil (lines 6390-6422)
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
void wtrc_vap_distil(const ScalarT& alpha,
                     const ScalarT& vtotold,
                     const ScalarT& vtotnew,
                     const ScalarT& visoold,
                     ScalarT& visonew,
                     ScalarT& dvapiso)
{
  // Purpose: Apply Rayleigh distillation to vapor phase
  //
  // Uses classic Rayleigh distillation formula:
  //   R_new = R_old * (f^alpha)
  // where f is the fraction of vapor remaining
  //
  // Args:
  //   alpha - fractionation factor (dimensionless)
  //   vtotold - initial total vapor mass (kg/kg)
  //   vtotnew - final total vapor mass (kg/kg)
  //   visoold - initial isotope vapor mass (kg/kg)
  //   visonew - final isotope vapor mass (kg/kg) [OUT]
  //   dvapiso - change in vapor isotope mass (kg/kg) [OUT]
  
  using ekat::impl::pow;
  using ekat::impl::max;
  using ekat::impl::min;
  
  constexpr ScalarT zero = 0.0;
  constexpr ScalarT one = 1.0;
  constexpr ScalarT tiny = 1.e-20;
  
  // Calculate fraction of vapor remaining
  ScalarT frac = one;
  if (vtotold > tiny) {
    frac = vtotnew / vtotold;
  }
  
  // Constrain fraction to physical range [0, 1]
  frac = min(one, max(zero, frac));
  
  // Apply Rayleigh fractionation
  visonew = visoold * pow(frac, alpha);
  
  // Calculate change
  dvapiso = visonew - visoold;
}

// Equilibrium implicit factor
// Fortran: wtrc_efac (lines 6721-6748)
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT wtrc_efac(const ScalarT& alpha,
                  const ScalarT& vapnew,
                  const ScalarT& liqnew)
{
  // Purpose: Calculate equilibrium factor for implicit time integration
  //
  // This factor represents the effective equilibration timescale
  // for isotopic exchange between liquid and vapor phases.
  //
  // Args:
  //   alpha - equilibrium fractionation factor
  //   vapnew - vapor mass (kg/kg)
  //   liqnew - liquid mass (kg/kg)
  //
  // Returns:
  //   Equilibrium factor (dimensionless)
  
  using ekat::impl::abs;
  
  constexpr ScalarT zero = 0.0;
  constexpr ScalarT tiny = 1.e-20;
  
  ScalarT denom = liqnew + alpha * vapnew;
  
  if (abs(denom) > tiny) {
    return (liqnew * alpha * vapnew) / denom;
  } else {
    return zero;
  }
}

// Calculate equilibrium change for distillation
// Fortran: wtrc_dqequil (lines 6344-6387)
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT wtrc_dqequil(const ScalarT& alpha,
                     const ScalarT& feq0,
                     const ScalarT& vtotnew,
                     const ScalarT& ltotnew,
                     const ScalarT& visoold,
                     const ScalarT& lisoold)
{
  // Purpose: Calculate change in isotope distribution for equilibration
  //
  // Args:
  //   alpha - equilibrium fractionation factor
  //   feq0 - fraction of equilibration (0-1)
  //   vtotnew - new total vapor mass (kg/kg)
  //   ltotnew - new total liquid mass (kg/kg)
  //   visoold - old isotope vapor mass (kg/kg)
  //   lisoold - old isotope liquid mass (kg/kg)
  //
  // Returns:
  //   Change in liquid isotope mass (kg/kg)
  
  using ekat::impl::abs;
  using ekat::impl::max;
  using ekat::impl::min;
  
  constexpr ScalarT zero = 0.0;
  constexpr ScalarT tiny = 1.e-20;
  
  // Total isotope mass (conserved)
  ScalarT isotot = visoold + lisoold;
  
  // Calculate equilibrated distribution
  ScalarT denom = ltotnew + alpha * vtotnew;
  ScalarT liqeq = zero;
  
  if (abs(denom) > tiny) {
    // Calculate equilibrium factor
    ScalarT feq = wtrc_efac(alpha, vtotnew, ltotnew);
    feq = feq0 * feq;
    
    // Calculate equilibrated liquid isotope mass
    liqeq = (ltotnew * isotot + feq * alpha * (ltotnew * visoold - vtotnew * lisoold)) / denom;
    
    // Limit to physical range
    liqeq = max(zero, min(isotot, liqeq));
  }
  
  // Return change in liquid isotope
  return liqeq - lisoold;
}

//======================================================================================
// GROUP 3: RATE CALCULATIONS
//======================================================================================

// Add single process rate to matrix
// Fortran: wtrc_add_rate (lines 901-943)
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
void wtrc_add_rate(View5D<ScalarT>& process_rates,
                   int icol, int iz,
                   int isrctype, int idsttype, int rtype,
                   const ScalarT& rate,
                   bool do_reverse)
{
  // Purpose: Add process rate to matrix
  //
  // Method: Add a process rate for the production of idsttype from isrctype
  //         to a matrix of process rates. This routine will add the reverse
  //         rate for the loss term for isrctype to idsttype if do_reverse is true.
  //
  // Args:
  //   process_rates - 5D array of process rates [col, lev, dst, src, rtype] (kg/kg/s)
  //   icol - column index
  //   iz - vertical level index
  //   isrctype - source water type index
  //   idsttype - destination water type index
  //   rtype - water type used to calculate ratio
  //   rate - process rate for src -> dst (kg/kg/s)
  //   do_reverse - apply reverse rates too (default=true)
  
  // Add forward rate: source -> destination
  process_rates(icol, iz, idsttype, isrctype, rtype) += rate;
  
  // Add reverse rate: destination <- source (loss from source)
  if ((isrctype != idsttype) && do_reverse) {
    process_rates(icol, iz, isrctype, idsttype, rtype) -= rate;
  }
}

// Add process rates array to matrix
// Fortran: wtrc_add_rates (lines 947-991)
template<typename ScalarT>
void wtrc_add_rates(View5D<ScalarT>& process_rates,
                    int ncol, int top_lev,
                    int isrctype, int idsttype, int rtype,
                    const View2D<ScalarT>& rate,
                    bool do_reverse,
                    const View2D<ScalarT>* wtfri)
{
  // Purpose: Add process rates array to matrix
  //
  // Method: Add process rates for the production of idsttype from isrctype
  //         to a matrix of process rates. Loops over columns and levels,
  //         calling wtrc_add_rate for each grid point. Optionally handles
  //         special case where freezing rain goes to cloud ice instead of snow.
  //
  // Args:
  //   process_rates - 5D array of process rates [col, lev, dst, src, rtype] (kg/kg/s)
  //   ncol - number of columns
  //   top_lev - top vertical level
  //   isrctype - source water type index
  //   idsttype - destination water type index
  //   rtype - water type used to calculate ratio
  //   rate - 2D array of process rates [col, lev] (kg/kg/s) for src -> dst
  //   do_reverse - apply reverse rates too (default=true)
  //   wtfri - optional 2D array indicating freezing rain to cloud ice [col, lev]
  
  using namespace water_types;
  constexpr ScalarT zero = 0.0;
  
  // Handle special case for rain freezing to cloud ice (MG2 microphysics)
  if (wtfri != nullptr) {
    // Loop with special freezing rain logic
    auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
      {0, top_lev}, {ncol, rate.extent(1)}
    );
    
    Kokkos::parallel_for("wtrc_add_rates_with_wtfri", policy,
      KOKKOS_LAMBDA(int icol, int k) {
        if ((*wtfri)(icol, k) > zero) {
          // Freezing rain goes to cloud ice instead of snow
          wtrc_add_rate(process_rates, icol, k, isrctype, 
                        static_cast<int>(WaterType::CloudIce), rtype, 
                        rate(icol, k), do_reverse);
        } else {
          // Normal destination
          wtrc_add_rate(process_rates, icol, k, isrctype, idsttype, rtype,
                        rate(icol, k), do_reverse);
        }
      }
    );
  } else {
    // Standard loop without special freezing logic
    auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
      {0, top_lev}, {ncol, rate.extent(1)}
    );
    
    Kokkos::parallel_for("wtrc_add_rates", policy,
      KOKKOS_LAMBDA(int icol, int k) {
        wtrc_add_rate(process_rates, icol, k, isrctype, idsttype, rtype,
                      rate(icol, k), do_reverse);
      }
    );
  }
}

// Initialize process rate matrix
// Fortran: wtrc_init_rates (lines 995-1016)
template<typename ScalarT>
void wtrc_init_rates(int top_lev, View5D<ScalarT>& process_rates)
{
  // Purpose: Initialize process rate matrix
  //
  // Method: Initializes the matrix of process rates to zero
  //
  // Args:
  //   top_lev - top vertical level to initialize from
  //   process_rates - 5D array of process rates (kg/kg/sec)
  
  constexpr ScalarT zero = 0.0;
  
  // Zero out rates from top_lev downward
  // Note: Full implementation would use Kokkos parallel_for
  // For now, provide the structure
  
  auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<5>>(
    {0, top_lev, 0, 0, 0},
    {process_rates.extent(0), process_rates.extent(1), 
     process_rates.extent(2), process_rates.extent(3), process_rates.extent(4)}
  );
  
  Kokkos::parallel_for("wtrc_init_rates", policy,
    KOKKOS_LAMBDA(int icol, int iz, int idst, int isrc, int irtype) {
      process_rates(icol, iz, idst, isrc, irtype) = zero;
    }
  );
}

//======================================================================================
// GROUP 11: MASS FIXING & UTILITIES
//======================================================================================

// Calculate fraction of rain that has experienced equilibration
// Fortran: wtrc_equil_time (lines 6548-6718)
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
void wtrc_equil_time(int ispec,
                     const ScalarT& temp,
                     const ScalarT& pres,
                     const ScalarT& rdrop,
                     const ScalarT& zdel,
                     const ScalarT& alpha,
                     const ScalarT& difrm,
                     ScalarT& fequil)
{
  // Purpose: Calculate the fraction of rain that has experienced equilibration
  //          with surrounding vapor during its fall through a model layer.
  //
  // Method: Based on Stewart (1975) equation for isotopic exchange during
  //         raindrop fall. Accounts for drop size, fall distance, and diffusion.
  //
  // Args:
  //   ispec - isotope species index
  //   temp - temperature (K)
  //   pres - pressure (Pa)
  //   rdrop - raindrop radius (m)
  //   zdel - layer thickness (m)
  //   alpha - equilibrium fractionation factor
  //   difrm - isotopic diffusivity (m^2/s)
  //   fequil - fraction equilibrated (0-1) [OUT]
  //
  // Reference: Stewart, M. K. (1975). Stable isotope fractionation due to
  //            evaporation and isotopic exchange of falling waterdrops.
  //            J. Geophys. Res., 80(9), 1133-1146.
  
  using ekat::impl::sqrt;
  using ekat::impl::exp;
  using ekat::impl::min;
  using ekat::impl::max;
  
  constexpr ScalarT zero = 0.0;
  constexpr ScalarT one = 1.0;
  constexpr ScalarT pi = 3.14159265358979323846;
  constexpr ScalarT tiny = 1.e-20;
  
  // Fall velocity for raindrop (m/s)
  // Using empirical formula: v = 9.58 * (1 - exp(-1180*r^1.147))
  ScalarT vfall = 9.58 * (one - exp(-1180.0 * pow(rdrop, 1.147)));
  
  // Time for drop to fall through layer (s)
  ScalarT tfall = zero;
  if (vfall > tiny) {
    tfall = zdel / vfall;
  }
  
  // Characteristic equilibration time (s)
  // tau = (r^2) / (3 * D * alpha)
  ScalarT tau = zero;
  if (difrm > tiny && alpha > tiny) {
    tau = (rdrop * rdrop) / (3.0 * difrm * alpha);
  }
  
  // Fraction equilibrated = 1 - exp(-tfall/tau)
  fequil = zero;
  if (tau > tiny) {
    fequil = one - exp(-tfall / tau);
  }
  
  // Limit to physical range [0, 1]
  fequil = min(one, max(zero, fequil));
}

// Get R_alpha for O2 fractionation (used in chemistry)
// Fortran: wtrc_get_rao2 (lines 7242-7267)
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT wtrc_get_rao2(int isp)
{
  // Purpose: Get fractionation factor for O2 reactions
  //          Used in methane oxidation chemistry
  //
  // Args:
  //   isp - isotope species index
  //
  // Returns:
  //   O2 fractionation factor (dimensionless)
  
  using namespace water_isotopes;
  
  constexpr ScalarT one = 1.0;
  
  // Values from Luz and Barkan (2010)
  // Different for different oxygen isotopes
  
  if (isp == static_cast<int>(IsoSpecies::H218O) - 1) {
    // 18O fractionation
    return 0.9967;  // Approximately 1 - 3.3 permil
  } else if (isp == static_cast<int>(IsoSpecies::H217O) - 1) {
    // 17O fractionation (mass-dependent scaling)
    return 0.9984;  // Approximately 1 - 1.6 permil
  } else {
    // No fractionation for other species (H2O, HDO, HTO)
    return one;
  }
}

} // namespace water_tracers
} // namespace scream

#endif // SCREAM_WATER_TRACERS_IMPL_HPP
