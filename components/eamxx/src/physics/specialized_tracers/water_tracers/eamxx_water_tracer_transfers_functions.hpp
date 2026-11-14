#ifndef EAMXX_WATER_TRACER_TRANSFERS_FUNCTIONS_HPP
#define EAMXX_WATER_TRACER_TRANSFERS_FUNCTIONS_HPP

#include "eamxx_water_tracers_functions.hpp"  // tracer_ratio

namespace scream {
namespace water_tracers {

/*
 * Device-callable transfers of the water_tracer_transfers process, one
 * function per class of bulk water transfer, so that each class can get its
 * own fractionation factor later.
 *
 * Each function takes the donor phase's tracer and bulk water, and the bulk
 * mass moved out of the donor phase, and returns the tracer mass moved
 * [kg/kg]. Same ScalarT convention as eamxx_water_tracers_functions.hpp (an
 * ekat::Pack; tracer_ratio masks dry donors, so these are FPE-safe).
 */

// -----------------------------------------------------------------------
// SHOC's net condensation `amount`: > 0 moves vapor to cloud liquid (the
// donor is vapor), < 0 moves cloud liquid to vapor (the donor is cloud
// liquid, and the result has the sign of amount too).
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT shoc_condensation (const ScalarT& tr_donor, const ScalarT& q_donor, const ScalarT& amount)
{
  // alpha = 1: the tracer moves at the donor's ratio. This is where the
  // liquid/vapor fractionation factor of SHOC's phase change goes.
  return amount*tracer_ratio(tr_donor, q_donor);
}

// -----------------------------------------------------------------------
// P3 process rates (the Part2* components of p3_water_transfers). The donor
// is the mass available to the donor phase over the P3 call (its post-part1
// reservoir plus its part2 inflows; see replay_p3), and amount = rate*dt >= 0
// is the bulk mass moved [kg/kg].
// -----------------------------------------------------------------------

// Vapor deposition and deposition/condensation-freezing nucleation (v -> i)
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT ice_deposition (const ScalarT& tr_v, const ScalarT& q_v, const ScalarT& amount)
{
  // alpha = 1. This is where the vapor/ice fractionation factor goes, with
  // its kinetic effect in supersaturation (Jouzel & Merlivat 1984).
  return amount*tracer_ratio(tr_v, q_v);
}

// Sublimation (i -> v)
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT ice_sublimation (const ScalarT& tr_i, const ScalarT& q_i, const ScalarT& amount)
{
  // alpha = 1. This is where a sublimation factor goes, if any.
  return amount*tracer_ratio(tr_i, q_i);
}

// Rain evaporation (r -> v)
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT rain_evaporation (const ScalarT& tr_r, const ScalarT& q_r, const ScalarT& amount)
{
  // alpha = 1. This is where the liquid/vapor fractionation factor goes, with
  // the kinetic effect of evaporation into subsaturated air (Stewart 1975).
  return amount*tracer_ratio(tr_r, q_r);
}

// Accretion, autoconversion and shedding (c -> r)
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT liquid_collection (const ScalarT& tr_c, const ScalarT& q_c, const ScalarT& amount)
{
  // alpha = 1. No phase change, so no factor is expected here.
  return amount*tracer_ratio(tr_c, q_c);
}

// Bergeron process (c -> i): cloud liquid evaporates and deposits on ice
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT bergeron (const ScalarT& tr_c, const ScalarT& q_c, const ScalarT& amount)
{
  // alpha = 1. The mass passes through vapor, so this is where the factor of
  // that path goes (likely the deposition one).
  return amount*tracer_ratio(tr_c, q_c);
}

// Riming, collection of rain by ice, and immersion freezing (c -> i, r -> i)
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT liquid_freezing (const ScalarT& tr_donor, const ScalarT& q_donor, const ScalarT& amount)
{
  // alpha = 1. This is where a freezing factor goes, if any.
  return amount*tracer_ratio(tr_donor, q_donor);
}

// Melting (i -> r)
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT ice_melting (const ScalarT& tr_i, const ScalarT& q_i, const ScalarT& amount)
{
  // alpha = 1. This is where a melting factor goes, if any.
  return amount*tracer_ratio(tr_i, q_i);
}

// -----------------------------------------------------------------------
// Clip-type stages: P3's clips of small masses and its homogeneous
// freezing move a whole bulk reservoir, `amount` (which can be negative,
// where P3 clips a small negative mass). The tracer's whole donor reservoir
// moves with it, wherever amount != 0. Where amount == 0 nothing moves: P3
// did not run that stage there, or the bulk reservoir was exactly 0 (the
// tracer's is then 0 too, up to roundoff of earlier transfers).
//
// Moving a whole reservoir cannot fractionate.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT whole_reservoir (const ScalarT& tr_donor, const ScalarT& amount)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  ScalarT moved(0);
  moved.set(amount != RealT(0), tr_donor);
  return moved;
}

// -----------------------------------------------------------------------
// Sedimentation: the tracer flux out of the bottom of a level [kg/m2/s].
//
// The bulk flux_out leaves at the ratio of the mass that was in the level
// during the call: its own tr and q plus what fell in from above (tr_flux_in,
// flux_in), with a = g*dt/dp converting a flux to a mixing ratio change
// [kg/kg per kg/m2/s]. Where the level held no water the flux out is what
// fell in, at the incoming ratio. Exact for a uniform ratio; P3 substeps its
// sedimentation, so with vertically varying ratios this is an approximation.
//
// No phase change, so this cannot fractionate.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT sedimentation_flux (const ScalarT& tr, const ScalarT& q,
                            const ScalarT& tr_flux_in, const ScalarT& flux_in,
                            const ScalarT& flux_out, const ScalarT& a)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const ScalarT q_avail = q + flux_in*a;
  ScalarT ratio = tracer_ratio(tr + tr_flux_in*a, q_avail);
  const auto empty = ekat::abs(q_avail) < RealT(water_tracers_qmin);
  if (empty.any()) {
    ratio.set(empty, tracer_ratio(tr_flux_in, flux_in));
  }
  return flux_out*ratio;
}

// -----------------------------------------------------------------------
// A change of bulk vapor that no export accounts for, `amount` [kg/kg]: a
// property-check repair (SHOC's or P3's qv bounds), or roundoff of P3's
// wet/dry conversions. The tracer follows at its ratio to q_v, the bulk vapor
// before the change. Not a physical process, so it does not fractionate.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT vapor_repair (const ScalarT& tr_v, const ScalarT& q_v, const ScalarT& amount)
{
  return amount*tracer_ratio(tr_v, q_v);
}

} // namespace water_tracers
} // namespace scream

#endif // EAMXX_WATER_TRACER_TRANSFERS_FUNCTIONS_HPP
