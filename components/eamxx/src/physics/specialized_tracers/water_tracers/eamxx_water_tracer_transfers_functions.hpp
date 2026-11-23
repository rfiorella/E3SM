#ifndef EAMXX_WATER_TRACER_TRANSFERS_FUNCTIONS_HPP
#define EAMXX_WATER_TRACER_TRANSFERS_FUNCTIONS_HPP

#include "eamxx_water_tracers_functions.hpp"  // tracer_ratio

#include "physics/p3/p3_water_transfers.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/util/eamxx_universal_constants.hpp"

namespace scream {
namespace water_tracers {

/*
 * Components of the water_tracer_transfers field water_transfer_stage_T: the
 * temperature at the start of each replayed stage [K] (see
 * stage_temperatures). A stage that is not replayed holds stage_T_fill.
 */
struct TransferStage {
enum : int {
  ShocAfterTransport,  // SHOC after its turbulent transport, before its condensation
  P3Entry,             // P3's entry state
  P3Part2,             // P3 after its part1 clips: the state its part2 rates are computed from
  Count
};
};

// water_transfer_stage_T of a stage that is not replayed [K]: EAMxx's fill
// value, which output files mark as missing (_FillValue)
static constexpr Real stage_T_fill = constants::fill_value<Real>;

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

// -----------------------------------------------------------------------
// The temperature at the start of each replayed stage [K] (TransferStage
// order), rebuilt from the temperature now and the exports. SHOC and P3
// export no temperature, but both heat linearly in the water they move, so
// their heating can be undone. Stages that are not replayed get stage_T_fill.
//
//   T_now        T_mid now [K]: after P3 when P3 is replayed, otherwise after
//                SHOC (nothing between them and this process may change T_mid;
//                the driver checks it)
//   shoc_cond    SHOC's net condensation, LiquidAfterCondensation minus
//                LiquidAfterTransport of shoc_water_stages [kg/kg, wet]
//   phases       P3's heating table (p3::water_transfer_phases)
//   p3_amount    p3_amount(cmp): mass moved by p3_water_transfers component
//                cmp over the P3 call, rate*dt [kg/kg, wet]
//   wet_to_dry   pseudo_density/pseudo_density_dry [-]
//
// P3. Each transfer heats theta by inv_exner*(H(to) - H(from))*m_dry/cp,
// with m_dry the mass moved on P3's dry basis (see p3::WaterTransferPhases).
// The P3 interface then rescales the call's dT by dp_dry/dp
// (eamxx_p3_process_interface.hpp, p3_postamble), the factor that also turns
// m_dry into the exported wet amount m. Hence
//   T_P3Entry = T_now - sum_k (H(to_k) - H(from_k))*m_k/cp
// over all components. P3 computes its part2 rates from its temperature after
// the part1 clips, whose heating is not rescaled (it is on the dry basis):
//   T_P3Part2 = T_P3Entry + sum_{part1 clips} (H(to_k) - H(from_k))*m_k*wet_to_dry/cp
//
// SHOC. shoc_main sets dse from T = thl/inv_exner + (LatVap/cp)*ql
// (update_host_dse), and the SHOC interface sets T_mid from that dse; thl is
// not changed by SHOC's condensation (shoc_assumed_pdf only reads it). The
// temperature with the liquid after transport is therefore
//   T_ShocAfterTransport = T_after_shoc - (LatVap/cp)*shoc_cond
// with T_after_shoc = T_P3Entry when P3 is replayed, else T_now. It includes
// SHOC's energy fixer, which corrects dse after update_host_dse and is not
// exported separately.
//
// cp, LatVap and LatIce are the constants P3 and SHOC use. In the
// water_tracer_transfers tests (dt = 300 s), |T_now - T_P3Entry| reaches
// 0.8 K, |T_P3Part2 - T_P3Entry| 3e-5 K and SHOC's condensation heating 1.5 K.
// -----------------------------------------------------------------------
template <typename ScalarT, typename AmountFn>
KOKKOS_INLINE_FUNCTION
void stage_temperatures (const ScalarT& T_now,
                         const bool replay_shoc, const ScalarT& shoc_cond,
                         const bool replay_p3, const p3::WaterTransferPhases& phases,
                         const AmountFn& p3_amount, const ScalarT& wet_to_dry,
                         ScalarT (&T)[TransferStage::Count])
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
  using C  = physics::Constants<RealT>;
  using TS = TransferStage;
  using WT = p3::WaterTransfer;
  using WE = p3::WaterEntry;
  static_assert(WT::Part1ClipC2V == 0 && WT::Part1ClipI2R == 3,
                "The part1 clips must be the first four WaterTransfer components.");

  constexpr RealT latvap = C::LatVap.value;  // [J/kg]
  constexpr RealT latice = C::LatIce.value;  // [J/kg]
  constexpr RealT cp     = C::CP.value;      // [J/kg/K]
  constexpr RealT inv_cp = C::INV_CP.value;  // [kg K/J]
  // Local copy: on CUDA, device code may not bind a reference to a
  // namespace-scope variable (assigning it to a Pack would)
  constexpr RealT fill   = stage_T_fill;

  // Latent enthalpy of a phase relative to vapor [J/kg]
  const auto H = [&] (const int phase) -> RealT {
    return phase==WE::Vapor ? RealT(0) : (phase==WE::Ice ? latvap + latice : latvap);
  };

  T[TS::P3Entry] = fill;
  T[TS::P3Part2] = fill;
  ScalarT T_after_shoc = T_now;
  if (replay_p3) {
    // P3's heating over the call, and over its part1 clips [K]
    ScalarT heating(0), part1_heating(0);
    for (int cmp=0; cmp<WT::Count; ++cmp) {
      const ScalarT dT = (H(phases.to[cmp]) - H(phases.from[cmp]))*p3_amount(cmp)*inv_cp;
      heating += dT;
      if (cmp <= WT::Part1ClipI2R) {
        part1_heating += dT;
      }
    }
    T[TS::P3Entry] = T_now - heating;
    T[TS::P3Part2] = T[TS::P3Entry] + part1_heating*wet_to_dry;
    T_after_shoc   = T[TS::P3Entry];
  }

  T[TS::ShocAfterTransport] = fill;
  if (replay_shoc) {
    T[TS::ShocAfterTransport] = T_after_shoc - (latvap/cp)*shoc_cond;
  }
}

} // namespace water_tracers
} // namespace scream

#endif // EAMXX_WATER_TRACER_TRANSFERS_FUNCTIONS_HPP
