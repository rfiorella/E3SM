#ifndef EAMXX_WATER_TRACER_TRANSFERS_FUNCTIONS_HPP
#define EAMXX_WATER_TRACER_TRANSFERS_FUNCTIONS_HPP

#include "eamxx_water_tracers_functions.hpp"  // tracer_ratio

#include "physics/p3/p3_water_transfers.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_pack.hpp>
#include <ekat_pack_math.hpp>  // ekat::expm1, ekat::abs; Kokkos::log1p

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
 * Device-callable transfers of the water_tracer_transfers process. Same
 * ScalarT convention as eamxx_water_tracers_functions.hpp (an ekat::Pack;
 * every division, log and comparison that could see a zero is masked, so
 * these are FPE-safe).
 *
 * Every replayed phase change moves tracer out of a donor phase through
 * donor_transfers, with one fractionation factor per sink,
 *   beta = R_sink/R_donor,
 * the ratio of what the sink removes to the donor's own tracer ratio. Each
 * class of bulk water transfer gets its beta from one function below
 * (shoc_condensation_beta, part2_beta). Only isotope tracers fractionate;
 * every other tracer has beta = 1 for every class, which donor_transfers
 * turns into the linear transfer it has always had, bit for bit.
 *
 * The betas come from
 *   alpha_l = R_liquid/R_vapor and alpha_i = R_ice/R_vapor
 * (water_isotopes::WisoAlphaDir CondensedOverVapor) at the temperature at the
 * start of the stage (water_transfer_stage_T): the equilibrium factors, or,
 * with the process option kinetic_fractionation, the effective (kinetic)
 * factors of the transfers that have one (water_isotopes::WaterIsotopeKinetic
 * alpha_liquid_vapor_kinetic for net evaporation of liquid,
 * alpha_ice_vapor_kinetic for deposition on ice).
 */

// Bounds of beta in donor_transfers. The equilibrium alphas of the water
// isotopologues are within [0.5, 2] at any temperature of the atmosphere, but
// the fits are extrapolated outside their fitted range (to 50 K at worst,
// where water_isotopes::check_temperature stops the run) and can return values
// far from that there; the bounds only keep the Rayleigh step finite and
// positive for whatever they return.
static constexpr Real fractionation_beta_min = 1.e-6;
static constexpr Real fractionation_beta_max = 1.e6;

// -----------------------------------------------------------------------
// The tracer moved out of one donor phase by its sinks [kg/kg].
//
//   tr_d, q_d  the donor's tracer and bulk water available to the sinks
//              (its reservoir, plus its inflows in P3's part2 stage)
//   amount[k]  bulk mass moved by sink k [kg/kg], >= 0; 0 for an idle sink,
//              or for an entry that is not a sink of this donor
//   beta[k]    sink k's fractionation factor, R_sink/R_donor
//   moved[k]   on return, the tracer moved by sink k [kg/kg]
//
// Lanes where every beta is exactly 1 (every lane of a non-fractionating
// tracer) get the linear transfer,
//   moved[k] = amount[k]*tracer_ratio(tr_d, q_d),
// bit for bit as before fractionation.
//
// Elsewhere the sinks act together as a Rayleigh distillation of the donor.
// While the sinks, in fixed proportions, remove dq from the donor's bulk
// water q, sink k removes tracer at beta_k times the donor's ratio at that
// moment, so
//   d(tr) = beta_eff*(tr/q)*dq,   beta_eff = sum_k amount_k*beta_k / A,
//   A = sum_k amount_k,
// and integrating from q_d to q_d - A (the Rayleigh law R/R_0 =
// f^(beta-1); e.g. Gat 1996, Annu. Rev. Earth Planet. Sci. 24, 225-262;
// Mook 2000, Environmental Isotopes in the Hydrological Cycle, vol. I)
//   tr_left = tr_d*f^beta_eff,   f = (q_d - A)/q_d.
// The tracer moved, tr_d*(1 - f^beta_eff), is split among the sinks in
// proportion to amount_k*beta_k. This
//  - is the linear transfer at beta = 1, to roundoff;
//  - cannot overdraw the donor or strand tracer in it: f is clamped to [0, 1],
//    so the total moved is between 0 and tr_d, and all of tr_d leaves when
//    A >= q_d (and, if A > 0, when q_d < water_tracers_qmin: a donor without
//    bulk water keeps no tracer);
//  - does not depend on the order of the sinks, and moves the same tracer in
//    one step as in several (f multiplies).
// 1 - f^beta_eff is formed as -expm1(beta_eff*log1p(-A/q_d)), which keeps full
// precision when A << q_d. Each beta is clamped to [fractionation_beta_min,
// fractionation_beta_max] first.
//
// Returns the lanes that took the Rayleigh path.
// -----------------------------------------------------------------------
template <typename ScalarT, int N>
KOKKOS_INLINE_FUNCTION
ekat::Mask<ScalarT::n> donor_transfers (const ScalarT& tr_d, const ScalarT& q_d,
                                        const ScalarT (&amount)[N], const ScalarT (&beta)[N],
                                        ScalarT (&moved)[N])
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
  using Mask  = ekat::Mask<ScalarT::n>;

  // Local copies: on CUDA, device code may not bind a reference to a
  // namespace-scope variable (passing one to a Pack function would)
  constexpr RealT qmin  = water_tracers_qmin;
  constexpr RealT b_min = fractionation_beta_min;
  constexpr RealT b_max = fractionation_beta_max;
  constexpr RealT zero  = 0;
  constexpr RealT one   = 1;

  // The linear transfer, at the donor's ratio
  const ScalarT ratio = tracer_ratio(tr_d, q_d);
  Mask rayleigh(false);
  for (int k=0; k<N; ++k) {
    moved[k] = amount[k]*ratio;
    rayleigh = rayleigh || (beta[k] != one);
  }
  if (not rayleigh.any()) {
    return rayleigh;
  }

  // Total sink A and its beta-weighted sum S = beta_eff*A [kg/kg]
  ScalarT A(0), S(0), weight[N];
  for (int k=0; k<N; ++k) {
    weight[k] = amount[k]*ekat::min(ekat::max(beta[k], b_min), b_max);
    A += amount[k];
    S += weight[k];
  }
  const Mask sinks = A > zero;  // S > 0 there too
  const Mask wet   = q_d >= qmin;

  // Fraction of the donor's bulk water removed, x = 1 - f = A/q_d in [0, 1]
  // (1 for a dry donor)
  ScalarT q_safe = q_d;
  q_safe.set(!wet, one);
  ScalarT x = ekat::min(A/q_safe, one);
  x.set(!wet, one);

  // Fraction of the donor's tracer removed, F = 1 - f^beta_eff, for 0 <= x < 1
  const Mask partial = sinks && x < one;
  ScalarT A_safe = A, S_safe = S, x_safe = x;
  A_safe.set(!sinks, one);
  S_safe.set(!sinks, one);
  x_safe.set(!partial, zero);
  ScalarT log_f;
  for (int s=0; s<ScalarT::n; ++s) {
    log_f[s] = Kokkos::log1p(-x_safe[s]);
  }
  ScalarT F = -ekat::expm1((S_safe/A_safe)*log_f);
  F.set(sinks && !partial, one);  // the sinks take all of the donor
  F.set(!sinks, zero);

  const ScalarT total = tr_d*F;
  for (int k=0; k<N; ++k) {
    moved[k].set(rayleigh, total*(weight[k]/S_safe));
  }
  return rayleigh;
}

// -----------------------------------------------------------------------
// SHOC's net condensation `amount` [kg/kg]: > 0 moves vapor to cloud liquid
// (the donor is vapor), < 0 moves cloud liquid to vapor (the donor is cloud
// liquid). Returns the tracer moved, with the sign of amount.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT shoc_condensation (const ScalarT& tr_donor, const ScalarT& q_donor, const ScalarT& amount,
                           const ScalarT& beta)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  const ScalarT a[1] = { ekat::abs(amount) };
  const ScalarT b[1] = { beta };
  ScalarT moved[1];
  donor_transfers(tr_donor, q_donor, a, b, moved);
  moved[0].set(amount < RealT(0), -moved[0]);
  return moved[0];
}

// The fractionation factor of SHOC's net condensation `amount`, given the
// liquid/vapor factors at the temperature after SHOC's transport:
// condensation (amount > 0) leaves the new liquid in equilibrium with the
// vapor, beta = alpha_cond (the equilibrium alpha_l); evaporation (amount < 0)
// gives the vapor the liquid's ratio over alpha_evap, beta = 1/alpha_evap.
// alpha_evap is alpha_l too, or with kinetic fractionation the effective
// factor of evaporation into the level's subsaturated air (akel). The net
// transfer of a level is all that SHOC exports, so its sign stands in for
// "unpaired" evaporation (see WaterTracerTransfers::replay_shoc).
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT shoc_condensation_beta (const ScalarT& alpha_cond, const ScalarT& alpha_evap,
                                const ScalarT& amount)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

  ScalarT beta = alpha_cond;
  beta.set(amount < RealT(0), RealT(1)/alpha_evap);
  return beta;
}

// Equilibrium only: alpha_l for condensation, 1/alpha_l for evaporation
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT shoc_condensation_beta (const ScalarT& alpha_l, const ScalarT& amount)
{
  return shoc_condensation_beta(alpha_l, alpha_l, amount);
}

// -----------------------------------------------------------------------
// The fractionation factor of P3 part2 process rate `cmp` (a Part2*
// component of p3_water_transfers), given alpha_l and alpha_i at P3's part2
// temperature. alpha_l is only used for liquid that evaporates and alpha_i
// only for vapor that deposits, so with kinetic fractionation the caller
// passes the effective factors (akel, akci) and the kinetic Bergeron factor
// is akci/akel:
//  - vapor deposition and deposition/condensation-freezing nucleation
//    (v -> i): alpha_i;
//  - rain evaporation (r -> v): 1/alpha_l;
//  - Bergeron (c -> i): alpha_i/alpha_l. The mass passes through the vapor:
//    the liquid evaporates (1/alpha_l) and the vapor deposits on the ice
//    (alpha_i), so the ice gets the liquid's ratio times alpha_i/alpha_l.
//    A deliberate departure from the EAM Fortran (wtype_get_alpha in
//    share/util/water_types.F90), which has no Bergeron branch and so uses 1;
//  - sublimation (i -> v): 1, as in the Fortran (ice does not exchange
//    isotopes with the vapor; it sublimates layer by layer);
//  - accretion, autoconversion, shedding, riming, collection, freezing and
//    melting: 1 (no phase change through the vapor).
// For a non-fractionating tracer, alpha_l = alpha_i = 1 gives 1 everywhere.
// -----------------------------------------------------------------------
template <typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT part2_beta (const int cmp, const ScalarT& alpha_l, const ScalarT& alpha_i)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
  using WT = p3::WaterTransfer;

  switch (cmp) {
    case WT::Part2VapdepV2I:
    case WT::Part2NucleatV2I: return alpha_i;
    case WT::Part2EvapR2V:    return RealT(1)/alpha_l;
    case WT::Part2BergC2I:    return alpha_i/alpha_l;
    default:                  return ScalarT(1);
  }
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
