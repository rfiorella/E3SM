#ifndef EAMXX_WATER_TRACER_TRANSFERS_HPP
#define EAMXX_WATER_TRACER_TRANSFERS_HPP

#include "eamxx_water_tracers_registry.hpp"
#include "../water_isotopes/eamxx_water_isotopes.hpp"

#include "share/atm_process/atmosphere_process.hpp"
#include <ekat_parameter_list.hpp>
#include <ekat_pack.hpp>

#include <string>
#include <vector>

namespace scream
{
/*
 * Water tracer transfers: moves the water tracers between phases the way
 * SHOC and P3 moved bulk water, by replaying the transfers they export (they
 * always do in builds with EAMXX_ENABLE_WATER_TRACERS). Bulk water is only read.
 *
 * Each transfer moves the tracer at the ratio R = tracer/bulk of its donor
 * phase, fractionated for isotope tracers (see Fractionation below). The
 * tracers are those of the water_tracers process, read from the same
 * water_tracers_file.
 *
 * SHOC replay (replay_shoc, default true). The process must run after SHOC,
 * in the same subcycle: shoc_water_stages holds only SHOC's last call.
 * The tracers are diffused in SHOC's own solve, so only SHOC's condensation is
 * replayed (see shoc_water_transfers.hpp). Afterwards <tracer>_residual_shoc
 * is the bulk water SHOC returned minus the tracer, per phase.
 *
 * P3 replay (replay_p3, default true). The process must run right after P3,
 * in the same subcycle and with the same dt: the P3 exports hold only P3's
 * last call, as rates over it. Starting from p3_water_entry, the replay
 * applies P3's stages in P3's order (see p3_water_transfers.hpp):
 *   1. checkpoint 1 (below);
 *   2. the part1 clips;
 *   3. the part2 process rates, applied together from the post-part1 state,
 *      as P3 applies them. Each donor's ratio is that of its reservoir plus
 *      its part2 inflows: P3 lets a sink exceed the reservoir by the
 *      inflows, and the reservoir's own ratio would amplify its roundoff by
 *      that factor at every call (see replay_p3);
 *   4. the part2 clips;
 *   5. sedimentation of cloud, rain and ice: one top-down sweep per column
 *      (see water_tracers::sedimentation_flux);
 *   6. homogeneous freezing and the part3 clips;
 *   7. checkpoint 2 (below).
 * Clips and homogeneous freezing move the tracer's whole donor reservoir
 * wherever P3's amount is nonzero (see water_tracers::whole_reservoir).
 * Alongside each tracer the replay rebuilds the bulk water from the same
 * exports, which gives the donor ratios at every stage.
 *
 * Bulk vapor also changes where no export accounts for it: property-check
 * repairs of qv, and roundoff of P3's wet/dry conversions. The checkpoints
 * apply those changes to the tracer vapor, at its ratio to the bulk vapor
 * before the change:
 *  - checkpoint 1: p3_water_entry minus SHOC's vapor (from shoc_water_stages):
 *    SHOC's qv postcondition repair and P3's qv precondition repair. It needs
 *    replay_shoc; without it there is no record of the bulk vapor before P3,
 *    and nothing is applied (a P3 precondition repair goes unfollowed).
 *  - checkpoint 2: the bulk qv P3 returned minus the rebuilt bulk vapor: P3's
 *    qv postcondition repair and wet/dry roundoff.
 * Condensate mismatches are never corrected, so a missing export still shows
 * in qc, qr or qi. <tracer>_residual_p3 is, per phase, the bulk water P3
 * returned minus the tracer after the replay, before checkpoint 2: for a
 * tracer with R = 1 it is the closure error (in vapor, plus any repair of
 * P3's qv postcondition); for other tracers only the bulk-minus-tracer
 * difference.
 *
 * <tracer>_precip_{liq,ice}_surf_mass accumulate the tracer's sedimentation
 * out of the lowest level (cloud + rain, and ice) [kg/m2], like P3's
 * precip_{liq,ice}_surf_mass; the driver zeroes them at each atm step.
 *
 * water_transfer_stage_T is the temperature at the start of each replayed
 * stage [K], in water_tracers::TransferStage order, rebuilt once per call
 * (not per tracer) from T_mid and the exports (see
 * water_tracers::stage_temperatures). Nothing between SHOC and this process
 * other than P3 may change T_mid; the driver checks it. Stages that are not
 * replayed hold water_tracers::stage_T_fill.
 *
 * Fractionation. Only isotope tracers fractionate; every other tracer moves
 * at its donor's ratio. For an isotope, each phase change through
 * the vapor moves the tracer with an equilibrium fractionation factor beta =
 * R_sink/R_donor (water_tracers::shoc_condensation_beta, part2_beta), from
 * alpha_l and alpha_i evaluated once per level and tracer at the stage
 * temperature (ShocAfterTransport for SHOC, P3Part2 for P3's part2 rates).
 * Each donor's sinks then act as a Rayleigh distillation of the donor
 * (water_tracers::donor_transfers), which cannot overdraw it. Clips,
 * homogeneous freezing, sedimentation, the checkpoints and the surface
 * fluxes never fractionate. The isotope options (formulations of alpha) are
 * those of water_tracers; the driver copies them. The test-only parameter
 * alpha_override > 0 replaces alpha_l and alpha_i by that value.
 * water_transfer_alpha_extrapolations counts, per column and atm step, the
 * alpha evaluations at a temperature outside the fit's accepted range (the
 * fit is extrapolated); the run's total is logged once, at finalize.
 */

class WaterTracerTransfers : public AtmosphereProcess
{
public:
  using KT   = ekat::KokkosTypes<DefaultDevice>;
  using Pack = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  WaterTracerTransfers (const ekat::Comm& comm, const ekat::ParameterList& params);

  AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }

  std::string name () const override { return "water_tracer_transfers"; }

  void create_requests () override;

protected:

#ifdef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
public:
#endif

  void initialize_impl (const RunType /* run_type */) override;
  void run_impl        (const double dt) override;
  void finalize_impl   () override;

  // Steps of run_impl
  void rebuild_stage_temperatures (const Real dt);
  void replay_shoc ();
  void replay_p3 (const Real dt);

protected:

  int m_num_cols;
  int m_num_levs;

  std::shared_ptr<const AbstractGrid> m_grid;

  // Which exported transfers to replay
  bool m_replay_shoc;
  bool m_replay_p3;

  // One entry per tracer in water_tracers_file, sorted by name. Parsed in
  // create_requests, so the parameter can be set after construction.
  std::vector<water_tracers::TracerSpec> m_specs;

  // The isotope options (configured in create_requests, from the parameters
  // the driver copies from water_tracers): the formulations of alpha_l and
  // alpha_i, and the isotope specs' resolution
  water_isotopes::WaterIsotopes m_isotopes;

  // Test only: if > 0, alpha_l = alpha_i = m_alpha_override for every
  // isotope tracer, instead of the equilibrium alphas (0: off)
  Real m_alpha_override;

  // P3 replay: the bulk water rebuilt from P3's exports, per column, phase
  // (WaterPhase order) and level pack [kg/kg, wet]. Scratch, set by each
  // tracer's replay.
  KT::view_3d<Pack> m_bulk;

  // Per column: the alpha evaluations that extrapolated a fit, summed over the
  // run (see water_transfer_alpha_extrapolations). Read once, at finalize.
  KT::view_1d<Real> m_alpha_extrapolations_run;

}; // class WaterTracerTransfers

} // namespace scream

#endif // EAMXX_WATER_TRACER_TRANSFERS_HPP
