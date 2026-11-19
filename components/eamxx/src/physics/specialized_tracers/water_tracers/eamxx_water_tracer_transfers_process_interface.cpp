#include "eamxx_water_tracer_transfers_process_interface.hpp"
#include "eamxx_water_tracer_transfers_functions.hpp"
#include "../water_isotopes/eamxx_water_isotopes_fractionation.hpp"
#include "../water_isotopes/eamxx_water_isotopes_kinetic.hpp"

#include "physics/p3/p3_water_transfers.hpp"
#include "physics/shoc/shoc_water_transfers.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/physics/physics_functions.hpp"  // qv_sat_wet (ETI on CPU builds)

#include <ekat_assert.hpp>
#include <ekat_pack_kokkos.hpp>
#include <ekat_subview_utils.hpp>
#include <ekat_team_policy_utils.hpp>
#include <ekat_units.hpp>

namespace scream
{

namespace {

namespace wt = water_tracers;

using Pack    = WaterTracerTransfers::Pack;
using IntPack = ekat::Pack<int, Pack::n>;

constexpr int iv = etoi(wt::WaterPhase::Vapor);
constexpr int ic = etoi(wt::WaterPhase::CloudLiquid);
constexpr int ir = etoi(wt::WaterPhase::Rain);
constexpr int ii = etoi(wt::WaterPhase::Ice);
constexpr int nphases = etoi(wt::WaterPhase::Count);

// A clip-type P3 stage (a clip or homogeneous freezing): the bulk moved its
// whole `from` reservoir, `amount` [kg/kg], to `to`. tr and q are the
// tracer's and the rebuilt bulk water's phases. Like P3, the bulk donor is
// set to 0 where the stage moved something.
KOKKOS_INLINE_FUNCTION
void clip_stage (Pack (&tr)[nphases], Pack (&q)[nphases], const int from, const int to,
                 const Pack& amount)
{
  const Pack moved = wt::whole_reservoir(tr[from], amount);
  tr[from] -= moved;
  tr[to]   += moved;
  q[to]    += amount;
  q[from].set(amount != Real(0), Real(0));
}

// The saturation ratio qv/qv_sat of vapor qv [kg/kg, wet] over liquid (ice =
// false) or ice, at temperature T [K], dry pressure p_dry [Pa] and wet and dry
// pseudo densities dp, dp_dry [Pa]: the relative humidity of EAMxx's RH
// diagnostic (share/diagnostics/relative_humidity.cpp), Murphy & Koop (2005)
// saturation vapor pressure, qv_sat on the wet basis of qv. It is also P3's
// own ratio, which P3 forms on its dry basis (qv*dp/dp_dry over qv_sat_dry at
// p_dry_mid). Lanes outside mask return 1, and get in-range inputs before the
// saturation function (which evaluates every lane), so padding cannot trap.
KOKKOS_INLINE_FUNCTION
Pack saturation_ratio (const Pack& qv, const Pack& T, const Pack& p_dry, const Pack& dp,
                       const Pack& dp_dry, const bool ice, const ekat::Mask<Pack::n>& mask)
{
  using PF = physics::Functions<Real, DefaultDevice>;

  Pack qv_live(0), T_live(273.15), p_live(1.e5), dp_live(1), dp_dry_live(1);
  qv_live.set(mask, qv);
  T_live.set(mask, T);
  p_live.set(mask, p_dry);
  dp_live.set(mask, dp);
  dp_dry_live.set(mask, dp_dry);
  const Pack qv_sat = PF::qv_sat_wet(T_live, p_live, ice, mask, dp_live, dp_dry_live, PF::MurphyKoop,
                                     "water_tracer_transfers (saturation ratio)");
  Pack s(1);
  s.set(mask, qv_live/qv_sat);
  return s;
}

// P3's part2 process rates, in WaterTransfer order from Part2VapdepV2I
constexpr int first_part2 = p3::WaterTransfer::Part2VapdepV2I;
constexpr int num_part2   = p3::WaterTransfer::Part2MeltI2R - first_part2 + 1;
static_assert(num_part2 == 13, "Update the part2 tables below to match p3_water_transfers.hpp.");

// Fixed-point passes that refine the part2 donor ratios (see replay_p3). The
// longest chain of donors without a cycle is c -> r -> i -> v, so three
// passes carry the ratio of a well-filled reservoir to the end of any chain.
// (In the isotope tests, 2, 3 and 6 passes give the same result bit for bit;
// 1 pass differs by up to 1e-4 of a tracer phase.)
constexpr int num_ratio_passes = 3;

// One top-down sedimentation sweep of one species in one column (level 0 is
// the top). flux is P3's flux out of the bottom of each level [kg/m2/s]; tr
// and q are the tracer and the rebuilt bulk water [kg/kg], updated in place;
// dp is the wet pseudo density [Pa]. Returns the tracer's flux out of the
// lowest level [kg/m2/s]. Serial in k: each level's outgoing flux depends on
// what falls in from above.
template <typename TrView, typename QView, typename FluxView, typename DpView>
KOKKOS_INLINE_FUNCTION
Real sedimentation_sweep (const TrView& tr, const QView& q, const FluxView& flux,
                          const DpView& dp, const int nlev, const Real dt)
{
  using SPack = ekat::Pack<Real,1>;
  constexpr Real g = physics::Constants<Real>::gravit.value;  // [m/s2]

  SPack flux_in(0), tr_flux_in(0);
  for (int k=0; k<nlev; ++k) {
    // Mixing ratio change per unit flux over dt: g*dt/dp [kg/kg per kg/m2/s]
    const SPack a = g*dt/dp(k);
    const SPack flux_out = flux(k);
    const SPack tr_flux_out = wt::sedimentation_flux(SPack(tr(k)), SPack(q(k)),
                                                     tr_flux_in, flux_in, flux_out, a);
    tr(k) += ((tr_flux_in - tr_flux_out)*a)[0];
    q(k)  += ((flux_in - flux_out)*a)[0];
    flux_in    = flux_out;
    tr_flux_in = tr_flux_out;
  }
  return tr_flux_in[0];
}

} // anonymous namespace

// =========================================================================================
WaterTracerTransfers::WaterTracerTransfers (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  m_replay_shoc = m_params.get<bool>("replay_shoc", true);
  m_replay_p3   = m_params.get<bool>("replay_p3",   true);

  // Test only (not in the namelist defaults): see m_alpha_override
  m_alpha_override = m_params.get<double>("alpha_override", 0.0);
  EKAT_REQUIRE_MSG(m_alpha_override >= 0,
      "Error! 'water_tracer_transfers::alpha_override' must be > 0 (or 0, off), got " +
      std::to_string(m_alpha_override) + ".\n");

  // Kinetic fractionation of isotopes in the replay (see the class comment)
  m_kinetic = m_params.get<bool>("kinetic_fractionation", false);
  const auto ice_S = m_params.get<std::string>("kinetic_ice_supersaturation", "empirical");
  EKAT_REQUIRE_MSG(ice_S=="empirical" or ice_S=="model",
      "Error! Invalid 'water_tracer_transfers::kinetic_ice_supersaturation': '" + ice_S +
      "'. Valid options: empirical, model.\n");
  m_kinetic_ice_model_S = ice_S=="model";
}

// =========================================================================================
void WaterTracerTransfers::create_requests ()
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // Same tracers file, and the same isotope support, as the water_tracers process
  wt::RegistryOptions opts;
#ifdef EAMXX_HAS_WATER_ISOTOPES
  water_isotopes::WaterIsotopes::enable(opts);
#endif

  const std::string key = "water_tracers_file";
  EKAT_REQUIRE_MSG(m_params.isParameter(key) and m_params.isType<std::string>(key) and
                   m_params.get<std::string>(key)!="",
      "Error! Process '" + m_params.name() + "' is in atm_procs_list, but its parameter\n"
      "  'water_tracers_file' is missing or empty. Set it to the water_tracers process' file.\n");

  m_specs = wt::parse_tracers_file(m_params.get<std::string>(key), opts);

  // The isotope formulation options (the driver copies water_tracers' ones,
  // so both processes use the same alphas). Read even with no isotope
  // tracers, like water_tracers does; only isotopes evaluate alpha.
  m_isotopes.configure(m_params);
  m_isotopes.resolve(m_specs);

  m_grid = m_grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs();
  m_num_levs = m_grid->get_num_vertical_levels();

  constexpr int ps = SCREAM_PACK_SIZE;
  constexpr auto m2 = pow(m, 2);
  const auto scalar2d     = m_grid->get_2d_scalar_layout();
  const auto scalar3d_mid = m_grid->get_3d_scalar_layout(LEV);
  const auto phases3d_mid = m_grid->get_3d_vector_layout(LEV, etoi(wt::WaterPhase::Count), "water_phase");

  // Bulk water is only read: the tracers follow it, never the other way around
  for (const auto phase : wt::water_phases) {
    add_tracer<Required>(wt::bulk_field_name(phase), m_grid, kg/kg, ps);
  }
  add_field<Required>("pseudo_density", scalar3d_mid, Pa, grid_name, ps);

  // The temperature at the start of each replayed stage, rebuilt from T_mid
  // (see water_tracers::stage_temperatures). P3's part2 stage needs P3's dry basis.
  add_field<Required>("T_mid", scalar3d_mid, K, grid_name, ps);
  if (m_replay_p3 or m_kinetic) {
    add_field<Required>("pseudo_density_dry", scalar3d_mid, Pa, grid_name, ps);
  }
  // The saturation ratios of the kinetic factors (see saturation_ratio)
  if (m_kinetic) {
    add_field<Required>("p_dry_mid", scalar3d_mid, Pa, grid_name, ps);
  }
  const auto stage_T_layout = m_grid->get_3d_vector_layout(LEV, wt::TransferStage::Count, "water_transfer_stage");
  add_field<Computed>("water_transfer_stage_T", stage_T_layout, K, grid_name, ps);

  // Per column, the isotope alpha evaluations (levels x tracers x phases)
  // whose temperature was outside the fit's accepted range, so that the fit
  // was extrapolated, summed over the atm step like the precipitation [count]
  add_field<Computed>("water_transfer_alpha_extrapolations", scalar2d, none, grid_name, "ACCUMULATED");
  m_alpha_extrapolations_run = KT::view_1d<Real>("water_tracer_transfers_alpha_extrapolations_run", m_num_cols);

  // Bulk water after each SHOC stage (exported by SHOC in builds with water tracers)
  if (m_replay_shoc) {
    const auto stages_layout = m_grid->get_3d_vector_layout(LEV, shoc::WaterStage::Count, "shoc_water_stage");
    add_field<Required>("shoc_water_stages", stages_layout, kg/kg, grid_name, ps);
  }

  // P3's entry state, transfers and sedimentation fluxes (exported by P3 in
  // builds with water tracers; same layouts as P3 registers them)
  if (m_replay_p3) {
    const auto entry_layout     = m_grid->get_3d_vector_layout(LEV, p3::WaterEntry::Count, "p3_water_species");
    const auto transfers_layout = m_grid->get_3d_vector_layout(LEV, p3::WaterTransfer::Count, "p3_water_transfer");
    const auto sed_layout       = m_grid->get_3d_vector_layout(LEV, p3::SedFlux::Count, "p3_sed_species");
    add_field<Required>("p3_water_entry",     entry_layout,     kg/kg,   grid_name, ps);
    add_field<Required>("p3_water_transfers", transfers_layout, kg/kg/s, grid_name, ps);
    add_field<Required>("p3_sed_flux",        sed_layout,       kg/m2/s, grid_name, ps);

    m_bulk = KT::view_3d<Pack>("water_tracer_transfers_bulk", m_num_cols,
                               etoi(wt::WaterPhase::Count), ekat::npack<Pack>(m_num_levs));
  }

  // The tracers, one field per phase. Each is requested by name, never through
  // the tracers or turbulence_advected_tracers group: a group is allocated as
  // one monolithic field that also holds bulk qv, qc, qr and qi, which would
  // make bulk water writable here.
  for (const auto& spec : m_specs) {
    for (const auto phase : wt::water_phases) {
      add_tracer<Updated>(wt::phase_field_name(spec.name, phase), m_grid, kg/kg, ps,
                          TracerAdvection::DynamicsAndTurbulence);
    }
    // Per phase (WaterPhase order): the bulk water SHOC returned minus the
    // tracer after the SHOC replay [kg/kg]. This is the closure residual of
    // a tracer that has R = 1 in every phase; for other tracers it is only
    // the bulk-minus-tracer difference.
    if (m_replay_shoc) {
      add_field<Computed>(spec.name + "_residual_shoc", phases3d_mid, kg/kg, grid_name, ps);
    }
    if (m_replay_p3) {
      // Per phase: the bulk water P3 returned minus the tracer after the P3
      // replay, before checkpoint 2 [kg/kg] (see the class comment)
      add_field<Computed>(spec.name + "_residual_p3", phases3d_mid, kg/kg, grid_name, ps);
      // The tracer's surface precipitation, accumulated over the atm step like
      // P3's precip_{liq,ice}_surf_mass [kg/m2]. Computed, not Updated: the
      // driver zeroes ACCUMULATED fields at the start of each atm step, so
      // there is no previous value to require.
      add_field<Computed>(spec.name + "_precip_liq_surf_mass", scalar2d, kg/m2, grid_name, "ACCUMULATED");
      add_field<Computed>(spec.name + "_precip_ice_surf_mass", scalar2d, kg/m2, grid_name, "ACCUMULATED");
    }
  }
}

// =========================================================================================
void WaterTracerTransfers::initialize_impl (const RunType /* run_type */)
{
  // Stages that are not replayed hold the fill value, so output marks them as missing
  if (not (m_replay_shoc and m_replay_p3)) {
    get_field_out("water_transfer_stage_T", m_grid->name()).get_header().set_may_be_filled(true);
  }
}

// =========================================================================================
void WaterTracerTransfers::run_impl (const double dt)
{
  rebuild_stage_temperatures(dt);
  if (m_replay_shoc) {
    replay_shoc();
  }
  if (m_replay_p3) {
    replay_p3(dt);
  }
}

// =========================================================================================
void WaterTracerTransfers::finalize_impl ()
{
  // The alpha evaluations that extrapolated a fit, over the run and all
  // ranks: one log message per run, rather than one per pack, level and
  // step (the per-column, per-step count is water_transfer_alpha_extrapolations).
  // Reading it here keeps the run itself free of device-to-host copies.
  // (The view is unallocated if the run stopped before create_requests.)
  const auto run_count = m_alpha_extrapolations_run;
  if (not run_count.is_allocated()) {
    return;
  }
  Real num = 0;
  Kokkos::parallel_reduce("water_tracer_transfers_alpha_extrapolations_sum",
                          KT::RangePolicy(0, run_count.extent(0)),
                          KOKKOS_LAMBDA(const int icol, Real& sum) {
    sum += run_count(icol);
  }, num);
  m_comm.all_reduce(&num, 1, MPI_SUM);
  if (num > 0 and m_comm.am_i_root()) {
    m_atm_logger->warn("[water_tracer_transfers] " + std::to_string(static_cast<long long>(num)) +
                       " isotope fractionation factor evaluations were at a temperature outside the "
                       "fit's accepted range (see TemperatureBounds in eamxx_water_isotopes_parameters.hpp) "
                       "and extrapolated the fit. Per column and atm step: water_transfer_alpha_extrapolations.");
  }
}

// =========================================================================================
void WaterTracerTransfers::rebuild_stage_temperatures (const Real dt)
{
  using MDPolicy = Kokkos::MDRangePolicy<KT::ExeSpace, Kokkos::Rank<2>>;
  using WS = shoc::WaterStage;
  using TS = wt::TransferStage;

  const auto& gn = m_grid->name();
  const int nlev = m_num_levs;
  const int nlev_packs = ekat::npack<Pack>(nlev);
  const MDPolicy policy({0, 0}, {m_num_cols, nlev_packs});

  const auto T_mid   = get_field_in("T_mid", gn).get_view<const Pack**>();
  const auto stage_T = get_field_out("water_transfer_stage_T", gn).get_view<Pack***>();

  const bool replay_shoc = m_replay_shoc;
  const bool replay_p3   = m_replay_p3;
  Field::view_dev_t<const Pack***> stages, xfer;
  Field::view_dev_t<const Pack**>  dp, dp_dry;
  if (replay_shoc) {
    stages = get_field_in("shoc_water_stages", gn).get_view<const Pack***>();
  }
  if (replay_p3) {
    xfer   = get_field_in("p3_water_transfers", gn).get_view<const Pack***>();
    dp     = get_field_in("pseudo_density", gn).get_view<const Pack**>();
    dp_dry = get_field_in("pseudo_density_dry", gn).get_view<const Pack**>();
  }

  // P3's heating table, captured by value
  const auto phases = p3::water_transfer_phases();

  Kokkos::parallel_for("water_tracer_transfers_stage_T", policy,
                       KOKKOS_LAMBDA(const int icol, const int k) {
    // Padded lanes past the last level are read as `pad` (never what the
    // fields hold there), so the math cannot trap
    const auto in_range = ekat::range<IntPack>(k*Pack::n) < nlev;
    const auto get = [&] (const auto& x, const Real pad) {
      Pack y(pad);
      y.set(in_range, x);
      return y;
    };

    Pack shoc_cond(0);
    if (replay_shoc) {
      shoc_cond = get(stages(icol, WS::LiquidAfterCondensation, k), 0)
                - get(stages(icol, WS::LiquidAfterTransport,    k), 0);
    }
    Pack wet_to_dry(1);
    if (replay_p3) {
      wet_to_dry = get(dp(icol, k), 1)/get(dp_dry(icol, k), 1);
    }
    // Mass moved by a P3 component over the call [kg/kg] (the exports are rates)
    const auto p3_amount = [&] (const int cmp) {
      return get(xfer(icol, cmp, k), 0)*dt;
    };

    Pack T[TS::Count];
    wt::stage_temperatures(get(T_mid(icol, k), 0), replay_shoc, shoc_cond, replay_p3, phases,
                           p3_amount, wet_to_dry, T);
    for (int s=0; s<TS::Count; ++s) {
      stage_T(icol, s, k) = T[s];
    }
  });
}

// =========================================================================================
void WaterTracerTransfers::replay_shoc ()
{
  using MDPolicy = Kokkos::MDRangePolicy<KT::ExeSpace, Kokkos::Rank<2>>;
  using WS  = shoc::WaterStage;
  using TS  = wt::TransferStage;
  using WIF = water_isotopes::WaterIsotopeFractionation;
  using WIK = water_isotopes::WaterIsotopeKinetic;

  const auto& gn = m_grid->name();
  const int nlev = m_num_levs;
  const int nlev_packs = ekat::npack<Pack>(nlev);
  const MDPolicy policy({0, 0}, {m_num_cols, nlev_packs});

  const auto stages = get_field_in("shoc_water_stages", gn).get_view<const Pack***>();
  const auto qr     = get_field_in("qr", gn).get_view<const Pack**>();
  const auto qi     = get_field_in("qi", gn).get_view<const Pack**>();

  // Fractionation of isotope tracers (see the class comment)
  const auto stage_T     = get_field_out("water_transfer_stage_T", gn).get_view<const Pack***>();
  const auto num_extrap  = get_field_out("water_transfer_alpha_extrapolations", gn).get_view<Real*>();
  const auto extrap_run  = m_alpha_extrapolations_run;
  const auto iso_params  = m_isotopes.parameters();
  const Real alpha_fixed = m_alpha_override;

  // Kinetic fractionation: the state for the relative humidity
  const bool kinetic = m_kinetic and m_alpha_override <= 0;
  Field::view_dev_t<const Pack**> p_dry, dp, dp_dry;
  if (kinetic) {
    p_dry  = get_field_in("p_dry_mid", gn).get_view<const Pack**>();
    dp     = get_field_in("pseudo_density", gn).get_view<const Pack**>();
    dp_dry = get_field_in("pseudo_density_dry", gn).get_view<const Pack**>();
  }

  // With the P3 replay, P3 ran between SHOC and this process, so the bulk qr
  // and qi SHOC returned are P3's entry values, not the bulk fields
  const bool after_p3 = m_replay_p3;
  using WE = p3::WaterEntry;
  Field::view_dev_t<const Pack***> entry;
  if (after_p3) {
    entry = get_field_in("p3_water_entry", gn).get_view<const Pack***>();
  }

  // One kernel per tracer: tracers are few, their number is only known at
  // runtime, and each reads and writes its own fields.
  for (const auto& spec : m_specs) {
    const auto tr_qv = get_field_out(wt::phase_field_name(spec.name, wt::WaterPhase::Vapor), gn)
                         .get_view<Pack**>();
    const auto tr_qc = get_field_out(wt::phase_field_name(spec.name, wt::WaterPhase::CloudLiquid), gn)
                         .get_view<Pack**>();
    const auto tr_qr = get_field_in(wt::phase_field_name(spec.name, wt::WaterPhase::Rain), gn)
                         .get_view<const Pack**>();
    const auto tr_qi = get_field_in(wt::phase_field_name(spec.name, wt::WaterPhase::Ice), gn)
                         .get_view<const Pack**>();
    const auto residual = get_field_out(spec.name + "_residual_shoc", gn).get_view<Pack***>();

    // Only isotopes fractionate. The species is unused unless fractionates is true.
    const bool fractionates = spec.type==wt::TracerType::Isotope;
    const auto species = fractionates ? water_isotopes::WaterIsotopes::species(spec)
                                      : water_isotopes::WaterIsotopologues::H216O;

    Kokkos::parallel_for("water_tracer_transfers_shoc_" + spec.name, policy,
                         KOKKOS_LAMBDA(const int icol, const int k) {
      // Bulk water after SHOC's turbulent transport, and the cloud liquid
      // after its condensation [kg/kg]. The tracers come in already
      // transported (SHOC diffuses them in the same solve).
      const Pack qv_trans = stages(icol, WS::VaporAfterTransport,     k);
      const Pack qc_trans = stages(icol, WS::LiquidAfterTransport,    k);
      const Pack qc_cond  = stages(icol, WS::LiquidAfterCondensation, k);

      // Net condensation [kg/kg]: > 0 moves vapor to liquid (the donor is
      // the vapor), < 0 moves liquid to vapor (the donor is the liquid)
      const Pack cond = qc_cond - qc_trans;
      const auto evap = cond < Real(0);
      Pack tr_donor = tr_qv(icol, k);
      Pack q_donor  = qv_trans;
      tr_donor.set(evap, tr_qc(icol, k));
      q_donor.set(evap, qc_trans);

      // Its fractionation factor: 1, or for an isotope from alpha_l at the
      // temperature after SHOC's transport, evaluated where SHOC moved water
      // (elsewhere alpha is not used); padded lanes are masked out
      Pack beta(1);
      if (fractionates) {
        Pack alpha_l(alpha_fixed);
        Pack alpha_evap(alpha_fixed);  // the liquid/vapor factor of net evaporation
        if (alpha_fixed <= 0) {
          const auto in_range = ekat::range<IntPack>(k*Pack::n) < nlev;
          Pack T(0);
          T.set(in_range, stage_T(icol, TS::ShocAfterTransport, k));
          int nx = 0;
          alpha_l = WIF::alpha_liquid_vapor(T, species, water_isotopes::CondensedOverVapor, iso_params,
                                            in_range && cond != Real(0), &nx);
          if (nx > 0) {
            Kokkos::atomic_add(&num_extrap(icol), Real(nx));
            Kokkos::atomic_add(&extrap_run(icol), Real(nx));
          }
          alpha_evap = alpha_l;

          // Kinetic fractionation (D4 of the design): net evaporation is
          // taken as evaporation into the level's subsaturated air, with the
          // effective factor akel at the relative humidity over liquid of
          // SHOC's vapor after transport, at the same temperature. Net
          // condensation, and the cloud SHOC keeps, stay at equilibrium.
          // Two proxies: the sign of the net transfer stands in for
          // evaporation that no condensation in the same level pairs with
          // (SHOC exports only the net per level), and the grid-mean
          // relative humidity for that of the sub-grid part where the cloud
          // evaporates (SHOC's PDF gives none).
          const auto evap_lanes = in_range && evap;
          if (kinetic && evap_lanes.any()) {
            const Pack rh_l = saturation_ratio(qv_trans, T, p_dry(icol, k), dp(icol, k), dp_dry(icol, k),
                                               false, evap_lanes);
            alpha_evap = WIK::alpha_liquid_vapor_kinetic(alpha_l, rh_l, species, iso_params, evap_lanes);
          }
        }
        beta = wt::shoc_condensation_beta(alpha_l, alpha_evap, cond);
      }

      const Pack dtr = wt::shoc_condensation(tr_donor, q_donor, cond, beta);
      tr_qv(icol, k) -= dtr;
      tr_qc(icol, k) += dtr;

      // Checkpoint: bulk water as SHOC returned it, minus the tracer. SHOC
      // returns qv = qw - qc with qc = qc_cond, i.e. the stage sum below,
      // unless its qv postcondition repair (clip to [0, 0.2]) fired; the
      // stages do not see that repair. SHOC does not change qr and qi, so
      // P3's entry values, or else the bulk fields, stand in for them: valid
      // only while nothing else between SHOC and here changes them.
      residual(icol, iv, k) = (qv_trans + qc_trans - qc_cond) - tr_qv(icol, k);
      residual(icol, ic, k) = qc_cond - tr_qc(icol, k);
      residual(icol, ir, k) = (after_p3 ? entry(icol, WE::Rain, k) : qr(icol, k)) - tr_qr(icol, k);
      residual(icol, ii, k) = (after_p3 ? entry(icol, WE::Ice,  k) : qi(icol, k)) - tr_qi(icol, k);
    });
  }
}

// =========================================================================================
void WaterTracerTransfers::replay_p3 (const Real dt)
{
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  using WS  = shoc::WaterStage;
  using WE  = p3::WaterEntry;
  using WT  = p3::WaterTransfer;
  using SF  = p3::SedFlux;
  using TS  = wt::TransferStage;
  using WIF = water_isotopes::WaterIsotopeFractionation;
  using WIK = water_isotopes::WaterIsotopeKinetic;
  using Mask = ekat::Mask<Pack::n>;

  const auto& gn = m_grid->name();
  const int nlev = m_num_levs;
  const int nlev_packs = ekat::npack<Pack>(nlev);
  const auto policy = TPF::get_default_team_policy(m_num_cols, nlev_packs);

  const auto qv    = get_field_in("qv", gn).get_view<const Pack**>();
  const auto qc    = get_field_in("qc", gn).get_view<const Pack**>();
  const auto qr    = get_field_in("qr", gn).get_view<const Pack**>();
  const auto qi    = get_field_in("qi", gn).get_view<const Pack**>();
  const auto dp    = get_field_in("pseudo_density", gn).get_view<const Pack**>();
  const auto entry = get_field_in("p3_water_entry", gn).get_view<const Pack***>();
  const auto xfer  = get_field_in("p3_water_transfers", gn).get_view<const Pack***>();
  const auto sed   = get_field_in("p3_sed_flux", gn).get_view<const Pack***>();

  // Checkpoint 1 needs the bulk vapor before P3 (see the class comment)
  const bool has_shoc = m_replay_shoc;
  Field::view_dev_t<const Pack***> stages;
  if (has_shoc) {
    stages = get_field_in("shoc_water_stages", gn).get_view<const Pack***>();
  }

  const auto bulk = m_bulk;

  // Fractionation of isotope tracers (see the class comment)
  const auto stage_T     = get_field_out("water_transfer_stage_T", gn).get_view<const Pack***>();
  const auto num_extrap  = get_field_out("water_transfer_alpha_extrapolations", gn).get_view<Real*>();
  const auto extrap_run  = m_alpha_extrapolations_run;
  const auto iso_params  = m_isotopes.parameters();
  const Real alpha_fixed = m_alpha_override;

  // Kinetic fractionation: the state for the saturation ratios
  const bool kinetic     = m_kinetic and m_alpha_override <= 0;
  const bool ice_model_S = m_kinetic_ice_model_S;
  Field::view_dev_t<const Pack**> p_dry, dp_dry;
  if (kinetic) {
    p_dry  = get_field_in("p_dry_mid", gn).get_view<const Pack**>();
    dp_dry = get_field_in("pseudo_density_dry", gn).get_view<const Pack**>();
  }

  for (const auto& spec : m_specs) {
    const auto tr_qv = get_field_out(wt::phase_field_name(spec.name, wt::WaterPhase::Vapor), gn)
                         .get_view<Pack**>();
    const auto tr_qc = get_field_out(wt::phase_field_name(spec.name, wt::WaterPhase::CloudLiquid), gn)
                         .get_view<Pack**>();
    const auto tr_qr = get_field_out(wt::phase_field_name(spec.name, wt::WaterPhase::Rain), gn)
                         .get_view<Pack**>();
    const auto tr_qi = get_field_out(wt::phase_field_name(spec.name, wt::WaterPhase::Ice), gn)
                         .get_view<Pack**>();
    const auto residual   = get_field_out(spec.name + "_residual_p3", gn).get_view<Pack***>();
    const auto precip_liq = get_field_out(spec.name + "_precip_liq_surf_mass", gn).get_view<Real*>();
    const auto precip_ice = get_field_out(spec.name + "_precip_ice_surf_mass", gn).get_view<Real*>();

    // Only isotopes fractionate. The species is unused unless fractionates is true.
    const bool fractionates = spec.type==wt::TracerType::Isotope;
    const auto species = fractionates ? water_isotopes::WaterIsotopes::species(spec)
                                      : water_isotopes::WaterIsotopologues::H216O;

    // One team per column: the stages before and after sedimentation are
    // local to a level, the sedimentation sweep is serial in k.
    Kokkos::parallel_for("water_tracer_transfers_p3_" + spec.name, policy,
                         KOKKOS_LAMBDA(const KT::MemberType& team) {
      const int icol = team.league_rank();

      // Padded lanes past the last level are read as 0 and never written, so
      // the math below cannot meet whatever the fields hold there.
      const auto get = [&] (const auto& v, const int k, const auto& in_range) {
        Pack x(0);
        x.set(in_range, v(icol, k));
        return x;
      };
      const auto get_cmp = [&] (const auto& v, const int cmp, const int k, const auto& in_range) {
        Pack x(0);
        x.set(in_range, v(icol, cmp, k));
        return x;
      };
      // Mass a P3 stage moved over this call [kg/kg] (the exports are rates)
      const auto moved = [&] (const int cmp, const int k, const auto& in_range) {
        return get_cmp(xfer, cmp, k, in_range)*dt;
      };

      // ---- Checkpoint 1, part1 clips, part2 rates, part2 clips ----
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&] (const int k) {
        const auto in_range = ekat::range<IntPack>(k*Pack::n) < nlev;

        Pack tr[nphases] = { get(tr_qv, k, in_range), get(tr_qc, k, in_range),
                             get(tr_qr, k, in_range), get(tr_qi, k, in_range) };
        Pack q[nphases]  = { get_cmp(entry, WE::Vapor,       k, in_range),
                             get_cmp(entry, WE::CloudLiquid, k, in_range),
                             get_cmp(entry, WE::Rain,        k, in_range),
                             get_cmp(entry, WE::Ice,         k, in_range) };

        // Checkpoint 1: bulk vapor changes between SHOC and P3 entry (qv
        // repairs), from SHOC's returned vapor, as replay_shoc computed it
        if (has_shoc) {
          const Pack qv_shoc = get_cmp(stages, WS::VaporAfterTransport,     k, in_range)
                             + get_cmp(stages, WS::LiquidAfterTransport,    k, in_range)
                             - get_cmp(stages, WS::LiquidAfterCondensation, k, in_range);
          tr[iv] += wt::vapor_repair(tr[iv], qv_shoc, q[iv] - qv_shoc);
        }

        // part1 clips, in P3's order (I2V and I2R exclude each other)
        clip_stage(tr, q, ic, iv, moved(WT::Part1ClipC2V, k, in_range));
        clip_stage(tr, q, ir, iv, moved(WT::Part1ClipR2V, k, in_range));
        clip_stage(tr, q, ii, iv, moved(WT::Part1ClipI2V, k, in_range));
        clip_stage(tr, q, ii, ir, moved(WT::Part1ClipI2R, k, in_range));

        // part2 process rates, all applied together from the post-part1
        // state, as P3 applies them. Each moves the tracer at its donor's
        // ratio (times its fractionation factor, for isotopes). P3 limits a
        // species' sinks to its reservoir plus its sources in the same call
        // (p3_conservation_impl.hpp), so a sink can far exceed the reservoir
        // (a tiny rain reservoir fed by autoconversion can be collected by ice
        // many times over); the reservoir's own ratio would then amplify its
        // roundoff by sink/reservoir, every call. The donor is therefore the
        // mass available to it over the call, its reservoir plus its inflows:
        //   R_x = (tr_x + sum of tracer inflows)/(q_x + sum of bulk inflows).
        // The tracer inflows depend on the upstream donors' ratios, so the
        // ratios are refined by fixed-point passes starting from the
        // reservoirs' (num_ratio_passes). Exact for a uniform ratio.
        //   Donor and receiver of each component, in WaterTransfer order:
        //   VapdepV2I, NucleatV2I, SublimI2V, EvapR2V, AccretC2R, AutoconvC2R,
        //   ShedC2R, BergC2I, CollectC2I, FreezeC2I, CollectR2I, FreezeR2I,
        //   MeltI2R
        const int from[num_part2] = {iv, iv, ii, ir, ic, ic, ic, ic, ic, ic, ir, ir, ii};
        const int to[num_part2]   = {ii, ii, iv, iv, ir, ir, ir, ii, ii, ii, ii, ii, ir};

        Pack amount[num_part2];
        Pack q_in[nphases];  // bulk inflows [kg/kg] (Packs are zero-initialized)
        for (int j=0; j<num_part2; ++j) {
          amount[j] = moved(first_part2 + j, k, in_range);
          q_in[to[j]] += amount[j];
        }

        // The fractionation factor of each rate (see wt::part2_beta): 1 for
        // a tracer that does not fractionate. For an isotope, alpha_l and
        // alpha_i at P3's part2 temperature, evaluated on the lanes where a
        // rate that uses them moved water (elsewhere they are not used).
        Pack beta[num_part2];
        for (int j=0; j<num_part2; ++j) {
          beta[j] = Pack(1);
        }
        if (fractionates) {
          Pack alpha_l(alpha_fixed), alpha_i(alpha_fixed);
          if (alpha_fixed <= 0) {
            const Pack T = get_cmp(stage_T, TS::P3Part2, k, in_range);
            const auto uses = [&] (const int cmp) { return amount[cmp - first_part2] > Real(0); };
            const auto need_l = in_range && (uses(WT::Part2EvapR2V) || uses(WT::Part2BergC2I));
            const auto need_i = in_range && (uses(WT::Part2VapdepV2I) || uses(WT::Part2NucleatV2I) ||
                                             uses(WT::Part2BergC2I));
            int nx = 0;
            alpha_l = WIF::alpha_liquid_vapor(T, species, water_isotopes::CondensedOverVapor, iso_params,
                                              need_l, &nx);
            alpha_i = WIF::alpha_ice_vapor(T, species, water_isotopes::CondensedOverVapor, iso_params,
                                           need_i, &nx);
            if (nx > 0) {
              Kokkos::atomic_add(&num_extrap(icol), Real(nx));
              Kokkos::atomic_add(&extrap_run(icol), Real(nx));
            }

            // Kinetic fractionation, from P3's state before its part2 rates:
            // the bulk vapor after the part1 clips (q[iv]) at P3Part2. (P3
            // takes its qv_sat from its entry temperature, before the part1
            // clips, which differs from P3Part2 by < 1e-4 K in the tests.)
            // alpha_l only serves liquid that evaporates (rain evaporation,
            // and Bergeron's liquid) and alpha_i only vapor that deposits, so
            // both become the effective factors (see wt::part2_beta).
            if (kinetic) {
              constexpr Real tkini = WIK::tkini;
              if (need_l.any()) {
                const Pack rh_l = saturation_ratio(q[iv], T, get(p_dry, k, in_range), get(dp, k, in_range),
                                                   get(dp_dry, k, in_range), false, need_l);
                alpha_l = WIK::alpha_liquid_vapor_kinetic(alpha_l, rh_l, species, iso_params, need_l);
              }
              const auto cold_i = need_i && T < tkini;
              if (cold_i.any()) {
                const Pack S_i = ice_model_S
                    ? saturation_ratio(q[iv], T, get(p_dry, k, in_range), get(dp, k, in_range),
                                       get(dp_dry, k, in_range), true, cold_i)
                    : WIK::ice_supersaturation_empirical(T);
                alpha_i = WIK::alpha_ice_vapor_kinetic(alpha_i, T, S_i, species, iso_params, cold_i);
              }
            }
          }
          for (int j=0; j<num_part2; ++j) {
            beta[j] = wt::part2_beta(first_part2 + j, alpha_l, alpha_i);
          }
        }

        // The donors' tracer and bulk water: first the reservoirs, then the
        // reservoirs plus the inflows of the previous pass. Each donor's
        // sinks move tracer through wt::donor_transfers (the linear transfer
        // where all of the donor's betas are 1, else a Rayleigh step); the
        // components that are not its sinks enter as amount 0, beta 1.
        Pack tr_d[nphases], q_d[nphases], dtr[num_part2];
        Mask rayleigh[nphases];  // lanes where a donor took the Rayleigh step
        for (int p=0; p<nphases; ++p) {
          tr_d[p] = tr[p];
          q_d[p]  = q[p];
          rayleigh[p] = Mask(false);
        }
        for (int pass=0; pass<=num_ratio_passes; ++pass) {
          for (int p=0; p<nphases; ++p) {
            Pack a[num_part2], b[num_part2], out[num_part2];
            for (int j=0; j<num_part2; ++j) {
              a[j] = from[j]==p ? amount[j] : Pack(0);
              b[j] = from[j]==p ? beta[j]   : Pack(1);
            }
            rayleigh[p] = wt::donor_transfers(tr_d[p], q_d[p], a, b, out);
            for (int j=0; j<num_part2; ++j) {
              if (from[j]==p) dtr[j] = out[j];
            }
          }
          if (pass < num_ratio_passes) {
            Pack tr_in[nphases];
            for (int j=0; j<num_part2; ++j) {
              tr_in[to[j]] += dtr[j];
            }
            for (int p=0; p<nphases; ++p) {
              tr_d[p] = tr[p] + tr_in[p];
              q_d[p]  = q[p]  + q_in[p];
            }
          }
        }

        // Guard of the Rayleigh donors. A Rayleigh step never moves more than
        // the donor's tracer plus its inflows, but the inflows are those of
        // the previous pass: where they changed in the last pass (a cycle of
        // donors that the passes have not converged, or roundoff) the donor's
        // outflows can exceed its tracer plus its final inflows. Scale them
        // down to that, which conserves the tracer, one donor at a time along
        // the longest chain (c, r, i, v: a donor's inflows are then final
        // unless they come from a later donor that is scaled too, i.e. from
        // vapor, which does not run dry). Linear donors are left as they
        // are, as for non-fractionating tracers.
        constexpr int guard_order[nphases] = {ic, ir, ii, iv};
        for (const int p : guard_order) {
          if (not rayleigh[p].any()) continue;
          Pack outflow(0), inflow(0);
          for (int j=0; j<num_part2; ++j) {
            if (from[j]==p) outflow += dtr[j];
            if (to[j]==p)   inflow  += dtr[j];
          }
          const Pack avail = ekat::max(tr[p] + inflow, Pack(0));
          const auto over = rayleigh[p] && outflow > avail;
          if (over.any()) {
            Pack outflow_safe = outflow;
            outflow_safe.set(!over, Real(1));
            const Pack scale = avail/outflow_safe;
            for (int j=0; j<num_part2; ++j) {
              if (from[j]==p) dtr[j].set(over, dtr[j]*scale);
            }
          }
        }

        for (int j=0; j<num_part2; ++j) {
          tr[from[j]] -= dtr[j];
          tr[to[j]]   += dtr[j];
          q[from[j]]  -= amount[j];
          q[to[j]]    += amount[j];
        }
        // Clamp what is left below 0 in a Rayleigh donor (roundoff of the
        // guard, or of a donor that came in at -roundoff): not conservative,
        // but only of roundoff size
        for (int p=0; p<nphases; ++p) {
          tr[p].set(rayleigh[p] && tr[p] < Real(0), Real(0));
        }

        // part2 clips, from the post-update state
        clip_stage(tr, q, ic, iv, moved(WT::Part2ClipC2V, k, in_range));
        clip_stage(tr, q, ir, iv, moved(WT::Part2ClipR2V, k, in_range));
        clip_stage(tr, q, ii, iv, moved(WT::Part2ClipI2V, k, in_range));

        tr_qv(icol, k).set(in_range, tr[iv]);
        tr_qc(icol, k).set(in_range, tr[ic]);
        tr_qr(icol, k).set(in_range, tr[ir]);
        tr_qi(icol, k).set(in_range, tr[ii]);
        for (int p=0; p<nphases; ++p) {
          bulk(icol, p, k) = q[p];
        }
      });
      team.team_barrier();

      // ---- Sedimentation of cloud, rain and ice, and the tracer's surface
      // precipitation (cloud + rain is liquid, as in P3) ----
      Kokkos::single(Kokkos::PerTeam(team), [&] () {
        const auto s_dp = ekat::scalarize(ekat::subview(dp, icol));
        const auto sweep = [&] (const auto& tr_view, const int phase, const int species) {
          return sedimentation_sweep(ekat::scalarize(ekat::subview(tr_view, icol)),
                                     ekat::scalarize(ekat::subview(bulk, icol, phase)),
                                     ekat::scalarize(ekat::subview(sed, icol, species)),
                                     s_dp, nlev, dt);
        };
        const Real liq = sweep(tr_qc, ic, SF::CloudLiquid) + sweep(tr_qr, ir, SF::Rain);
        const Real ice = sweep(tr_qi, ii, SF::Ice);
        // [kg/m2]: the surface fluxes [kg/m2/s] over this call
        precip_liq(icol) += liq*dt;
        precip_ice(icol) += ice*dt;
      });
      team.team_barrier();

      // ---- Homogeneous freezing, part3 clips, checkpoint 2 ----
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&] (const int k) {
        const auto in_range = ekat::range<IntPack>(k*Pack::n) < nlev;

        Pack tr[nphases] = { get(tr_qv, k, in_range), get(tr_qc, k, in_range),
                             get(tr_qr, k, in_range), get(tr_qi, k, in_range) };
        Pack q[nphases];
        for (int p=0; p<nphases; ++p) {
          q[p] = bulk(icol, p, k);
        }

        clip_stage(tr, q, ic, ii, moved(WT::HomFreezeC2I, k, in_range));
        clip_stage(tr, q, ir, ii, moved(WT::HomFreezeR2I, k, in_range));

        clip_stage(tr, q, ic, iv, moved(WT::Part3ClipC2V, k, in_range));
        clip_stage(tr, q, ir, iv, moved(WT::Part3ClipR2V, k, in_range));
        clip_stage(tr, q, ii, iv, moved(WT::Part3ClipI2V, k, in_range));

        // The bulk water P3 returned, minus the tracer
        const Pack q_out[nphases] = { get(qv, k, in_range), get(qc, k, in_range),
                                      get(qr, k, in_range), get(qi, k, in_range) };
        for (int p=0; p<nphases; ++p) {
          residual(icol, p, k) = q_out[p] - tr[p];
        }

        // Checkpoint 2: bulk vapor changes no export accounts for (P3's qv
        // postcondition repair, wet/dry roundoff)
        tr[iv] += wt::vapor_repair(tr[iv], q[iv], q_out[iv] - q[iv]);

        tr_qv(icol, k).set(in_range, tr[iv]);
        tr_qc(icol, k).set(in_range, tr[ic]);
        tr_qr(icol, k).set(in_range, tr[ir]);
        tr_qi(icol, k).set(in_range, tr[ii]);
      });
    });
  }
}

} // namespace scream
