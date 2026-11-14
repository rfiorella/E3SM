#include "eamxx_water_tracer_transfers_process_interface.hpp"
#include "eamxx_water_tracer_transfers_functions.hpp"
#include "../water_isotopes/eamxx_water_isotopes.hpp"

#include "physics/p3/p3_water_transfers.hpp"
#include "physics/shoc/shoc_water_transfers.hpp"
#include "share/physics/physics_constants.hpp"

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

// P3's part2 process rates, in WaterTransfer order from Part2VapdepV2I
constexpr int first_part2 = p3::WaterTransfer::Part2VapdepV2I;
constexpr int num_part2   = p3::WaterTransfer::Part2MeltI2R - first_part2 + 1;
static_assert(num_part2 == 13, "Update the part2 tables below to match p3_water_transfers.hpp.");

// The tracer moved by part2 component `cmp` [kg/kg], given its donor's tracer
// and bulk water, and the bulk mass moved, through its transfer class
KOKKOS_INLINE_FUNCTION
Pack part2_transfer (const int cmp, const Pack& tr_d, const Pack& q_d, const Pack& amount)
{
  using WT = p3::WaterTransfer;
  switch (cmp) {
    case WT::Part2VapdepV2I:
    case WT::Part2NucleatV2I:  return wt::ice_deposition(tr_d, q_d, amount);
    case WT::Part2SublimI2V:   return wt::ice_sublimation(tr_d, q_d, amount);
    case WT::Part2EvapR2V:     return wt::rain_evaporation(tr_d, q_d, amount);
    case WT::Part2AccretC2R:
    case WT::Part2AutoconvC2R:
    case WT::Part2ShedC2R:     return wt::liquid_collection(tr_d, q_d, amount);
    case WT::Part2BergC2I:     return wt::bergeron(tr_d, q_d, amount);
    case WT::Part2CollectC2I:
    case WT::Part2FreezeC2I:
    case WT::Part2CollectR2I:
    case WT::Part2FreezeR2I:   return wt::liquid_freezing(tr_d, q_d, amount);
    default:                   return wt::ice_melting(tr_d, q_d, amount);  // Part2MeltI2R
  }
}

// Fixed-point passes that refine the part2 donor ratios (see replay_p3). The
// longest chain of donors without a cycle is c -> r -> i -> v, so three
// passes carry the ratio of a well-filled reservoir to the end of any chain.
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

  // Bulk water after each SHOC stage (SHOC with export_water_transfers=true)
  if (m_replay_shoc) {
    const auto stages_layout = m_grid->get_3d_vector_layout(LEV, shoc::WaterStage::Count, "shoc_water_stage");
    add_field<Required>("shoc_water_stages", stages_layout, kg/kg, grid_name, ps);
  }

  // P3's entry state, transfers and sedimentation fluxes (P3 with
  // export_water_transfers=true; same layouts as P3 registers them)
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
void WaterTracerTransfers::run_impl (const double dt)
{
  if (m_replay_shoc) {
    replay_shoc();
  }
  if (m_replay_p3) {
    replay_p3(dt);
  }
}

// =========================================================================================
void WaterTracerTransfers::replay_shoc ()
{
  using MDPolicy = Kokkos::MDRangePolicy<KT::ExeSpace, Kokkos::Rank<2>>;
  using WS = shoc::WaterStage;

  const auto& gn = m_grid->name();
  const int nlev_packs = ekat::npack<Pack>(m_num_levs);
  const MDPolicy policy({0, 0}, {m_num_cols, nlev_packs});

  const auto stages = get_field_in("shoc_water_stages", gn).get_view<const Pack***>();
  const auto qr     = get_field_in("qr", gn).get_view<const Pack**>();
  const auto qi     = get_field_in("qi", gn).get_view<const Pack**>();

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

    Kokkos::parallel_for("water_tracer_transfers_shoc_" + spec.name, policy,
                         KOKKOS_LAMBDA(const int icol, const int k) {
      // Bulk water after SHOC's turbulent transport, and the cloud liquid
      // after its condensation [kg/kg]. The tracers come in already
      // transported (SHOC diffuses them in the same solve).
      const Pack qv_trans = stages(icol, WS::VaporAfterTransport,     k);
      const Pack qc_trans = stages(icol, WS::LiquidAfterTransport,    k);
      const Pack qc_cond  = stages(icol, WS::LiquidAfterCondensation, k);

      // Net condensation [kg/kg]: > 0 moves vapor to liquid at the vapor's
      // ratio, < 0 moves liquid to vapor at the liquid's ratio
      const Pack cond = qc_cond - qc_trans;
      const auto evap = cond < Real(0);
      Pack tr_donor = tr_qv(icol, k);
      Pack q_donor  = qv_trans;
      tr_donor.set(evap, tr_qc(icol, k));
      q_donor.set(evap, qc_trans);

      const Pack dtr = wt::shoc_condensation(tr_donor, q_donor, cond);
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
        // ratio. P3 limits a species' sinks to its reservoir plus its sources
        // in the same call (p3_conservation_impl.hpp), so a sink can far
        // exceed the reservoir (a tiny rain reservoir fed by autoconversion
        // can be collected by ice many times over); the reservoir's own ratio
        // would then amplify its roundoff by sink/reservoir, every call. The
        // donor ratio is therefore that of the mass available to the donor
        // over the call, its reservoir plus its inflows:
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

        // The donors' tracer and bulk water: first the reservoirs, then the
        // reservoirs plus the inflows of the previous pass
        Pack tr_d[nphases], q_d[nphases], dtr[num_part2];
        for (int p=0; p<nphases; ++p) {
          tr_d[p] = tr[p];
          q_d[p]  = q[p];
        }
        for (int pass=0; pass<=num_ratio_passes; ++pass) {
          for (int j=0; j<num_part2; ++j) {
            dtr[j] = part2_transfer(first_part2 + j, tr_d[from[j]], q_d[from[j]], amount[j]);
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
        for (int j=0; j<num_part2; ++j) {
          tr[from[j]] -= dtr[j];
          tr[to[j]]   += dtr[j];
          q[from[j]]  -= amount[j];
          q[to[j]]    += amount[j];
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
