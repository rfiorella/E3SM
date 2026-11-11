#include "eamxx_water_tracer_transfers_process_interface.hpp"
#include "eamxx_water_tracer_transfers_functions.hpp"
#include "../water_isotopes/eamxx_water_isotopes.hpp"

#include "physics/shoc/shoc_water_transfers.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>

namespace scream
{

namespace {

namespace wt = water_tracers;

} // anonymous namespace

// =========================================================================================
WaterTracerTransfers::WaterTracerTransfers (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  m_replay_shoc = m_params.get<bool>("replay_shoc", true);
  m_replay_p3   = m_params.get<bool>("replay_p3",   true);

  EKAT_REQUIRE_MSG(not m_replay_p3,
      "Error! Process '" + m_params.name() + "': the P3 replay is not implemented yet.\n"
      "  Set replay_p3 to false.\n");
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
  }
}

// =========================================================================================
void WaterTracerTransfers::run_impl (const double /* dt */)
{
  if (m_replay_shoc) {
    replay_shoc();
  }
}

// =========================================================================================
void WaterTracerTransfers::replay_shoc ()
{
  using MDPolicy = Kokkos::MDRangePolicy<KT::ExeSpace, Kokkos::Rank<2>>;
  using WS = shoc::WaterStage;

  constexpr int iv = etoi(wt::WaterPhase::Vapor);
  constexpr int ic = etoi(wt::WaterPhase::CloudLiquid);
  constexpr int ir = etoi(wt::WaterPhase::Rain);
  constexpr int ii = etoi(wt::WaterPhase::Ice);

  const auto& gn = m_grid->name();
  const int nlev_packs = ekat::npack<Pack>(m_num_levs);
  const MDPolicy policy({0, 0}, {m_num_cols, nlev_packs});

  const auto stages = get_field_in("shoc_water_stages", gn).get_view<const Pack***>();
  const auto qr     = get_field_in("qr", gn).get_view<const Pack**>();
  const auto qi     = get_field_in("qi", gn).get_view<const Pack**>();

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
      // the bulk fields stand in for them: valid only while nothing between
      // SHOC and this process changes them.
      residual(icol, iv, k) = (qv_trans + qc_trans - qc_cond) - tr_qv(icol, k);
      residual(icol, ic, k) = qc_cond - tr_qc(icol, k);
      residual(icol, ir, k) = qr(icol, k) - tr_qr(icol, k);
      residual(icol, ii, k) = qi(icol, k) - tr_qi(icol, k);
    });
  }
}

} // namespace scream
