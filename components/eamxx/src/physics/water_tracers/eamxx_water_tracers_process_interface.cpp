#include "eamxx_water_tracers_process_interface.hpp"

#include <ekat_assert.hpp>

namespace scream
{

// =========================================================================================
WaterTracers::WaterTracers(const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
  , m_tracer_count(0)
  , m_qv_iso_group("qv_iso")
  , m_qc_iso_group("qc_iso")
  , m_qi_iso_group("qi_iso")
  , m_qr_iso_group("qr_iso")
{
  // Read tracer count from parameter list (default 0)
  m_tracer_count = m_params.get<int>("tracer_count", 0);

  // Allocate storage for tracer names
  if (m_tracer_count > 0) {
    m_tracer_names.resize(m_tracer_count);
  }
}

// =========================================================================================
void WaterTracers::create_requests()
{
  // Get the grid for this process
  m_grid = m_grids_manager->get_grid("physics");
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Skip field registration if tracer_count is 0 (disabled)
  if (m_tracer_count == 0) {
    return;
  }

  // Register tracer fields for each component
  // Each isotope component is registered as a separate scalar tracer
  // that joins the "tracers" group (for dycore advection) and a
  // family-specific group (for multi-component indexing)
  const int pack_size = 1;  // Default pack size
  register_tracer_fields(m_tracer_count, m_grid, pack_size);
}

// =========================================================================================
void WaterTracers::initialize_impl(const RunType /* run_type */)
{
  // Skip initialization if tracer_count is 0
  if (m_tracer_count == 0) {
    return;
  }

  // Attach metadata to tracer fields
  attach_tracer_metadata();

  // Retrieve field groups for multi-component access
  retrieve_tracer_groups();

  // Initialize all tracer fields to zero
  for (const auto& names_vec : m_tracer_names) {
    for (const auto& name : names_vec) {
      auto field = get_field_out(name);
      field.deep_copy(0.0);
      field.sync_to_host();
    }
  }
}

// =========================================================================================
void WaterTracers::run_impl(const double /* dt */)
{
  // No-op: tracers are advected by the dynamics and turbulence processes
  // via their membership in the "tracers" group.
  // Fractionation physics will be added in subsequent specs of the water isotope campaign.
  // For now, tracers maintain their zero-initialized values, ensuring no impact
  // on existing model physics.
}

// =========================================================================================
void WaterTracers::finalize_impl()
{
  // Nothing to finalize at this stage
}

} // namespace scream
