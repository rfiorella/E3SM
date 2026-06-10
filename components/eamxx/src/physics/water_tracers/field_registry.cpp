#include "physics/water_tracers/eamxx_water_tracers_process_interface.hpp"

namespace scream {

void WaterTracers::register_tracer_fields(
    const int tracer_count,
    const std::shared_ptr<const AbstractGrid>& grid,
    const int pack_size)
{
  using namespace ekat::units;

  const auto kg_per_kg = kg/kg;

  // Field families to register: qv, qc, qi, qr
  // Each component is a separate scalar tracer that joins the "tracers" group
  // and a family-specific group for multi-component access

  for (int k = 0; k < tracer_count; ++k) {
    // Vapor isotope tracers
    const std::string qv_name = "qv_iso_" + std::to_string(k);
    add_tracer<Updated>(qv_name, grid, kg_per_kg, pack_size);

    // Cloud liquid isotope tracers
    const std::string qc_name = "qc_iso_" + std::to_string(k);
    add_tracer<Updated>(qc_name, grid, kg_per_kg, pack_size);

    // Cloud ice isotope tracers
    const std::string qi_name = "qi_iso_" + std::to_string(k);
    add_tracer<Updated>(qi_name, grid, kg_per_kg, pack_size);

    // Rain isotope tracers
    const std::string qr_name = "qr_iso_" + std::to_string(k);
    add_tracer<Updated>(qr_name, grid, kg_per_kg, pack_size);

    // Store names for later retrieval
    m_tracer_names[k] = {qv_name, qc_name, qi_name, qr_name};
  }

  // Request family groups for multi-component access
  // These groups will contain all components of each phase
  add_group<Updated>("qv_iso", grid->name(), pack_size);
  add_group<Updated>("qc_iso", grid->name(), pack_size);
  add_group<Updated>("qi_iso", grid->name(), pack_size);
  add_group<Updated>("qr_iso", grid->name(), pack_size);
}

void WaterTracers::attach_tracer_metadata()
{
  // Attach long_name metadata to each tracer field
  for (const auto& names_vec : m_tracer_names) {
    for (size_t phase_idx = 0; phase_idx < names_vec.size(); ++phase_idx) {
      const auto& name = names_vec[phase_idx];
      auto field = m_fields_mgr->get_field(name);
      auto& header = field.get_header();

      // Set long_name based on phase
      std::string long_name;
      if (phase_idx == 0) {
        long_name = "Water vapor isotope tracer mixing ratio";
      } else if (phase_idx == 1) {
        long_name = "Cloud liquid water isotope tracer mixing ratio";
      } else if (phase_idx == 2) {
        long_name = "Cloud ice isotope tracer mixing ratio";
      } else if (phase_idx == 3) {
        long_name = "Rain isotope tracer mixing ratio";
      }

      header.get_tracking().get_metadata().set_string("long_name", long_name);
    }
  }
}

void WaterTracers::retrieve_tracer_groups()
{
  // Retrieve field groups for multi-component access
  m_qv_iso_group = m_fields_mgr->get_field_group("qv_iso", m_grid->name());
  m_qc_iso_group = m_fields_mgr->get_field_group("qc_iso", m_grid->name());
  m_qi_iso_group = m_fields_mgr->get_field_group("qi_iso", m_grid->name());
  m_qr_iso_group = m_fields_mgr->get_field_group("qr_iso", m_grid->name());

  // Add individual tracer fields to their family groups
  for (const auto& names_vec : m_tracer_names) {
    m_fields_mgr->add_to_group(names_vec[0], "qv_iso");  // vapor
    m_fields_mgr->add_to_group(names_vec[1], "qc_iso");  // cloud liquid
    m_fields_mgr->add_to_group(names_vec[2], "qi_iso");  // cloud ice
    m_fields_mgr->add_to_group(names_vec[3], "qr_iso");  // rain
  }
}

} // namespace scream
