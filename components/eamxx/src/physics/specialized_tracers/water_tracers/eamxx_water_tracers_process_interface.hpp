#ifndef EAMXX_WATER_TRACERS_HPP
#define EAMXX_WATER_TRACERS_HPP

#include "eamxx_water_tracers_registry.hpp"
#include "../water_isotopes/eamxx_water_isotopes.hpp"

#include "share/atm_process/atmosphere_process.hpp"
#include <ekat_parameter_list.hpp>

#include <string>
#include <vector>

namespace scream
{
/*
 * Water tracers: extra water species that follow the bulk water (qv, qc, qr,
 * qi) through the model.
 *
 * The tracers are defined in a YAML file (the water_tracers_file parameter; see
 * eamxx_water_tracers_registry.hpp). Each one is registered as four advected
 * tracers <name>_qv, <name>_qc, <name>_qr, <name>_qi [kg/kg], one per bulk
 * water phase, plus the surface flux of its vapor, <name>_qv_surf_flux
 * [kg/m2/s, upward positive], in the group "tracer_surface_fluxes". The flux
 * is the tracer's share of the surface evaporation: the ocean part is computed
 * every step from the ocean bulk state imported from the coupler (ocn_ssq,
 * ocn_re, ocn_ustar) and the bottom-level atmosphere, and the rest of the
 * merged surf_evap (land and sea ice) evaporates at the tracer's r_land (see
 * surface_tracer_flux).
 *
 * This process only computes the fluxes; it never changes bulk water or a
 * tracer after initialization. SHOC applies them (in the same implicit solve
 * as the bulk surface evaporation, when the driver sets its parameter
 * apply_tracer_surface_fluxes), so water_tracers must come before shoc in the
 * atm_procs_list: it reads surf_evap at the start of the step.
 *
 * Tracer kinds: generic (bulk_proportional or bulk_formula) and isotopes.
 * Isotope entries are only accepted when EAMXX_HAS_WATER_ISOTOPES is
 * defined; isotope-specific state and math live in the composed
 * water_isotopes::WaterIsotopes member, m_isotopes.
 */

class WaterTracers : public AtmosphereProcess
{
public:
  using KT = ekat::KokkosTypes<DefaultDevice>;

  WaterTracers (const ekat::Comm& comm, const ekat::ParameterList& params);

  AtmosphereProcessType type () const override { return AtmosphereProcessType::Physics; }

  std::string name () const override { return "water_tracers"; }

  void create_requests () override;

protected:

#ifdef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
public:
#endif

  void initialize_impl (const RunType run_type) override;
  void run_impl        (const double dt) override;
  void finalize_impl   () override {}

  // Steps of run_impl
  void compute_bulk_state ();
  void compute_tracer_fluxes ();

  // Sets every phase of every init: ratio tracer to r_init times its bulk
  // phase (e.g. <name>_qc = r_init*qc) on the given grid
  void apply_ratio_initial_conditions (const std::string& grid_name);

protected:

  int m_num_cols;
  int m_num_levs;

  std::shared_ptr<const AbstractGrid> m_grid;

  // One entry per tracer in water_tracers_file, sorted by name. r_ocn/r_land/r_init and
  // ocean_flux for isotope specs are resolved by m_isotopes.resolve() in the
  // constructor, before initialize/run time reads them.
  std::vector<water_tracers::TracerSpec> m_specs;

  // Isotope-specific state (species table, runtime formulation options).
  // Always constructed; inert unless EAMXX_HAS_WATER_ISOTOPES is defined and
  // the tracers file defines isotope entries. See compute_tracer_fluxes().
  water_isotopes::WaterIsotopes m_isotopes;

  // Per-column bottom-level bulk state, filled at the start of run_impl:
  //  - m_sfc_f: ocnfrac, set to 0 where the ocean state is not valid
  //  - m_sfc_ustar, m_sfc_re, m_sfc_ssq: ocn_ustar, ocn_re, ocn_ssq, set to
  //    0 where the ocean state is not valid (so they never hold fill values)
  //  - m_sfc_rho: bottom-level air density [kg/m3]
  //  - m_sfc_evap: E_ocn = bulk_ocean_flux(f,rho,u,r,s,q) [kg/m2/s], upward
  //    positive, exactly 0 where the ocean state is not valid
  // compute_tracer_fluxes reads qv, <name>_qv and sst from their fields directly.
  KT::view_1d<Real> m_sfc_f, m_sfc_rho, m_sfc_ustar, m_sfc_re, m_sfc_ssq,
                    m_sfc_evap;

private:

  // The grid on which the driver sets initial conditions, if it differs from
  // m_grid (physics_gll in a pgN run), else null. See create_requests.
  std::shared_ptr<const AbstractGrid> m_ic_grid;

}; // class WaterTracers

} // namespace scream

#endif // EAMXX_WATER_TRACERS_HPP
