#ifndef EAMXX_WATER_TRACERS_REGISTRY_HPP
#define EAMXX_WATER_TRACERS_REGISTRY_HPP

#include "eamxx_water_tracers_functions.hpp"  // OceanFlux

#include "share/core/eamxx_types.hpp"
#include "share/util/eamxx_utils.hpp"  // etoi

#include <ekat_parameter_list.hpp>

#include <iterator>
#include <string>
#include <vector>

namespace scream {
namespace water_tracers {

/*
 * Host-only parser for the water tracer YAML file (the water_tracers_file parameter
 * of the water_tracers process).
 *
 * The file has one top-level key, `tracers`, a map keyed by tracer name, and
 * freely mixes generic and isotope entries. Flow style and block style are
 * both accepted:
 *
 *   tracers:
 *     wtrc_bp:  {type: generic, ocean_flux: bulk_proportional, r_ocn: 1.0}
 *     wtrc_bf:  {type: generic, ocean_flux: bulk_formula, r_ocn: 1.0, init: ratio}
 *     wtrc_hdo: {type: isotope, species: HDO}
 *
 * Each tracer follows all of the bulk water, one field per phase (see
 * WaterPhase): tracer wtrc_hdo has fields wtrc_hdo_qv, wtrc_hdo_qc,
 * wtrc_hdo_qr and wtrc_hdo_qi. Tracer names must start with `wtrc_` (so
 * none of these can alias any other field) and must not end in
 * `_surf_flux` (each tracer also registers <name>_surf_flux).
 *
 * Keys by type (anything else is an error, so a misspelled key cannot be
 * silently ignored):
 *   generic:  type, ocean_flux (bulk_proportional [default] | bulk_formula),
 *             r_ocn (default 1), init (ratio [default] | zero),
 *             r_init (default r_ocn)
 *   isotope:  type, species (required), init (default ratio), r_init
 *             (default: the species' ocean ratio, resolved by the process at
 *             construction). Only accepted when the process is built with
 *             EAMXX_ENABLE_WATER_ISOTOPES=ON.
 *
 * The returned specs are sorted by tracer name: ekat::ParameterList stores
 * sublists in a std::map. Nothing downstream depends on the order.
 */

// What the calling process can handle.
struct RegistryOptions {
  bool accept_isotopes = false;
  // Valid `species:` values for isotope entries (used only if accept_isotopes)
  std::vector<std::string> isotope_species;
};

enum class TracerType { Generic, Isotope };

// The bulk water species a water tracer follows. Each tracer has one field
// per phase, named <tracer name>_<bulk field name>, e.g. wtrc_hdo_qc.
enum class WaterPhase { Vapor, CloudLiquid, Rain, Ice, Count };

constexpr WaterPhase water_phases[] = {
  WaterPhase::Vapor, WaterPhase::CloudLiquid, WaterPhase::Rain, WaterPhase::Ice
};
static_assert(std::size(water_phases) == static_cast<std::size_t>(etoi(WaterPhase::Count)),
              "water_phases must list every WaterPhase.");

// The bulk field a phase follows: qv, qc, qr or qi.
std::string bulk_field_name (WaterPhase phase);

// The field of tracer `name` for a phase, e.g. wtrc_hdo_qc.
std::string phase_field_name (const std::string& name, WaterPhase phase);

// How a tracer is initialized on an Initial run:
//  - Ratio: every phase is r_init times its bulk phase (<name>_qv = r_init*qv,
//    <name>_qc = r_init*qc, ...).
//  - Zero:  the tracer keeps the initial conditions the driver gave it (buildnml
//    sets initial_conditions::<name>_<phase> to 0.0 unless the user sets it).
enum class InitType { Ratio, Zero };

struct TracerSpec {
  std::string name;
  TracerType  type = TracerType::Generic;

  // Isotope tracers only: one of RegistryOptions::isotope_species
  // (water_isotopes::WaterIsotopes::species() maps it to the species enum)
  std::string species;

  // Ocean ratio R_ocn. Generic: from the file (default 1). Isotope: set by the
  // water_isotopes process from its parameter tables.
  Real r_ocn = 1;

  // Initialization (see InitType). r_init_given records whether the file set
  // r_init, so a subclass that resolves r_ocn later knows whether to default
  // r_init to it.
  InitType init = InitType::Ratio;
  bool r_init_given = false;
  Real r_init = 1;

  // Generic: from the file (default bulk_proportional). Isotope: set to
  // BulkFormula by the water_isotopes process, like r_ocn.
  OceanFlux ocean_flux = OceanFlux::BulkProportional;

  // The tracer's raw entry from the file. Extension point for keys added by
  // later changes; nothing reads it yet.
  ekat::ParameterList extra;
};

// Parse and validate a tracer file. Errors (std::runtime_error, via
// EKAT_REQUIRE_MSG) name the file, the tracer and the offending key.
std::vector<TracerSpec> parse_tracers_file (const std::string& path,
                                            const RegistryOptions& opts);

} // namespace water_tracers
} // namespace scream

#endif // EAMXX_WATER_TRACERS_REGISTRY_HPP
