#include "physics/p3/eamxx_p3_process_interface.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"
// Needed for p3_init, the only F90 code still used.
#include "physics/p3/p3_functions.hpp"
#include "physics/share/physics_constants.hpp"
#include "physics/p3/p3_f90.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{

// =========================================================================================
WaterTracers::WaterTracers (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}



}
