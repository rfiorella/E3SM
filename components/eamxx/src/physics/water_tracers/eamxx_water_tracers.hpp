#ifndef SCREAM_WATER_TRACERS_HPP
#define SCREAM_WATER_TRACERS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "physics/water_tracers/water_tracers.hpp"
#include "share/util/scream_common_physics_functions.hpp"

#include <string>

namespace scream
{
// =========================================================================================
WaterTracers::WaterTracers (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}

// Set initial Variable Declarations
    bool trace_water;
    bool wisotope;
    bool wtrc_alpha_kinetic;
    bool wtrc_use_ice_supsat;

    Real wtrc_fixed_alpha = 1.
    Real wtrc_fixed_rstd  = 1.
    Real wtrc_qchkmin = 1e-14
    Real wtrc_qmin    = 1e-18

    int wtrc_ncnst = 0 // number of water tracer sets


}

#endif // SCREAM_WATER_TRACERS_HPP