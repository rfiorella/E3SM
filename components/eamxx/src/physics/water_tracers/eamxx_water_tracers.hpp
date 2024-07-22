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

    static const Real wtrc_fixed_alpha = 1.; // default standard fractionation factor
    static const Real wtrc_fixed_rstd  = 1.; // default standard water isotope ratio (1 chosen for numerics)
    static const Real wtrc_qchkmin = 1e-14;
    static const Real wtrc_qmin    = 1e-18;

    int wtrc_ncnst = 0; // number of water tracer sets
    int wtrc_niter = 10; // number of iterations when applying process rates


    // Tunable parameters related to water isotopic fraciontation schemes.
    static const Real dkfac   = 0.58; // diffusive evaporative kinetic power law exponent (Stewart 1975)
    static const Real ssatmx  = 2.0;  // maximum value of supersaturation parameters
    static const Real fsata   = 1.;   // calculate supersaturation parameter s = fsata + fsatb
    static const Real fsatb   = -0.002; // supersaturation parameter s = fsata + fsatb 

    static const double hlhto = 3.888e8; // seconds

}

#endif // SCREAM_WATER_TRACERS_HPP