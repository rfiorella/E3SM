#include "water_tracers.hpp"
// #include "water_isotopes.hpp"  // Not needed for passive-copy hooks
#include "water_tracer_hooks.hpp"

#include <vector>
#include <string>

// Note: These includes are not currently needed for the passive-copy implementation
// #include "share/core/eamxx_types.hpp"
// #include "share/physics/physics_constants.hpp"
// #include <ekat/util/ekat_string_utils.hpp>
// #include <ekat/ekat_scalar_traits.hpp>
// #include <ekat/logging/ekat_logger.hpp>

namespace scream { 
namespace WaterTracers { 
// =========================================================================================
// WaterTracers::WaterTracers (const ekat::Comm& comm, const ekat::ParameterList& params)
//   : AtmosphereProcess(comm, params)
// {
//   // Nothing to do here
// }

// void WtypeGetItype(name) {
//   // Retrieve WaterType based on type name
// };


void WtypeGetAlpha() {
  // retrieve fractionation factor for process that converts
  // source water type (isrctype) to desination water type (idsttype)

};

// // Set of functions to determine the type of water tracer
// bool trc_is_wtrc(tracer_num) {
//   // returns true if tracer is a water tracer
// };




void OceanTracerFlux() {
  // compute water tracer exchange from ocean
  // see L. 531-635 of wiso_flxoce in CAM
};

} // namespace WaterTracers

#ifdef SCREAM_TRACE_WATER

// Hook implementations for SCREAM_TRACE_WATER=ON

namespace WaterTracerHooks {

// Passive-copy implementation (for N=2 test and N=1 no-op)
// Copies tendency from bulk tracer (CMP=0) to all other tracers
void passive_copy_condensation(const view_3d& qv, const view_3d& qc,
                                const view_2d& T, const view_2d& p,
                                const view_2d& dqv_tend, int ncol, int nlev) {
  // WTRC_MAX_CNST is a preprocessor define, not a namespace member
  if (WTRC_MAX_CNST == 1) return;  // No work for single tracer

  // For each additional tracer, apply same tendency as bulk
  for (int i = 1; i < WTRC_MAX_CNST; ++i) {
    for (int icol = 0; icol < ncol; ++icol) {
      for (int ilev = 0; ilev < nlev; ++ilev) {
        qv(icol, i, ilev) += dqv_tend(icol, ilev);
        qc(icol, i, ilev) -= dqv_tend(icol, ilev);
      }
    }
  }
}

// Similar passive-copy implementations for other hooks
void passive_copy_evaporation(const view_3d& qv, const view_3d& qc,
                               const view_2d& T, const view_2d& p,
                               const view_2d& dqc_tend, int ncol, int nlev) {
  // WTRC_MAX_CNST is a preprocessor define, not a namespace member
  if (WTRC_MAX_CNST == 1) return;

  for (int i = 1; i < WTRC_MAX_CNST; ++i) {
    for (int icol = 0; icol < ncol; ++icol) {
      for (int ilev = 0; ilev < nlev; ++ilev) {
        qc(icol, i, ilev) += dqc_tend(icol, ilev);
        qv(icol, i, ilev) -= dqc_tend(icol, ilev);
      }
    }
  }
}

void passive_copy_deposition(const view_3d& qv, const view_3d& qi,
                              const view_2d& T, const view_2d& p,
                              const view_2d& dqv_tend, int ncol, int nlev) {
  // WTRC_MAX_CNST is a preprocessor define, not a namespace member
  if (WTRC_MAX_CNST == 1) return;

  for (int i = 1; i < WTRC_MAX_CNST; ++i) {
    for (int icol = 0; icol < ncol; ++icol) {
      for (int ilev = 0; ilev < nlev; ++ilev) {
        qv(icol, i, ilev) += dqv_tend(icol, ilev);
        qi(icol, i, ilev) -= dqv_tend(icol, ilev);
      }
    }
  }
}

void passive_copy_sublimation(const view_3d& qi, const view_3d& qv,
                               const view_2d& T, const view_2d& p,
                               const view_2d& dqi_tend, int ncol, int nlev) {
  // WTRC_MAX_CNST is a preprocessor define, not a namespace member
  if (WTRC_MAX_CNST == 1) return;

  for (int i = 1; i < WTRC_MAX_CNST; ++i) {
    for (int icol = 0; icol < ncol; ++icol) {
      for (int ilev = 0; ilev < nlev; ++ilev) {
        qi(icol, i, ilev) += dqi_tend(icol, ilev);
        qv(icol, i, ilev) -= dqi_tend(icol, ilev);
      }
    }
  }
}

void passive_copy_freezing(const view_3d& qc, const view_3d& qi,
                            const view_2d& T, const view_2d& p,
                            const view_2d& dqc_tend, int ncol, int nlev) {
  // WTRC_MAX_CNST is a preprocessor define, not a namespace member
  if (WTRC_MAX_CNST == 1) return;

  for (int i = 1; i < WTRC_MAX_CNST; ++i) {
    for (int icol = 0; icol < ncol; ++icol) {
      for (int ilev = 0; ilev < nlev; ++ilev) {
        qc(icol, i, ilev) += dqc_tend(icol, ilev);
        qi(icol, i, ilev) -= dqc_tend(icol, ilev);
      }
    }
  }
}

void passive_copy_melting(const view_3d& qi, const view_3d& qc,
                           const view_2d& T, const view_2d& p,
                           const view_2d& dqi_tend, int ncol, int nlev) {
  // WTRC_MAX_CNST is a preprocessor define, not a namespace member
  if (WTRC_MAX_CNST == 1) return;

  for (int i = 1; i < WTRC_MAX_CNST; ++i) {
    for (int icol = 0; icol < ncol; ++icol) {
      for (int ilev = 0; ilev < nlev; ++ilev) {
        qi(icol, i, ilev) += dqi_tend(icol, ilev);
        qc(icol, i, ilev) -= dqi_tend(icol, ilev);
      }
    }
  }
}

// Initialize hooks to passive-copy implementations
void initialize_water_tracer_hooks() {
  condensation_hook  = passive_copy_condensation;
  evaporation_hook   = passive_copy_evaporation;
  deposition_hook    = passive_copy_deposition;
  sublimation_hook   = passive_copy_sublimation;
  freezing_hook      = passive_copy_freezing;
  melting_hook       = passive_copy_melting;
}

} // namespace WaterTracerHooks

#endif // SCREAM_TRACE_WATER

} // namespace scream
