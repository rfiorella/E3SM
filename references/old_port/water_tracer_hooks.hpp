#ifndef WATER_TRACER_HOOKS_HPP
#define WATER_TRACER_HOOKS_HPP

#include "share/field/field.hpp"
#include <Kokkos_Core.hpp>

namespace scream {
namespace WaterTracerHooks {

// Type aliases for Kokkos views
using view_2d = Field::view_host_t<Real**>;
using view_3d = Field::view_host_t<Real***>;

// Process-specific hook function signatures
// Each hook receives the relevant water fields (rank-3 with CMP dim),
// the phase-change tendency, and thermodynamic state (T, p)

// Condensation: vapor → liquid (qv decreases, qc increases)
using CondensationHookFn = void(*)(
  const view_3d& qv,        // (COL, CMP, LEV) - vapor
  const view_3d& qc,        // (COL, CMP, LEV) - cloud liquid
  const view_2d& T,         // (COL, LEV) - temperature
  const view_2d& p,         // (COL, LEV) - pressure
  const view_2d& dqv_tend,  // (COL, LEV) - bulk qv tendency (negative for condensation)
  int ncol, int nlev
);

// Evaporation: liquid → vapor (qc decreases, qv increases)
using EvaporationHookFn = void(*)(
  const view_3d& qv,
  const view_3d& qc,
  const view_2d& T,
  const view_2d& p,
  const view_2d& dqc_tend,  // bulk qc tendency (negative for evaporation)
  int ncol, int nlev
);

// Deposition: vapor → ice (qv decreases, qi increases)
using DepositionHookFn = void(*)(
  const view_3d& qv,
  const view_3d& qi,
  const view_2d& T,
  const view_2d& p,
  const view_2d& dqv_tend,
  int ncol, int nlev
);

// Sublimation: ice → vapor (qi decreases, qv increases)
using SublimationHookFn = void(*)(
  const view_3d& qi,
  const view_3d& qv,
  const view_2d& T,
  const view_2d& p,
  const view_2d& dqi_tend,
  int ncol, int nlev
);

// Freezing: liquid → ice (qc decreases, qi increases)
using FreezingHookFn = void(*)(
  const view_3d& qc,
  const view_3d& qi,
  const view_2d& T,
  const view_2d& p,
  const view_2d& dqc_tend,
  int ncol, int nlev
);

// Melting: ice → liquid (qi decreases, qc increases)
using MeltingHookFn = void(*)(
  const view_3d& qi,
  const view_3d& qc,
  const view_2d& T,
  const view_2d& p,
  const view_2d& dqi_tend,
  int ncol, int nlev
);

// Default no-op implementations (inline, eliminated by optimizer)
// These are compiled in ALL builds and have zero runtime cost
inline void default_noop_condensation(const view_3d&, const view_3d&,
                                      const view_2d&, const view_2d&,
                                      const view_2d&, int, int) {}

inline void default_noop_evaporation(const view_3d&, const view_3d&,
                                     const view_2d&, const view_2d&,
                                     const view_2d&, int, int) {}

inline void default_noop_deposition(const view_3d&, const view_3d&,
                                    const view_2d&, const view_2d&,
                                    const view_2d&, int, int) {}

inline void default_noop_sublimation(const view_3d&, const view_3d&,
                                     const view_2d&, const view_2d&,
                                     const view_2d&, int, int) {}

inline void default_noop_freezing(const view_3d&, const view_3d&,
                                  const view_2d&, const view_2d&,
                                  const view_2d&, int, int) {}

inline void default_noop_melting(const view_3d&, const view_3d&,
                                 const view_2d&, const view_2d&,
                                 const view_2d&, int, int) {}

// Function pointer table (initialized to no-ops in all builds)
inline CondensationHookFn  condensation_hook  = default_noop_condensation;
inline EvaporationHookFn   evaporation_hook   = default_noop_evaporation;
inline DepositionHookFn    deposition_hook    = default_noop_deposition;
inline SublimationHookFn   sublimation_hook   = default_noop_sublimation;
inline FreezingHookFn      freezing_hook      = default_noop_freezing;
inline MeltingHookFn       melting_hook       = default_noop_melting;

// Hook initialization (called only when SCREAM_TRACE_WATER=ON)
#ifdef SCREAM_TRACE_WATER
void initialize_water_tracer_hooks();
#endif

} // namespace WaterTracerHooks
} // namespace scream

#endif // WATER_TRACER_HOOKS_HPP
