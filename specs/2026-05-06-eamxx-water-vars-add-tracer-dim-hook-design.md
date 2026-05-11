# WaterTracerHook Interface Design

## Design Decision: Process-Specific Function Pointer Table

Based on user feedback, the hook interface uses **process-specific function pointers** rather than a single generic hook. This provides:

1. **Clarity at call sites** - obvious what physics is happening
2. **Right data for each process** - tailored signatures per pathway
3. **Different fractionation physics** - each pathway has distinct fractionation factors
4. **Zero-cost in default builds** - inline no-ops eliminated by optimizer

## Hook Table Structure

```cpp
// water_tracers/water_tracer_hooks.hpp

namespace scream {
namespace WaterTracerHooks {

// Process-specific hook signatures
// Each hook receives the relevant water fields (rank-3 with CMP dim),
// the phase-change tendency, and thermodynamic state (T, p)

// Condensation: vapor → liquid (qv decreases, qc increases)
using CondensationHookFn = void(*)(
  const view_3d& qv,        // (COL, CMP, LEV) - vapor
  const view_3d& qc,        // (COL, CMP, LEV) - cloud liquid
  const view_2d& T,         // (COL, LEV) - temperature
  const view_2d& p,         // (COL, LEV) - pressure
  const view_2d& dqv_tend,  // (COL, LEV) - bulk qv tendency (negative)
  int ncol, int nlev
);

// Evaporation: liquid → vapor (qc decreases, qv increases)
using EvaporationHookFn = void(*)(
  const view_3d& qv,
  const view_3d& qc,
  const view_2d& T,
  const view_2d& p,
  const view_2d& dqc_tend,  // bulk qc tendency (negative)
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

// Function pointer table (initialized to no-ops in all builds)
inline CondensationHookFn  condensation_hook  = default_noop_condensation;
inline EvaporationHookFn   evaporation_hook   = default_noop_evaporation;
inline DepositionHookFn    deposition_hook    = default_noop_deposition;
inline SublimationHookFn   sublimation_hook   = default_noop_sublimation;
inline FreezingHookFn      freezing_hook      = default_noop_freezing;
inline MeltingHookFn       melting_hook       = default_noop_melting;

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

// Hook initialization (called only when SCREAM_TRACE_WATER=ON)
#ifdef SCREAM_TRACE_WATER
void initialize_water_tracer_hooks();
#endif

} // namespace WaterTracerHooks
} // namespace scream
```

## Call Site Examples

### In P3 microphysics (condensation)

```cpp
// After computing condensation tendency
auto dqv_cond = compute_condensation(...);
qv -= dqv_cond;  // bulk water updated
qc += dqv_cond;

// Hook call (unconditional, inlines to nothing in default builds)
WaterTracerHooks::condensation_hook(qv, qc, T, p, dqv_cond, ncol, nlev);
```

### In P3 microphysics (evaporation)

```cpp
// After computing evaporation tendency  
auto dqc_evap = compute_evaporation(...);
qc -= dqc_evap;  // bulk water updated
qv += dqc_evap;

// Hook call
WaterTracerHooks::evaporation_hook(qv, qc, T, p, dqc_evap, ncol, nlev);
```

### In SHOC (turbulent mixing with condensation)

```cpp
// After turbulent mixing produces condensation
auto dqv_turb_cond = ...;
qv -= dqv_turb_cond;
qc += dqv_turb_cond;

// Hook call
WaterTracerHooks::condensation_hook(qv, qc, T, p, dqv_turb_cond, ncol, nlev);
```

## Implementation Under SCREAM_TRACE_WATER=ON

```cpp
// water_tracers/eamxx_water_tracers.cpp (compiled only when flag is ON)

#ifdef SCREAM_TRACE_WATER

#include "water_tracer_hooks.hpp"
#include "water_isotopes.hpp"

namespace scream {
namespace WaterTracerHooks {

// Real isotope-aware implementation for N=1 (still no-op, no isotopes)
void isotope_condensation_n1(const view_3d& qv, const view_3d& qc,
                              const view_2d& T, const view_2d& p,
                              const view_2d& dqv_tend, int ncol, int nlev) {
  // With WTRC_MAX_CNST=1, CMP=0 is bulk only, nothing to do
  return;
}

// Real isotope-aware implementation for N>=2
void isotope_condensation(const view_3d& qv, const view_3d& qc,
                          const view_2d& T, const view_2d& p,
                          const view_2d& dqv_tend, int ncol, int nlev) {
  // For each tracer index i > 0:
  //   auto cat_idx = tracer_catalog_index(i);  // from registry (future spec)
  //   auto alpha = AlphaEqLiquidVapor(cat_idx, T);
  //   Apply fractionation: qv(i) and qc(i) evolve with alpha factor
  // (Full implementation in fractionation-physics spec)
}

void initialize_water_tracer_hooks() {
  if (WTRC_MAX_CNST == 1) {
    // N=1: register no-op hooks (BFB with default build)
    condensation_hook  = isotope_condensation_n1;
    evaporation_hook   = ..._n1;
    // etc
  } else {
    // N>=2: register real isotope-aware hooks
    condensation_hook  = isotope_condensation;
    evaporation_hook   = isotope_evaporation;
    deposition_hook    = isotope_deposition;
    sublimation_hook   = isotope_sublimation;
    freezing_hook      = isotope_freezing;
    melting_hook       = isotope_melting;
  }
}

} // namespace WaterTracerHooks
} // namespace scream

#endif // SCREAM_TRACE_WATER
```

## Notes for This Spec

For the **qv-only spec** (this one), the hooks are:

1. **Defined** in water_tracer_hooks.hpp (table + no-op defaults)
2. **Called unconditionally** from P3/SHOC at identified phase-change sites
3. **Implemented** as no-ops even under SCREAM_TRACE_WATER=ON when N=1 (BFB check)
4. **Implemented** as passive-copy for N=2 test (tracer[1] follows tracer[0])

Full fractionation physics (using AlphaEq* from water_isotopes.hpp) is deferred to follow-up spec.

## Additional Hooks Needed (Future)

As qc, qi, qr, qm get lifted to rank-3:
- Autoconversion (qc → qr)
- Accretion (qc + qr → qr)
- Bergeron process (qc → qi)
- Ice shedding (qi → qr)
- Riming (qc + qi → qi)

Each will follow the same pattern: process-specific signature, inline no-op default.
