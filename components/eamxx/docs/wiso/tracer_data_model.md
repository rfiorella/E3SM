# EAMxx Water Tracer Data Model

## Overview

This document describes the water tracer data model in EAMxx, including field
layouts, tracer ratio computation, and validation methodology.

## Field Layout

All water mass fields in EAMxx support multiple tracer constituents via a leading
tracer dimension:

- **3D fields**: `(tracer, col, lev)` - qv, qc, qi, qr, qm
- **2D surface fields**: `(tracer, col)` - precip_liq_surf_mass, precip_ice_surf_mass

The number of tracers is set at compile time via `SCREAM_NUM_TRACERS` CMake variable.

### Slot 0: Bulk Water

Slot 0 (tracer index 0) always contains the bulk water mass. This is the canonical
reservoir used for:
- Physical parameterizations (phase changes, microphysics)
- Conservation checking
- Backward compatibility (BFB with pre-tracer code when `SCREAM_NUM_TRACERS=1`)

### Slot 1+: Tagged Tracers

Slots 1 and higher contain tagged tracer constituents:
- Water isotopologs (H₂¹⁶O, H₂¹⁸O, HDO)
- Source-tagged water (ocean, land, ice)
- Age-tagged water
- Custom tracers defined via `add_water_tracer()` CMake function

Tracer slots follow the same physical processes as bulk water (advection, mixing,
phase changes, precipitation) but with optional modifications to surface fluxes
for isotope fractionation or tagging purposes.

## Tracer Ratio Computation

### Definition

The tracer ratio is defined as:

```
ratio(N) = tracer(N) / tracer(0)
```

where:
- `tracer(N)` is the amount of tracer N at a grid point
- `tracer(0)` is the bulk water amount at the same point
- Division by zero is handled by returning 0.0 when bulk < threshold (default: 1e-20)

### Use Cases

1. **Passive tracer validation**: For tracers with scaled surface fluxes and no
   fractionation, ratios should be preserved to machine precision throughout
   model integration. This validates correct tracer transport.

2. **Isotope ratio tracking**: For water isotopologs, the ratio represents the
   isotopic composition (δ¹⁸O, δD). Ratios change due to fractionation during
   phase transitions.

3. **Source attribution**: For tagged tracers, the ratio represents the fraction
   of water from a specific source.

### Implementation

The `compute_tracer_ratio` utility (defined in `water_tracer_ratio.hpp`) provides:

- **Scalar function**: Compute ratio at a single point (inline, device-compatible)
- **Field function**: Compute ratio field in parallel (Kokkos MDRangePolicy)
- **Mean function**: Compute global mean ratio with masking

Example usage:

```cpp
#include "water_tracer_ratio.hpp"

// Compute ratio at a point
Real ratio = compute_tracer_ratio(qv(1, icol, ilev), qv(0, icol, ilev));

// Compute ratio field in parallel
Kokkos::View<Real**> ratio_field("ratio", ncols, nlevs);
compute_tracer_ratio_field(qv, 1, ratio_field);

// Compute global mean
Real mean_ratio = compute_mean_tracer_ratio(qv, 1);
```

## Validation Methodology

### Tracer Scaling Test

The tracer scaling test validates passive tracer transport by verifying ratio
preservation to machine precision.

**Test procedure**:
1. Build EAMxx with `SCREAM_NUM_TRACERS=2`
2. Initialize tracer 1 = 0.5 × tracer 0 for all reservoirs (qv, qc, qi, qr, surface)
3. Scale surface evaporation flux for tracer 1 by exactly 0.5
4. Run model for sufficient timesteps (typically 100 steps at 30-min timestep)
5. Compute tracer ratios for all reservoirs at final timestep
6. Verify: `|ratio - 0.5| < atol + rtol × 0.5` where rtol=1e-12, atol=1e-15

**Expected result**: All water reservoirs maintain ratio of 0.5 to within
numerical precision (relative error < 1e-12 for double precision).

**Interpretation**:
- PASS: Tracer transport is correct, no spurious sources/sinks
- FAIL: Bug in tracer handling (see debugging guide below)

### Running the Test

```bash
cd components/eamxx/tests/water_tracers
./run_tracer_scaling_test.sh
```

Output:

```
========================================================================
EAMxx Water Tracer Ratio Validation
========================================================================
qv                            : mean ratio = 0.500000000000000, max rel error = 4.4e-16, PASS
qc                            : mean ratio = 0.500000000000000, max rel error = 2.2e-16, PASS
qi                            : mean ratio = 0.500000000000000, max rel error = 3.3e-16, PASS
qr                            : mean ratio = 0.500000000000000, max rel error = 5.5e-16, PASS
precip_liq_surf_mass          : mean ratio = 0.500000000000000, max rel error = 1.1e-16, PASS
precip_ice_surf_mass          : mean ratio = 0.500000000000000, max rel error = 2.2e-16, PASS
========================================================================
RESULT: All reservoirs passed tracer ratio validation
========================================================================
```

### Debugging Ratio Failures

If validation fails (ratio error > 1e-12), investigate:

1. **Initialization**: Verify tracer 1 initialized to exactly 0.5 × tracer 0
2. **Surface fluxes**: Confirm tracer 1 fluxes scaled by 0.5
3. **Process tracer loops**: Check all physics processes loop over tracer dimension
4. **Hardcoded slot-0**: Look for `qv(0, icol, ilev)` without surrounding tracer loop
5. **Field operations**: Verify copies, subviews, reshapes preserve tracer dimension
6. **Conditional logic**: Check for conditions that skip tracer slots
7. **Sedimentation**: Verify vertical transport handles tracer dimension

Common bug patterns:

```cpp
// BAD: Hardcoded slot 0, should loop over tracers
qv_new(0, icol, ilev) = qv_old(0, icol, ilev) + delta;

// BAD: Loop over tracers but use slot 0 inside
for (int itracer = 0; itracer < ntracers; ++itracer) {
  qv_new(itracer, icol, ilev) = qv_old(0, icol, ilev) + delta;  // Wrong!
}

// GOOD: Loop with correct indexing
for (int itracer = 0; itracer < ntracers; ++itracer) {
  qv_new(itracer, icol, ilev) = qv_old(itracer, icol, ilev) + delta;
}

// BAD: Conditional that excludes tracers
if (itracer == 0) {  // Processes only slot 0!
  // ... physics ...
}

// GOOD: Process all tracers
for (int itracer = 0; itracer < ntracers; ++itracer) {
  // ... physics ...
}
```

### Tolerance Justification

- **rtol=1e-12**: Near machine precision for double precision (ε ≈ 2.2e-16)
- **atol=1e-15**: Handles cases where expected ratio is exactly 0.0
- Any error larger than this indicates incorrect tracer handling, not numerical roundoff

For single precision builds, use rtol=1e-6, atol=1e-8.

## BFB Requirement

When `SCREAM_NUM_TRACERS=1`, the model must be bit-for-bit identical to pre-tracer
baseline. This ensures:
- Zero performance overhead for non-tracer configurations
- Backward compatibility with existing runs
- No unintended changes to bulk water physics

## Performance Considerations

The tracer dimension adds memory overhead and computational cost:

- **Memory**: Linear in number of tracers (2× for 2 tracers, 3× for 3, etc.)
- **Compute**: Overhead depends on loop structure:
  - Explicit loops over tracers: proportional to tracer count
  - Kokkos parallel over tracers: same but better vectorization
  - Subviews with index 0: minimal overhead for slot-0 access

**Performance targets** (Group 1 completion):
- SCREAM_NUM_TRACERS=1: < 2% overhead vs. pre-tracer baseline
- SCREAM_NUM_TRACERS=2: < 5% overhead vs. SCREAM_NUM_TRACERS=1

## References

- Campaign: `campaigns/wiso-group1-revised.campaign.md`
- Spec: `specs/2026-05-28-tracer-validation.md`
- Utility: `src/physics/water_tracers/water_tracer_ratio.hpp`
- Test: `tests/water_tracers/test_tracer_scaling.cpp`
- Validation: `tests/water_tracers/validate_tracer_scaling.py`
