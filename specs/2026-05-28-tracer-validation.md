---
spec_id: 2026-05-28-tracer-validation
spec_type: model-e3sm
spec_version: 1
title: "Tracer ratio utility and validation test"
created: 2026-05-28T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

depends_on:
  - specs/2026-05-28-extend-precip-tracer.md

inputs:
  source_files:
    - wiso_group1_campaign_revision.md
    - specs/2026-05-28-extend-qv-tracer.md
    - specs/2026-05-28-extend-cloud-tracer.md
    - specs/2026-05-28-extend-precip-tracer.md
  data: []
  baseline: null

deliverables:
  - components/eamxx/src/physics/water_tracers/water_tracer_ratio.hpp
  - components/eamxx/src/physics/water_tracers/water_tracer_ratio.cpp
  - tests/water_tracers/test_tracer_scaling.cpp
  - tests/water_tracers/validate_tracer_scaling.py
  - tests/water_tracers/run_tracer_scaling_test.sh
  - tests/water_tracers/CMakeLists.txt
  - docs/wiso/tracer_data_model.md

success_criteria:
  - id: compile-ratio-utility
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr5 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=2 && cmake --build build/pr5 --target water_tracer_ratio -j"
    expect: exit_zero
    phase: implementation
    verifies:
      - deliverable: components/eamxx/src/physics/water_tracers/water_tracer_ratio.hpp
      - deliverable: components/eamxx/src/physics/water_tracers/water_tracer_ratio.cpp
      - claim: "Water tracer ratio utility compiles without errors"

  - id: unit-test-ratio-function
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/pr5 -R test_tracer_ratio --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - claim: "compute_tracer_ratio handles division by zero, produces correct ratios"

  - id: compile-test-case-n2
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr5-test -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=2 && cmake --build build/pr5-test -j"
    expect: exit_zero
    phase: implementation
    verifies:
      - deliverable: tests/water_tracers/test_tracer_scaling.cpp

  - id: test-case-runs
    type: shell
    cmd: "cd tests/water_tracers && bash run_tracer_scaling_test.sh"
    expect: exit_zero
    phase: testing
    verifies:
      - claim: "Doubly-periodic test with SCREAM_NUM_TRACERS=2 runs to completion"

  - id: tracer-ratio-validation
    type: comparison
    cmd: "cd tests/water_tracers && python3 validate_tracer_scaling.py output.nc 1e-12 1e-15"
    expect: all_pass
    rtol: 1e-12
    atol: 1e-15
    gate: blocking
    phase: testing
    verifies:
      - claim: "All water reservoirs maintain tracer ratio of 0.5 to within numerical precision"
      - claim: "Maximum relative error < 1e-12 for double precision"
    on_fail: halt_after_investigation
    resolution_notes: "If validation fails: ratio errors >1e-12 indicate bug in tracer handling. Bisect to find which process breaks ratio preservation. Common issues: missing tracer loops, hardcoded slot-0, incorrect subview. See planning doc section 'Action on tracer ratio failure' for debugging steps."

  - id: performance-baseline-n1
    type: shell
    cmd: "cd components/eamxx && ./scripts/benchmark-eamxx.sh --tracers 1 > bench_n1.txt"
    expect: exit_zero
    phase: testing
    verifies:
      - claim: "Baseline performance benchmark with SCREAM_NUM_TRACERS=1 completes"

  - id: performance-scream-num-tracers-2
    type: shell
    cmd: "cd components/eamxx && ./scripts/benchmark-eamxx.sh --tracers 2 > bench_n2.txt && python3 scripts/compare_benchmarks.py bench_n1.txt bench_n2.txt --max-overhead 0.05"
    expect: exit_zero
    gate: advisory
    phase: testing
    depends_on: [performance-baseline-n1]
    verifies:
      - claim: "SCREAM_NUM_TRACERS=2 overhead < 5% vs SCREAM_NUM_TRACERS=1"

  - id: documentation-updated
    type: shell
    cmd: "grep -q 'Tracer ratio computation' docs/wiso/tracer_data_model.md && grep -q 'Validation methodology' docs/wiso/tracer_data_model.md"
    expect: exit_zero
    phase: implementation
    verifies:
      - deliverable: docs/wiso/tracer_data_model.md

out_of_scope:
  - Isotope fractionation - future campaigns
  - Multiple tracers beyond test_tracer - future campaigns
  - Tagged tracers - future campaigns
  - Production surface flux infrastructure - basic test harness only

resolved_decisions:
  test_design: "Scale surface evaporation flux for tracer 1 by exactly 0.5. Initialize all tracer 1 fields to 0.5 * tracer 0. Run until precipitation develops. Verify ratios preserved to rtol=1e-12."
  tolerance: "rtol=1e-12 for double precision, rtol=1e-6 for single precision builds. This is numerical precision - any larger error indicates bug."
  test_harness: "Use doubly-periodic configuration (F2010-SCREAMv1-DP, ne4pg2). Simpler than full-grid, sufficient for passive transport validation."

ask_before:
  - Relaxing validation tolerance beyond 1e-12
  - Skipping validation test
  - Modifying test to use different scaling factor

model_specific:
  validation_tier: 1
  compset: F2010-SCREAMv1-DP
  resolution: ne4pg2
  baseline_tag: null
  dp_proxy: true
---

# Tracer ratio utility and validation test

## Purpose

Add `compute_tracer_ratio` utility and validation test that proves passive tracer
transport correctness. This is the critical test that Group 1 is correct.

Test scales surface evaporation flux by 0.5 for tracer 1, verifies precipitation
ratios converge to 0.5 within numerical precision (rtol=1e-12). If ratios drift,
there's a bug in tracer handling.

## Approach

Reference skills: `eamxx-cpp-conventions`, `tracer-conservation`, `e3sm-build-and-test`.

1. Implement `compute_tracer_ratio` utility (Kokkos parallel, handles division by zero)
2. Create doubly-periodic test case with `SCREAM_NUM_TRACERS=2`
3. Initialize tracer 1 to exactly 0.5 * tracer 0 for all reservoirs
4. Scale surface evaporation flux for tracer 1 by exactly 0.5
5. Run until precipitation develops (sufficient timesteps for quasi-steady state)
6. Python validation script computes ratios, verifies within rtol=1e-12
7. Update documentation with tracer ratio methodology

## Tracer ratio utility

```cpp
namespace scream {
namespace water_tracers {

// Compute ratio of tracer N to tracer 0 (bulk)
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT compute_tracer_ratio(
  const ScalarT& tracer_amount,
  const ScalarT& bulk_amount,
  const ScalarT& min_bulk = 1e-20
) {
  if (bulk_amount < min_bulk) {
    return 0.0;
  }
  return tracer_amount / bulk_amount;
}

// Kokkos-parallel field version
template<typename ViewT>
void compute_tracer_ratio_field(
  const ViewT& tracer_field,  // (tracer, col, lev)
  int tracer_idx,
  ViewT& ratio_field          // (col, lev)
);

} // namespace water_tracers
} // namespace scream
```

## Test configuration

```cpp
// Initialize tracer 1 = 0.5 * tracer 0
Kokkos::parallel_for(ncols, KOKKOS_LAMBDA(int icol) {
  for (int ilev = 0; ilev < nlevs; ++ilev) {
    qv(1, icol, ilev) = 0.5 * qv(0, icol, ilev);
    qc(1, icol, ilev) = 0.5 * qc(0, icol, ilev);
    // ... same for qi, qr, qs
  }
});

// Scale surface flux
Kokkos::parallel_for(ncols, KOKKOS_LAMBDA(int icol) {
  qv_surf_flux(1, icol) = 0.5 * qv_surf_flux(0, icol);
});
```

## Validation script

```python
#!/usr/bin/env python3
import netCDF4 as nc
import numpy as np

def validate_tracer_scaling(output_file, rtol=1e-12, atol=1e-15):
    reservoirs = ['qv', 'qc', 'qi', 'qr', 'qs', 'rain', 'snow']
    expected_ratio = 0.5
    
    for var in reservoirs:
        bulk = data[-1, 0, ...]   # Final timestep, slot 0
        tracer = data[-1, 1, ...]  # Final timestep, slot 1
        
        mask = bulk > 1e-20
        ratio = tracer[mask] / bulk[mask]
        
        max_error = np.max(np.abs(ratio - expected_ratio))
        assert max_error <= atol + rtol * expected_ratio, \
            f"{var} failed: max error {max_error:.3e}"
```

## Test execution

```bash
cd tests/water_tracers
./run_tracer_scaling_test.sh

# Expected output:
# qv: mean ratio = 0.500000000000, max rel error = 4.4e-16, PASS
# qc: mean ratio = 0.500000000000, max rel error = 2.2e-16, PASS
# ...
# All reservoirs passed tracer ratio validation
```

## Failure diagnosis

If validation fails (ratio error > 1e-12), indicates bug:

1. Check initialization: tracer 1 = 0.5 * tracer 0
2. Check surface flux scaling applied correctly
3. Check all physics processes apply same rates to slots 0 and 1
4. Look for inadvertent resets of slot 1
5. Check sedimentation handles tracer dimension
6. Verify no hardcoded slot-0 access without tracer loop

Common bugs:
- Missing `qv(0, icol, ilev)` → using `qv(icol, ilev)` (no slot index)
- Loop over tracer but hardcoded index 0 inside
- Conditional that only processes slot 0
- Field copy not preserving tracer dimension

## Critical requirement

This test must pass for Group 1 to be complete. Ratios preserved to machine
precision proves passive tracer transport is correct - no spurious numerical
diffusion, no tracer corruption.
