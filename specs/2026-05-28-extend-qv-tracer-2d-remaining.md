---
spec_id: 2026-05-28-extend-qv-tracer-2d-remaining
spec_type: model-e3sm
spec_version: 1
title: "Extend qv to tracer dimension - Part 2d: Remaining Processes"
created: 2026-05-29T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

dependencies:
  - 2026-05-28-water-tracer-metadata-and-gate
  - 2026-05-28-extend-qv-tracer-2a-infrastructure
  - 2026-05-28-extend-qv-tracer-2b-shoc
  - 2026-05-28-extend-qv-tracer-2c-p3

inputs:
  source_files:
    - specs/2026-05-28-extend-qv-tracer-2a-infrastructure.md
    - specs/2026-05-28-extend-qv-tracer-2b-shoc.md
    - specs/2026-05-28-extend-qv-tracer-2c-p3.md
    - components/eamxx/src/physics/rrtmgp/eamxx_rrtmgp_process_interface.cpp
    - components/eamxx/src/physics/zm/eamxx_zm_process_interface.cpp
    - components/eamxx/src/dynamics/homme/eamxx_homme_process_interface.cpp
    - components/eamxx/src/physics/surface_coupling/eamxx_surface_coupling_exporter.cpp
  data: []
  baseline: /data/baselines/pre-group1

deliverables:
  - components/eamxx/src/physics/rrtmgp/eamxx_rrtmgp_process_interface.cpp
  - components/eamxx/src/physics/rrtmgp/impl/rrtmgp_utils.hpp
  - components/eamxx/src/physics/zm/eamxx_zm_process_interface.cpp
  - components/eamxx/src/physics/zm/impl/zm_conv_evap_impl.hpp
  - components/eamxx/src/dynamics/homme/eamxx_homme_process_interface.cpp
  - components/eamxx/src/dynamics/homme/interface/scream_homme_interface.cpp
  - components/eamxx/src/physics/surface_coupling/eamxx_surface_coupling_exporter.cpp
  - tests/water_tracers/test_qv_rrtmgp.cpp
  - tests/water_tracers/test_qv_zm.cpp
  - tests/water_tracers/test_qv_homme.cpp
  - tests/water_tracers/test_qv_surface.cpp
  - tests/water_tracers/CMakeLists.txt

success_criteria:
  - id: compile-all-processes-n1
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr2d-n1 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=1 -DSCREAM_TRACER_ACCESS=SUBVIEW && cmake --build build/pr2d-n1 -j"
    expect: exit_zero
    phase: implementation
    verifies:
      - claim: "All remaining processes compile with tracer-aware qv"

  - id: unit-tests-all-processes
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/pr2d-n1 -R 'test_qv_(rrtmgp|zm|homme|surface)' --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - claim: "All remaining processes correctly access qv slot-0"

  - id: rrtmgp-standalone-bfb
    type: comparison
    cmd: "cd components/eamxx && ctest --test-dir build/pr2d-n1 -R rrtmgp_stand_alone --output-on-failure"
    expect: all_bfb
    rtol: 0
    atol: 0
    gate: blocking
    phase: testing
    verifies:
      - claim: "RRTMGP standalone tests pass BFB"
    on_fail: "Halt. Investigate RRTMGP-specific issue."

  - id: full-atmosphere-bfb
    type: comparison
    cmd: "cd components/eamxx && ./scripts/test-all-eamxx -m docker-local --baseline-dir /data/baselines/pre-group1 -t sp -c"
    expect: all_bfb
    rtol: 0
    atol: 0
    gate: blocking
    phase: testing
    verifies:
      - claim: "Full atmosphere configuration (all processes) passes BFB with SCREAM_NUM_TRACERS=1"
    on_fail: "Halt. Bisect to find which process or process interaction broke BFB."
    resolution_notes: "If full-atmosphere fails but all standalone tests passed, issue is in process coupling or driver. Check field manager, boundary exchange, restart I/O."

  - id: compile-all-processes-n3
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr2d-n3 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=3 && cmake --build build/pr2d-n3 -j"
    expect: exit_zero
    phase: testing
    verifies:
      - claim: "All processes compile with multiple tracers"

  - id: performance-overhead
    type: shell
    cmd: "cd components/eamxx && ./scripts/benchmark-eamxx.sh --tracers 1 > bench_pr2d.txt && python3 scripts/compare_benchmarks.py /data/baselines/pre-group1/benchmark.txt bench_pr2d.txt --max-overhead 0.02"
    expect: exit_zero
    gate: advisory
    phase: testing
    verifies:
      - claim: "Full atmosphere performance overhead < 2% vs pre-campaign baseline"

model_specific:
  validation_tier: 2
  compset: F2010-SCREAMv1
  resolution: ne4pg2_ne4pg2
  machine: docker-local
  baseline_tag: pre-group1
  dp_proxy: false

---

# Purpose

Extend the remaining EAMxx processes to use tracer-aware qv field:
- **RRTMGP**: Radiation (qv needed for optical properties)
- **ZM**: Deep convection (qv for convective evaporation)
- **HOMME**: Dynamics (qv advection)
- **Surface Coupling**: Export qv to surface models

This is the final spec for qv tracer extension. After this, all EAMxx processes handle tracer-aware qv.

# Approach

Follow the pattern from specs 2b (SHOC) and 2c (P3). Each process:
1. Update field registration
2. Apply accessor pattern
3. Update kernel calls
4. Unit test + BFB validation

## Process-Specific Notes

### RRTMGP (Radiation)
- Uses qv for calculating water vapor optical properties
- Relatively simple qv usage (diagnostic, not prognostic)
- Small scope: 1-2 kernel files

### ZM (Deep Convection)
- Uses qv in convective evaporation and detrainment
- Moderate scope: 2-3 kernel files
- Coupled with P3 (convective precipitation)

### HOMME (Dynamics)
- Advects qv in spectral element grid
- Interface between SCREAM and HOMME dynamical core
- Moderate scope: 1-2 interface files

### Surface Coupling
- Exports qv to coupler/surface models
- Simple scope: 1 file
- Important for boundary conditions

# Implementation Pattern

Same as specs 2b and 2c:
1. Field registration with tracer layout
2. Accessor pattern (EXPLICIT or SUBVIEW from spec 2a)
3. Kernel updates per process
4. Process-by-process BFB validation
5. Full-atmosphere integration test

The key difference: This spec combines multiple processes because each is smaller scope than SHOC or P3.

# Resolved Decisions

## Grouping strategy
Group these four processes into one spec because:
1. Each has smaller scope than SHOC or P3
2. Pattern is well-established by specs 2b and 2c
3. Risk is low - if BFB fails, can bisect to specific process
4. Faster campaign completion

## Accessor pattern
Use `{{PATTERN}}` from spec 2a, same as specs 2b and 2c.

## Testing strategy
- Unit tests per process (smoke tests for field access)
- Standalone tests per process (where available)
- Full-atmosphere integration test (critical gate)

The full-atmosphere test is the most important - it validates all processes work together.

# Ask Before

1. If any process standalone test fails BFB (suggests process-specific issue)
2. If full-atmosphere test fails BFB but all standalone tests passed (suggests coupling issue)
3. If implementation reveals one of these processes is more complex than estimated

# Out of Scope

- Other water species (qc, qi, qr, qs, qm, rain, snow) - specs 3 and 4
- Multi-tracer validation - spec 5
- Field manager changes - handled in spec 2a
- Restart/history I/O changes - TBD if needed

# Notes

## Why Group These Processes?

After specs 2b (SHOC) and 2c (P3), the pattern is proven. These four processes are:
- Smaller scope individually
- Lower BFB risk (simpler qv usage)
- Natural to group by function: radiation, convection, dynamics, coupling

Grouping speeds up the campaign while maintaining safety via BFB gates.

## Full-Atmosphere Test as Critical Gate

The `full-atmosphere-bfb` criterion is the most important in this spec because:
1. It's the first time all processes run together with tracer-aware qv
2. Tests process coupling and field manager integration
3. Tests boundary exchange and halo updates
4. Tests restart I/O (if restart tests included)

If this fails, the investigation may reveal issues not visible in standalone tests.

## Process Update Order

Within this spec, update processes in this order:
1. RRTMGP (simplest, diagnostic use)
2. ZM (moderate, coupled with P3)
3. HOMME (moderate, dynamics interface)
4. Surface coupling (simplest, just export)

This order minimizes risk - start with simpler processes.

## Deliverables Uncertainty

The deliverables list is an estimate. Actual kernel file count may differ. This is acceptable - the pattern applies uniformly.

## Success After This Spec

After spec 2d passes:
- All EAMxx processes handle tracer-aware qv
- Full atmosphere runs BFB with SCREAM_NUM_TRACERS=1
- Pattern is validated across all process types
- Ready to extend other water species (specs 3-4) using same pattern
