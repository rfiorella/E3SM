---
spec_id: 2026-05-28-extend-qv-tracer-2c-p3
spec_type: model-e3sm
spec_version: 1
title: "Extend qv to tracer dimension - Part 2c: P3 Microphysics"
created: 2026-05-29T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

dependencies:
  - 2026-05-28-water-tracer-metadata-and-gate
  - 2026-05-28-extend-qv-tracer-2a-infrastructure
  - 2026-05-28-extend-qv-tracer-2b-shoc

inputs:
  source_files:
    - specs/2026-05-28-extend-qv-tracer-2a-infrastructure.md
    - specs/2026-05-28-extend-qv-tracer-2b-shoc.md
    - components/eamxx/src/physics/p3/eamxx_p3_process_interface.hpp
    - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
    - components/eamxx/src/physics/p3/impl/p3_main_impl.hpp
  data: []
  baseline: /data/baselines/pre-group1

deliverables:
  - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
  - components/eamxx/src/physics/p3/impl/p3_main_impl.hpp
  - components/eamxx/src/physics/p3/impl/p3_evaporation_sublimation_impl.hpp
  - components/eamxx/src/physics/p3/impl/p3_cloud_water_conservation_impl.hpp
  - tests/water_tracers/test_qv_p3.cpp
  - tests/water_tracers/CMakeLists.txt

success_criteria:
  - id: compile-p3-n1
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr2c-n1 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=1 -DSCREAM_TRACER_ACCESS={{PATTERN}} && cmake --build build/pr2c-n1 --target p3 -j"
    expect: exit_zero
    phase: implementation
    verifies:
      - deliverable: components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
      - claim: "P3 process compiles with tracer-aware qv field"

  - id: unit-test-p3-qv-access
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/pr2c-n1 -R test_qv_p3 --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: tests/water_tracers/test_qv_p3.cpp
      - claim: "P3 kernels correctly access qv slot-0"

  - id: p3-standalone-bfb
    type: comparison
    cmd: "cd components/eamxx && ctest --test-dir build/pr2c-n1 -R p3_stand_alone --output-on-failure"
    expect: all_bfb
    rtol: 0
    atol: 0
    gate: blocking
    phase: testing
    verifies:
      - claim: "P3 standalone tests pass BFB vs pre-campaign baseline"
    on_fail: "Halt. Investigate which P3 kernel broke BFB."
    resolution_notes: "P3 has many kernels. Bisect to find which phase process (vapor→cloud, cloud→rain, etc.) broke BFB."

  - id: p3-integration-bfb
    type: comparison
    cmd: "cd components/eamxx && ctest --test-dir build/pr2c-n1 -R dp_.*p3 --output-on-failure"
    expect: all_bfb
    rtol: 0
    atol: 0
    gate: blocking
    phase: testing
    verifies:
      - claim: "P3 integration tests pass BFB"
    on_fail: "Halt. Check field registration and coupling with SHOC."

  - id: compile-p3-n3
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr2c-n3 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=3 && cmake --build build/pr2c-n3 --target p3 -j"
    expect: exit_zero
    phase: testing
    verifies:
      - claim: "P3 compiles with multiple tracers"

model_specific:
  validation_tier: 2
  compset: F2010-SCREAMv1
  resolution: ne4pg2_ne4pg2
  machine: docker-local
  baseline_tag: pre-group1
  dp_proxy: true

---

# Purpose

Extend P3 (Predicted Particle Properties microphysics) to use tracer-aware qv field. P3 is the most complex process in this campaign because it handles all phase transitions involving water vapor (evaporation, sublimation, condensation).

P3 is second after SHOC because:
1. SHOC establishes the pattern
2. P3 is the heaviest qv consumer (all microphysics processes)
3. P3's BFB validation is critical - any microphysics error compounds quickly
4. Success here proves the pattern works for complex multi-kernel processes

# Approach

Follow the pattern established in spec 2b (SHOC), using the accessor pattern from spec 2a.

## Step 1: Update Field Registration

In `eamxx_p3_process_interface.cpp`:
```cpp
auto qv_layout = grid->get_3d_tracer_layout(SCREAM_NUM_TRACERS);
add_field<Required>("qv", qv_layout, kg/kg, grid);
```

## Step 2: Update Kernel Calls

P3 has more kernel call sites than SHOC, but the pattern is identical:
- Update process interface to pass tracer-aware qv
- Apply accessor pattern (EXPLICIT or SUBVIEW from spec 2a)
- Verify each kernel implementation handles qv correctly

## Step 3: Update P3 Kernels

Key P3 kernels touching qv:
- `p3_main_impl.hpp` - main driver
- `p3_evaporation_sublimation_impl.hpp` - vapor production from rain/ice
- `p3_cloud_water_conservation_impl.hpp` - vapor<->cloud phase change
- Additional impl files as discovered during implementation

Each follows the same pattern as SHOC (spec 2b).

## Step 4: Testing

Same strategy as SHOC:
- Unit test for qv field access
- Standalone P3 tests with BFB requirement
- Integration tests (doubly-periodic with P3+SHOC)

# Implementation Pattern

Identical to spec 2b (SHOC):
1. Field registration with tracer layout
2. Accessor pattern application
3. Kernel updates
4. BFB validation via existing test suite

The only difference is scope: P3 has more kernels than SHOC.

# Resolved Decisions

## Follow SHOC pattern
Use exact same approach as spec 2b. Any deviation should be justified and documented.

## Accessor pattern
Use `{{PATTERN}}` from spec 2a (EXPLICIT or SUBVIEW).

## Kernel update order
Update P3 kernels in phase-process order:
1. Main driver
2. Evaporation/sublimation (vapor production)
3. Cloud water conservation (vapor consumption)
4. Remaining kernels as needed

This order minimizes risk - get the main driver working first, then phase processes.

# Ask Before

1. If P3 standalone tests fail BFB and investigation reveals issue with accessor pattern
2. If P3 requires special handling different from SHOC (suggests pattern isn't general)
3. If number of kernel files significantly exceeds estimate (~3-5 files)

# Out of Scope

- Other water species (qc, qi, qr, qs) - those are extended in specs 3 and 4
- RRTMGP, ZM, dynamics, surface coupling - specs 2d-2g
- Multi-tracer microphysics - future campaigns

# Notes

## Why P3 Second?

P3 is the critical test of the pattern:
1. **Complexity**: Many more kernels than SHOC
2. **Impact**: Microphysics errors compound quickly, so BFB is essential
3. **Coupling**: P3 receives qv from SHOC, so this tests the chain
4. **Confidence**: Success here means the pattern is robust

## Kernel Count Uncertainty

P3 has many kernel impl files. The deliverables list (~3 files) is an estimate. Actual count may be higher. This is acceptable - the pattern applies uniformly.

## BFB Risk

P3 has the highest BFB risk of all processes because:
1. Most kernel call sites touching qv
2. Complex phase transition logic with many floating-point operations
3. Any evaporation/condensation error shows up in precipitation

This is why P3 comes early in the sequence - fail fast if the pattern doesn't work.

## Success Criteria Placeholders

The `{{PATTERN}}` placeholder in commands should be replaced with the result from spec 2a (EXPLICIT or SUBVIEW). This carries forward from spec 2b.
