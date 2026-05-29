---
spec_id: 2026-05-28-extend-qv-tracer
spec_type: model-e3sm
spec_version: 1
title: "Extend qv to tracer dimension"
created: 2026-05-28T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

inputs:
  source_files:
    - wiso_group1_campaign_revision.md
    - components/eamxx/src/share/field/field_layout.hpp
    - components/eamxx/src/share/field/field_manager.hpp
    - components/eamxx/src/share/field/field_manager.cpp
    - components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
    - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
    - components/eamxx/src/physics/rrtmgp/eamxx_rrtmgp_process_interface.cpp
    - components/eamxx/src/physics/zm/eamxx_zm_process_interface.cpp
    - components/eamxx/src/dynamics/homme/eamxx_homme_process_interface.cpp
    - components/eamxx/src/physics/surface_coupling/eamxx_surface_coupling_exporter.cpp
  data: []
  baseline: /data/baselines/pre-group1

deliverables:
  - components/eamxx/src/share/field/field_layout.hpp
  - components/eamxx/src/share/field/field_layout.cpp
  - components/eamxx/src/share/field/field_manager.cpp
  - components/eamxx/src/share/field/field_manager.hpp
  - components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
  - components/eamxx/src/physics/shoc/impl/shoc_main_impl.hpp
  - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
  - components/eamxx/src/physics/p3/impl/p3_main_impl.hpp
  - components/eamxx/src/physics/rrtmgp/eamxx_rrtmgp_process_interface.cpp
  - components/eamxx/src/physics/zm/eamxx_zm_process_interface.cpp
  - components/eamxx/src/dynamics/homme/eamxx_homme_process_interface.cpp
  - components/eamxx/src/physics/surface_coupling/eamxx_surface_coupling_exporter.cpp
  - tests/water_tracers/test_qv_tracer_access.cpp
  - tests/water_tracers/test_qv_transport.cpp
  - tests/water_tracers/CMakeLists.txt

success_criteria:
  - id: compile-scream-num-tracers-1
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr2-n1 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=1 && cmake --build build/pr2-n1 -j"
    expect: exit_zero
    phase: implementation

  - id: compile-scream-num-tracers-3
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr2-n3 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=3 && cmake --build build/pr2-n3 -j"
    expect: exit_zero
    phase: implementation

  - id: unit-test-qv-access
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/pr2-n1 -R test_qv_tracer_access --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: tests/water_tracers/test_qv_tracer_access.cpp
      - claim: "Slot-0 access equivalent to old scalar access, tracer dimension correctly sized"

  - id: component-test-qv-transport
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/pr2-n3 -R test_qv_transport --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: tests/water_tracers/test_qv_transport.cpp
      - claim: "Independent tracer slots remain independent through physics timestep"

  - id: bfb-vs-pre-campaign
    type: comparison
    cmd: "cd components/eamxx && ./scripts/test-all-eamxx -m docker-local --baseline-dir /data/baselines/pre-group1 -c"
    expect: all_bfb
    rtol: 0
    atol: 0
    gate: blocking
    phase: testing
    verifies:
      - claim: "SCREAM_NUM_TRACERS=1 produces BFB-identical results vs pre-campaign baseline"
    on_fail: halt_after_investigation
    resolution_notes: "If BFB fails: bisect to find which kernel broke BFB. Check subview usage, slot-0 indexing. If architectural issue confirmed by PR 1 gate, apply documented fallback strategy."

  - id: performance-regression
    type: shell
    cmd: "cd components/eamxx && ./scripts/benchmark-eamxx.sh --tracers 1 > bench_pr2.txt && python3 scripts/compare_benchmarks.py /data/baselines/pre-group1/benchmark.txt bench_pr2.txt --max-overhead 0.02"
    expect: exit_zero
    gate: advisory
    phase: testing
    verifies:
      - claim: "Performance overhead < 2% vs pre-campaign baseline"
    on_fail: skip

  - id: no-segfaults
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/pr2-n1 --output-on-failure"
    expect: exit_zero
    phase: testing

out_of_scope:
  - Other water species (qc, qi, qr, qs, qm, rain, snow) - PRs 3-4
  - Multi-tracer validation test - PR 5
  - Isotope-specific kernels - future campaigns

resolved_decisions:
  field_layout: "Leading tracer dimension: qv(tracer, col, lev). Kokkos vectorization favors column-innermost."
  subview_pattern: "Use Kokkos::subview(qv, 0, Kokkos::ALL, Kokkos::ALL) to pass bulk slice to existing kernels."
  slot_indexing: "All existing physics accesses qv(0, icol, ilev) explicitly for slot-0 bulk water."

ask_before:
  - Modifying existing test expectations
  - Changing field registration API beyond adding tracer dimension
  - Force-pushing if BFB fails

model_specific:
  validation_tier: 2
  compset: F2010-SCREAMv1
  resolution: ne4pg2_ne4pg2
  baseline_tag: pre-group1
  dp_proxy: false
---

# Extend qv to tracer dimension

## Purpose

Extend water vapor field `qv` from shape `(col, lev)` to `(tracer, col, lev)`.
Establishes pattern for PRs 3-4. Must preserve BFB for `SCREAM_NUM_TRACERS=1`.

## Approach

Reference skills: `eamxx-cpp-conventions`, `regression-baseline`, `e3sm-build-and-test`.

1. Add `TRACER` tag to `FieldLayout` enum
2. Update `FieldManager` to allocate tracer dimension based on `SCREAM_NUM_TRACERS`
3. Modify field registration in all processes touching `qv`: SHOC, P3, RRTMGP, ZM, dynamics, surface coupling
4. Update kernels: `qv(icol, ilev)` → `qv(0, icol, ilev)` or use subview pattern
5. Update restart/history I/O to handle tracer dimension
6. Write unit and component tests

## Implementation pattern

For each process interface file:

```cpp
// Old
add_field<Required>("qv", scalar3d_mid, kg/kg, grid);

// New
const int num_tracers = SCREAM_NUM_TRACERS;
auto layout = grid.get_3d_vector_layout(true, num_tracers, "tracer");
add_field<Required>("qv", layout, kg/kg, grid);
```

For kernels:

```cpp
// Old
qv(icol, ilev) = ...

// New - explicit slot 0
qv(0, icol, ilev) = ...

// Or use subview for entire kernel
auto qv_bulk = Kokkos::subview(qv, 0, Kokkos::ALL, Kokkos::ALL);
// pass qv_bulk to existing kernel
```

## Test execution

```bash
cd components/eamxx
./scripts/test-all-eamxx -m docker-local --baseline-dir /data/baselines/pre-group1 -c
```

Expected: All tests PASS BFB with rtol=0, atol=0 when `SCREAM_NUM_TRACERS=1`.

## Critical BFB requirement

This is a hard gate. If BFB fails:
1. Use debugger to find where bit pattern diverges
2. Check Kokkos subview doesn't change memory layout
3. Verify all slot-0 accesses use explicit index 0
4. If architectural issue (confirmed by PR 1 gate results), apply documented fallback
