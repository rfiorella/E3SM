---
spec_id: 2026-05-28-water-tracer-metadata-and-gate
spec_type: model-e3sm
spec_version: 1
title: "Water-tracer metadata, types, and BFB feasibility gate"
created: 2026-05-28T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

inputs:
  source_files:
    - wiso_group1_campaign_revision.md
    - components/eamxx/src/share/field/field_layout.hpp
    - components/eamxx/src/share/field/field_manager.hpp
  data: []
  baseline: null

deliverables:
  - components/eamxx/src/physics/water_tracers/water_tracer_metadata.hpp
  - components/eamxx/src/physics/water_tracers/water_tracer_metadata.cpp
  - components/eamxx/src/physics/water_tracers/water_tracer_registry.hpp
  - components/eamxx/src/physics/water_tracers/water_tracer_registry.cpp
  - components/eamxx/src/physics/water_tracers/CMakeLists.txt
  - components/eamxx/src/physics/water_tracers/water_tracer_config.hpp.in
  - components/eamxx/src/physics/water_tracers/prototype/qv_extension_test.cpp
  - cmake/add_water_tracer.cmake
  - docs/wiso/tracer_data_model.md
  - tests/water_tracers/CMakeLists.txt

success_criteria:
  - id: compile-metadata-headers
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr1 -DCMAKE_BUILD_TYPE=Debug && cmake --build build/pr1 --target water_tracers -j"
    expect: exit_zero
    phase: implementation
    verifies:
      - deliverable: components/eamxx/src/physics/water_tracers/water_tracer_metadata.hpp

  - id: cmake-add-tracer-valid
    type: shell
    cmd: "cd components/eamxx && rm -rf build/pr1-tracer && cmake -S . -B build/pr1-tracer -DCMAKE_BUILD_TYPE=Debug && cmake --build build/pr1-tracer -j && grep -q 'SCREAM_NUM_TRACERS=2' build/pr1-tracer/CMakeCache.txt"
    expect: exit_zero
    phase: implementation
    verifies:
      - deliverable: cmake/add_water_tracer.cmake
      - claim: "add_water_tracer increments SCREAM_NUM_TRACERS"

  - id: prototype-compiles
    type: shell
    cmd: "cd components/eamxx/src/physics/water_tracers/prototype && g++ -std=c++17 -I../../../.. qv_extension_test.cpp -o qv_extension_test"
    expect: exit_zero
    phase: implementation
    verifies:
      - deliverable: components/eamxx/src/physics/water_tracers/prototype/qv_extension_test.cpp

  - id: prototype-bfb-gate
    type: shell
    cmd: "cd components/eamxx/src/physics/water_tracers/prototype && ./qv_extension_test && grep -E '(BFB: PASS|Max diff: 0\\.0)' qv_test_output.txt"
    expect: exit_zero
    gate: blocking
    phase: testing
    verifies:
      - claim: "SCREAM_NUM_TRACERS=1 produces BFB-identical results vs scalar baseline OR rtol < 1e-12"
    on_fail: halt_after_investigation
    resolution_notes: "If BFB fails: update docs/wiso/tracer_data_model.md with fallback strategy (template specialization or relaxed tolerance rtol=1e-12). Get confirmation before proceeding to PRs 2-5."

  - id: prototype-performance-gate
    type: shell
    cmd: "cd components/eamxx/src/physics/water_tracers/prototype && ./qv_extension_test --benchmark && awk '$0 ~ /Overhead:/ && $2 < 2.0 {exit 0} $0 ~ /Overhead:/ && $2 >= 2.0 {exit 1}' qv_test_output.txt"
    expect: exit_zero
    gate: blocking
    phase: testing
    verifies:
      - claim: "Performance overhead < 2% for SCREAM_NUM_TRACERS=1 vs scalar"
    on_fail: halt_after_investigation
    resolution_notes: "If >2% overhead: document fallback strategy in docs/wiso/tracer_data_model.md (e.g., template specialization for scalar path). Get confirmation before proceeding."

  - id: design-doc-exists
    type: shell
    cmd: "test -f docs/wiso/tracer_data_model.md && grep -q 'Slot-0 semantics' docs/wiso/tracer_data_model.md"
    expect: exit_zero
    phase: implementation
    verifies:
      - deliverable: docs/wiso/tracer_data_model.md
      - claim: "Design document defines slot-0 semantics, BFB strategy, field layout"

out_of_scope:
  - Actual field registration changes (PR 2)
  - Kernel modifications (PRs 2-4)
  - Multi-tracer validation test (PR 5)
  - Isotope-specific logic (future campaigns)
  - Fractionation functions (future campaigns)

resolved_decisions:
  bfb_gate_policy: "If prototype achieves exact BFB, PRs 2-4 enforce BFB as hard gate. If prototype fails BFB but within rtol=1e-12, document fallback in tracer_data_model.md and proceed with relaxed tolerance. If >1e-12, halt campaign for architectural review."
  approval_gate: "Skip human approval - proceed automatically if both gates pass"
  performance_fallback: "If overhead >2%, document template specialization strategy in design doc, proceed with implementation in PR 2"
  platform: "Run prototype on Charliecloud container (rfiorella/model-containers:e3sm-openmpi-dev-latest) on x86 host for consistent environment and performance measurement"

ask_before:
  - Proceeding to PR 2 if BFB gate fails
  - Proceeding to PR 2 if performance gate fails
  - Modifying any existing EAMxx test files

model_specific:
  validation_tier: 0
  compset: null
  resolution: null
  baseline_tag: null
  dp_proxy: false
---

# Water-tracer metadata, types, and BFB feasibility gate

## Purpose

Establish water-tracer metadata system, CMake registration function, and run
feasibility gates for BFB preservation and performance. This is a blocking PR
that must complete before PRs 2-5 can start.

## Approach

Reference skills: `eamxx-cpp-conventions`, `e3sm-build-and-test`, `scientific-modeling-conventions`.

1. Define `WaterTracerKind` enum and `WaterTracerMetadata` struct in namespace `scream::water_tracers`
2. Implement `WaterTracerRegistry` singleton for compile-time registration
3. Create `add_water_tracer()` CMake function that increments `SCREAM_NUM_TRACERS` and registers metadata
4. Write prototype test extending `qv(col,lev)` → `qv(tracer,col,lev)` to measure BFB and performance
5. Document results in `docs/wiso/tracer_data_model.md`

## Success criteria

All shell commands in `success_criteria` block must pass. BFB and performance
gates are blocking - if either fails, campaign pauses for design doc update and
confirmation before proceeding.

## Test execution

```bash
cd components/eamxx/src/physics/water_tracers/prototype
./qv_extension_test
./qv_extension_test --benchmark
```

Prototype must demonstrate:
- BFB preservation OR max difference < rtol=1e-12
- Performance overhead < 2%
