---
spec_id: 2026-05-28-extend-cloud-tracer
spec_type: model-e3sm
spec_version: 1
title: "Extend qc, qi, qm to tracer dimension"
created: 2026-05-28T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

inputs:
  source_files:
    - wiso_group1_campaign_revision.md
    - specs/2026-05-28-extend-qv-tracer.md
    - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
    - components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
    - components/eamxx/src/physics/cld_fraction/eamxx_cld_fraction_process_interface.cpp
    - components/eamxx/src/physics/rrtmgp/eamxx_rrtmgp_process_interface.cpp
    - components/eamxx/src/physics/mam/eamxx_mam_process_interface.cpp
  data: []
  baseline: /data/baselines/pre-group1

deliverables:
  - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
  - components/eamxx/src/physics/p3/impl/p3_main_impl.hpp
  - components/eamxx/src/physics/p3/impl/p3_ice_sed_impl.hpp
  - components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
  - components/eamxx/src/physics/shoc/impl/shoc_main_impl.hpp
  - components/eamxx/src/physics/cld_fraction/eamxx_cld_fraction_process_interface.cpp
  - components/eamxx/src/physics/rrtmgp/eamxx_rrtmgp_process_interface.cpp
  - components/eamxx/src/physics/mam/eamxx_mam_process_interface.cpp
  - tests/water_tracers/test_cloud_tracer_access.cpp
  - tests/water_tracers/test_cloud_transport.cpp

success_criteria:
  - id: compile-n1
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr3-n1 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=1 && cmake --build build/pr3-n1 -j"
    expect: exit_zero
    phase: implementation

  - id: compile-n3
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr3-n3 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=3 && cmake --build build/pr3-n3 -j"
    expect: exit_zero
    phase: implementation

  - id: unit-test-cloud-access
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/pr3-n1 -R test_cloud_tracer_access --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: tests/water_tracers/test_cloud_tracer_access.cpp

  - id: component-test-cloud-transport
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/pr3-n3 -R test_cloud_transport --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: tests/water_tracers/test_cloud_transport.cpp

  - id: p3-standalone-bfb
    type: comparison
    cmd: "cd components/eamxx && ctest --test-dir build/pr3-n1 -R p3_standalone --output-on-failure"
    expect: all_bfb
    rtol: 0
    atol: 0
    phase: testing
    verifies:
      - claim: "P3 microphysics passes BFB vs baseline"

  - id: shoc-standalone-bfb
    type: comparison
    cmd: "cd components/eamxx && ctest --test-dir build/pr3-n1 -R shoc_standalone --output-on-failure"
    expect: all_bfb
    rtol: 0
    atol: 0
    phase: testing
    verifies:
      - claim: "SHOC turbulence passes BFB vs baseline"

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

  - id: cumulative-performance
    type: shell
    cmd: "cd components/eamxx && ./scripts/benchmark-eamxx.sh --tracers 1 > bench_pr3.txt && python3 scripts/compare_benchmarks.py /data/baselines/pre-group1/benchmark.txt bench_pr3.txt --max-overhead 0.02"
    expect: exit_zero
    gate: advisory
    phase: testing
    verifies:
      - claim: "Cumulative overhead (PR 2 + PR 3) < 2%"

out_of_scope:
  - Precipitation species (qr, qs, rain, snow) - PR 4
  - Phase-change fractionation logic - future campaigns
  - Multi-tracer validation test - PR 5

resolved_decisions:
  pattern_reuse: "Follow exact pattern from PR 2 (extend qv). Field registration, kernel updates, subview usage all identical."
  derived_fields: "Fields like inv_qc_relvar, eff_radius_qc, eff_radius_qi remain scalar if computed solely from slot-0 bulk values."

ask_before:
  - Modifying P3 or SHOC test expectations
  - Changes to cloud fraction calculation
  - Modifications to radiation interface

model_specific:
  validation_tier: 2
  compset: F2010-SCREAMv1
  resolution: ne4pg2_ne4pg2
  baseline_tag: pre-group1
  dp_proxy: false
---

# Extend qc, qi, qm to tracer dimension

## Purpose

Extend cloud water fields `qc` (cloud liquid), `qi` (cloud ice), `qm` (mixed-phase)
from shape `(col, lev)` to `(tracer, col, lev)`. Follows PR 2 pattern.

## Approach

Reference skills: `eamxx-cpp-conventions`, `regression-baseline`, `e3sm-build-and-test`.

1. Apply PR 2 pattern to cloud water fields
2. Update field registration in P3, SHOC, cloud fraction, RRTMGP, MAM
3. Update kernels: `qc(icol, ilev)` → `qc(0, icol, ilev)`, same for qi, qm
4. Verify derived diagnostic fields (effective radius, cloud variance) remain scalar
5. Write unit and component tests

## Implementation pattern

Identical to PR 2, but for qc/qi/qm:

```cpp
// Field registration
const int num_tracers = SCREAM_NUM_TRACERS;
auto layout = grid.get_3d_vector_layout(true, num_tracers, "tracer");
add_field<Required>("qc", layout, kg/kg, grid);
add_field<Required>("qi", layout, kg/kg, grid);
add_field<Required>("qm", layout, kg/kg, grid);

// Kernel access
qc(0, icol, ilev) = ...
qi(0, icol, ilev) = ...
qm(0, icol, ilev) = ...
```

## Test execution

```bash
cd components/eamxx
./scripts/test-all-eamxx -m docker-local --baseline-dir /data/baselines/pre-group1 -c
```

Expected: All tests PASS BFB with rtol=0, atol=0 when `SCREAM_NUM_TRACERS=1`.

Priority single-process tests:
- `p3_standalone` - primary microphysics
- `shoc_standalone` - turbulence and cloud fraction

## Notes

More processes touch cloud water than vapor (radiation, aerosols, cloud fraction),
so file count ~40 vs PR 2's ~35. Pattern remains mechanical: add tracer dimension,
update slot-0 indexing.
