---
spec_id: 2026-05-28-extend-precip-tracer
spec_type: model-e3sm
spec_version: 1
title: "Extend qr, qs, rain, snow to tracer dimension"
created: 2026-05-28T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

inputs:
  source_files:
    - wiso_group1_campaign_revision.md
    - specs/2026-05-28-extend-qv-tracer.md
    - specs/2026-05-28-extend-cloud-tracer.md
    - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
    - components/eamxx/src/physics/p3/impl/p3_rain_sed_impl.hpp
    - components/eamxx/src/physics/p3/impl/p3_ice_sed_impl.hpp
  data: []
  baseline: /data/baselines/pre-group1

deliverables:
  - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
  - components/eamxx/src/physics/p3/impl/p3_main_impl.hpp
  - components/eamxx/src/physics/p3/impl/p3_rain_sed_impl.hpp
  - components/eamxx/src/physics/p3/impl/p3_ice_sed_impl.hpp
  - components/eamxx/src/physics/p3/impl/p3_autoconversion_impl.hpp
  - components/eamxx/src/physics/p3/impl/p3_collection_impl.hpp
  - tests/water_tracers/test_precip_tracer_access.cpp
  - tests/water_tracers/test_precip_transport.cpp

success_criteria:
  - id: compile-n1
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr4-n1 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=1 && cmake --build build/pr4-n1 -j"
    expect: exit_zero
    phase: implementation

  - id: compile-n3
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/pr4-n3 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_NUM_TRACERS=3 && cmake --build build/pr4-n3 -j"
    expect: exit_zero
    phase: implementation

  - id: unit-test-precip-access
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/pr4-n1 -R test_precip_tracer_access --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: tests/water_tracers/test_precip_tracer_access.cpp

  - id: component-test-precip-transport
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/pr4-n3 -R test_precip_transport --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: tests/water_tracers/test_precip_transport.cpp

  - id: p3-with-precip-bfb
    type: comparison
    cmd: "cd components/eamxx && ctest --test-dir build/pr4-n1 -R p3_standalone --output-on-failure"
    expect: all_bfb
    rtol: 0
    atol: 0
    phase: testing
    verifies:
      - claim: "P3 microphysics with precipitation passes BFB"

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

  - id: full-atm-multiprocess-bfb
    type: comparison
    cmd: "cd components/eamxx && ctest --test-dir build/pr4-n1 -R 'homme.*p3.*rrtmgp' --output-on-failure"
    expect: all_bfb
    rtol: 0
    atol: 0
    phase: testing
    verifies:
      - claim: "Full atmosphere multi-process integration passes BFB"

  - id: cumulative-performance
    type: shell
    cmd: "cd components/eamxx && ./scripts/benchmark-eamxx.sh --tracers 1 > bench_pr4.txt && python3 scripts/compare_benchmarks.py /data/baselines/pre-group1/benchmark.txt bench_pr4.txt --max-overhead 0.02"
    expect: exit_zero
    gate: advisory
    phase: testing
    verifies:
      - claim: "Cumulative overhead (PRs 2+3+4) < 2%"

out_of_scope:
  - Isotope-aware sedimentation - future campaigns
  - Multi-tracer validation test - PR 5
  - Convective vs stratiform precipitation separation

resolved_decisions:
  pattern_reuse: "Follow exact pattern from PR 2-3. Field registration, kernel updates, subview usage identical."
  sedimentation: "Sedimentation handles tracer dimension. For slots 1+, apply same sedimentation velocity as slot 0 (passive advection)."
  microphysics_processes: "Autoconversion, collection, etc. remain bulk-only (slot 0). Future campaigns add fractionation."

ask_before:
  - Modifying P3 microphysics test expectations
  - Changes to sedimentation algorithm
  - Modifications to precipitation rate calculations

model_specific:
  validation_tier: 2
  compset: F2010-SCREAMv1
  resolution: ne4pg2_ne4pg2
  baseline_tag: pre-group1
  dp_proxy: false
---

# Extend qr, qs, rain, snow to tracer dimension

## Purpose

Extend precipitation fields `qr` (rain mixing ratio), `qs` (snow mixing ratio),
`rain` (rain rate), `snow` (snow rate) from shape `(col, lev)` to `(tracer, col, lev)`.
Completes Group 1 array extension. Follows PR 2-3 pattern.

## Approach

Reference skills: `eamxx-cpp-conventions`, `regression-baseline`, `e3sm-build-and-test`.

1. Apply PR 2-3 pattern to precipitation fields
2. Update field registration in P3 microphysics
3. Update kernels: `qr(icol, ilev)` → `qr(0, icol, ilev)`, same for qs, rain, snow
4. Update sedimentation to handle tracer dimension (passive advection for slots 1+)
5. Verify autoconversion, collection remain bulk-only (slot 0)
6. Write unit and component tests

## Implementation pattern

Identical to PR 2-3:

```cpp
// Field registration
const int num_tracers = SCREAM_NUM_TRACERS;
auto layout = grid.get_3d_vector_layout(true, num_tracers, "tracer");
add_field<Required>("qr", layout, kg/kg, grid);
add_field<Required>("qs", layout, kg/kg, grid);
add_field<Required>("rain", layout, kg/m^2/s, grid);
add_field<Required>("snow", layout, kg/m^2/s, grid);

// Kernel access - bulk only for microphysics processes
qr(0, icol, ilev) = ...
qs(0, icol, ilev) = ...

// Sedimentation - loop over tracers for passive advection
for (int itracer = 0; itracer < num_tracers; ++itracer) {
  // Apply same sedimentation velocity to all tracers
  qr(itracer, icol, ilev) -= sed_flux / dz;
}
```

## Test execution

```bash
cd components/eamxx
./scripts/test-all-eamxx -m docker-local --baseline-dir /data/baselines/pre-group1 -c
```

Expected: All tests PASS BFB with rtol=0, atol=0 when `SCREAM_NUM_TRACERS=1`.

Critical tests:
- `p3_standalone` - microphysics with precipitation
- Full atmosphere multi-process integration

## Notes

Primarily P3 microphysics changes. Sedimentation must handle tracer dimension
properly - for slots 1+, apply same velocity as slot 0 (passive tracer behavior).

After PR 4 merges, all water reservoirs have tracer dimension. Ready for PR 5
validation test.
