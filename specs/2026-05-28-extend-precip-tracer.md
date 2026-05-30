---
spec_id: 2026-05-28-extend-precip-tracer
spec_type: model-e3sm
spec_version: 1
title: "Extend qr and surface precipitation fluxes to tracer dimension"
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
  - tests/water_tracers/test_precip_tracer_access.cpp

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

# Extend qr and surface precipitation fluxes to tracer dimension

## Purpose

Extend P3 precipitation fields to support tracers:
- `qr` (rain mixing ratio): `(col, lev)` → `(tracer, col, lev)`
- `precip_liq_surf_mass` (liquid precip flux): `(col)` → `(tracer, col)`
- `precip_ice_surf_mass` (ice precip flux): `(col)` → `(tracer, col)`

Completes Group 1 array extension. Follows PR 2-6 SUBVIEW pattern.

**Note**: P3 does not use separate `qs`, `rain`, `snow` 3D fields. All ice mixing 
ratio is in `qi` (extended in PR 3). Surface precipitation is tracked via 2D 
accumulated fluxes.

## Approach

Reference skills: `eamxx-cpp-conventions`, `regression-baseline`, `e3sm-build-and-test`.

1. Update field registration in P3 process interface:
   - Convert `qr` from 2D to 3D tracer layout
   - Convert `precip_liq_surf_mass` from scalar2d to 2D tracer layout
   - Convert `precip_ice_surf_mass` from scalar2d to 2D tracer layout
2. Update P3 kernels using SUBVIEW pattern (established in PRs 2b-6)
3. Update rain sedimentation to handle tracer dimension (passive for slots 1+)
4. Verify microphysics processes remain bulk-only (slot 0)
5. Write unit test for field access patterns

## Implementation pattern

Follow SUBVIEW pattern from PRs 2b-6:

```cpp
// Field registration in eamxx_p3_process_interface.cpp
const int num_tracers = SCREAM_NUM_TRACERS;

// 3D rain mixing ratio
auto qr_layout = m_grid->get_3d_tracer_layout(num_tracers);
add_field<Updated>("qr", qr_layout, kg/kg, grid_name);

// 2D surface precipitation fluxes
auto precip_layout = m_grid->get_2d_tracer_layout(num_tracers);
add_field<Updated>("precip_liq_surf_mass", precip_layout, kg/m2, grid_name);
add_field<Updated>("precip_ice_surf_mass", precip_layout, kg/m2, grid_name);

// Kernel access - use subview for slot 0
auto qr_bulk = Kokkos::subview(qr, 0, Kokkos::ALL, Kokkos::ALL);
auto precip_liq_bulk = Kokkos::subview(precip_liq_surf_mass, 0, Kokkos::ALL);
auto precip_ice_bulk = Kokkos::subview(precip_ice_surf_mass, 0, Kokkos::ALL);

// Sedimentation - loop over tracers for passive advection
for (int itracer = 0; itracer < num_tracers; ++itracer) {
  auto qr_t = Kokkos::subview(qr, itracer, Kokkos::ALL, Kokkos::ALL);
  // Apply same sedimentation velocity to all tracers
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

Primarily P3 microphysics changes:
- `qr` is the only 3D precipitation field in P3 (rain mixing ratio)
- P3 does NOT use `qs` - all ice is tracked via `qi` (extended in PR 3)
- P3 does NOT use 3D `rain`/`snow` rates - uses 2D surface flux accumulators
- Rain sedimentation handles tracer dimension - slots 1+ use same velocity as slot 0 (passive)
- Microphysics processes (autoconversion, collection) remain bulk-only (slot 0)

After PR 4 merges, all water reservoir fields have tracer dimension. Ready for PR 5
validation test.
