---
spec_id: 2026-05-26-extend-qv-multi-tracer
spec_type: model-eamxx
spec_version: 1
title: "Extend water-vapor array (qv) to multi-tracer-species dimension"
created: 2026-05-26
author: rfiorella
project: EAMxx-wiso

inputs:
  source_files:
    - components/eamxx/src/share/atm_process/atmosphere_process.hpp
    - components/eamxx/src/share/field/field_manager.hpp
    - components/eamxx/src/physics/water_tracers/        # created by prior spec (PR 1 in chain)
  data: []
  baseline:
    path: /data/baselines/wiso-dev-precampaign/<dp_compset>-<dp_resolution>
    type: bfb-comparison
  # Chain context: this spec's branch (wiso-02-extend-qv-multi-tracer) is
  # based on the prior spec's branch (wiso-01-water-species-concept), NOT
  # on wiso-dev directly. The species enum + add_tracer wiring from PR 1
  # must be present in the working tree for this spec to build.

deliverables:
  - components/eamxx/src/physics/water_tracers/qv_multi_tracer.hpp
  - components/eamxx/src/physics/water_tracers/qv_multi_tracer.cpp
  - components/eamxx/tests/water_tracers/test_qv_slice0_bfb.cpp
  - components/eamxx/tests/water_tracers/test_qv_multi_tracer.cpp
  - components/eamxx/tests/water_tracers/CMakeLists.txt

success_criteria:
  - id: compile-clean
    type: shell
    cmd: "docker exec e3sm-orch-2026-05-26-extend-qv-multi-tracer bash -lc './case.build'"
    expect: exit_zero
    phase: implementation

  - id: lint-new-lines-cpp
    # Diff against the immediate parent branch in the chain (PR 1), not wiso-dev,
    # so this spec only lints lines it actually added.
    type: shell
    cmd: "docker exec e3sm-orch-2026-05-26-extend-qv-multi-tracer bash -lc 'git clang-format --diff $(git merge-base HEAD wiso-01-water-species-concept)..HEAD -- components/eamxx/src/physics/water_tracers/ | tee /tmp/clang-fmt.diff; test ! -s /tmp/clang-fmt.diff'"
    expect: exit_zero
    phase: implementation

  - id: unit-multi-tracer
    type: shell
    cmd: "docker exec e3sm-orch-2026-05-26-extend-qv-multi-tracer bash -lc 'ctest -R water_tracers_qv_multi_tracer --output-on-failure'"
    expect: exit_zero
    phase: testing

  - id: dp-run-no-nan
    type: shell
    cmd: "docker exec e3sm-orch-2026-05-26-extend-qv-multi-tracer bash -lc './case.submit && ./check_for_NaN.sh'"
    expect: exit_zero
    phase: testing

  - id: slice0-bfb-regression
    type: comparison
    cmd: "docker exec e3sm-orch-2026-05-26-extend-qv-multi-tracer bash -lc 'ncks -O -d tracer_species,0,0 -v qv ./run/output/*.h0.*.nc /tmp/slice0.nc && cprnc /home/e3smuser/baselines/wiso-dev-precampaign/<dp_compset>-<dp_resolution>/h0.final.nc /tmp/slice0.nc'"
    expect: "BFB"
    phase: testing

  - id: water-mass-closure
    type: tolerance
    cmd: "docker exec e3sm-orch-2026-05-26-extend-qv-multi-tracer bash -lc 'python scripts/check_water_balance.py ./run/output/*.h0.*.nc'"
    metric: "fractional_imbalance"
    rtol: 1e-12
    atol: 1e-15
    phase: testing

  - id: multi-tracer-config-n2
    type: shell
    cmd: "docker exec e3sm-orch-2026-05-26-extend-qv-multi-tracer bash -lc 'ctest -R water_tracers_qv_slice0_bfb --output-on-failure'"
    expect: exit_zero
    phase: testing

  - id: physical-reasonableness
    type: human-review
    description: "qv field at final timestep shows no unexpected spatial structure; slice-0 visually identical to baseline; slice-1 (HDO) field has plausible magnitude and pattern when tracer_count=2."
    phase: integration

out_of_scope:
  - Modifying transport routines (HOMME or any advection driver)
  - Extending qc, qi, qr, qs (those are PRs 3 and 4)
  - Implementing any fractionation arithmetic (groups 2 and 3)
  - Changing the field-manager API beyond adding species-dim metadata
  - Touching any code outside components/eamxx/src/physics/water_tracers/
    and the field-manager headers listed in inputs.source_files

resolved_decisions:
  - "Tracer-species dimension is the leading dimension of multi-tracer arrays: A[species, col, lev] in C++ Kokkos view layout."
  - "Species index 0 = bulk H2O. This is the slice that must be BFB against the pre-campaign baseline."
  - "Multi-tracer arrays are registered through FieldManager with a `tracer_species` dimension; tracer indices are never hard-coded."
  - "CMake variable `SCREAM_TRACER_COUNT` controls compile-time tracer dimension; default 1."
  - "Single-tracer build (`SCREAM_TRACER_COUNT=1`) must be slice-0 BFB; this is the regression criterion."

ask_before:
  - Modifying any file outside components/eamxx/src/physics/water_tracers/ and the field-manager headers listed in inputs.source_files
  - Adding new namelist parameters or YAML config keys outside the wiso block
  - Bumping `SCREAM_TRACER_COUNT` default away from 1
  - Submitting full-grid Tier-2 runs (this spec uses DP only)

execution:
  mode: checkpoint
  checkpoints:
    - after: planning
      requires: user-confirmation
    - after: implementation
      requires: tests-pass
    - after: testing
      requires: user-confirmation
    - after: integration
      requires: user-confirmation
  max_iterations_per_phase: 5
  parallelization:
    allowed: false
    max_subagents: 1
    plan: null

post_completion_review:
  enabled: true
  reviewers:
    - performance-specialist
    - code-reviewer
  blocking: false

model_specific:
  validation_tier: 2
  target_compset: <dp_compset>               # inherits from project-config dp_compset (DP harness used for this spec)
  target_resolution: <dp_resolution>
  stop_n: 5
  stop_option: ndays
  platform: docker-local
  container_image: rfiorella/model-containers:e3sm-openmpi-dev-latest
  tracer_count: 1                            # default build under regression; multi-tracer tested via unit test
  dp_proxy: true                             # DP is primary validation surface for this spec
  cmake_option: SCREAM_TRACER_COUNT
---

# Extend water-vapor array (qv) to multi-tracer-species dimension

## Context

This spec extends EAMxx's water-vapor field `qv` from a rank-2 array
(column, level) to a rank-3 array (tracer_species, column, level). The
species dimension is the leading dimension to match the convention used
throughout `scream::water_tracers` and to make slice operations
cache-coherent in the layout Kokkos uses on this platform. The spec
follows immediately after the species-enum spec
(`2026-05-26-water-species-concept.md`) and is the foundational
structural change for the wiso campaign.

The core challenge is proving that this is a pure refactor: when the
build is configured with `SCREAM_TRACER_COUNT=1`, every value at
`qv[0, :, :]` must equal the value the previous scalar code produced
at `qv[:, :]`. This is the `slice-1-unchanged` invariant from the
`tracer-conservation` skill. Failure to maintain it means existing
physics has been silently perturbed, which is not acceptable.

Group-1 PRs 3 and 4 (cloud and precipitation arrays) mirror the
pattern set here. The work done in this spec defines the contract those
specs follow.

## Approach

Reference skills:
- `eamxx-cpp-conventions` — Kokkos view layout, EKAT Pack interaction
  with the new leading dimension, FieldManager registration pattern
- `tracer-conservation` — invariants, tolerance defaults
- `regression-baseline` — BFB proof procedure
- `e3sm-build-and-test` — Tier 2 mechanics
- `e3sm-platforms` (DP reference doc) — DP harness invocation
- `scientific-modeling-conventions` — comments, units

Decomposition:

1. **Field-manager metadata.** Add `tracer_species` dimension metadata
   to FieldManager and a registration helper in
   `qv_multi_tracer.hpp`. Bulk H2O is registered at species index 0.
2. **CMake plumbing.** Introduce `SCREAM_TRACER_COUNT` (default 1).
   Compile-time size of the species dim.
3. **Field rewiring.** Change every read/write of `qv` to use the
   species-aware accessor. For `SCREAM_TRACER_COUNT=1`, the accessor
   must compile down to the previous scalar access pattern (verified
   by the slice-0 BFB criterion).
4. **Tests.** Two Catch2 unit tests under
   `components/eamxx/tests/water_tracers/`:
   - `test_qv_slice0_bfb.cpp` — instantiates the multi-tracer array
     with N=1, runs a pinned mini-driver against pinned inputs,
     compares against scalar-path output bit-for-bit.
   - `test_qv_multi_tracer.cpp` — instantiates with N=2 and N=4, checks
     that slice-0 remains identical when the additional species are
     zero-initialized, and that the species-aware accessor returns the
     expected pointer arithmetic.
5. **DP run.** With `SCREAM_TRACER_COUNT=1`, run the DP compset for
   5 days. Compare final-timestep history against the pre-campaign
   baseline. Must be BFB.

Risks:

- Kokkos layout: adding a leading dim can change memory layout and
  hence FP operation order. If full-grid BFB later turns out to be
  unachievable for this reason, the DP-harness BFB stands in as the
  proxy (see `regression-baseline` skill, non-BFB fallback section).
- EKAT Pack interaction: if `qv` is currently packed in the level
  dimension, the new leading species dim must not break pack
  alignment. Pack size must remain the same as in the rest of the file
  (see eamxx-cpp-conventions).
- FieldManager API change visibility: this PR touches field-manager
  headers used elsewhere in EAMxx. Confined to additive changes; no
  existing call sites should require modification.

## References

- `wiso_campaign_plan.md` — group 1 section
- `2026-05-26-water-species-concept.md` — prior spec defining the
  species enum
- CAM water-tracer implementation in `references/CAM/wiso/` for layout
  precedent (note: CAM uses Fortran column-major; EAMxx uses Kokkos
  views, layout chosen at view construction)
- Bony et al. 2008 — fractionation philosophy (relevant from PR 6
  onward; cited here only for design context)
