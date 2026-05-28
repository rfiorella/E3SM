---
spec_id: 2026-05-26-extend-qv-multi-tracer
spec_type: model-eamxx
spec_version: 1
title: "EAMxx: extend qv to (COL, species, LEV) via water_tracers helper"
created: 2026-05-26T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

inputs:
  source_files:
    - components/eamxx/src/physics/water_tracers/water_tracers.hpp
    - components/eamxx/src/physics/water_tracers/water_species.hpp
    - components/eamxx/src/physics/water_tracers/water_tracer_registration.hpp
    - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
    - components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
    - components/eamxx/src/share/atm_process/atmosphere_process.hpp
    - components/eamxx/src/share/grid/point_grid.cpp
    - components/eamxx/src/share/field/field_group.hpp
    - specs/2026-05-26-water-species-concept.md
    - references/old_specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md
    - wiso_campaign_plan.md
  data: []
  baseline:
    path: /data/baselines/wiso-dev-F2010-SCREAMv1-ne4pg2_ne4pg2/h0.0001-01-06.nc
    type: bfb-comparison

deliverables:
  - "Edit components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp line 78: replace `add_tracer<Updated>(\"qv\", m_grid, kg/kg, ps)` with a call routed through scream::water_tracers::add_tracer_multi (delivered by PR 1) that produces a (COL, CMP=WTRC_MAX_CNST, LEV) layout with CMP dim name \"water_tracer\". qc/qi/qr/qm registrations at lines 79-82 stay scalar (deferred to PRs 3-4)."
  - "Edit components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp line 65: same change to the qv add_tracer<Updated> call. qc registration at line 84 stays scalar (PR 3)."
  - "In each consumer kernel that previously read qv as a (COL, LEV) view, wire the read through get_bulk_water_subview (delivered by PR 1) so the kernel body is unchanged. No `#ifdef SCREAM_TRACE_WATER` in P3 or SHOC source files. Touched files: P3 process interface .cpp + any P3 impl headers that take a qv view by name; SHOC process interface .cpp + any SHOC impl headers that take a qv view by name. Precise list lands in the planning checkpoint."
  - "components/eamxx/src/physics/water_tracers/tests/water_tracers_qv_passive_copy_test.cpp — Catch2 test, gated SCREAM_TRACE_WATER=ON AND SCREAM_NUM_WATER_TRACERS>=2. Allocates qv with N=2, copies the CMP=0 state into CMP=1, advances through the P3 + SHOC qv-touching paths for a small number of steps with the standalone test harness, and asserts qv(:,1,:) == qv(:,0,:) to machine epsilon at every output step."
  - "scripts/wiso/extract_slice0.sh — wrapper around ncks that extracts CMP=0 of every variable carrying the water_tracer dim from a history file, writing a sibling .slice0.nc. Used by the BFB comparison criterion."
  - "Updated comment in components/eamxx/src/physics/water_tracers/water_tracers.hpp listing remaining species to lift (qc, qi for PR 3; qr, qm for PR 4) and pointing at this spec for the established pattern."

success_criteria:
  - id: compile-default-off
    type: shell
    cmd: "cd components/eamxx && rm -rf build/wiso-default && cmake -S . -B build/wiso-default -DCMAKE_BUILD_TYPE=Debug && cmake --build build/wiso-default -j"
    expect: exit_zero
    phase: implementation

  - id: compile-on-n1
    type: shell
    cmd: "cd components/eamxx && rm -rf build/wiso-on-n1 && cmake -S . -B build/wiso-on-n1 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_TRACE_WATER=ON -DSCREAM_NUM_WATER_TRACERS=1 && cmake --build build/wiso-on-n1 -j"
    expect: exit_zero
    phase: implementation

  - id: compile-on-n2
    type: shell
    cmd: "cd components/eamxx && rm -rf build/wiso-on-n2 && cmake -S . -B build/wiso-on-n2 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_TRACE_WATER=ON -DSCREAM_NUM_WATER_TRACERS=2 && cmake --build build/wiso-on-n2 -j"
    expect: exit_zero
    phase: implementation

  - id: clang-format-check
    type: shell
    cmd: "git diff --name-only $(git merge-base HEAD wiso-dev)..HEAD -- 'components/eamxx/**/*.hpp' 'components/eamxx/**/*.cpp' | xargs -r clang-format --dry-run --Werror"
    expect: exit_zero
    gate: advisory
    on_fail: skip
    phase: implementation

  - id: existing-tests-pass-default
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/wiso-default --output-on-failure"
    expect: exit_zero
    phase: testing

  - id: existing-tests-pass-on-n1
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/wiso-on-n1 --output-on-failure -E 'water_tracers_qv_passive_copy'"
    expect: exit_zero
    phase: testing

  - id: passive-copy-qv-n2
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/wiso-on-n2 -R '^water_tracers_qv_passive_copy$' --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: components/eamxx/src/physics/water_tracers/tests/water_tracers_qv_passive_copy_test.cpp
      - claim: "qv(:,1,:) tracks qv(:,0,:) bit-for-bit through P3 + SHOC qv-touching paths"

  - id: tier2-short-run-stable
    type: shell
    cmd: "cd cime/scripts && ./case.submit --no-batch --case /code/E3SM/EAMXX-wiso/cases/wiso-02-qv && /code/E3SM/EAMXX-wiso/cases/wiso-02-qv/check_for_NaN.sh"
    expect: exit_zero
    phase: testing

  - id: bfb-slice0-vs-pre-campaign-baseline
    type: comparison
    cmd: "bash scripts/wiso/extract_slice0.sh /code/E3SM/EAMXX-wiso/cases/wiso-02-qv/run/h0.0001-01-06.nc && cprnc /data/baselines/wiso-dev-F2010-SCREAMv1-ne4pg2_ne4pg2/h0.0001-01-06.nc /code/E3SM/EAMXX-wiso/cases/wiso-02-qv/run/h0.0001-01-06.slice0.nc"
    expect: "BFB"
    phase: testing
    verifies:
      - claim: "slice-1-unchanged: extending qv to (COL, CMP, LEV) does not change physics at species index 0"
      - deliverable: scripts/wiso/extract_slice0.sh

  - id: water-mass-closure
    type: tolerance
    cmd: "python scripts/check_water_balance.py /code/E3SM/EAMXX-wiso/cases/wiso-02-qv/run/h0.0001-01-06.nc"
    metric: "fractional_imbalance"
    rtol: 1e-12
    atol: 1e-15
    phase: testing

  - id: physical-reasonableness-qv
    type: human-review
    description: "Inspect qv at slice 0 against the pre-campaign baseline qv: spatial pattern, vertical structure, and time-evolution must look identical to eye. Inspect qv at slice 1 (N=2 build only): for a passive duplicate initialization it should also track slice 0 visually."
    phase: integration

  - id: architectural-readiness-for-pr3
    type: human-review
    description: "Confirm the pattern is mechanical enough that PR 3 can apply the same edits to qc, qi (and PR 4 to qr, qm) without redesign: the add_tracer_multi call replaces add_tracer<Updated>; consumer reads route through get_bulk_water_subview; the passive-copy test extends to additional species."
    phase: integration

out_of_scope:
  - "qc, qi, qr, qm extension (PR 3-4 territory)."
  - "nc, nr, ni, bm number-concentration tracers."
  - "Isotope physics, fractionation hooks, WaterTracerHook abstraction."
  - "Surface fluxes / boundary conditions for water isotopes."
  - "IO schema rewrite beyond accepting the new water_tracer dim."
  - "Restart / checkpoint handling for the new dim."
  - "Modification of atmosphere_process.hpp add_tracer template (the multi-tracer path layers above it via add_tracer_multi)."
  - "Any changes under components/eam/ (legacy Fortran EAM) or other components/."
  - "Performance tuning for large N."
  - "Changes to the SCREAM_TRACE_WATER / SCREAM_NUM_WATER_TRACERS CMake mechanism (inherited from PR 1; do not re-litigate)."

resolved_decisions:
  - "Bulk-H2O index is 0 (per spec PR 1, scream::water_tracers::is_bulk(0) == true). The slice-1-unchanged invariant from the tracer-conservation skill is the load-bearing test: bfb-slice0-vs-pre-campaign-baseline is the proof."
  - "Default-OFF build keeps qv at (COL, CMP=1, LEV) — not (COL, LEV). The unit CMP dim is the trivial FieldGroup case (field_group.hpp:34-37, m_subview_dim = 1). This unifies the data path across OFF / ON builds; consumer kernels see a single (COL, LEV) view via get_bulk_water_subview regardless."
  - "Pre-campaign baseline lives at /data/baselines/wiso-dev-F2010-SCREAMv1-ne4pg2_ne4pg2/. The orchestrator runs the regression-baseline workflow against this directory before declaring the spec done. If the baseline does not exist, the orchestrator pauses at the planning checkpoint and asks the user to generate it per the regression-baseline skill."
  - "Compset/resolution for the Tier-2 BFB run match the pre-campaign baseline: F2010-SCREAMv1 at ne4pg2_ne4pg2, 5-day run (stop_n=5, stop_option=ndays). These mirror EAMxx-developer compset conventions in AGENTS.md and are what wiso_campaign_plan.md group 1 boundary regenerates."
  - "IO: the qv field in NetCDF output gains a water_tracer dim of size 1 in default builds, exactly as in the rejected predecessor spec (references/old_specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md, resolved_decisions). Downstream tooling that hard-codes rank-2 qv may break; the slice0 extraction script in deliverables is the mitigation for baseline comparison."
  - "P3 + SHOC are the only call sites touched in this PR. Dynamics, MAM, RRTMGP, surface coupler — none read qv directly with rank assumptions today; they read via field manager lookup which returns the registered layout. If grep during planning surfaces a hard-coded rank-2 qv consumer outside P3 / SHOC, planning escalates rather than silently extending scope."

ask_before:
  - Modifying any file outside components/eamxx/.
  - Modifying atmosphere_process.hpp add_tracer template.
  - Modifying FieldTag enum or PointGrid layout-factory signatures.
  - Modifying the SCREAM_TRACE_WATER / SCREAM_NUM_WATER_TRACERS CMake blocks (PR 1 owns these).
  - Modifying water_species.hpp / water_tracer_registration.hpp delivered by PR 1.
  - Generating, replacing, or deleting any baseline file under /data/baselines/ — the pre-campaign baseline is consumed read-only.
  - Adding qc/qi/qr/qm to the rank-lift here (deferred to PRs 3-4).
  - Wiring any WaterTracerHook calls into P3 / SHOC (deferred to group 3).
  - Touching components/eam/ (legacy Fortran EAM).

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
  target_compset: F2010-SCREAMv1
  target_resolution: ne4pg2_ne4pg2
  stop_n: 5
  stop_option: ndays
  platform: docker-local
  container_image: rfiorella/model-containers:e3sm-openmpi-dev-latest
  build_mode: eamxx-standalone-cmake
  test_driver: components/eamxx/scripts/test-all-eamxx
---

# EAMxx: extend qv to (COL, species, LEV) via water_tracers helper

## Context

This is PR 2 of the group-1 chain in the EAMxx-wiso campaign
(`campaigns/wiso.campaign.md`). PR 1 landed the foundation —
`scream::water_tracers` namespace, the species enum, the
`add_tracer_multi` helper, `WTRC_MAX_CNST`, and `get_bulk_water_subview`.
This PR is the first one to actually change a production field's
allocation: `qv` is lifted from `(COL, LEV)` to `(COL, CMP, LEV)` in
both P3 and SHOC, in all builds.

The reason `qv` is its own PR rather than bundled with the
condensates is risk surface. `qv` is read or written by every
parameterization that touches water: P3, SHOC, dynamics, the surface
coupler, RRTMGP via specific humidity. Most of those reads go through
field-manager lookup and are layout-agnostic, but any hard-coded
rank-2 assumption surfaces as a compile error or a slice-0 BFB
failure. Isolating `qv` lets that risk land cleanly before the
condensates compound it.

The success criterion that matters is `bfb-slice0-vs-pre-campaign-baseline`.
The `slice-1-unchanged` invariant (`tracer-conservation` skill) says
that after extending an array by a leading or trailing species
dimension with N=1, the bulk-water slice must reproduce the pre-
extension scalar values bit-for-bit. With `WTRC_MAX_CNST=1` the qv
View has shape `(COL, 1, LEV)`; the CMP=0 subview is byte-identical
in layout to the previous `(COL, LEV)` qv. If physics changes, the
diff catches it.

The N=2 passive-copy test is the second proof. A passive duplicate of
qv at species index 1 must track species index 0 bit-for-bit through
P3 and SHOC. That demonstrates the CMP dim is plumbed end-to-end
through the data path, not merely allocated.

## Approach

Files in parentheses are likely edit sites; precise diffs land in the
planning checkpoint.

1. **Prerequisite check.** PR 1's deliverables must be on this
   branch: `scream::water_tracers::add_tracer_multi`,
   `get_bulk_water_subview`, the species enum,
   `SCREAM_TRACE_WATER`/`SCREAM_NUM_WATER_TRACERS` CMake options.
   Confirmed by `compile-default-off` and `compile-on-n1` succeeding
   in the implementation phase.

2. **P3 qv registration.** In
   `eamxx_p3_process_interface.cpp:78`, replace
   `add_tracer<Updated>("qv", m_grid, kg/kg, ps)` with a call
   through `scream::water_tracers::add_tracer_multi("qv", m_grid,
   kg/kg, ps)` that forwards the returned FieldRequest into
   `add_field<Updated>`. The lines 79-82 (qc, qr, qi, qm) are not
   touched in this PR.

3. **SHOC qv registration.** Same change at
   `eamxx_shoc_process_interface.cpp:65`. Line 84 (qc) is not
   touched in this PR.

4. **Consumer-side reads.** Audit every TU in P3 and SHOC that
   takes `qv` by name from the field manager and assumes
   `(COL, LEV)`. Each such site routes its qv view through
   `get_bulk_water_subview(qv_rank3)` so the kernel body is
   unchanged. No `#ifdef SCREAM_TRACE_WATER` in P3/SHOC sources.

5. **Audit beyond P3/SHOC.** Grep for `"qv"` reads across EAMxx
   physics + share trees. Any consumer outside P3/SHOC that assumes
   rank-2 surfaces here. If found, escalate at planning checkpoint
   — do not silently expand scope.

6. **Passive-copy unit test.**
   `tests/water_tracers_qv_passive_copy_test.cpp`. Allocates qv with
   N=2 through a standalone test harness, copies the CMP=0 state
   into CMP=1, runs the qv-touching physics for a few steps, asserts
   `qv(:,1,:) == qv(:,0,:)` to machine epsilon at every step.
   Registered in `tests/CMakeLists.txt` gated on
   `SCREAM_TRACE_WATER=ON AND SCREAM_NUM_WATER_TRACERS GREATER_EQUAL 2`.

7. **slice0 extraction script.**
   `scripts/wiso/extract_slice0.sh`: wraps `ncks` to extract the
   CMP=0 slice of every variable carrying the `water_tracer` dim
   from an EAMxx history file, writing a `<file>.slice0.nc` sibling.
   Consumed by the BFB criterion. The cprnc comparison is then
   `cprnc <pre-campaign baseline> <run>.slice0.nc`.

8. **Tier-2 5-day run + slice-0 BFB.** Build a CIME case (F2010-
   SCREAMv1 at ne4pg2_ne4pg2, 5 days) on the working branch.
   Confirm no NaN/Inf. Run the slice0-extraction + cprnc against
   `/data/baselines/wiso-dev-F2010-SCREAMv1-ne4pg2_ne4pg2/h0.0001-01-06.nc`.
   Must be BFB.

9. **Water-mass closure.** The standard column water-mass closure
   diagnostic (qv tendencies + surface fluxes integrated over the
   column) holds within `rtol=1e-12, atol=1e-15` per the
   `tracer-conservation` default. This is invariant against the rank
   lift; failure would indicate a tendency-write site missed the
   subview accessor.

Risks:

- **Hidden rank-2 qv consumers outside P3/SHOC.** The audit in step
  5 is the mitigation. If the grep surfaces unexpected callers,
  planning halts and the user decides whether to expand scope or
  patch with a temporary squeeze.
- **Pack alignment.** `LEV` is innermost so packs vectorize over
  levels. The `(COL, CMP, LEV)` layout preserves that. A unit CMP
  dim adds one extra outer stride and does not break Pack alignment.
  Verified by the slice-0 BFB.
- **IO schema change.** Default-OFF builds gain a `water_tracer` dim
  of size 1 on qv in NetCDF output. Downstream tools that hard-code
  rank-2 qv may break. The slice0 extraction script is the
  workaround for baseline comparison; broader downstream impact is
  documented but not fixed here.
- **Baseline staleness.** If `wiso-dev` has advanced and the
  pre-campaign baseline no longer reflects parent SHA at branch
  point, the regression-baseline skill's refresh decision tree
  applies. Planning checkpoint surfaces this if `BASELINE.txt` in
  the baseline directory does not match `git merge-base HEAD master`.

Skills relied on: `regression-baseline` (slice-0 BFB workflow,
baseline directory layout, `BASELINE.txt` discipline),
`tracer-conservation` (slice-1-unchanged invariant + water-mass
closure tolerances), `eamxx-cpp-conventions` (subview discipline,
Kokkos device callability, EKAT idioms), `e3sm-build-and-test`
(Tier-2 conventions, cprnc usage), `e3sm-platforms` (docker-local
Tier-2 run mechanics).

## References

- `wiso_campaign_plan.md` — group 1 definition and Tier-2 boundary.
- `specs/2026-05-26-water-species-concept.md` — PR 1; this spec
  consumes its deliverables.
- `references/old_specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md`
  — prior attempt; informs the always-rank-3 / unit-CMP-OFF design.
- `components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp:78`
  — qv registration site.
- `components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp:65`
  — qv registration site.
- `components/eamxx/src/share/grid/point_grid.cpp:93-104` —
  `get_3d_vector_layout` factory.
- `components/eamxx/src/share/field/field_group.hpp:34-37` —
  `m_subview_dim = 1` invariant.
- `claude-workflows/skills/tracer-conservation/SKILL.md` —
  invariants and tolerances.
- `claude-workflows/skills/regression-baseline/SKILL.md` — slice-0
  BFB workflow.
