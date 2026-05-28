---
spec_id: 2026-05-26-extend-cloud-multi-tracer
spec_type: model-eamxx
spec_version: 1
title: "EAMxx: extend qc, qi to (COL, species, LEV)"
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
    - specs/2026-05-26-water-species-concept.md
    - specs/2026-05-26-extend-qv-multi-tracer.md
    - wiso_campaign_plan.md
  data: []
  baseline:
    path: /data/baselines/wiso-dev-F2010-SCREAMv1-ne4pg2_ne4pg2/h0.0001-01-06.nc
    type: bfb-comparison

deliverables:
  - "Edit components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp lines 79 (qc) and 81 (qi): route both through scream::water_tracers::add_tracer_multi following the pattern PR 2 established for qv. Lines 80 (qr) and 82 (qm) stay scalar (PR 4)."
  - "Edit components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp line 84 (qc): same routing. SHOC does not register qi."
  - "Audit and update every P3 / SHOC consumer that reads qc or qi as a (COL, LEV) view; route through get_bulk_water_subview. Includes diagnostic fields that derive from qc/qi (eff_radius_qc, eff_radius_qi, qc_sed, qi_sed) — these are macrostep diagnostics; verify whether their internal qc/qi reads need the subview wrapper. Precise file list in planning checkpoint."
  - "components/eamxx/src/physics/water_tracers/tests/water_tracers_cloud_passive_copy_test.cpp — Catch2 test, gated SCREAM_TRACE_WATER=ON AND SCREAM_NUM_WATER_TRACERS>=2. Allocates qv + qc + qi with N=2, copies CMP=0 into CMP=1 for each, advances through the P3 + SHOC paths, asserts species[k](:,1,:) == species[k](:,0,:) for k in {qv, qc, qi} to machine epsilon. Exercises cross-species accretion / autoconversion / deposition / sublimation paths in P3's level loop."
  - "Updated comment in components/eamxx/src/physics/water_tracers/water_tracers.hpp crossing qc, qi off the deferred-species list and pointing at this spec."

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
    cmd: "cd components/eamxx && ctest --test-dir build/wiso-on-n1 --output-on-failure -E 'water_tracers_cloud_passive_copy'"
    expect: exit_zero
    phase: testing

  - id: passive-copy-cloud-n2
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/wiso-on-n2 -R '^water_tracers_cloud_passive_copy$' --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: components/eamxx/src/physics/water_tracers/tests/water_tracers_cloud_passive_copy_test.cpp
      - claim: "qc(:,1,:), qi(:,1,:) track CMP=0 bit-for-bit through P3 + SHOC cross-species paths"

  - id: tier2-short-run-stable
    type: shell
    cmd: "cd cime/scripts && ./case.submit --no-batch --case /code/E3SM/EAMXX-wiso/cases/wiso-03-cloud && /code/E3SM/EAMXX-wiso/cases/wiso-03-cloud/check_for_NaN.sh"
    expect: exit_zero
    phase: testing

  - id: bfb-slice0-vs-pre-campaign-baseline
    type: comparison
    cmd: "bash scripts/wiso/extract_slice0.sh /code/E3SM/EAMXX-wiso/cases/wiso-03-cloud/run/h0.0001-01-06.nc && cprnc /data/baselines/wiso-dev-F2010-SCREAMv1-ne4pg2_ne4pg2/h0.0001-01-06.nc /code/E3SM/EAMXX-wiso/cases/wiso-03-cloud/run/h0.0001-01-06.slice0.nc"
    expect: "BFB"
    phase: testing
    verifies:
      - claim: "slice-1-unchanged: extending qc, qi to (COL, CMP, LEV) does not change physics at species index 0"

  - id: water-mass-closure
    type: tolerance
    cmd: "python scripts/check_water_balance.py /code/E3SM/EAMXX-wiso/cases/wiso-03-cloud/run/h0.0001-01-06.nc"
    metric: "fractional_imbalance"
    rtol: 1e-12
    atol: 1e-15
    phase: testing

  - id: physical-reasonableness-cloud
    type: human-review
    description: "Inspect qc and qi at slice 0 against the pre-campaign baseline: vertical profiles in cloud-bearing columns, time evolution, and global means must look identical to eye. Inspect slice 1 (N=2 build) — passive duplicate must track slice 0."
    phase: integration

  - id: architectural-readiness-for-pr4
    type: human-review
    description: "Confirm the pattern is mechanical enough that PR 4 can apply the same edits to qr, qm. Specifically: the cross-species level-loop interactions in P3 (autoconversion, accretion, sedimentation between cloud and precip) compiled and ran correctly under N=2 here, so PR 4's edits to the precip side should be a straightforward extension."
    phase: integration

out_of_scope:
  - "qv extension (PR 2 territory; this PR consumes its result)."
  - "qr, qm extension (PR 4 territory)."
  - "nc, nr, ni, bm number-concentration tracers."
  - "Isotope physics, fractionation hooks, WaterTracerHook."
  - "Surface fluxes / boundary conditions for water isotopes."
  - "IO schema rewrite beyond accepting the new dim."
  - "Restart / checkpoint handling for the new dim."
  - "Modification of atmosphere_process.hpp add_tracer template."
  - "Changes to the SCREAM_TRACE_WATER / SCREAM_NUM_WATER_TRACERS CMake mechanism (inherited from PR 1)."
  - "Changes to water_species.hpp, water_tracer_registration.hpp, or get_bulk_water_subview (inherited from PR 1)."
  - "Any changes under components/eam/ (legacy Fortran EAM) or other components/."
  - "Performance tuning for large N."

resolved_decisions:
  - "Strict downstream of PR 2 (extend-qv-multi-tracer). The qv path is already lifted; this PR applies the identical pattern to qc and qi. No mechanism decisions are re-litigated."
  - "Species set: qc (cloud liquid) and qi (cloud ice) only. Precipitation species qr (rain) and qm (rimed-ice mass) are deferred to PR 4. EAMxx P3 has no separate qs; qm represents combined rimed-ice mass."
  - "Combined passive-copy test exercises qv + qc + qi together rather than qc and qi in isolation. Reason: P3's level loop interleaves accretion, autoconversion, deposition, sublimation across species; bugs in the rank lift that only manifest in cross-species interactions would be missed by per-species tests."
  - "Same pre-campaign baseline as PR 2 — the slice-0 BFB target is unchanged. The new field reshape is qc / qi only; qr / qm still scalar at this PR. Baseline at /data/baselines/wiso-dev-F2010-SCREAMv1-ne4pg2_ne4pg2/."
  - "Diagnostic fields (eff_radius_qc, eff_radius_qi, qc_sed, qi_sed) are *not* lifted to the multi-tracer dim. They are scalar macrostep diagnostics of the bulk tracer. If their internal qc/qi reads bypass the field-manager subview path, planning escalates."

ask_before:
  - Modifying any file outside components/eamxx/.
  - Modifying atmosphere_process.hpp, FieldTag, PointGrid signatures.
  - Modifying SCREAM_TRACE_WATER / SCREAM_NUM_WATER_TRACERS CMake (PR 1).
  - Modifying water_species.hpp / water_tracer_registration.hpp (PR 1).
  - Modifying the qv lift in P3 / SHOC (PR 2).
  - Adding qr / qm to the rank-lift here (PR 4).
  - Wiring any WaterTracerHook calls (group 3).
  - Generating, replacing, or deleting any baseline file under /data/baselines/.
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

# EAMxx: extend qc, qi to (COL, species, LEV)

## Context

PR 3 in the group-1 chain. PR 2 lifted `qv`; this PR extends the
same pattern to `qc` (cloud liquid) and `qi` (cloud ice). PR 4 will
follow with `qr` (rain) and `qm` (rimed-ice mass). The species split
across three PRs (not four) is because P3's cross-species level loop
groups vapor/condensates/precip into natural review units —
condensates and precip have distinct phase-change call sets and
distinct cross-species accretion paths.

The risk profile is similar to PR 2 with one addition: P3's level
loop interleaves autoconversion (qc → qr), accretion (qc + qr → qr,
qi + qr → qi), and riming (qc + qi → qm) within a single column
sweep. A passive-copy test that runs each species in isolation would
miss bugs that only manifest when qc → qr while qi → qr happens in
the same level. The combined `qv + qc + qi` N=2 test deliberately
exercises those interactions; cross-precip interactions (qc → qr,
qi → qr) wait until PR 4 because qr is still scalar here.

The slice-0 BFB criterion is the same shape as PR 2 against the
same pre-campaign baseline, but it now also covers qc / qi paths.
If the rank lift quietly perturbs cloud condensation or
autoconversion arithmetic, the cprnc diff catches it.

## Approach

Files in parentheses are likely edit sites; precise diffs land in
the planning checkpoint.

1. **Prerequisite check.** PR 2's deliverables on this branch: qv
   registered through `add_tracer_multi`, slice0 extraction script
   under `scripts/wiso/`, passive-copy harness pattern. Confirmed
   by the implementation-phase compile criteria.

2. **P3 qc, qi registration.**
   `eamxx_p3_process_interface.cpp:79, 81` — route both through
   `scream::water_tracers::add_tracer_multi`. Lines 80 (qr), 82
   (qm) untouched.

3. **SHOC qc registration.**
   `eamxx_shoc_process_interface.cpp:84` — same routing. SHOC has
   no qi registration.

4. **Consumer reads.** Every P3/SHOC TU that pulls `qc` or `qi` by
   name from the field manager and assumes `(COL, LEV)` is wrapped
   through `get_bulk_water_subview`. No `#ifdef SCREAM_TRACE_WATER`
   in process sources. Diagnostic fields (`eff_radius_qc`,
   `eff_radius_qi`, `qc_sed`, `qi_sed`) audited during planning —
   they should already go through the FM subview path, but if any
   read qc/qi from a raw device pointer that bypasses field
   manager, escalate.

5. **Combined passive-copy test.**
   `tests/water_tracers_cloud_passive_copy_test.cpp`. Allocates qv
   + qc + qi with N=2, copies CMP=0 → CMP=1 for each, advances
   through P3 + SHOC, asserts species[k](:,1,:) == species[k](:,0,:)
   for k in {qv, qc, qi}. Gated on `SCREAM_TRACE_WATER=ON AND
   SCREAM_NUM_WATER_TRACERS GREATER_EQUAL 2`.

6. **Tier-2 5-day run + slice-0 BFB.** Same compset/resolution/length
   as PR 2 (F2010-SCREAMv1, ne4pg2_ne4pg2, 5 days). Same baseline
   at `/data/baselines/wiso-dev-F2010-SCREAMv1-ne4pg2_ne4pg2/`. The
   slice0 extraction script (delivered by PR 2) handles the cprnc
   prep.

7. **Water-mass closure.** Per `tracer-conservation`,
   `rtol=1e-12, atol=1e-15`. The closure equation now sums qv + qc
   + qi tendencies; the qr + qm terms are still scalar (single
   slice) and contribute as before.

Risks:

- **Cross-species level-loop interactions.** Autoconversion qc → qr,
  riming qc + qi → qm, deposition qv → qi, sublimation qi → qv all
  share the same level loop. The combined N=2 test catches bugs
  here; the slice-0 BFB confirms physics at the bulk slice is
  unchanged.
- **Diagnostic outputs.** `eff_radius_qc/qi`, `qc_sed/qi_sed`
  derive from qc/qi internally. If they take raw views that bypass
  the FM subview path, this PR must wrap them. Audit during
  planning.
- **Cross-process subview alias.** If P3 returns a `qc` View and
  SHOC reads it, both must see the same CMP=0 slice. The FM
  subview accessor is single-source-of-truth; verified by the
  combined test.

Skills relied on: `regression-baseline`, `tracer-conservation`,
`eamxx-cpp-conventions`, `e3sm-build-and-test`, `e3sm-platforms`.

## References

- `wiso_campaign_plan.md` — group 1 chain.
- `specs/2026-05-26-water-species-concept.md` — PR 1 (foundation).
- `specs/2026-05-26-extend-qv-multi-tracer.md` — PR 2; pattern
  inherited.
- `references/old_specs/2026-05-06-eamxx-water-condensates-add-tracer-dim.md`
  — prior attempt covering qc/qi/qr/qm in one PR; lessons absorbed
  (split per-class for review tractability, combined cross-species
  test).
- `components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp:79, 81`
  — qc, qi registration sites.
- `components/eamxx/src/physics/p3/impl/p3_main_impl_part2.hpp`
  (around lines 410-450) — `update_prognostic_ice` /
  `update_prognostic_liquid`; cited as context for why hook
  wiring is deferred (group 3).
- `claude-workflows/skills/tracer-conservation/SKILL.md`.
- `claude-workflows/skills/regression-baseline/SKILL.md`.
