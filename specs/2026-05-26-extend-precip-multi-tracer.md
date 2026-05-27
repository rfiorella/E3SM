---
spec_id: 2026-05-26-extend-precip-multi-tracer
spec_type: model-e3sm
spec_version: 1
title: "EAMxx: extend qr, qm to (COL, species, LEV) — group 1 closer"
created: 2026-05-26T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

inputs:
  source_files:
    - components/eamxx/src/physics/water_tracers/water_tracers.hpp
    - components/eamxx/src/physics/water_tracers/water_species.hpp
    - components/eamxx/src/physics/water_tracers/water_tracer_registration.hpp
    - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
    - specs/2026-05-26-water-species-concept.md
    - specs/2026-05-26-extend-qv-multi-tracer.md
    - specs/2026-05-26-extend-cloud-multi-tracer.md
    - wiso_campaign_plan.md
  data: []
  baseline:
    path: /data/baselines/wiso-dev-F2010-SCREAMv1-ne4pg2_ne4pg2/h0.0001-01-06.nc
    type: bfb-comparison

deliverables:
  - "Edit components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp lines 80 (qr) and 82 (qm): route both through scream::water_tracers::add_tracer_multi following the pattern PRs 2-3 established. After this PR, all five P3 mass-mixing-ratio tracers (qv, qc, qr, qi, qm) carry the (COL, CMP, LEV) layout."
  - "Audit and update every P3 consumer that reads qr or qm as a (COL, LEV) view; route through get_bulk_water_subview. Sedimentation diagnostics (qr_sed) and any kernel that reads qr in cross-species accretion paths (qc → qr autoconv, qi → qr melt, qr → qv evap) are the highest-risk audit sites."
  - "components/eamxx/src/physics/water_tracers/tests/water_tracers_all_mass_passive_copy_test.cpp — Catch2 test, gated SCREAM_TRACE_WATER=ON AND SCREAM_NUM_WATER_TRACERS>=2. Allocates qv + qc + qi + qr + qm with N=2, copies CMP=0 into CMP=1 for each, advances through the full P3 + SHOC path, asserts species[k](:,1,:) == species[k](:,0,:) for k in {qv, qc, qi, qr, qm} to machine epsilon. This is the strongest cross-species integration check in the campaign — it exercises every P3 phase-change path under the multi-tracer layout."
  - "Updated comment in components/eamxx/src/physics/water_tracers/water_tracers.hpp marking the group-1 mass-tracer lift complete (qv, qc, qi, qr, qm all carry the species dim). Remaining deferred items: nc, nr, ni, bm (number concentrations) and WaterTracerHook wiring (group 3)."
  - "Regenerated post-group-1 Tier-2 baseline at /data/baselines/wiso-dev-postG1-F2010-SCREAMv1-ne4pg2_ne4pg2/ following the regression-baseline skill. This is the new baseline that groups 2-3 will compare against. BASELINE.txt records parent SHA = HEAD-after-merge of this PR onto wiso-dev, generation date, compset/resolution/stop options, and a note that the post-G1 baseline supersedes the pre-campaign baseline for downstream specs."

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
    cmd: "cd components/eamxx && ctest --test-dir build/wiso-on-n1 --output-on-failure -E 'water_tracers_all_mass_passive_copy'"
    expect: exit_zero
    phase: testing

  - id: passive-copy-all-mass-n2
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/wiso-on-n2 -R '^water_tracers_all_mass_passive_copy$' --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: components/eamxx/src/physics/water_tracers/tests/water_tracers_all_mass_passive_copy_test.cpp
      - claim: "All five P3 mass tracers (qv, qc, qi, qr, qm) track CMP=0 bit-for-bit through full P3 + SHOC cross-species level loop"

  - id: tier2-short-run-stable
    type: shell
    cmd: "cd cime/scripts && ./case.submit --no-batch --case /code/E3SM/EAMXX-wiso/cases/wiso-04-precip && /code/E3SM/EAMXX-wiso/cases/wiso-04-precip/check_for_NaN.sh"
    expect: exit_zero
    phase: testing

  - id: bfb-slice0-vs-pre-campaign-baseline
    type: comparison
    cmd: "bash scripts/wiso/extract_slice0.sh /code/E3SM/EAMXX-wiso/cases/wiso-04-precip/run/h0.0001-01-06.nc && cprnc /data/baselines/wiso-dev-F2010-SCREAMv1-ne4pg2_ne4pg2/h0.0001-01-06.nc /code/E3SM/EAMXX-wiso/cases/wiso-04-precip/run/h0.0001-01-06.slice0.nc"
    expect: "BFB"
    phase: testing
    verifies:
      - claim: "slice-1-unchanged: extending qr, qm to (COL, CMP, LEV) — closing group 1 — does not change physics at species index 0"

  - id: water-mass-closure
    type: tolerance
    cmd: "python scripts/check_water_balance.py /code/E3SM/EAMXX-wiso/cases/wiso-04-precip/run/h0.0001-01-06.nc"
    metric: "fractional_imbalance"
    rtol: 1e-12
    atol: 1e-15
    phase: testing

  - id: post-group1-baseline-generated
    type: shell
    cmd: "test -f /data/baselines/wiso-dev-postG1-F2010-SCREAMv1-ne4pg2_ne4pg2/BASELINE.txt && test -f /data/baselines/wiso-dev-postG1-F2010-SCREAMv1-ne4pg2_ne4pg2/h0.0001-01-06.nc && cprnc /data/baselines/wiso-dev-F2010-SCREAMv1-ne4pg2_ne4pg2/h0.0001-01-06.nc /data/baselines/wiso-dev-postG1-F2010-SCREAMv1-ne4pg2_ne4pg2/h0.0001-01-06.nc"
    expect: "BFB"
    phase: integration
    verifies:
      - deliverable: /data/baselines/wiso-dev-postG1-F2010-SCREAMv1-ne4pg2_ne4pg2/

  - id: physical-reasonableness-precip
    type: human-review
    description: "Inspect qr and qm at slice 0 against the pre-campaign baseline: precipitation rates, vertical profiles in active-precip columns, time evolution must look identical. Inspect slice 1 (N=2) — passive duplicate tracks slice 0."
    phase: integration

  - id: architectural-readiness-for-group2
    type: human-review
    description: "Confirm the group-1 structural extension is complete and ready for group 2 (fractionation primitives) to build on top: (a) qv, qc, qi, qr, qm all carry the (COL, CMP, LEV) layout in all builds; (b) consumer kernels are unchanged for the bulk path via get_bulk_water_subview; (c) the post-group-1 Tier-2 baseline is in place at /data/baselines/wiso-dev-postG1-...; (d) the passive-copy test pattern is generic enough that group-3 hook PRs can extend it without redesign."
    phase: integration

out_of_scope:
  - "qv, qc, qi extension (PRs 2-3; consumed here)."
  - "nc, nr, ni, bm number-concentration tracers (post-campaign cleanup)."
  - "Isotope physics, fractionation hooks, WaterTracerHook abstraction (group 2-3)."
  - "Surface fluxes / boundary conditions for water isotopes (group 3)."
  - "IO schema rewrite beyond accepting the new dim."
  - "Restart / checkpoint handling for the new dim."
  - "Modification of atmosphere_process.hpp add_tracer template."
  - "Changes to the SCREAM_TRACE_WATER / SCREAM_NUM_WATER_TRACERS CMake mechanism (PR 1)."
  - "Changes to water_species.hpp / water_tracer_registration.hpp / get_bulk_water_subview (PR 1)."
  - "Any changes under components/eam/ (legacy Fortran EAM) or other components/."
  - "Performance tuning for large N."

resolved_decisions:
  - "Strict downstream of PRs 2 and 3. The qv, qc, qi paths are already lifted; this PR applies the identical pattern to the precip species. No mechanism decisions are re-litigated."
  - "Species set: qr (rain) and qm (rimed-ice mass). EAMxx P3 has no separate qs; qm represents combined rimed-ice mass and serves as P3's combined snow/graupel proxy (per references/old_specs/2026-05-06-eamxx-water-condensates-add-tracer-dim.md resolved_decisions). Number-concentration tracers nc/nr/ni/bm are deferred."
  - "Combined N=2 passive-copy test now spans all five mass species. This is the strongest integration check in the campaign — it exercises P3's full level loop including precip-side autoconversion (qc → qr), riming (qc + qi → qm), melting (qi → qr), and rain evaporation (qr → qv) under the multi-tracer layout."
  - "Post-group-1 Tier-2 baseline at /data/baselines/wiso-dev-postG1-F2010-SCREAMv1-ne4pg2_ne4pg2/. This becomes the new BFB target for groups 2-3 (per wiso_campaign_plan.md group 1 boundary: Tier-2 full-grid baseline regenerated with the array machinery in place but tracer_count=1). The pre-campaign baseline is retained read-only as historical reference; BASELINE.txt in the new directory records the relationship."
  - "Baseline equivalence check: the new post-G1 baseline must be BFB against the pre-campaign baseline when reduced to slice-0 (since group 1 is by design a non-functional refactor). post-group1-baseline-generated criterion enforces this."

ask_before:
  - Modifying any file outside components/eamxx/.
  - Modifying atmosphere_process.hpp, FieldTag, PointGrid signatures.
  - Modifying SCREAM_TRACE_WATER / SCREAM_NUM_WATER_TRACERS CMake (PR 1).
  - Modifying water_species.hpp / water_tracer_registration.hpp (PR 1).
  - Modifying the qv / qc / qi lifts (PRs 2-3).
  - Wiring any WaterTracerHook calls (group 3).
  - Generating, replacing, or deleting any existing baseline file under /data/baselines/ — the new wiso-dev-postG1- baseline is created fresh here; the pre-campaign baseline is read-only.
  - Touching components/eam/ (legacy Fortran EAM).
  - Adding nc / nr / ni / bm to the rank-lift (post-campaign).

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

# EAMxx: extend qr, qm to (COL, species, LEV) — group 1 closer

## Context

PR 4 — the final PR of the group-1 chain. PR 2 lifted `qv`; PR 3
lifted `qc` and `qi`; this PR lifts `qr` (rain) and `qm` (rimed-ice
mass), completing the structural extension of every bulk water mass
mixing ratio in EAMxx P3. After this PR merges, the campaign's group
boundary triggers regeneration of the Tier-2 full-grid baseline at
`tracer_count = 1`, which becomes the BFB target for groups 2-3.

The combined passive-copy test in this PR is the strongest cross-
species integration check in the entire campaign. With all five
mass species at N=2, every phase-change path P3 supports —
autoconversion `qc → qr`, accretion `qc + qr → qr`, riming `qc + qi
→ qm`, melting `qi → qr`, sedimentation `qr ↓`, deposition `qv →
qi`, sublimation `qi → qv`, freezing `qr → qm` — runs simultaneously
across two species slices. A bug in the rank lift that only
manifests when, say, rain evaporation interacts with snow
sublimation would surface here but nowhere else in the chain.

The post-group-1 baseline generation is part of this PR's
deliverables. The `regression-baseline` skill defines the workflow:
build at the post-merge SHA, run the same compset/resolution/length,
copy the final-timestep history to a new baseline directory with
`BASELINE.txt`. The baseline equivalence check (slice-0 BFB against
the pre-campaign baseline) confirms group 1 is non-functional by
construction.

## Approach

Files in parentheses are likely edit sites; precise diffs land in
the planning checkpoint.

1. **Prerequisite check.** PR 2 and PR 3 deliverables on this
   branch: qv, qc, qi all routed through `add_tracer_multi`; the
   passive-copy test pattern; the slice0 extraction script; the
   pre-campaign baseline still readable. Implementation-phase
   compiles are the gate.

2. **P3 qr, qm registration.**
   `eamxx_p3_process_interface.cpp:80, 82` — route both through
   `scream::water_tracers::add_tracer_multi`. SHOC does not
   register qr or qm.

3. **Consumer reads.** Every P3 TU that pulls `qr` or `qm` by name
   from the field manager and assumes `(COL, LEV)` is wrapped
   through `get_bulk_water_subview`. Sedimentation paths and
   cross-species precip generation (autoconv, accretion, melt,
   freeze, rain evap) are the highest-risk audit sites; cite the
   list in the planning checkpoint.

4. **Combined passive-copy test.**
   `tests/water_tracers_all_mass_passive_copy_test.cpp`. Allocates
   qv + qc + qi + qr + qm with N=2, copies CMP=0 → CMP=1 for each,
   advances through the full P3 + SHOC level loop, asserts every
   species[k](:,1,:) == species[k](:,0,:) at every output step.
   Gated on `SCREAM_TRACE_WATER=ON AND SCREAM_NUM_WATER_TRACERS
   GREATER_EQUAL 2`.

5. **Tier-2 5-day run + slice-0 BFB.** Same compset/resolution/length
   as PRs 2-3. Same pre-campaign baseline. Slice0 extraction +
   cprnc must report BFB.

6. **Post-group-1 baseline regeneration.** Per the
   `regression-baseline` skill: with this PR merged onto
   `wiso-dev`, build at the new HEAD, run F2010-SCREAMv1 at
   ne4pg2_ne4pg2 for 5 days, confirm no NaN/Inf, copy the final-
   timestep history file to
   `/data/baselines/wiso-dev-postG1-F2010-SCREAMv1-ne4pg2_ne4pg2/`,
   write `BASELINE.txt` recording parent SHA, date, compset,
   resolution, stop options, and a pointer back to the pre-campaign
   baseline noting that this baseline supersedes it for groups
   2-3.

7. **Baseline equivalence check.** Run slice0 extraction on the new
   baseline and cprnc against the pre-campaign baseline. Must be
   BFB — group 1 is non-functional. This is the
   `post-group1-baseline-generated` criterion.

8. **Update the deferred-species note.** Mark qv, qc, qi, qr, qm
   complete. Remaining: nc, nr, ni, bm (number concentrations,
   post-campaign cleanup) and WaterTracerHook wiring (group 3).

Risks:

- **Cross-precip interactions hidden until now.** Rain evaporation
  back to qv, melting qi → qr, freezing qr → qm — these are the
  paths that go untested in PRs 2-3 because at least one of the
  participating species was still scalar. The all-mass N=2 test is
  the catcher.
- **Baseline generation timing.** The post-group-1 baseline is
  generated at the SHA the orchestrator computes after this PR
  merges into `wiso-dev`. If the orchestrator runs baseline
  generation pre-merge, the SHA in `BASELINE.txt` will not match.
  The integration checkpoint covers this — the user confirms
  generation happens at the right SHA.
- **Baseline equivalence failure.** If slice-0 BFB against the
  pre-campaign baseline fails on the regenerated post-G1
  baseline, group 1 has accidentally become functional and the
  diff must be root-caused before the chain proceeds.

Skills relied on: `regression-baseline` (slice-0 BFB workflow +
post-group baseline regeneration procedure), `tracer-conservation`
(slice-1-unchanged + water-mass closure), `eamxx-cpp-conventions`,
`e3sm-build-and-test` (Tier-2 mechanics), `e3sm-platforms`.

## References

- `wiso_campaign_plan.md` — group 1 boundary description and
  Tier-2 baseline regeneration policy.
- `specs/2026-05-26-water-species-concept.md` — PR 1 foundation.
- `specs/2026-05-26-extend-qv-multi-tracer.md` — PR 2 pattern.
- `specs/2026-05-26-extend-cloud-multi-tracer.md` — PR 3
  pattern; this PR extends the test from {qv, qc, qi} to {qv, qc,
  qi, qr, qm}.
- `references/old_specs/2026-05-06-eamxx-water-condensates-add-tracer-dim.md`
  — prior attempt; resolved_decisions on qm semantics inherited
  (qm is rimed-ice mass, P3's combined snow/graupel proxy).
- `components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp:80, 82`
  — qr, qm registration sites.
- `components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp:129-143`
  — P3 phase-change tendency registration list; informs hook
  placement in group 3 (out of scope here).
- `claude-workflows/skills/tracer-conservation/SKILL.md`.
- `claude-workflows/skills/regression-baseline/SKILL.md`.
