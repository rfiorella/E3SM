---
spec_id: 2026-05-06-eamxx-water-condensates-add-tracer-dim
spec_type: model-e3sm
spec_version: 1
title: "EAMxx: extend qc, qi, qr, qm with the tracer (CMP) dimension"
created: 2026-05-06T00:00:00-06:00
author: rfiorella
project: EAMXX-wiso

inputs:
  source_files:
    - components/eamxx/src/physics/water_tracers/water_tracers.hpp
    - components/eamxx/src/physics/water_tracers/water_types.hpp
    - components/eamxx/src/physics/water_tracers/water_isotopes.hpp
    - components/eamxx/src/physics/water_tracers/eamxx_water_tracers.cpp
    - components/eamxx/src/physics/water_tracers/CMakeLists.txt
    - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
    - components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
    - components/eamxx/src/share/atm_process/atmosphere_process.hpp
    - specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md
  data: []
  baseline: null

deliverables:
  - "qc, qi, qr, and qm fields each allocated as a 3D vector layout (COL, CMP, LEV) with CMP dim named \"water_tracer\" and CMP size = WTRC_MAX_CNST, in *all* builds. The mechanism (whatever was settled in the prior spec — wrapper helper around add_tracer, modified add_tracer, or per-call edits) is reused without redesign. With SCREAM_TRACE_WATER=OFF the CMP dim has size 1 and the layout is memory-equivalent to the prior scalar layouts."
  - "Subview-along-CMP-at-0 accessor extended (or trivially reused if the qv-spec accessor was generic) to qc, qi, qr, qm so that P3 (eamxx_p3_process_interface.cpp:79-82), SHOC (eamxx_shoc_process_interface.cpp), and any other bulk-water consumer see the same (COL, LEV) view per species as today. Consumer source files contain *no* `#ifdef SCREAM_TRACE_WATER`."
  - "Existing WaterTracerHook calls (established by the prior spec at qv-touching sites only) are *not* extended to additional phase-change sites in this spec — see resolved_decisions for why this is deferred. The hook interface itself is unchanged."
  - "New combined unit test under components/eamxx/src/physics/water_tracers/tests/, only registered when SCREAM_TRACE_WATER=ON and SCREAM_NUM_WATER_TRACERS >= 2, that allocates qv + qc + qi + qr + qm with N=2, copies the CMP=0 state into CMP=1 for each species, advances all five through the existing qv/qc/qi/qr/qm-touching physics for a small number of steps, and asserts species[i](:, 1, :) == species[i](:, 0, :) to machine epsilon at every step for every species. Cross-species accretion / autoconversion / evaporation paths are exercised because P3 mixes them within the same level loop."
  - "Updated follow-up note (block comment in water_tracers.hpp or sibling follow-up-specs.md) crossing qc, qi, qr, qm off the deferred-species list and adding a pointer to the in-substep-hook spec that will wire WaterTracerHook calls into P3 / SHOC at the per-substep phase-change sites. nc, nr, ni, bm remain on the deferred list."

success_criteria:
  - id: prerequisite-prior-spec-merged
    type: human-review
    description: "Confirm specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md has been completed and merged. Specifically: SCREAM_TRACE_WATER is a real CMake option (default OFF) defined per a real `option(...)` block at water_tracers/CMakeLists.txt:2 (not the unconditional definition that exists at start of this spec); WTRC_MAX_CNST is a build-time constant with the chosen mechanism wired; qv is allocated as (COL, CMP, LEV) in all builds; the WaterTracerHook interface and default no-op implementation exist; the qv-only passive-copy and OFF-vs-ON-N=1 BFB tests pass. This spec assumes those conditions and does not re-establish them."
    phase: planning

  - id: scope-confirmation
    type: human-review
    description: "Confirm the species set (qc, qi, qr, qm — the EAMxx P3 mass-mixing-ratio tracers excluding qv which the prior spec covered) and confirm the deferred items: nc, nr, ni, bm (number concentrations) and the wiring of WaterTracerHook calls to P3's per-substep phase-change emission sites (qr2qv_evap, qi2qv_sublim, qc2qr_accret, qc2qr_autoconv, qv2qi_vapdep, qc2qi_berg, qc2qr_ice_shed, qc2qi_collect, qr2qi_collect, qc2qi_hetero_freeze, qr2qi_immers_freeze, qi2qr_melt, sedimentation) and the equivalent SHOC sites. Both deferrals are intentional given that P3's update_prognostic_* calls update th_atm mid-loop, making post-hoc reading of registered tendency fields physically wrong for isotope fractionation."
    phase: planning

  - id: compile-clean-default-trace-water-off
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/default -DCMAKE_BUILD_TYPE=Debug && cmake --build build/default -j"
    expect: exit_zero
    phase: implementation

  - id: compile-clean-trace-water-on-n1
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/twon1 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_TRACE_WATER=ON -DSCREAM_NUM_WATER_TRACERS=1 && cmake --build build/twon1 -j"
    expect: exit_zero
    phase: implementation

  - id: compile-clean-trace-water-on-n2
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/twon2 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_TRACE_WATER=ON -DSCREAM_NUM_WATER_TRACERS=2 && cmake --build build/twon2 -j"
    expect: exit_zero
    phase: implementation

  - id: clang-format-check
    type: shell
    cmd: "git diff --name-only origin/master...HEAD -- 'components/eamxx/**/*.hpp' 'components/eamxx/**/*.cpp' | xargs -r clang-format --dry-run --Werror"
    expect: exit_zero
    gate: advisory
    on_fail: skip
    phase: implementation

  - id: existing-tests-pass-default-build
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/default --output-on-failure"
    expect: exit_zero
    phase: testing

  - id: existing-tests-pass-trace-water-on-n1
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/twon1 --output-on-failure -E 'water_mass_n2_passive_copy|water_tracers_qv_n2_passive_copy'"
    expect: exit_zero
    phase: testing

  - id: bfb-condensates-trace-water-on-n1-vs-default
    type: tolerance
    cmd: "cd components/eamxx && ctest --test-dir build/twon1 -R '^water_tracers_condensates_trace_water_on_n1_bfb$' --output-on-failure"
    metric: "max_abs_diff"
    rtol: 0
    atol: 0
    phase: testing

  - id: passive-copy-n2-all-mass-species
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/twon2 -R '^water_mass_n2_passive_copy$' --output-on-failure"
    expect: exit_zero
    phase: testing

  - id: architectural-readiness-for-followup
    type: human-review
    description: "Confirm the (COL, CMP, LEV) plumbing now established for qv + qc + qi + qr + qm is generic enough that (a) nc, nr, ni, bm follow the same pattern in their own follow-up spec (they are number concentrations with units 1/kg rather than kg/kg but the rank-lift is identical) and (b) the in-substep hook spec can wire WaterTracerHook calls into P3 and SHOC without further changes to field allocation or the subview accessor. If either is not true, capture the gap before declaring this spec done."
    phase: integration

out_of_scope:
  - "nc, nr, ni, bm — the number-concentration tracers — are deferred to a follow-up spec"
  - "Wiring WaterTracerHook calls into P3's per-substep phase-change emission sites (qr2qv_evap, qi2qv_sublim, qc2qr_accret, qc2qr_autoconv, qv2qi_vapdep, qc2qi_berg, qc2qr_ice_shed, qc2qi_collect, qr2qi_collect, qc2qi_hetero_freeze, qr2qi_immers_freeze, qi2qr_melt, plus sedimentation) — deferred to a separate in-substep-hook spec"
  - "Wiring WaterTracerHook calls into SHOC's vapor↔liquid phase-change sites — deferred to the same in-substep-hook spec"
  - "Isotope physics (equilibrium fractionation, kinetic fractionation, Rayleigh) — separate downstream specs"
  - "Surface fluxes / boundary conditions for water isotopes (ocean, land, sea ice)"
  - "IO writing of new tracer fields beyond what is needed for the in-test assertion"
  - "Restart / checkpoint of the new tracer dimension"
  - "Any changes under components/eam/ (legacy EAM Fortran) or other components/"
  - "Performance tuning of the (COL, CMP, LEV) layout for large N"
  - "Modifying the CMakeLists.txt option block for SCREAM_TRACE_WATER (the prior spec converts the unconditional definition into a real option; this spec inherits that)"
  - "Modifying or extending the WaterTracerHook interface (the prior spec settled it; only inherited call sites are touched here)"

resolved_decisions:
  - "This spec is a strict downstream of specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md. The prior spec must be completed and merged before this one starts. The prerequisite-prior-spec-merged success criterion is the explicit gate. All compile-time-mechanism decisions, hook-form decisions, IO-schema decisions, and the SCREAM_TRACE_WATER CMake-option conversion are inherited from the prior spec without re-litigation."
  - "The species in scope are qc, qi, qr, qm — the four bulk water mass mixing ratios in EAMxx P3 (eamxx_p3_process_interface.cpp:79-82) other than qv. There is no separate qs in EAMxx P3; qm represents rimed-ice mass and serves as P3's combined snow/graupel proxy. The wtrc_bulk_names array in water_tracers.hpp:58 (which still mentions RAINQM/SNOWQM/RAINQC/SNOWQC) is a holdover from EAM Fortran and is not the source of truth for which species exist in EAMxx; this spec aligns with the EAMxx P3 reality."
  - "Wiring WaterTracerHook calls to P3's additional phase-change emission sites is *deferred* to a separate spec, not absorbed here. Reason: P3 calls update_prognostic_ice and update_prognostic_liquid inside its tendency-computation level loop (p3_main_impl_part2.hpp around lines 410-450), and these update th_atm — and therefore T_atm — mid-step. The registered tendency fields (qr2qv_evap, qi2qv_sublim, qc2qr_accret, qc2qr_autoconv, qv2qi_vapdep, qc2qi_berg, qc2qr_ice_shed, qc2qi_collect, qr2qi_collect, qc2qi_hetero_freeze, qr2qi_immers_freeze, qi2qr_melt) are macrostep aggregates that have already lost the temperature trajectory needed for isotope fractionation. Adding a hook that reads those registered fields after-the-fact would produce physically incorrect fractionation. The correct hook placement is inside the substep loop adjacent to the per-substep tendency compute, which is non-trivial design work that deserves its own spec rather than being bundled with a structural rank-lift. Same logic applies to SHOC. Cited evidence: p3_main_impl.hpp:295,303,310 (adaptive sedimentation substepping) and p3_main_impl_part2.hpp `update_prognostic_*` calls inside the level loop."
  - "qc, qi, qr, qm are *always* allocated as (COL, CMP, LEV), inheriting the rank-lift discipline established for qv. With SCREAM_TRACE_WATER=OFF the CMP dim has size 1 and the layout is memory-equivalent to the prior scalar layouts; SIMD over LEV is unaffected. Default-build IO output for qc/qi/qr/qm gains a water_tracer dim of size 1, exactly as for qv (out-of-scope per the prior spec)."
  - "The new passive-copy unit test runs qv + qc + qi + qr + qm together with N=2, rather than per-species in isolation. This deliberately exercises P3's cross-species level-loop interactions (autoconversion, accretion, evaporation, deposition, sublimation, riming) and asserts each species' CMP=1 slice matches its CMP=0 slice through all those interactions. It is a stronger integration check than the prior qv-only test."
  - "All other deferred items from the prior spec (isotope physics, surface fluxes, IO/restart, EAM legacy code) remain deferred — this spec does not relax any of those exclusions."

ask_before:
  - Modifying any file outside components/eamxx/
  - Touching components/eam/ (legacy Fortran EAM)
  - Modifying cime_config/ or any compset / coupler XML
  - Changing the FieldTag enum or short-name aliases in field_tag.hpp
  - Modifying PointGrid::get_*_vector_layout signatures
  - Modifying the SCREAM_TRACE_WATER option block in water_tracers/CMakeLists.txt (inherited from the prior spec; do not re-litigate the form here)
  - Modifying the WaterTracerHook interface or its default no-op implementation (inherited from the prior spec)
  - Adding a new namelist parameter
  - Generating, replacing, or deleting any baseline file under /data/baselines/
  - Wiring WaterTracerHook calls to P3 or SHOC phase-change sites beyond what the prior spec already did (this is the deferred work; reject if encountered)

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
    max_subagents: 3
    plan: null

post_completion_review:
  enabled: false
  reviewers: []
  blocking: false

model_specific:
  validation_tier: 0
  target_compset: null
  target_resolution: ne4
  stop_n: 0
  stop_option: ndays
  platform: docker-local
  container_image: rfiorella/model-containers:e3sm-openmpi-dev-latest
  build_mode: eamxx-standalone-cmake
  test_driver: components/eamxx/scripts/test-all-eamxx
---

# EAMxx: extend qc, qi, qr, qm with the tracer (CMP) dimension

## Context

The prior spec, `specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md`, lifted
`qv` from `(COL, LEV)` to `(COL, CMP, LEV)` in all builds, with default builds
carrying a unit `CMP` dim. It established: a real `SCREAM_TRACE_WATER` CMake
option (default OFF), a build-time `WTRC_MAX_CNST`, a single subview-along-
`CMP`-at-0 accessor used by all bulk-water consumers, and a `WaterTracerHook`
interface with a no-op default that is called from the small set of qv-touching
phase-change emission sites the prior spec identified. This spec extends the
exact same rank-lift discipline to the four other bulk water mass mixing
ratios that EAMxx P3 carries: `qc` (cloud liquid), `qi` (cloud ice),
`qr` (rain), and `qm` (rimed-ice mass). After this spec, every bulk water mass
field in EAMxx physics has the same `(COL, CMP, LEV)` layout in all builds and
the same subview-at-0 access pattern.

This spec is **strictly downstream** of the prior spec. It assumes the
`SCREAM_TRACE_WATER` option block, the `WTRC_MAX_CNST` mechanism, the subview
accessor, and the `WaterTracerHook` interface all already exist and are
working; the `prerequisite-prior-spec-merged` planning-phase criterion is the
gate that confirms this. None of those decisions are re-litigated here. The
work in this spec is mechanical application of the prior spec's pattern to
four more fields, plus a stronger integration-style passive-copy unit test
that runs all five mass species (`qv, qc, qi, qr, qm`) together so cross-
species interactions in P3 (autoconversion, accretion, evaporation,
deposition, sublimation, riming) are exercised end-to-end on the new layout.

What is **not** in this spec — and the reasoning matters because it shapes
follow-up work — is the wiring of `WaterTracerHook` calls into the additional
phase-change emission sites that come with `qc, qi, qr, qm`. P3 already
exposes a long catalogue of registered phase-change tendency fields
(`qr2qv_evap, qi2qv_sublim, qc2qr_accret, qc2qr_autoconv, qv2qi_vapdep,
qc2qi_berg, qc2qr_ice_shed, qc2qi_collect, qr2qi_collect, qc2qi_hetero_freeze,
qr2qi_immers_freeze, qi2qr_melt`, plus sedimentation diagnostics). It is
tempting to wire a hook that reads each of those fields after P3 returns and
applies fractionation, completing the isotope-physics story in one go. That
approach is **physically wrong** and so it is deliberately deferred. The
reason: P3 calls `update_prognostic_ice` and `update_prognostic_liquid` inside
its level loop (`p3_main_impl_part2.hpp` around lines 410-450), and those
calls update `th_atm` — hence `T_atm` — within the same step. Sedimentation
is also adaptively substepped (`p3_main_impl.hpp:295, 303, 310`). The
registered tendency fields are macrostep aggregates that have already lost the
temperature trajectory across the substeps. Equilibrium and kinetic
fractionation factors depend on `T`, so applying fractionation to an
aggregated tendency at an aggregated temperature gives the wrong isotopic
delta. The correct hook placement is inside the per-substep tendency compute,
adjacent to where the bulk tendency is calculated and just before the
prognostic update. That is non-trivial design work — choosing where to
intercept, what state to pass to the hook (per-substep `T`, `p`, the local
tendency, the previous-state qv/qc/qi/qr/qm slices), how to handle the
ordering of multiple phase-change paths in the same level loop, and how to
keep the hook's runtime cost negligible when the default no-op is in force.
That work is scoped into a separate spec and explicitly excluded here. The
same reasoning applies to SHOC's vapor↔liquid path.

The science motivation behind this whole arc — the reason `WaterTypes`,
`WaterIsotopologue`, and an `OceanTracerFlux` stub already exist on this
branch — is water-isotope-enabled SCREAM (EAMxx-wiso). This spec, like the
qv-only one before it, is structural plumbing for that science work, not the
science itself.

## Approach

The decomposition is narrower than the prior spec because the patterns are
already established. Files in parentheses are likely edit sites; precise
diffs land in the planning checkpoint.

1. **Confirm prerequisite is in place.** The `prerequisite-prior-spec-merged`
   planning criterion checks the prior spec's deliverables landed: real
   `option(SCREAM_TRACE_WATER ... OFF)` block in
   `water_tracers/CMakeLists.txt`, `WTRC_MAX_CNST` driven by the chosen
   mechanism, `qv` allocated rank-3 in all builds, `WaterTracerHook`
   interface and default no-op present, qv passive-copy and ON-N=1-vs-OFF
   BFB tests passing.
2. **Apply the same registration pattern to qc, qi, qr, qm.** In
   `eamxx_p3_process_interface.cpp` (lines 79-82, the four
   `add_tracer<Updated>("qc"|"qr"|"qi"|"qm", ...)` calls) and any sibling
   registration sites (e.g., SHOC), wire each through whatever helper the
   prior spec produced — most likely a wrapper around `add_tracer` that
   allocates a rank-3 vector field with the `water_tracer` CMP dim. If the
   prior spec instead produced per-call edits, replicate them here for the
   four additional species.
3. **Reuse or extend the subview accessor.** If the prior spec's accessor is
   generic (takes any rank-3 water field and returns the CMP=0 slice), no
   new accessor is needed and consumer sites for qc, qi, qr, qm are updated
   to use it. If the accessor is qv-specific, this spec generalizes it. The
   bulk consumer sites are in `eamxx_p3_process_interface.cpp` and
   `eamxx_shoc_process_interface.cpp`; they should not gain any
   `#ifdef SCREAM_TRACE_WATER`.
4. **Combined passive-copy unit test.** Add
   `components/eamxx/src/physics/water_tracers/tests/water_mass_n2_passive_copy_test.cpp`
   (and register it in the local CMakeLists, gated on
   `SCREAM_TRACE_WATER=ON AND SCREAM_NUM_WATER_TRACERS GREATER_EQUAL 2`).
   The test allocates qv, qc, qi, qr, qm with N=2; copies CMP=0 to CMP=1
   for each; advances all five species through the qv/qc/qi/qr/qm-touching
   physics for several steps; and asserts species[k](:, 1, :) ==
   species[k](:, 0, :) for every species k at every output step. This test
   deliberately exercises P3 cross-species interactions (autoconversion,
   accretion, evaporation, deposition, sublimation, riming) — those are
   precisely the level-loop paths where the rank lift could go wrong.
5. **Condensates ON-vs-OFF BFB harness.** Add
   `water_tracers_condensates_trace_water_on_n1_bfb` ctest target — analogous
   to the prior spec's qv-only `water_tracers_qv_trace_water_on_n1_bfb` —
   that runs an identical minimal qc/qi/qr/qm-evolution case in both builds
   and bit-compares.
6. **Update the follow-up note.** In whichever file holds the deferred-
   species list (water_tracers.hpp comment block or follow-up-specs.md),
   strike qc, qi, qr, qm. Add a pointer to this spec. Add a new entry
   pointing at the to-be-written in-substep-hook spec, with a one-line
   explanation that links back to the temperature-trajectory argument.

Risks:

- **Hidden assumptions of scalar `qc/qi/qr/qm`.** Some kernels in P3, SHOC,
  or microphysics-adjacent code may take rank-2 device pointers or assume
  rank-2 layout for these species. Same risk pattern as the qv lift; the
  default build is the regression detector. Address as found.
- **Cross-species level-loop interactions.** The combined passive-copy test
  is deliberately stronger than five separate per-species tests because
  P3's level loop interleaves accretion, autoconversion, evaporation, and
  riming across species. A bug in the rank lift that only manifests when
  e.g. `qc → qr` interacts with `qi → qr` would not be caught by per-species
  tests. The combined test catches it.
- **Diagnostic outputs.** Some diagnostic fields (e.g., `eff_radius_qc`,
  `qc_sed`) reference qc/qi/qr/qm. Their layouts are unchanged by this spec
  (they are macrostep diagnostics, not the tracer fields themselves), but
  any internal copies that take a rank-2 view of qc into a diagnostic need
  the subview accessor. Audit during planning.
- **Discipline against scope creep.** The temptation to wire hooks at the
  many qc/qi/qr/qm phase-change sites in P3 is real. The
  `architectural-readiness-for-followup` integration criterion is the
  guardrail; if it is being relaxed to absorb hook wiring "since it's
  almost ready," stop and cut a new spec instead.

Relevant skills: **scientific-modeling-conventions** (rank-change discipline,
BFB requirements, ask-before for boundary-crossing edits),
**e3sm-platforms** (docker-local path mapping for builds),
**e3sm-build-and-test** (validation-tier discipline; this is Tier 0 — compile
+ unit tests, no coupled-model run).

## References

- `specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md` — prior spec; this
  spec is its strict downstream and assumes its deliverables are merged.
- `components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp:78-86` —
  the `add_tracer<Updated>` registrations for qv, qc, qr, qi, qm, nc, nr,
  ni, bm; lines 79-82 are the four mass-tracer registrations this spec
  modifies.
- `components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp:129-143` —
  the registered phase-change tendency fields (qr2qv_evap, qi2qv_sublim,
  qc2qr_accret, qc2qr_autoconv, qv2qi_vapdep, qc2qi_berg, qc2qr_ice_shed,
  qc2qi_collect, qr2qi_collect, qc2qi_hetero_freeze, qr2qi_immers_freeze,
  qi2qr_melt, qr_sed, qc_sed, qi_sed). Listed here only to document what
  is *deferred* from this spec.
- `components/eamxx/src/physics/p3/impl/p3_main_impl.hpp:295, 303, 310` —
  comments documenting adaptive substepping for cloud, rain, and ice
  sedimentation. Cited as evidence that P3 substeps internally.
- `components/eamxx/src/physics/p3/impl/p3_main_impl_part2.hpp` (around
  lines 410-450) — `update_prognostic_ice` and `update_prognostic_liquid`
  calls inside the level loop, updating `th_atm` mid-step. Cited as the
  load-bearing evidence for deferring hook wiring.
- `components/eamxx/src/physics/water_tracers/water_tracers.hpp:16, 58-60` —
  `WTRC_MAX_CNST` and the `wtrc_bulk_names` array (the latter still lists
  RAINQM/SNOWQM-style names from EAM that are not the source of truth for
  EAMxx P3).
- Project config: `.claude/project-config.yml` — standalone CMake build mode,
  ne4, docker-local container.
