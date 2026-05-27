---
spec_id: 2026-05-06-eamxx-water-vars-add-tracer-dim
spec_type: model-e3sm
spec_version: 1
title: "EAMxx: extend qv with a tracer (CMP) dimension"
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
    - components/eamxx/src/share/field/field_tag.hpp
    - components/eamxx/src/share/grid/abstract_grid.cpp
    - components/eamxx/src/share/grid/point_grid.cpp
    - components/eamxx/src/share/field/field_group.hpp
    - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
    - components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
  data: []
  baseline: null

deliverables:
  - "Convert the existing unconditional definition `target_compile_definitions(water_tracers PUBLIC SCREAM_TRACE_WATER)` (water_tracers/CMakeLists.txt:2) into a real CMake option `SCREAM_TRACE_WATER`, default OFF. The flag's job is *not* to gate qv's rank — qv is always rank-3 (see below). The flag instead gates (a) the upper bound on WTRC_MAX_CNST (must be 1 when OFF; user-settable via SCREAM_NUM_WATER_TRACERS when ON) and (b) the registration of a non-default `WaterTracerHook` implementation into the water_tracers library."
  - "Promoted WTRC_MAX_CNST (currently water_tracers.hpp:16) to a build-time constant. When SCREAM_TRACE_WATER=OFF, it is fixed at 1 (configure-time error if SCREAM_NUM_WATER_TRACERS != 1). When ON, it equals SCREAM_NUM_WATER_TRACERS."
  - "qv field is *always* allocated as a 3D vector layout (COL, CMP, LEV) with CMP dim named \"water_tracer\" and CMP size = WTRC_MAX_CNST, driven by the standard PointGrid::get_3d_vector_layout factory. There is no `#ifdef` around the registration. With SCREAM_TRACE_WATER=OFF the CMP dim has size 1 — the trivial case the FieldGroup machinery (field_group.hpp:34-37) is designed to handle. Pack/SIMD over LEV is unaffected; the unit CMP dim adds only an extra stride in the View."
  - "Subview-along-CMP-at-0 accessor compiled in *all* builds, returning a (COL, LEV) view of the bulk water from the rank-3 qv field. P3, SHOC, dynamics, and any other bulk-water consumer call this accessor and continue to operate on a (COL, LEV) view exactly as today. Consumer source files contain *no* `#ifdef SCREAM_TRACE_WATER`."
  - "WaterTracerHook abstraction defined in water_tracers and called unconditionally from the small set of sites that emit phase-change tendencies (e.g., inside P3 substep, inside SHOC's vapor↔liquid path). The default implementation is a zero-cost no-op compiled in all builds. Under SCREAM_TRACE_WATER=ON, the water_tracers library replaces the default registration with an isotope-aware implementation. Concrete form of the hook (functor / virtual base / function-pointer table) is gated by mechanism-decision-resolved in planning."
  - "New unit test under components/eamxx/src/physics/water_tracers/tests/ — only registered with the test suite when SCREAM_TRACE_WATER=ON and SCREAM_NUM_WATER_TRACERS >= 2 — that configures a passive duplicate of qv at CMP index 1, advances both through the existing qv-touching physics, and asserts tracer[1] == tracer[0] to machine epsilon."
  - "Code-level note (block comment in water_tracers.hpp or a follow-up-specs.md) listing the remaining bulk water species (qc, qi, qr, qs and number concentrations nc, ni, nr) that need the same treatment, plus a note that follow-ups will plug real isotope-aware logic into the WaterTracerHook interface established here."

success_criteria:
  - id: mechanism-decision-resolved
    type: human-review
    description: "Before writing implementation code, the user confirms two design choices: (1) the compile-time count mechanism — (a) a CMake preprocessor define SCREAM_NUM_WATER_TRACERS that resolves WTRC_MAX_CNST, (b) a class template parameter, or (c) a hybrid; and (2) the WaterTracerHook form — (a) a functor / std::function injected at process construction, (b) a virtual base class with a no-op default override, or (c) a free-function-pointer table. The hook interface lives in headers compiled in all builds; only the non-no-op implementation is gated by SCREAM_TRACE_WATER=ON."
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
    cmd: "cd components/eamxx && ctest --test-dir build/twon1 --output-on-failure -E 'water_tracers_qv_n2_passive_copy'"
    expect: exit_zero
    phase: testing

  - id: bfb-qv-trace-water-on-n1-vs-default
    type: tolerance
    cmd: "cd components/eamxx && ctest --test-dir build/twon1 -R '^water_tracers_qv_trace_water_on_n1_bfb$' --output-on-failure"
    metric: "max_abs_diff"
    rtol: 0
    atol: 0
    phase: testing

  - id: passive-copy-n2
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/twon2 -R '^water_tracers_qv_n2_passive_copy$' --output-on-failure"
    expect: exit_zero
    phase: testing

  - id: architectural-readiness-for-followup
    type: human-review
    description: "Confirm the (COL, CMP, LEV) plumbing established for qv (under SCREAM_TRACE_WATER) is generic enough that subsequent specs (qc, qi, qr, qs, and the number concentrations) can mechanically follow the same pattern without redesign, and that the SCREAM_TRACE_WATER=OFF default path remains the no-op fallback. If not, capture the gap before declaring this spec done."
    phase: integration

out_of_scope:
  - Isotope physics (equilibrium fractionation, kinetic fractionation, Rayleigh)
  - Surface fluxes / boundary conditions for water isotopes (ocean, land, sea ice)
  - IO writing of new tracer fields beyond what is needed for the in-test assertion (no NetCDF output of water_tracer dim in this spec)
  - Restart / checkpoint of the new tracer dimension
  - Extending qc, qi, qr, qs, nc, ni, nr to the new dimension (deferred to follow-up specs)
  - Any changes under components/eam/ (legacy EAM Fortran) or other components/
  - Performance tuning of the (COL, CMP, LEV) layout for large N

resolved_decisions:
  - "qv is *always* allocated as (COL, CMP, LEV) — there is no `#ifdef` on its rank. With SCREAM_TRACE_WATER=OFF the CMP dim has size 1 (the trivial case the FieldGroup machinery in field_group.hpp:34-37 was designed for); with SCREAM_TRACE_WATER=ON it has size SCREAM_NUM_WATER_TRACERS. This decision was made specifically so that physics processes that need isotope-aware hooks (P3, SHOC, dynamics — wherever phase changes happen, especially within subcycles) see a *single*, stable qv signature across all builds. Gating the rank by `#ifdef` would either force ifdef regions into those parameterizations or force a post-process-only fractionation pattern that is physically incorrect for sub-cycling schemes."
  - "SCREAM_TRACE_WATER's job is to gate (a) the upper bound on WTRC_MAX_CNST — 1 when OFF, user-set via SCREAM_NUM_WATER_TRACERS when ON — and (b) the *registration* of a non-default WaterTracerHook implementation. Process source files (P3, SHOC, dynamics) compile identically and call the hook unconditionally; the default no-op hook is compiled in all builds, so OFF builds incur no runtime cost beyond an empty function call (which the inliner removes)."
  - "Tracer dimension uses FieldTag::Component (CMP), giving 3D water-var layout (COL, CMP, LEV) and 2D layout (COL, CMP). This is dictated by PointGrid::get_3d_vector_layout (point_grid.cpp:93-104) and PointGrid::get_2d_vector_layout (point_grid.cpp:52-58); CMP is also the documented subview axis for monolithic groups (field_group.hpp:34-37, m_subview_dim is always 1). Putting the tracer first would fight LayoutRight + Pack vectorization (LEV must be innermost) and break per-column parallelism."
  - "User-facing dim string is 'water_tracer' (passed as vec_dim_name to the layout factory and surfaced in IO/diagnostics). FieldTag stays CMP; only the printed/IO name changes."
  - "Tracer index 0 is bulk water (the existing qv); higher indices are isotopologues / passive tracers added in future specs. The first non-zero index in this spec is a passive duplicate of qv used only by the new unit test."
  - "Tracer count N is a compile-time constant (per user direction). Exact mechanism (CMake preprocessor define vs. class template parameter) is the planning-phase open question gated by mechanism-decision-resolved. The constant is consumed in all builds (default OFF builds use N=1)."
  - "Bulk-water consumers (P3, SHOC, dynamics) access qv via a single subview-along-CMP-at-0 accessor that returns a (COL, LEV) view. The accessor is compiled in all builds and contains no `#ifdef`. The unit-CMP-stride layout (default builds) and the multi-CMP-stride layout (ON builds) produce identical (COL, LEV) views at index 0, so consumer code is bit-for-bit identical between flag settings."
  - "NetCDF output of qv in default builds will include the water_tracer dim of size 1 — i.e., the IO schema for qv changes for *all* EAMxx users, not just isotope users. This is the simplest implementation and avoids ON/OFF schema divergence. If the unit dim turns out to break downstream tooling, a follow-up spec adds an IO-layer unit-dim squeeze; that work is out of scope here."
  - "All other bulk water species (qc, qi, qr, qs, nc, ni, nr) are explicitly deferred to follow-up specs; this spec lands the qv plumbing and the WaterTracerHook interface only. Follow-ups will plug real isotope-aware logic into the same hook interface and apply the same always-rank-3 pattern to the other species."

ask_before:
  - Modifying any file outside components/eamxx/
  - Touching components/eam/ (legacy Fortran EAM)
  - Modifying cime_config/ or any compset / coupler XML
  - Changing the FieldTag enum or short-name aliases in field_tag.hpp
  - Modifying PointGrid::get_*_vector_layout signatures
  - Adding a new namelist parameter beyond a single integer for water-tracer count
  - Generating, replacing, or deleting any baseline file under /data/baselines/

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

# EAMxx: extend qv with a tracer (CMP) dimension

## Context

This is the first concrete plumbing step in templating water tracers in EAMxx
(SCREAM). Earlier work on this branch (commits a7ce65fc1f through b791d6ffb5)
laid down a `components/eamxx/src/physics/water_tracers/` subdirectory with
`water_tracers.hpp`, `water_types.hpp`, `water_isotopes.hpp`, and a CMakeLists,
plus a stub atmosphere process at `eamxx_water_tracers.cpp`. The scaffolding
defines a compile-time constant `WTRC_MAX_CNST` (currently `1`,
`water_tracers.hpp:16`) and a name array `wtrc_bulk_names` listing the seven
bulk water species (Q, CLDLIQ, CLDICE, RAINQM, SNOWQM, RAINQC, SNOWQC). What
that scaffolding does **not** yet do is change the rank of any water field
that the rest of EAMxx allocates and consumes: `qv` is still a scalar
`(COL, LEV)` field everywhere it appears in process interfaces.

This spec lands the smallest end-to-end change that makes the tracer-dim
plumbing real for one variable: extend `qv` from `(COL, LEV)` to
`(COL, CMP, LEV)`, with `CMP` being a compile-time count `N` of water
tracers (`WTRC_MAX_CNST`). `qv` takes this layout in *all* builds, including
the default `SCREAM_TRACE_WATER=OFF` build, where the CMP dim has size 1.
The reason for not gating the rank itself by `SCREAM_TRACE_WATER` is
forward-looking: subsequent specs will need to insert isotope-aware hooks
inside parameterizations that perform phase changes — P3, SHOC, dynamics —
and especially inside their subcycles, where state variables (T, p, q)
evolve mid-step and post-process fractionation cannot recover the correct
result. Those hooks need a single, stable `qv` signature across builds. A
rank that flips by `#ifdef` would force either `#ifdef` regions inside
those parameterizations or a post-process-only fractionation pattern that
is physically incorrect for subcycling schemes — neither is acceptable.
Always-rank-3 with a unit CMP dim in default builds is the trivial case the
FieldGroup machinery already supports (`field_group.hpp:34-37` documents
`m_subview_dim = 1` for exactly this layout family).

`SCREAM_TRACE_WATER` therefore changes role: it gates (a) the upper bound
on `WTRC_MAX_CNST` (1 when OFF, user-settable via
`SCREAM_NUM_WATER_TRACERS` when ON) and (b) the registration of a
non-default `WaterTracerHook` implementation that actually does isotope
work. The hook *interface* is compiled in all builds and called
unconditionally from the small set of phase-change emission sites; the
*default implementation* is a zero-cost no-op that the inliner removes.
The current line at `water_tracers/CMakeLists.txt:2`, which adds
`SCREAM_TRACE_WATER` as a public compile definition unconditionally, is
incorrect for this intent and is part of this spec's edit set: it becomes
a proper CMake `option(... OFF)` controlling only the upper bound on N
and the hook registration, never the qv rank.

With `SCREAM_TRACE_WATER=ON` and `N=1`, the new build must be bit-for-bit
identical to the default `SCREAM_TRACE_WATER=OFF` build on the same
workload. This is a stronger property than under a rank-gated design,
because both builds now use the same field shape and the same data path;
the *only* difference is whether the hook registration replaces the
default no-op. At `N=1` no isotope tracers exist, so the registered hook
must short-circuit to no-op as well — and the BFB check is what verifies
that short-circuit is working. With `SCREAM_TRACE_WATER=ON` and `N=2`,
tracer index 1 is set up as a passive duplicate of `qv`, and both indices
must remain identical to machine epsilon after every qv-touching physics
process. That second test proves the CMP dim is plumbed end-to-end through
the data path, not merely allocated.

The remaining bulk species (`qc, qi, qr, qs`, and the number concentrations
`nc, ni, nr`) need the same treatment, but are intentionally deferred to
follow-up specs. Once the qv path works and the architectural-readiness
checkpoint passes, those follow-ups should be mechanical applications of the
same pattern.

The wider goal — the reason `WaterTypes`, `WaterIsotopologue`, fractionation
factors, and an `OceanTracerFlux` stub already exist on this branch — is to
support water-isotope-enabled SCREAM (EAMxx-wiso). That science work depends
on this rank change being correct; this spec is the foundation, not the
science.

## Approach

The EAMxx field-allocation machinery already supports the desired layout
without modification:

- `FieldTag::Component` (`field_tag.hpp:36`, short alias `CMP` at line 54) is
  the established tag for a "thing this field is a vector of." There is no
  separate `Tracer` tag and there does not need to be — distinguishing
  *isotopes* from other component dims is done by the dimension's *string
  name*, not its `FieldTag`.
- `PointGrid::get_3d_vector_layout` (`point_grid.cpp:93-104`) returns
  `FieldLayout({COL, CMP, vtag}, {ncols, vector_dim, nlevs})` and lets the
  caller pass a string name for the CMP dim. That is exactly the layout we
  want for `qv`. The 2D analogue at lines 52-58 (`{COL, CMP}`) covers any
  surface-level companion fields that arise later.
- `FieldGroup` documentation (`field_group.hpp:34-37`) confirms that
  monolithic-group subviews are taken along `m_subview_dim = 1` — i.e., along
  `CMP`. Existing per-process consumers can therefore receive a
  subview-along-CMP-at-index-0 and continue to operate on a `(COL, LEV)`
  view, with no churn in their kernels.

The work decomposes as follows. Files in parentheses are the most likely
edit sites; precise diffs land in the planning checkpoint.

1. **Make `SCREAM_TRACE_WATER` a real option.** Replace the unconditional
   `target_compile_definitions(water_tracers PUBLIC SCREAM_TRACE_WATER)` at
   `water_tracers/CMakeLists.txt:2` with an `option(SCREAM_TRACE_WATER "..."
   OFF)` and a guarded `if (SCREAM_TRACE_WATER) target_compile_definitions(...
   PUBLIC SCREAM_TRACE_WATER) endif()` (or equivalent — exact form to be
   confirmed in planning). The flag must propagate from any caller who
   consumes the `water_tracers` interface library, so the public-scope
   propagation is preserved. Default OFF. CMake also rejects, at configure
   time, the combination `SCREAM_TRACE_WATER=OFF` with
   `SCREAM_NUM_WATER_TRACERS != 1`.
2. **Compile-time count.** Promote `WTRC_MAX_CNST` from a literal `1`
   (water_tracers.hpp:16) to a value driven by the chosen mechanism (CMake
   preprocessor define `SCREAM_NUM_WATER_TRACERS` or class template
   parameter — see `mechanism-decision-resolved`). The constant exists in
   *all* builds: it is `1` when OFF, and `SCREAM_NUM_WATER_TRACERS` when ON.
   It is *not* wrapped in `#ifdef SCREAM_TRACE_WATER`.
3. **qv allocation site (always rank-3).** In the upstream registration /
   field-manager wiring for `qv`, replace the existing scalar
   `get_3d_scalar_layout(LEV)` call with an unconditional
   `get_3d_vector_layout(LEV, WTRC_MAX_CNST, "water_tracer")`. There is no
   `#ifdef` and no `#else` branch. With `WTRC_MAX_CNST=1` (default OFF
   builds), the resulting layout is `(COL, CMP=1, LEV)` — a unit CMP dim,
   memory-equivalent to the prior scalar layout because LEV is innermost
   contiguous and a unit-stride extra outer dim is free at runtime. (The
   exact location of the qv registration is identified during planning and
   added to `inputs.source_files` at that time.)
4. **Subview accessor (unconditional).** Provide a helper that returns the
   bulk-water subview at `CMP=0` of any rank-3 water field. Used in *all*
   builds wherever `qv` is consumed — `eamxx_p3_process_interface.cpp`,
   `eamxx_shoc_process_interface.cpp`, dynamics. P3/SHOC source files
   contain no `#ifdef SCREAM_TRACE_WATER`. With `WTRC_MAX_CNST=1` the
   subview is byte-identical in layout to the prior scalar `qv` view; with
   `WTRC_MAX_CNST>1` it is the bulk-tracer slice with the same shape.
5. **`WaterTracerHook` interface (unconditional).** Define an abstract hook
   interface in the `water_tracers` library — exact form (functor /
   virtual base / function-pointer table) gated by
   `mechanism-decision-resolved`. The interface lives in a header compiled
   into all builds. Identify the small set of phase-change emission sites
   inside P3 (post-microphysics-step, post-condensation-substep) and SHOC
   (vapor↔liquid path) and call the hook unconditionally from each. The
   hook receives the rank-3 qv view, the phase-change tendency / state
   needed for fractionation, and the current substep T/p; the precise
   signature is captured during planning.
6. **Default no-op hook implementation (unconditional).** Compiled in all
   builds. Empty body. Verify (in the bfb check) that the inliner removes
   the call so the default-build runtime is identical to pre-spec.
7. **Hook registration gated by `SCREAM_TRACE_WATER`.** Under
   `SCREAM_TRACE_WATER=ON` the `water_tracers` library replaces the default
   no-op registration with a real `WaterTracerHook` implementation. For
   this spec's scope the "real" implementation still no-ops at `N=1`
   (there is no isotope work to do with only the bulk tracer); for `N=2`
   it implements the passive-copy logic the unit test exercises. Future
   specs swap in fractionation-aware implementations without touching the
   call sites.
8. **Passive-copy unit test (`ON, N=2` only).** Add
   `components/eamxx/src/physics/water_tracers/tests/qv_n2_passive_copy_test.cpp`
   (and register it in the local `CMakeLists.txt`, gated on
   `SCREAM_TRACE_WATER=ON AND SCREAM_NUM_WATER_TRACERS GREATER_EQUAL 2`).
   The test allocates `qv` with `N=2`, copies the initial state at
   `CMP=0` into `CMP=1`, runs the qv-touching physics for a small number
   of steps, and asserts `qv(:, 1, :) == qv(:, 0, :)` to machine epsilon at
   every step.
9. **`OFF`-vs-`ON, N=1` BFB harness.** Add a small driver — likely a
   ctest-driven shell harness in
   `components/eamxx/src/physics/water_tracers/tests/` — that runs an
   identical minimal qv-evolution case in both the `SCREAM_TRACE_WATER=OFF`
   build and the `SCREAM_TRACE_WATER=ON, N=1` build, captures qv at the
   same timesteps, and bit-compares. With always-rank-3 storage, both
   builds use the same FieldLayout, the same View shapes, and the same
   data path; the BFB check therefore validates that (a) the no-op hook
   default and the registered-but-no-isotope-work `N=1` hook are
   observationally identical, and (b) the inliner removes the empty
   default-hook call.
10. **Follow-up note.** Drop a short comment block in `water_tracers.hpp`
    (or a sibling `follow-up-specs.md`) listing `qc, qi, qr, qs, nc, ni,
    nr` as still needing the same lift, with a pointer back to this spec,
    and noting that follow-up specs plug fractionation-aware
    implementations into the `WaterTracerHook` interface established here.

Risks:

- **Hidden assumptions of scalar `qv`.** Some kernels may take the raw
  device pointer or assume rank-2. The subview accessor mitigates most
  cases; remaining ones surface as compile errors in the default build
  (since the rank changes there now too). This is the price of always-
  rank-3 — the default build is itself the regression detector. Address
  as found.
- **Hook-interface design.** Getting the `WaterTracerHook` signature wrong
  on the first pass means rework when isotope-aware follow-ups land.
  Mitigation: treat the *first* call site (the one this spec exercises in
  the unit test) as the design probe, and bias toward passing more state
  to the hook than seems strictly necessary for the passive-copy case —
  cheap to ignore an arg, expensive to add one later across N call sites.
- **IO output schema change for default builds.** With qv always
  `(COL, CMP, LEV)`, the NetCDF output gains a `water_tracer` dim of size
  1 in *all* builds, not just isotope ones. Downstream tools that read qv
  with rank-2 expectations may break. Resolved decision: accept the
  schema change and document it; if it bites in practice, a follow-up
  spec adds an IO-layer unit-dim squeeze. Out of scope here.
- **No-op-hook overhead in default builds.** Empty function calls at
  every phase-change site theoretically cost a function-call boundary.
  In practice the inliner removes them. Mitigate by marking the default
  implementation `inline` / providing it in the header and verifying with
  the bfb check that runtime is unchanged.
- **Pack alignment.** `LEV` is innermost specifically so packs vectorize
  over levels. The `(COL, CMP, LEV)` layout preserves that. A unit CMP
  dim adds one extra outer stride and does not break Pack alignment;
  verify with the bfb check.
- **Build-system mechanism.** Route (a) (CMake `-D` flag) couples N to
  configure-time and forces a clean rebuild on change. Route (b) (template
  parameter) avoids that but ripples a template parameter through every
  consumer. The planning checkpoint exists to settle this before either
  rabbit hole is committed to.

Relevant skills: **scientific-modeling-conventions** (rank-change discipline,
BFB requirements, ask-before for boundary-crossing edits),
**e3sm-platforms** (docker-local path mapping for builds),
**e3sm-build-and-test** (validation-tier discipline; this is Tier 0 — compile
+ unit tests, no coupled-model run).

## References

- `components/eamxx/src/share/field/field_tag.hpp:28-55` — `FieldTag` enum and
  `ShortFieldTagsNames` (`COL`, `CMP`, `LEV`, `ILEV`).
- `components/eamxx/src/share/grid/point_grid.cpp:52-58` and `93-104` —
  `get_2d_vector_layout` and `get_3d_vector_layout` factory definitions.
- `components/eamxx/src/share/field/field_group.hpp:25-37` — documentation of
  the `(COL, CMP, LEV)` group layout and the `m_subview_dim = 1` invariant.
- `components/eamxx/src/share/grid/abstract_grid.cpp:100-118` — convenience
  overloads that default the CMP dim name to `e2str(CMP)`; we override with
  `"water_tracer"`.
- `components/eamxx/src/physics/water_tracers/water_tracers.hpp:16` — current
  `WTRC_MAX_CNST = 1` literal that this spec promotes.
- `components/eamxx/src/physics/water_tracers/CMakeLists.txt:2` — current
  unconditional `target_compile_definitions(water_tracers PUBLIC
  SCREAM_TRACE_WATER)` that this spec converts into a real CMake option
  (default OFF).
- `components/eamxx/src/physics/water_tracers/water_types.hpp:22-44` —
  `WaterTypes` enum + `wtrc_bulk_names`/`wtrc_suffix` arrays informing the
  follow-up species list.
- Branch context: commits `1dd0f11ad2` (`add water_tracers subdir and
  SCREAM_TRACE_WATER`), `2471f72864` (`wip - definition of WaterIsotopologue
  and WaterTracer structs`), `e907f483d8` (`move Water* structs to header
  file`), `b791d6ffb5` (`add equilibrium fractionation unit tests`).
- Project config: `.claude/project-config.yml` — standalone CMake build mode,
  ne4, docker-local container.
