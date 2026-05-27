---
spec_id: 2026-05-06-eamxx-water-tracer-isotopologue-registry
spec_type: model-e3sm
spec_version: 1
title: "EAMxx: water-tracer ↔ isotopologue assignment registry"
created: 2026-05-06T00:00:00-06:00
author: rfiorella
project: EAMXX-wiso

inputs:
  source_files:
    - components/eamxx/src/physics/water_tracers/water_tracers.hpp
    - components/eamxx/src/physics/water_tracers/water_isotopes.hpp
    - components/eamxx/src/physics/water_tracers/water_types.hpp
    - components/eamxx/src/physics/water_tracers/eamxx_water_tracers.cpp
    - components/eamxx/src/physics/water_tracers/CMakeLists.txt
    - specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md
    - specs/2026-05-06-eamxx-water-condensates-add-tracer-dim.md
  data: []
  baseline: null

deliverables:
  - "A compile-time tracer registry under namespace scream::WaterTracers that exposes, for each tracer index `i` in `[0, WTRC_MAX_CNST)`: a per-tracer name (string-view, distinguishes multiple tracers that share an isotopologue), the isotopologue catalog index (the existing 0–5 indices into the `WaterIsotopologues<Scalar>::*` constant arrays in water_isotopes.hpp), and an is-tag flag (true for tracers that follow bulk H2O physics but track source attribution; false for true isotopologues). The registry is a single source of truth — the loose arrays currently declared but unpopulated in water_tracers.hpp:49-77 (`wtrc_names`, `wtrc_species_names`, `wtrc_species`, `wtrc_is_tag`, etc.) are either populated from the registry or replaced by it."
  - "Build-time configuration that produces this registry. The output of configuration is a compile-time-known fixed-length array of `(name, catalog_index, is_tag)` tuples. The exact CMake input format (single delimited list, sibling CMake file with `add_water_tracer(...)` calls, or parallel lists) is the planning-phase open question gated by `mechanism-decision-resolved`."
  - "Configure-time validation in CMake that rejects: an empty configured tracer list when `SCREAM_TRACE_WATER=ON`; a list whose entry 0 does not name an `H2O` tracer; any entry naming an unknown isotopologue (i.e., not in `WaterIsotopologueNames` from water_isotopes.hpp:60); duplicate per-tracer names. The error must surface at configure time, not link or run time."
  - "A query API in the registry, callable from device kernels and host code: `tracer_isotopologue(i) → catalog_index`, `tracer_name(i) → string_view`, `tracer_is_tag(i) → bool`, plus a host-only convenience `find_tracer_by_name(name) → std::optional<int>`. All compile-time-resolvable for `if constexpr` branching in templated downstream code."
  - "Deprecation of `SCREAM_NUM_WATER_TRACERS` as a directly user-settable CMake variable. `WTRC_MAX_CNST` is now derived from the configured tracer list length. CMake emits a warning (or hard error — settled in planning) if `SCREAM_NUM_WATER_TRACERS` is set on the command line, with a pointer to the new mechanism."
  - "New unit test under components/eamxx/src/physics/water_tracers/tests/ — only registered when `SCREAM_TRACE_WATER=ON` — that builds with a configured list of 4 tracers `[bulk H2O, passive-tag H2O, HDO, H218O]` and asserts: `WTRC_MAX_CNST == 4`; tracer 0 is H2O (catalog idx 0) and `is_tag=false`; tracer 1 is H2O (catalog idx 0) and `is_tag=true`; tracer 2 is HDO (catalog idx 2); tracer 3 is H218O (catalog idx 3); per-tracer names are unique and round-trip through `find_tracer_by_name`."
  - "Updated follow-up note (water_tracers.hpp comment block or follow-up-specs.md) noting that fractionation physics implementation now has the metadata it needs (registry → isotopologue → fractionation factor via `AlphaEqIceVapor`/`AlphaEqLiquidVapor`/etc.) and is unblocked, but is itself a separate downstream spec."

success_criteria:
  - id: prerequisite-prior-specs-merged
    type: human-review
    description: "Confirm both prior specs are completed and merged: specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md (qv lifted to (COL, CMP, LEV) in all builds; SCREAM_TRACE_WATER is a real CMake option; WTRC_MAX_CNST is a build-time constant; subview accessor and WaterTracerHook interface in place; passive-copy and ON-vs-OFF BFB tests pass) and specs/2026-05-06-eamxx-water-condensates-add-tracer-dim.md (same lift applied to qc, qi, qr, qm; combined passive-copy test passes). This spec assumes that infrastructure and does not re-establish it."
    phase: planning

  - id: mechanism-decision-resolved
    type: human-review
    description: "Before writing implementation code, the user confirms (1) the CMake input format for the tracer list — (a) a single delimited string like `-DSCREAM_WATER_TRACERS=\"bulk:H2O;passive:H2O:tag;hdo:HDO;h218o:H218O\"`, (b) a CMake config file consumed via `-DSCREAM_WATER_TRACERS_FILE=...` containing a sequence of `add_water_tracer(NAME ... ISOTOPOLOGUE ... [TAG])` invocations, or (c) parallel CMake lists; (2) whether `SCREAM_NUM_WATER_TRACERS` becomes a hard error or a deprecation warning when set; and (3) the precise C++-level form of the per-tracer entry — a POD struct `WaterTracerInfo { string_view name; int catalog_idx; bool is_tag; }` vs. a tuple, etc. The catalog-index naming and the registry's read-only constexpr API are already settled (see resolved_decisions)."
    phase: planning

  - id: compile-clean-default-trace-water-off
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/default -DCMAKE_BUILD_TYPE=Debug && cmake --build build/default -j"
    expect: exit_zero
    phase: implementation

  - id: compile-clean-trace-water-on-bulk-only
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/twoBulk -DCMAKE_BUILD_TYPE=Debug -DSCREAM_TRACE_WATER=ON $(./scripts/water-tracers-config-flags bulk_only) && cmake --build build/twoBulk -j"
    expect: exit_zero
    phase: implementation

  - id: compile-clean-trace-water-on-registry-n4
    type: shell
    cmd: "cd components/eamxx && cmake -S . -B build/registry_n4 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_TRACE_WATER=ON $(./scripts/water-tracers-config-flags registry_n4) && cmake --build build/registry_n4 -j"
    expect: exit_zero
    phase: implementation

  - id: configure-time-rejects-bad-config
    type: shell
    cmd: "cd components/eamxx && ./scripts/water-tracers-config-validation-tests"
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

  - id: existing-tests-pass-trace-water-on-bulk-only
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/twoBulk --output-on-failure -E 'water_tracer_registry_n4'"
    expect: exit_zero
    phase: testing

  - id: bfb-bulk-only-vs-default
    type: tolerance
    cmd: "cd components/eamxx && ctest --test-dir build/twoBulk -R '^water_tracer_registry_bulk_vs_default_bfb$' --output-on-failure"
    metric: "max_abs_diff"
    rtol: 0
    atol: 0
    phase: testing

  - id: registry-queries-correct-n4
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/registry_n4 -R '^water_tracer_registry_n4$' --output-on-failure"
    expect: exit_zero
    phase: testing

  - id: architectural-readiness-for-followup
    type: human-review
    description: "Confirm the registry's API surface (tracer_isotopologue, tracer_name, tracer_is_tag, find_tracer_by_name) is sufficient for the next downstream spec — fractionation-physics implementation — to proceed without further changes to the registry itself. Specifically: a fractionation-aware WaterTracerHook implementation should be able to do `auto cat_idx = tracer_isotopologue(i); auto alpha = AlphaEqLiquidVapor(WaterIsotopologueNames[cat_idx], T);` and similar without needing additional metadata. If a missing query is identified, capture the gap before declaring this spec done."
    phase: integration

out_of_scope:
  - "Actual fractionation physics. The registry merely provides the metadata; turning a tracer index into a fractionation factor and applying it to a tendency is a separate downstream spec."
  - "A non-trivial WaterTracerHook implementation. The hook continues to default to no-op as established by the prior specs; this spec does not register an isotope-aware hook."
  - "Wiring WaterTracerHook calls to additional in-substep phase-change sites in P3 or SHOC. Still deferred to the in-substep-hook spec — physically blocked on the temperature-trajectory issue documented in specs/2026-05-06-eamxx-water-condensates-add-tracer-dim.md."
  - "Number concentration tracers (nc, nr, ni, bm). Still on the deferred list; their rank-lift is a separate follow-up. (Note: when that spec lands, the registry established here applies to those fields equally — the CMP→isotopologue mapping is the same.)"
  - "IO output of per-tracer field names. The registry now provides the per-tracer name as a string-view, but plumbing those names into NetCDF output dimension labels or per-slice metadata is deferred."
  - "Restart / checkpoint of the new tracer dimension or registry."
  - "Tag-tracer semantics beyond the `is_tag` flag itself. Source-attribution physics (e.g., how tag tracers are initialized at boundaries, how they evolve through the convection scheme) is a separate spec."
  - "Removal or repurposing of the WaterTypes::WaterType enum (Vapor, CloudLiquid, etc.) in water_types.hpp. That enum encodes phase, which the registry deliberately does *not*. They coexist; this spec does not modify water_types.hpp."
  - "Any changes under components/eam/ (legacy EAM Fortran) or other components/."
  - "Performance tuning for large registries (N > ~10). The registry is constexpr, so larger N imposes no per-tracer-call cost, but the build-time CMake parsing of the configured list may need attention for very large N — out of scope here."

resolved_decisions:
  - "This spec is strictly downstream of both specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md and specs/2026-05-06-eamxx-water-condensates-add-tracer-dim.md. Both must be completed and merged before this one starts. The prerequisite-prior-specs-merged success criterion is the explicit gate. The CMake option block for SCREAM_TRACE_WATER, the WTRC_MAX_CNST machinery, the (COL, CMP, LEV) layout for qv/qc/qi/qr/qm, the subview accessor, and the WaterTracerHook interface are inherited and not re-litigated."
  - "Per-tracer metadata is `(name, catalog_index, is_tag)` only. It is *not* phase-aware — the same tracer index in qv, qc, qi, qr, qm refers to the same isotopologue, just in a different phase. Phase is determined by which field is being accessed, not by the registry. The `wtrc_types` array in water_tracers.hpp:68 (which currently sits as unused scaffolding) is dropped or repurposed; it is not part of the registry."
  - "Two distinct index spaces are formalized and named: (a) **catalog index** (0–5): the existing `IsotopologueToIndex` from water_isotopes.hpp:61-68. References the constant tables `WaterIsotopologues<Scalar>::mwiso, rnat, alpal, ...`. Fixed at compile time of the codebase, never changes per build. (b) **tracer index** (0..WTRC_MAX_CNST−1): the CMP slot in the `(COL, CMP, LEV)` field layout. Build-dependent; settled at configure time from the user's tracer list. The registry maps tracer index → catalog index. Multiple tracer indices may map to the same catalog index (e.g., two tracers carrying H2O, one bulk, one tag)."
  - "The configured tracer list must include H2O (catalog index 0) at tracer index 0, and configure-time CMake validation enforces this. The reason: the prior specs already commit to tracer index 0 being bulk water, and downstream physics (the subview-at-CMP-0 accessor used by every bulk-water consumer) hard-codes that. Allowing tracer 0 to be anything else would silently break every existing consumer."
  - "Multiple tracer indices may share the same catalog index. The most common case is `[bulk H2O, tag H2O, HDO, H218O, ...]` where the tag H2O behaves like bulk physically but tracks source attribution. The registry distinguishes them by per-tracer name and the `is_tag` flag. Fractionation hooks (in a future spec) will short-circuit to no-op when `tracer_isotopologue(i) == 0` (catalog index for H2O) regardless of `is_tag`, because H2O does not fractionate against itself."
  - "`SCREAM_NUM_WATER_TRACERS`, the user-settable CMake variable introduced by the qv spec, is deprecated as a direct user input. `WTRC_MAX_CNST` is derived from the configured tracer list length. The exact deprecation behavior (warn vs. hard-error when set on the command line) is gated by mechanism-decision-resolved."
  - "The registry's query API is constexpr and callable from both host and device contexts: `tracer_isotopologue(i)`, `tracer_name(i)`, `tracer_is_tag(i)`, plus a host-only `find_tracer_by_name(name)`. The exact C++-level form of the per-entry struct (POD vs. tuple) is gated by mechanism-decision-resolved, but the read-only constexpr API and the field set are settled."
  - "All other deferred items from the prior specs (isotope physics, surface fluxes, IO/restart, EAM legacy code, in-substep hooks, number concentrations) remain deferred."

ask_before:
  - Modifying any file outside components/eamxx/
  - Touching components/eam/ (legacy Fortran EAM)
  - Modifying cime_config/ or any compset / coupler XML
  - Modifying water_isotopes.hpp's `WaterIsotopologueNames` or `IsotopologueToIndex` ordering — those are the catalog-index source of truth and reordering them is a breaking change for any callsite that hard-codes a catalog index. If a new isotopologue needs to be appended (legitimate growth), do so explicitly with a planning-checkpoint discussion.
  - Modifying water_types.hpp's `WaterType` enum (phase). Out of scope per resolved_decisions.
  - Modifying the SCREAM_TRACE_WATER CMake option block (inherited from prior specs)
  - Modifying the WaterTracerHook interface (inherited from prior specs)
  - Removing or relocating the loose `wtrc_*` arrays in water_tracers.hpp without first checking whether anything currently reads them — they are scaffolding but downstream code may consume them.
  - Adding a new namelist parameter
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

# EAMxx: water-tracer ↔ isotopologue assignment registry

## Context

The two prior specs in this arc — `specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md`
and `specs/2026-05-06-eamxx-water-condensates-add-tracer-dim.md` — produce
the data path: `qv, qc, qi, qr, qm` are all allocated as `(COL, CMP, LEV)` in
all builds, with `WTRC_MAX_CNST` slots along `CMP`, accessed via a single
subview-along-`CMP`-at-`0` helper for the bulk path and via a stable
`WaterTracerHook` interface for any future per-tracer logic. What is
*missing*, after those specs land, is **identity**: nothing in the codebase
yet knows what `qv(:, 2, :)` *is*. Future fractionation physics needs to
look up `qv(:, n, :) → "this slot is HDO"` to choose the correct equilibrium
fractionation factor (`AlphaEqLiquidVapor("HDO", T)` from
`water_isotopes.hpp:96-121`), and there is currently no compile-time table
that answers that question.

Two existing pieces of scaffolding inform the registry's shape. First,
`water_isotopes.hpp` already provides a six-isotopologue catalog: the
`WaterIsotopologueNames` array at line 60 (`{"H2O","H216O","HDO","H218O",
"H217O","HTO"}`), the parallel `IsotopologueToIndex` map at lines 61-68, and
the `WaterIsotopologues<Scalar>` constexpr struct (lines 19-58) holding
fractionation polynomial coefficients, molecular weights, VSMOW ratios, and
kinetic-fractionation parameters indexed by that catalog. The catalog is the
existing source of truth for "what does `HDO` mean physically." Second,
`water_tracers.hpp:49-77` declares (but does not yet populate) a parallel
set of arrays — `wtrc_names`, `wtrc_species_names`, `wtrc_species`,
`wtrc_is_tag`, and others — sized `WTRC_MAX_CNST`. These are clearly
intended to hold per-CMP-slot metadata, but they currently sit as empty
`{}` placeholders. This spec connects the two: it adds a build-time-driven
mechanism that populates a per-tracer-slot table whose entries reference
the catalog by index, and it exposes a constexpr query API on top of that
table.

A subtlety that came out during interview must be documented up front:
**multiple tracer slots can share an isotopologue**. The simplest case is a
"tag" tracer — a passive scalar that follows bulk H2O physics but tracks
some source attribution (e.g., evaporation from a particular ocean basin).
Tag tracers carry catalog index `0` (H2O) just like the actual bulk slot,
but they are distinct tracer slots, with distinct names, and the registry
must keep them separable. The existing `wtrc_is_tag` array already
anticipates this; this spec wires it up. The implication for the public API
is that `tracer_isotopologue(i)` returns a catalog index, and multiple `i`
values may return the same catalog index — downstream code needing to
distinguish them must consult `tracer_name(i)` and/or `tracer_is_tag(i)`.

Two index spaces are formalized so future spec text and code stay
unambiguous. The **catalog index** (0–5) is the existing
`IsotopologueToIndex` and addresses the constant arrays in
`water_isotopes.hpp` — it is fixed at compile time of the codebase. The
**tracer index** (0..`WTRC_MAX_CNST`−1) is the `CMP` slot in the field
layout — it is build-dependent, set by the configured tracer list. Catalog
indices are universal; tracer indices are per-build. The registry is the
mapping `tracer_index → catalog_index`, plus per-tracer name and is-tag
flag.

`SCREAM_NUM_WATER_TRACERS`, introduced by the qv spec as the user-facing
count knob, is deprecated by this spec. With a configured list of tracers
the count is derived (= list length); having a parallel knob would invite
inconsistency. The deprecation is either a configure-time warning or a
hard error — the choice is the planning-phase open question.

## Approach

This spec is structurally smaller than the two prior ones — there are no
field-allocation changes and no consumer-side edits in P3 or SHOC. The
work is concentrated in `components/eamxx/src/physics/water_tracers/` and
in CMake. Files in parentheses are likely edit sites; precise diffs land
in the planning checkpoint.

1. **Settle the CMake input format.** The planning-checkpoint
   `mechanism-decision-resolved` criterion picks one of: (a) a single
   delimited list `-DSCREAM_WATER_TRACERS="bulk:H2O;passive:H2O:tag;hdo:HDO;h218o:H218O"`,
   (b) a CMake config file with `add_water_tracer(NAME ... ISOTOPOLOGUE ... [TAG])`
   calls consumed via `-DSCREAM_WATER_TRACERS_FILE=...`, or (c) parallel
   CMake lists `-DSCREAM_WATER_TRACER_NAMES=... -DSCREAM_WATER_TRACER_ISOTOPES=... -DSCREAM_WATER_TRACER_TAGS=...`.
   I lean (b) for scalability and explicitness, but (a) wins on ergonomics
   for small configurations. The choice does not affect the C++ side.
2. **Configure-time validation in CMake.** Whichever input format is chosen,
   add validation that emits a hard error at `cmake -S . -B ...` time on
   any of: empty list when `SCREAM_TRACE_WATER=ON`, entry 0 not naming an
   `H2O` isotopologue, any entry naming an isotopologue not in
   `WaterIsotopologueNames`, duplicate per-tracer names. Errors must be
   actionable (point at the offending entry). This is the
   `configure-time-rejects-bad-config` success criterion's surface.
3. **Generate the C++ registry header from the configured list.**
   `water_tracers/CMakeLists.txt` writes a header — likely
   `water_tracer_registry.gen.hpp` — that defines the constexpr
   per-tracer-slot table. Mechanism (`configure_file`, `file(WRITE)`, or a
   small helper script) is implementation detail. The generated header is
   re-built when the configured list changes.
4. **Define the public registry API.** In a non-generated header
   (`water_tracer_registry.hpp` or extension of `water_tracers.hpp`),
   declare the constexpr query functions: `tracer_isotopologue(i)`,
   `tracer_name(i)`, `tracer_is_tag(i)`, `find_tracer_by_name(name)`.
   Implementations consult the generated table. All four are constexpr;
   the first three are usable in device code and `if constexpr` template
   contexts. `find_tracer_by_name` is host-only convenience.
5. **Wire up (or remove) the loose `wtrc_*` arrays.** The arrays declared
   at `water_tracers.hpp:49-77` are either populated from the registry (if
   any callsite already consumes them) or removed (if dead). A grep pass
   during planning settles which.
6. **Deprecate `SCREAM_NUM_WATER_TRACERS`.** When a user passes
   `-DSCREAM_NUM_WATER_TRACERS=N` on the command line, CMake either emits
   a deprecation warning (and proceeds, ignoring the value) or hard-errors
   with a pointer to the new mechanism. Decision gated by
   `mechanism-decision-resolved`. `WTRC_MAX_CNST` itself becomes a derived
   constexpr equal to the configured list length.
7. **Unit test exercising the n=4 registry.** Add
   `components/eamxx/src/physics/water_tracers/tests/water_tracer_registry_n4_test.cpp`,
   registered in the local `CMakeLists` only when `SCREAM_TRACE_WATER=ON`
   and the configured list has length 4. The test asserts: `WTRC_MAX_CNST
   == 4`; `tracer_isotopologue(0) == 0` (H2O catalog idx) and
   `tracer_is_tag(0) == false`; `tracer_isotopologue(1) == 0` and
   `tracer_is_tag(1) == true`; `tracer_isotopologue(2) == 2` (HDO);
   `tracer_isotopologue(3) == 3` (H218O); per-tracer names match the
   configured list and round-trip through `find_tracer_by_name`.
8. **CMake-validation regression tests.** A small driver script
   `scripts/water-tracers-config-validation-tests` (or similar) runs
   `cmake` repeatedly with deliberately malformed configurations (empty,
   missing H2O at index 0, unknown isotopologue, duplicate name) and
   asserts each fails at configure time with a clear error message. This
   is the `configure-time-rejects-bad-config` criterion.
9. **Bulk-only-vs-default BFB regression.** Reuse the prior specs' BFB
   pattern: with the configured list `[bulk H2O]` (length 1) and
   `SCREAM_TRACE_WATER=ON`, qv/qc/qi/qr/qm evolution must remain
   bit-for-bit identical to the `SCREAM_TRACE_WATER=OFF` default build.
   This is essentially the prior specs' ON-N=1-vs-OFF check, rewritten
   against the new registry-driven configuration. The
   `bfb-bulk-only-vs-default` criterion's ctest target wraps it.
10. **Update the follow-up note.** Strike "registry" off the deferred list.
    Add an explicit pointer to the next spec (fractionation physics
    implementation) and to the in-substep-hook spec — both of which
    consume the registry.

Risks:

- **Catalog ordering churn.** Future additions to `WaterIsotopologueNames`
  must append, never reorder. If anyone reorders the catalog, every
  build's tracer-index→catalog-index mapping silently changes meaning.
  Mitigation: a comment block on `WaterIsotopologueNames` in
  `water_isotopes.hpp` warning against reorder, plus the `ask_before` entry
  on this spec.
- **CMake parser fragility.** Whichever input format wins, the CMake-level
  parsing has corner cases (empty entries, whitespace, quoting). The
  validation tests in step 8 are the regression net.
- **`is_tag` semantics drift.** This spec wires up the `is_tag` flag
  without committing to any *physical* semantics for tag tracers (no
  source-attribution boundary conditions, no tag-specific advection
  logic). Future specs will define what a tag tracer actually does. The
  flag is provided so those specs do not have to retroactively introduce
  it. Mitigation: explicit `out_of_scope` entry stating tag-tracer
  semantics are deferred.
- **Multiple H2O tracers and BFB.** With the registry, `[bulk H2O,
  passive H2O tag]` is a valid 2-tracer configuration. The qv-spec and
  condensates-spec passive-copy tests, framed as "two H2O tracers must
  evolve identically through physics," are essentially this case. Verify
  that those tests still pass against the registry-driven configuration
  with no reinterpretation needed.
- **Constexpr-on-device.** The query functions must be callable from
  device kernels (constexpr + `KOKKOS_INLINE_FUNCTION` or equivalent
  annotation). Verify with the n=4 test that includes a device-side
  registry query.

Relevant skills: **scientific-modeling-conventions** (registry as
single-source-of-truth discipline, ask-before for catalog changes),
**e3sm-platforms** (docker-local builds), **e3sm-build-and-test**
(validation tier 0 — compile + unit tests, no coupled-model run).

## References

- `specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md` — qv rank lift;
  prerequisite.
- `specs/2026-05-06-eamxx-water-condensates-add-tracer-dim.md` — qc/qi/qr/qm
  rank lift; prerequisite.
- `components/eamxx/src/physics/water_tracers/water_isotopes.hpp:60` —
  `WaterIsotopologueNames` catalog (`{"H2O","H216O","HDO","H218O","H217O","HTO"}`).
- `components/eamxx/src/physics/water_tracers/water_isotopes.hpp:61-68` —
  `IsotopologueToIndex` name→catalog-index map.
- `components/eamxx/src/physics/water_tracers/water_isotopes.hpp:19-58` —
  `WaterIsotopologues<Scalar>` constexpr struct (mwiso, rnat, alpal, alpbl,
  alpcl, alpdl, alpel, alpai, alpbi, alpci, etc.) — the per-catalog-index
  fractionation parameters that downstream specs will read via
  `tracer_isotopologue(i)`.
- `components/eamxx/src/physics/water_tracers/water_isotopes.hpp:73-121` —
  `AlphaEqIceVapor`, `AlphaEqLiquidVapor` template functions; the form of
  query the registry's clients will perform.
- `components/eamxx/src/physics/water_tracers/water_tracers.hpp:49-77` —
  unpopulated `wtrc_*` arrays this spec replaces or fills.
- `components/eamxx/src/physics/water_tracers/water_types.hpp:22-44` —
  `WaterType` (phase) enum; orthogonal to and not modified by this spec.
- Project config: `.claude/project-config.yml` — standalone CMake build
  mode, ne4, docker-local container.
