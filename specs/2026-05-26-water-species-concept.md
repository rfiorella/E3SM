---
spec_id: 2026-05-26-water-species-concept
spec_type: model-eamxx
spec_version: 1
title: "EAMxx: water-species enum + add_tracer helper in scream::water_tracers"
created: 2026-05-26T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

inputs:
  source_files:
    - components/eamxx/src/physics/water_tracers/water_tracers.hpp
    - components/eamxx/src/physics/water_tracers/water_isotopes.hpp
    - components/eamxx/src/physics/water_tracers/water_types.hpp
    - components/eamxx/src/physics/water_tracers/eamxx_water_tracers.cpp
    - components/eamxx/src/physics/water_tracers/CMakeLists.txt
    - components/eamxx/src/physics/CMakeLists.txt
    - components/eamxx/src/share/atm_process/atmosphere_process.hpp
    - components/eamxx/src/share/field/field_tag.hpp
    - components/eamxx/src/share/grid/point_grid.cpp
    - components/eamxx/src/share/field/field_group.hpp
    - references/old_port/water_tracers.hpp
    - references/old_port/water_isotopes.hpp
    - references/old_port/water_types.hpp
    - wiso_campaign_plan.md
    - campaigns/wiso.campaign.md
  data: []
  baseline: null

deliverables:
  - "components/eamxx/src/physics/water_tracers/water_species.hpp — header defining a compile-time water-species enum in namespace scream::water_tracers. Index 0 = bulk H2O (catalog index 0 in water_isotopes.hpp WaterIsotopologueNames). Subsequent indices reserved for isotopologues (H216O, HDO, H218O, H217O, HTO) per WaterIsotopologueNames ordering. Header is constexpr / KOKKOS_INLINE_FUNCTION-callable. Provides species_count(), species_name(int), species_catalog_index(int), and an is_bulk(int) predicate."
  - "components/eamxx/src/physics/water_tracers/water_tracer_registration.hpp — header declaring add_tracer_multi(name, grid, units, pack_size, tracer_advection) free function in namespace scream::water_tracers. Function returns the FieldRequest a process would submit via the existing add_field<RT> path; the rank-3 (COL, CMP, LEV) layout uses PointGrid::get_3d_vector_layout with CMP dim name \"water_tracer\" and CMP size = WTRC_MAX_CNST. The returned FieldRequest carries FieldLayout{COL, CMP, LEV} in that tag order; the underlying Kokkos View is LayoutRight (the EAMxx default), matching the MAM4 (COL, nmodes, LEV) precedent in components/eamxx/src/physics/mam/eamxx_mam_microphysics_process_interface.cpp:151. The companion .cpp implements it without touching atmosphere_process.hpp."
  - "Convert components/eamxx/src/physics/water_tracers/CMakeLists.txt line 2 from `target_compile_definitions(water_tracers PUBLIC SCREAM_TRACE_WATER)` (unconditional, defective) into a real option block: `option(SCREAM_TRACE_WATER \"Enable water tracer / isotope machinery\" OFF)` driven at the parent CMakeLists, plus an integer cache var `SCREAM_NUM_WATER_TRACERS` (default 1). When OFF, force `SCREAM_NUM_WATER_TRACERS=1`; when ON, accept user value ≥ 1. Both propagate via `target_compile_definitions(water_tracers PUBLIC SCREAM_TRACE_WATER=1 WTRC_MAX_CNST=${SCREAM_NUM_WATER_TRACERS})`. Configure-time error if SCREAM_TRACE_WATER=OFF and SCREAM_NUM_WATER_TRACERS≠1."
  - "Update components/eamxx/src/physics/water_tracers/water_tracers.hpp: replace literal `constexpr int WTRC_MAX_CNST = 1` (water_tracers.hpp:16) with `#ifndef WTRC_MAX_CNST #error \"WTRC_MAX_CNST must be defined by CMake\" #endif`, matching the convention in references/old_port/water_tracers.hpp:19-21. Add bulk-water subview accessor template (`get_bulk_water_subview`) following references/old_port/water_tracers.hpp:113-117. No other behavior changes."
  - "components/eamxx/src/physics/water_tracers/tests/water_species_test.cpp — Catch2 unit test asserting: species_count() ≥ 1; species_name(0) == \"H2O\"; species_catalog_index(0) == 0; is_bulk(0) == true; for N>1 builds (gated SCREAM_TRACE_WATER=ON, SCREAM_NUM_WATER_TRACERS>=2) is_bulk(1) == false. Tests for add_tracer_multi: given a mock PointGrid with ncols=4, nlevs=8, the returned FieldRequest's FieldIdentifier carries layout (COL, CMP, LEV) with extents (4, WTRC_MAX_CNST, 8) and CMP dim name == \"water_tracer\"."
  - "Local CMakeLists.txt entry in components/eamxx/src/physics/water_tracers/tests/ registering the new water_species_test, only when SCREAM_TRACE_WATER=ON (test exercises the multi-tracer path which is meaningful only in that build mode)."
  - "components/eamxx/src/physics/water_tracers/water_tracer_access.hpp — device-callable accessors that codify the layout decisions for every downstream PR. Three accessors, all KOKKOS_INLINE_FUNCTION and device-safe (no std::string, no std::array in device code): (a) get_species_subview(field, i, s) returns the (LEV/pack,) uview_1d<Pack> view for column i, species s — at N=1 this is shape- and stride-identical to the current scalar (COL, LEV) prognostic subview, which is what makes slice-0 BFB achievable in PR 2; (b) get_bulk_water_subview(field, i) is an alias for get_species_subview(field, i, 0); (c) get_all_species_at_level(field, i, k_pack) returns a value-type Pack[WTRC_MAX_CNST] (or equivalent fixed-size struct) holding all species at fixed (col, level-pack), used by fractionation kernels — value type, NOT a Kokkos subview, so stack-allocated and aliasing-free. The header includes an `if constexpr (WTRC_MAX_CNST == 1)` specialization that reduces (a) and (b) to a direct scalar reference."
  - "components/eamxx/src/physics/water_tracers/tests/water_tracer_access_test.cpp — Catch2 unit test asserting: (a) get_species_subview on a (ncol=4, ncmp=1, nlev=8) field returns a uview_1d<Pack> whose data() pointer equals the data() pointer of the equivalent (ncol=4, nlev=8) scalar field at column i (proves slice-0 stride identity at N=1); (b) on a (ncol=4, ncmp=4, nlev=8) field, get_species_subview(field, i, s) for s in [0, N) returns views whose data() pointers differ by exactly nlev_pack * sizeof(Pack) (proves LayoutRight CMP stride); (c) get_all_species_at_level returns Pack[WTRC_MAX_CNST] holding the correct values for a known test field."

success_criteria:
  - id: configure-default-off
    type: shell
    cmd: "cd components/eamxx && rm -rf build/wiso-default && cmake -S . -B build/wiso-default -DCMAKE_BUILD_TYPE=Debug 2>&1 | tee build/wiso-default/configure.log && grep -q 'SCREAM_TRACE_WATER.*OFF' build/wiso-default/configure.log"
    expect: exit_zero
    phase: implementation

  - id: configure-on-rejects-n0
    type: shell
    cmd: "cd components/eamxx && rm -rf build/wiso-bad && ! cmake -S . -B build/wiso-bad -DCMAKE_BUILD_TYPE=Debug -DSCREAM_TRACE_WATER=ON -DSCREAM_NUM_WATER_TRACERS=0 2>build/wiso-bad-err.log; grep -q 'SCREAM_NUM_WATER_TRACERS' build/wiso-bad-err.log"
    expect: exit_zero
    phase: implementation

  - id: configure-off-rejects-n2
    type: shell
    cmd: "cd components/eamxx && rm -rf build/wiso-off-n2 && ! cmake -S . -B build/wiso-off-n2 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_TRACE_WATER=OFF -DSCREAM_NUM_WATER_TRACERS=2 2>build/wiso-off-n2-err.log; grep -q 'SCREAM_NUM_WATER_TRACERS' build/wiso-off-n2-err.log"
    expect: exit_zero
    phase: implementation

  - id: compile-default-off
    type: shell
    cmd: "cd components/eamxx && cmake --build build/wiso-default -j --target water_tracers"
    expect: exit_zero
    phase: implementation

  - id: compile-on-n1
    type: shell
    cmd: "cd components/eamxx && rm -rf build/wiso-on-n1 && cmake -S . -B build/wiso-on-n1 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_TRACE_WATER=ON -DSCREAM_NUM_WATER_TRACERS=1 && cmake --build build/wiso-on-n1 -j --target water_tracers"
    expect: exit_zero
    phase: implementation

  - id: compile-on-n4
    type: shell
    cmd: "cd components/eamxx && rm -rf build/wiso-on-n4 && cmake -S . -B build/wiso-on-n4 -DCMAKE_BUILD_TYPE=Debug -DSCREAM_TRACE_WATER=ON -DSCREAM_NUM_WATER_TRACERS=4 && cmake --build build/wiso-on-n4 -j --target water_tracers"
    expect: exit_zero
    phase: implementation

  - id: clang-format-check
    type: shell
    cmd: "git diff --name-only $(git merge-base HEAD wiso-dev)..HEAD -- 'components/eamxx/**/*.hpp' 'components/eamxx/**/*.cpp' | xargs -r clang-format --dry-run --Werror"
    expect: exit_zero
    gate: advisory
    on_fail: skip
    phase: implementation
    verifies:
      - deliverable: components/eamxx/src/physics/water_tracers/water_species.hpp

  - id: unit-test-water-species-n4
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/wiso-on-n4 -R '^water_species_test$' --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: components/eamxx/src/physics/water_tracers/tests/water_species_test.cpp
      - claim: "species_count() == WTRC_MAX_CNST and add_tracer_multi produces (COL, CMP=N, LEV) layout"

  - id: existing-tests-pass-default
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/wiso-default --output-on-failure -E 'water_species_test'"
    expect: exit_zero
    phase: testing

  - id: layout-stride-test-n1
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/wiso-on-n1 -R '^water_tracer_access_test$' --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: components/eamxx/src/physics/water_tracers/water_tracer_access.hpp
      - claim: "Memory layout is LayoutRight matching MAM4 (COL, nmodes, LEV) precedent"
      - claim: "Pack dim stays over LEV; species subview is uview_1d<Pack> of shape (nlev_pack,)"
      - claim: "N=1 species-0 subview is data-pointer-identical to the equivalent (COL, LEV) scalar field subview"

  - id: layout-stride-test-n4
    type: shell
    cmd: "cd components/eamxx && ctest --test-dir build/wiso-on-n4 -R '^water_tracer_access_test$' --output-on-failure"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: components/eamxx/src/physics/water_tracers/tests/water_tracer_access_test.cpp
      - claim: "LayoutRight CMP stride equals nlev_pack * sizeof(Pack)"
      - claim: "Cross-species accessor get_all_species_at_level returns Pack[WTRC_MAX_CNST] value type, not a Kokkos subview"

  - id: scalar-codegen-equivalence
    type: shell
    cmd: "cd components/eamxx && bash scripts/wiso/check_n1_codegen.sh build/wiso-default build/wiso-on-n1"
    expect: exit_zero
    gate: advisory
    on_fail: skip
    phase: testing
    description: "Compares object-file disassembly of a small reference kernel compiled against the scalar (COL, LEV) prognostic vs the N=1 (COL, CMP=1, LEV) prognostic. Confirms the `if constexpr (WTRC_MAX_CNST == 1)` specialization collapses to the same instructions. Advisory because compilers may legitimately reorder; the binding contract is the slice-0-BFB run criterion in PR 2."
    verifies:
      - claim: "N=1 if-constexpr specialization collapses species loop to scalar code path"
      - claim: "FMA reordering at N=1 does not produce 1-ULP differences vs the pre-extension scalar code"
      - claim: "Slice-0 BFB at rtol=0/atol=0 is structurally achievable in PR 2"

  - id: arch-readiness-fractionation
    type: human-review
    description: "Confirm get_all_species_at_level returns a value type (Pack[WTRC_MAX_CNST] or equivalent fixed-size struct, NOT a Kokkos subview) so that fractionation kernels can stack-allocate species state without aliasing the underlying Field memory. Confirm the `if constexpr (WTRC_MAX_CNST == 1)` specialization for the species-subview accessors reduces the species loop to a single direct reference (verified by reading the generated header or by compiler explorer)."
    phase: integration
    verifies:
      - deliverable: components/eamxx/src/physics/water_tracers/water_tracer_access.hpp
      - claim: "Cross-species access API is bound and downstream-stable for PR 6-9 fractionation kernels"

  - id: architectural-readiness-for-pr2
    type: human-review
    description: "Confirm the add_tracer_multi helper and the (COL, CMP, LEV) layout path are generic enough that PR2 (extend qv) can wire P3 and SHOC's qv add_tracer<Updated> calls through them with no further changes to water_tracers headers. Specifically: an interface process can call scream::water_tracers::add_tracer_multi(\"qv\", grid, kg/kg, ps) and forward the returned FieldRequest through its existing add_field<RT> path. Additionally: PR 2 can pass the uview_1d<Pack> returned by get_species_subview directly into existing P3/SHOC kernels that accept uview_1d<Pack>, with no change to those kernels' signatures at N=1."
    phase: integration
    verifies:
      - deliverable: components/eamxx/src/physics/water_tracers/water_tracer_registration.hpp

out_of_scope:
  - Any change to qv, qc, qi, qr, qm allocation in P3 or SHOC (PR 2-4 territory).
  - Any change to atmosphere_process.hpp's add_tracer template (specs/2026-05-26-extend-qv-multi-tracer.md decides whether to modify it or layer above it).
  - WaterTracerHook abstraction, isotope physics, fractionation hooks, surface fluxes.
  - IO schema changes; restart/checkpoint handling.
  - Number-concentration tracers (nc, nr, ni, bm).
  - Modification of WaterIsotopologueNames or IsotopologueToIndex ordering in water_isotopes.hpp.
  - Any change under components/eam/ (legacy Fortran EAM) or other components/.
  - Performance tuning for large N.
  - Workspace policy for 2D (CMP × LEV) scratch buffers. PRs 2-4 may need species-aware scratch but take N separate 1D buffers via the existing WorkspaceManager API; a true 2D WorkspaceManager extension is out of campaign scope.
  - Changing the pack dim from LEV to CMP. Pack-over-CMP is rejected by this spec because it would force nk_pack=1 and destroy level-vector parallelism in every existing physics kernel; any future revisit is its own spec.

resolved_decisions:
  - "Namespace is scream::water_tracers (lowercase), per wiso_campaign_plan.md design principles. The existing scream::WaterTracers (PascalCase) namespace in water_tracers.hpp is retained for the prior arrays but new code lives in scream::water_tracers."
  - "Bulk H2O lives at species index 0 (C++ 0-based). This is the slice-1-unchanged invariant from tracer-conservation skill: PR2-4 must keep the species-0 slice bit-for-bit identical to the pre-extension scalar fields."
  - "Two index spaces are formalized: catalog index (0–5, fixed; IsotopologueToIndex in water_isotopes.hpp:60-67) vs. tracer/species index (0..WTRC_MAX_CNST−1, build-dependent CMP slot). This spec defines the tracer-index ↔ catalog-index mapping as identity-with-truncation: species index i maps to catalog index i for i < min(WTRC_MAX_CNST, 6). A future registry spec (post-campaign) may generalize this."
  - "SCREAM_TRACE_WATER is a real CMake option, default OFF. SCREAM_NUM_WATER_TRACERS is an integer cache var, default 1, forced to 1 when SCREAM_TRACE_WATER=OFF. WTRC_MAX_CNST is the C++ preprocessor define carrying that value into headers."
  - "add_tracer_multi returns a FieldRequest (or equivalent struct) rather than directly registering with the field manager. The caller (P3, SHOC) still owns the add_field<RT> invocation. This keeps water_tracers free of atmosphere_process.hpp dependencies and avoids re-litigating the AtmosphereProcess base-class API."
  - "Per regression-baseline skill: this spec is Tier 0 — no model run, no baseline comparison. PR 2 is the first spec carrying a slice-0 BFB criterion."
  - "Per eamxx-cpp-conventions: C++17, Kokkos device-callable functions tagged KOKKOS_INLINE_FUNCTION, no std::string in device code (species_name returns const char* or std::string_view)."
  - "Memory layout: LayoutRight (the EAMxx default). With FieldLayout tag order (COL, CMP, LEV), this gives unit-stride access along LEV, stride nlev along CMP, and stride ncmp*nlev along COL. Matches the MAM4 precedent of (COL, nmodes, LEV) used in components/eamxx/src/physics/mam/eamxx_mam_microphysics_process_interface.cpp:151. Cross-species access at fixed (col, lev) is non-unit-stride but remains L1-cache-resident for N ≤ 6."
  - "Pack dim stays over LEV. Pack-over-CMP is rejected: it would force nk_pack=1, destroy level-vector parallelism in every existing physics kernel, and require CMP ≤ SCREAM_PACK_SIZE. Consequence: a species slice subview(field, i, s) returns a uview_1d<Pack> of shape (nlev_pack,), identical in shape and stride to the current scalar (COL, LEV) prognostic subview. This identity is what makes slice-0 BFB achievable at PR 2."
  - "CMP iteration pattern in downstream kernels: serial-inside-thread by default. The species loop sits OUTSIDE the existing TeamVectorRange over levels. Sketch: `for (int s = 0; s < WTRC_MAX_CNST; ++s) { auto q_s = get_species_subview(qv_field, i, s); Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nk_pack), [&] (int k) { /* bulk-water physics on q_s(k) */ }); }`. This preserves level vectorization (the dominant axis on CPU pack lanes) and matches the slice-0 BFB requirement by construction at N=1. TeamThreadRange-over-CMP is allowed but not required; specs adopting it must document why (e.g., when CMP-parallel work fits a GPU wavefront better than LEV-parallel work for that kernel)."
  - "Cross-species access (for fractionation kernels) is provided via get_all_species_at_level, returning a value-type Pack[WTRC_MAX_CNST]. Non-unit-stride loads under the chosen LayoutRight (stride = nlev_pack * sizeof(Pack)); the cost is borne by the fractionation kernel where N ≤ 6 reads + N×N alpha ops + N writes fit in L1. Bulk microphysics never uses this access path."
  - "N=1 build path uses `if constexpr (WTRC_MAX_CNST == 1)` template specialization to collapse the species loop to a direct scalar reference to the s=0 slice. With pack-over-LEV stride identical to the current (COL, LEV) prognostic layout, the generated machine code at N=1 is BFB with the pre-extension scalar path. PR 2's slice-0 BFB criterion depends on this specialization; PR 1 ships the helper that downstream PRs use."
  - "Pack workspace policy: PR 1 does not change WorkspaceManager. PRs that need species-aware scratch buffers take N separate uview_1d<Pack> buffers via the existing take_many_and_reset API rather than introducing 2D (CMP × LEV) scratch. Rationale: the only physics scratch that actually needs species is the in-cloud condensate (qc_incld, qr_incld, qi_incld, qm_incld in p3_main_impl.hpp:131); the other ~39 P3 scratch buffers are bulk diagnostics. PRs 2-4 audit which scratch lifts and which stays bulk; PR 1 binds the policy."

ask_before:
  - Modifying any file outside components/eamxx/src/physics/water_tracers/ and its parent CMakeLists.
  - Modifying atmosphere_process.hpp add_tracer template (defer to PR 2 planning).
  - Modifying FieldTag enum or short-name aliases in field_tag.hpp.
  - Modifying PointGrid::get_*_vector_layout signatures.
  - Modifying WaterIsotopologueNames or IsotopologueToIndex ordering (catalog reorder is a breaking change for every downstream consumer).
  - Adding a new namelist parameter beyond the integer SCREAM_NUM_WATER_TRACERS CMake var.
  - Touching components/eam/ (legacy Fortran EAM).
  - Generating, replacing, or deleting any baseline file under /data/baselines/.
  - Changing the pack dim from LEV to CMP, or introducing a Kokkos memory layout other than LayoutRight, for any water-tracer field. These are binding architectural decisions of this spec and revisiting them invalidates the slice-0 BFB strategy in PR 2.
  - Replacing the if-constexpr N=1 specialization with a runtime branch. The compile-time specialization is what guarantees PR 2's rtol=0/atol=0 slice-0 BFB criterion.

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
  blocking: true                              # foundational spec — pre-merge review gate

model_specific:
  validation_tier: 0
  target_compset: null
  target_resolution: ne4pg2_ne4pg2
  stop_n: 0
  stop_option: ndays
  platform: docker-local
  container_image: rfiorella/model-containers:e3sm-openmpi-dev-latest
  build_mode: eamxx-standalone-cmake
  test_driver: components/eamxx/scripts/test-all-eamxx
  # model-eamxx-specific fields (per templates/spec-model-eamxx.md)
  tracer_count: 1                             # PR 1 builds at N=1 and N=4 but does not run; baseline carries N=1
  dp_proxy: false                             # no DP run; Tier-0 build-only spec
  cmake_option: SCREAM_TRACE_WATER            # this spec adds the option
  # Layout block populated for structural refactors (foundational spec contract)
  layout:
    kokkos_layout: LayoutRight                # EAMxx default; matches MAM4 (COL, nmodes, LEV) precedent
    pack_axis: LEV                            # pack stays over LEV; pack-over-CMP rejected (see resolved_decisions)
    cmp_dim_name: water_tracer
    n1_specialization: if-constexpr           # collapses species loop at WTRC_MAX_CNST=1 to scalar code path
---

# EAMxx: water-species enum + add_tracer helper in scream::water_tracers

## Context

This spec is PR 1 of the 4-PR group-1 chain in the EAMxx-wiso campaign
(`campaigns/wiso.campaign.md`). Group 1 lifts every bulk water mass
field — `qv`, `qc`, `qi`, `qr`, `qm` — from `(COL, LEV)` to
`(COL, CMP, LEV)` with the CMP dim carrying water-tracer species. This
PR lands the foundation: the species concept, the CMake plumbing that
sizes the CMP dim, the registration helper that callers will use, and
the bulk-slice subview accessor. No production field allocation
changes ship in this PR.

The reason for splitting the foundation out is mechanical: PRs 2-4
each touch one species class (vapor, cloud condensates, precip) across
P3 and SHOC, and each must demonstrate slice-0 BFB against the
pre-campaign baseline (the `slice-1-unchanged` invariant from the
`tracer-conservation` skill). Bundling the foundation with any of the
species-class lifts would conflate "did the helper compile" with "did
the rank-lift preserve bulk-water physics." Splitting also keeps PR 1
header-only + CMake (Tier 0, no model run), so the campaign begins
with a low-risk, fast-iterating PR rather than the higher-risk Tier-2
slice-0 BFB work that follows.

Two existing pieces of scaffolding inform the shape. First,
`water_isotopes.hpp` already defines a six-isotopologue catalog
(`WaterIsotopologueNames` line 60: `{H2O, H216O, HDO, H218O, H217O,
HTO}`) with constexpr fractionation coefficient tables. Second, the
prior `references/old_port/water_tracers.hpp` shows a working pattern
of `WTRC_MAX_CNST` as a CMake-define and a `get_bulk_water_subview`
template that returns the CMP=0 slice. Both patterns are inherited
here rather than redesigned.

## Approach

The work decomposes into four small pieces. Files in parentheses are
likely edit sites; precise diffs land in the planning checkpoint.

1. **CMake option + count.** Convert the unconditional
   `target_compile_definitions(... SCREAM_TRACE_WATER)` at
   `water_tracers/CMakeLists.txt:2` into a real
   `option(SCREAM_TRACE_WATER ... OFF)` declared in
   `components/eamxx/src/physics/CMakeLists.txt` (or higher). Add a
   cache var `SCREAM_NUM_WATER_TRACERS` (int, default 1). Configure-
   time validation: OFF + N≠1 is an error; ON + N<1 is an error.
   Both flags propagate to the `water_tracers` interface library via
   `target_compile_definitions(... PUBLIC SCREAM_TRACE_WATER=1
   WTRC_MAX_CNST=${SCREAM_NUM_WATER_TRACERS})`.

2. **Species enum + query API** (`water_species.hpp`). Define a tiny
   constexpr API in `namespace scream::water_tracers`:
   `species_count()`, `species_name(int)`, `species_catalog_index(int)`,
   `is_bulk(int)`. All `KOKKOS_INLINE_FUNCTION` and device-callable.
   Mapping is identity-with-truncation against
   `WaterIsotopologueNames`: species index `i` maps to catalog index
   `i` for `i < min(WTRC_MAX_CNST, 6)`. The map is settled at
   compile time, no per-build configuration.

3. **add_tracer_multi helper** (`water_tracer_registration.hpp` +
   `.cpp`). Free function in `namespace scream::water_tracers` that
   returns a `FieldRequest` carrying the rank-3 vector layout
   `(COL, CMP=WTRC_MAX_CNST, LEV)` with CMP dim name `"water_tracer"`,
   units, pack size, and the tracer-group memberships that
   `AtmosphereProcess::add_tracer` already puts a tracer into. Callers
   (PR 2's P3 / SHOC edits) feed the returned request through their
   existing `add_field<RT>` path. The helper does not depend on
   `atmosphere_process.hpp`; it depends only on grid + field-manager
   headers.

4. **Bulk-slice subview accessor + WTRC_MAX_CNST guard** (edits to
   `water_tracers.hpp`). Replace the literal `constexpr int
   WTRC_MAX_CNST = 1` at `water_tracers.hpp:16` with the CMake-driven
   `#ifndef WTRC_MAX_CNST #error ... #endif` block from
   `references/old_port/water_tracers.hpp:19-21`. Add a
   `get_bulk_water_subview` template matching
   `references/old_port/water_tracers.hpp:113-117` — returns the
   `(COL, LEV)` slice at CMP=0 of any rank-3 water field.

5. **Unit test** (`tests/water_species_test.cpp`). Catch2 test. The
   default-OFF build does not register the test (WTRC_MAX_CNST=1
   makes the multi-tracer assertions trivial). The ON, N=4 build
   asserts: enum entries 0-3 map to H2O, H216O, HDO, H218O;
   `is_bulk(0) == true`, `is_bulk(1) == false`; `add_tracer_multi`
   against a 4×8 PointGrid produces a FieldIdentifier whose layout
   has extents `(4, 4, 8)` and whose CMP dim name is `water_tracer`.

6. **Layout / pack / codegen accessors and test**
   (`water_tracer_access.hpp` + `tests/water_tracer_access_test.cpp`).
   These bind the four architectural decisions that every downstream
   PR inherits, and provide the device-callable surface PR 2-9 use
   to access species fields:

   - *Memory layout*: LayoutRight, FieldLayout tag order
     `(COL, CMP, LEV)`. Matches the MAM4 precedent at
     `components/eamxx/src/physics/mam/eamxx_mam_microphysics_process_interface.cpp:151`.
     Cross-species access at fixed `(col, lev)` is non-unit-stride
     (stride = `nlev_pack * sizeof(Pack)`) but L1-resident for N ≤ 6.
   - *Pack dim placement*: Pack stays over LEV. Pack-over-CMP is
     rejected — it would force `nk_pack=1` and destroy
     level-vectorization in every existing physics kernel.
     Consequence: `get_species_subview(field, i, s)` returns a
     `uview_1d<Pack>` of shape `(nlev_pack,)`, identical in shape
     and stride to the current scalar `(COL, LEV)` prognostic
     subview.
   - *Downstream iteration pattern*: species loop sits OUTSIDE the
     existing `TeamVectorRange` over levels (serial-inside-thread).
     This preserves level vectorization on CPU pack lanes and
     matches the slice-0 BFB requirement by construction at N=1.
     Kernels that want hierarchical species-parallelism via
     `TeamThreadRange`-over-CMP must justify it in their own spec.
   - *N=1 specialization*: `if constexpr (WTRC_MAX_CNST == 1)`
     reduces `get_species_subview` to a direct scalar reference.
     With pack-over-LEV stride identity, the generated machine
     code at N=1 matches the pre-extension scalar path. This is
     the structural precondition for PR 2's `rtol=0, atol=0`
     slice-0 BFB criterion; without it, FMA reordering produces
     1-ULP differences that fail the criterion on day 1.
   - *Cross-species access*: `get_all_species_at_level(field, i, k)`
     returns a value-type `Pack[WTRC_MAX_CNST]`, NOT a Kokkos
     subview. Stack-allocated, aliasing-free, used by fractionation
     kernels in PR 6-9 to load all species at a fixed
     `(col, level-pack)`.

   The Catch2 test asserts these properties directly: data-pointer
   identity at N=1 against a (COL, LEV) field, LayoutRight stride
   between species at N=4, and correct values from
   `get_all_species_at_level`. An advisory codegen-equivalence
   script compares object-file disassembly of a small reference
   kernel built against the scalar prognostic vs the N=1
   prognostic; differences flag the if-constexpr specialization is
   not collapsing cleanly.

Risks:

- **Pre-existing `scream::WaterTracers` namespace clash.** The
  current `water_tracers.hpp` uses PascalCase
  `scream::WaterTracers`; this spec adds lowercase
  `scream::water_tracers`. Both coexist during the transition. The
  PascalCase namespace and its arrays are not removed here — that is
  a cleanup that belongs after group 1 lands.
- **CMake-defined WTRC_MAX_CNST poisoning the default build.** With
  the `#error` guard in place, any TU that includes
  `water_tracers.hpp` without linking `water_tracers` will fail to
  compile. Mitigation: keep `water_tracers` an INTERFACE-link
  dependency for the small set of consumers; default-OFF builds
  still get `WTRC_MAX_CNST=1` via the unconditional propagation.
- **Tracer-index vs. catalog-index drift.** The identity-with-
  truncation mapping is convenient for an N=4 campaign but breaks
  the moment a user wants `[bulk H2O, HDO, H218O]` (N=3, skipping
  H216O). A registry spec (out of campaign scope) is the proper
  fix. Flagged in `out_of_scope`; documented in a follow-up note.
- **FMA reordering at N=1.** If the `if constexpr (WTRC_MAX_CNST == 1)`
  specialization is missing or the compiler does not fully collapse
  the species loop, fused-multiply-add reassociation around a
  one-element loop can produce 1-ULP differences vs the pre-
  extension scalar code, failing PR 2's `rtol=0, atol=0` slice-0
  BFB criterion. Mitigations: the explicit specialization in
  `water_tracer_access.hpp`; the data-pointer-identity unit test
  at N=1; the advisory codegen-equivalence check. All three are
  required because the failure mode is silent and only surfaces in
  PR 2.
- **Argument-count pressure in downstream kernels.** P3's
  `p3_main_impl_part2` already takes ~50 `uview_1d<Pack>` views as
  parameters. Adding species accessors raises the risk of hitting
  CUDA's 4 KB kernel-argument limit at N>1. Not a PR 1 risk per se,
  but the chosen accessor signature in PR 1 determines whether PR 2
  needs to bundle arguments into a struct. Mitigation: keep the
  accessors returning the same `uview_1d<Pack>` type the kernels
  already accept, so PR 2 can pass them through without changing
  signatures at N=1.

Skills relied on: `eamxx-cpp-conventions` (C++17, Kokkos device
callability, EKAT idioms), `scientific-modeling-conventions`
(provenance comments, named constants over magic numbers),
`spec-schema` (this spec's structure), `e3sm-build-and-test`
(validation tier 0, lint).

## References

- `wiso_campaign_plan.md` — campaign source of truth; group 1
  definition and slice-1-unchanged invariant.
- `campaigns/wiso.campaign.md` — manifest entry for this spec.
- `references/old_port/water_tracers.hpp:19-21, 113-117` — prior
  WTRC_MAX_CNST guard and `get_bulk_water_subview` accessor reused
  here.
- `references/old_port/water_isotopes.hpp:59-67` —
  `WaterIsotopologueNames` and `IsotopologueToIndex` catalog (source
  of truth for species ordering).
- `references/old_specs/2026-05-06-eamxx-water-vars-add-tracer-dim.md`
  — prior attempt at the broader rank-lift; lessons absorbed
  (always-rank-3, OFF=size-1 trivial case, separating the count
  mechanism from process-side edits).
- `components/eamxx/src/share/field/field_tag.hpp:36, 54` —
  `FieldTag::Component` / `CMP` short alias.
- `components/eamxx/src/share/grid/point_grid.cpp:93-104` —
  `get_3d_vector_layout` factory the helper calls.
- `components/eamxx/src/share/field/field_group.hpp:34-37` —
  `m_subview_dim = 1` invariant for monolithic groups, the documented
  CMP subview axis.
- `components/eamxx/src/share/atm_process/atmosphere_process.hpp:375-393`
  — current `add_tracer` template the helper sits beside (without
  modifying).
- `components/eamxx/src/physics/water_tracers/water_tracers.hpp:16` —
  literal `WTRC_MAX_CNST = 1` that this spec replaces with the CMake-
  driven `#error`-guarded define.
- `components/eamxx/src/physics/water_tracers/CMakeLists.txt:2` —
  unconditional `SCREAM_TRACE_WATER` definition that this spec
  converts into a real option.
