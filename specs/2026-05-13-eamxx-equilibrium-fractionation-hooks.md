---
spec_id: 2026-05-13-eamxx-equilibrium-fractionation-hooks
spec_type: model-e3sm
spec_version: 1.0
title: Apply equilibrium isotopic fractionation in the water-tracer phase-change hooks
created: 2026-05-13
author: rfiorella
project: EAMXX-wiso
work_summary: Apply equilibrium isotopic fractionation in the water-tracer phase-change hooks
priority: normal
estimated_effort_hours: 6

# Execution mode and checkpoints
execution:
  mode: checkpoint
  checkpoints:
    - after: planning
      requires: user-confirmation
    - after: implementation
      requires: tests-pass
    - after: testing
      requires: user-confirmation
  max_iterations_per_phase: 5
  parallelization:
    allowed: false
    max_subagents: 1
    plan: null

# Contact and ownership
primary_contact: rfiorella
reviewers: []

# Model-specific fields
model_specific:
  subsystem: eamxx
  build_mode: eamxx-standalone-cmake
  target_compset: null
  target_resolution: ne4
  platform: docker-local
  container_image: rfiorella/model-containers:e3sm-openmpi-dev-latest
  validation_tier: 0

# Inputs
inputs:
  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_isotopes.hpp
    description: Equilibrium fractionation functions AlphaEqLiquidVapor() and AlphaEqIceVapor()
    format: C++ header
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_tracer_hooks.hpp
    description: Process-specific hook function-pointer table (6 hooks). This spec wires real implementations for condensation, deposition, melting, freezing.
    format: C++ header
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_tracer_registry.hpp
    description: Compile-time tracer registry providing tracer_isotopologue(i) and tracer_is_tag(i)
    format: C++ header
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/eamxx_water_tracers.cpp
    description: Current passive-copy hook implementations and initialize_water_tracer_hooks(). Replaced by this spec.
    format: C++ source
    required: true

# Deliverables
deliverables:
  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/eamxx_water_tracers.cpp
    description: Replace passive_copy_* hook implementations with eq_* (equilibrium-fractionation, ratio-preserving) implementations for condensation, deposition, melting, freezing. initialize_water_tracer_hooks() registers eq_* hooks.
    format: C++ source
    validation_method: compiles under SCREAM_TRACE_WATER=ON for N=1 and N=4 configs; existing tests still pass

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/tests/test_equilibrium_fractionation_hooks.cpp
    description: Catch2 unit test exercising the four eq_* hooks with synthetic inputs. Covers tag ratio-preservation, isotopologue fractionation, hand-calculated reference values, mass conservation, q_min edge case, α=1 limit.
    format: C++ source
    validation_method: compiles and passes all REQUIRE statements

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/tests/CMakeLists.txt
    description: Register new unit test with CTest
    format: CMake
    validation_method: test appears in test-all-eamxx output

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/REGISTRY_README.md
    description: Update follow-up notes to reflect that equilibrium fractionation is now implemented; identify remaining work (kinetic, evaporation/sublimation, liquid-ice α)
    format: Markdown
    validation_method: manual review

# Success criteria
success_criteria:
  - id: SC1
    phase: implementation
    description: Default build (SCREAM_TRACE_WATER=OFF) compiles unchanged; no-op hook defaults in water_tracer_hooks.hpp are not modified
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | tee build_default.log && \
      ! grep -E "error:|Error " build_default.log
    expected_output: Build completes; no error lines from water_tracers translation unit
    blocking: true

  - id: SC2
    phase: implementation
    description: Build succeeds under SCREAM_TRACE_WATER=ON with N=4 isotopologue registry (bulk + HDO + H218O + one tag)
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      cmake -DSCREAM_TRACE_WATER=ON \
            -DSCREAM_WATER_TRACERS_FILE=cmake/water_tracers/registry_n4.cmake \
            -S . -B build_n4 && \
      cmake --build build_n4 --target water_tracers 2>&1 | tee build_n4.log
    expected_output: libwater_tracers.a builds; eq_condensation, eq_deposition, eq_melting, eq_freezing symbols present
    blocking: true

  - id: SC3
    phase: implementation
    description: New unit test compiles and is registered with CTest
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -q "test_equilibrium_fractionation_hooks"
    expected_output: Test name appears in CTest registration listing
    blocking: true

  - id: SC4
    phase: testing
    description: For tags (tracer_is_tag(i) == true), the ratio q_tag/q_bulk in each phase is unchanged across each of the four hooks (within FP tolerance)
    criterion_type: assertion
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | tee test.log && \
      grep -E "test_equilibrium_fractionation_hooks.*tag_ratio_preservation.*PASSED" test.log
    assertion: For each hook and each tag tracer, |R_after - R_before| < 1e-12 where R = q_tag / q_bulk in both source and destination phases
    blocking: true

  - id: SC5
    phase: testing
    description: For H2O (catalog_idx 0) and H216O (catalog_idx 1) routed through the eq_* hooks, the per-tracer tendency equals the ratio-preserving result with α=1 (i.e., these tracers behave as tags)
    criterion_type: assertion
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_equilibrium_fractionation_hooks.*h2o_h216o_alpha_one.*PASSED"
    assertion: Computed dq for catalog_idx in {0, 1} matches the α=1 ratio-preserving formula within 1e-14
    blocking: true

  - id: SC6
    phase: testing
    description: HDO condensation at T=273.15K and T=253.15K matches hand-calculated reference δ-values from the closed-system formula using AlphaEqLiquidVapor()
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_equilibrium_fractionation_hooks.*hdo_condensation_reference.*PASSED"
    tolerance: absolute 1e-10
    reference_values: Computed by hand from R_cond = α * R_vapor with α from Horita & Wesolowski (1994) at T=273.15K and T=253.15K
    blocking: true

  - id: SC7
    phase: testing
    description: HDO deposition at T=233K and T=253K matches hand-calculated reference δ-values using AlphaEqIceVapor()
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_equilibrium_fractionation_hooks.*hdo_deposition_reference.*PASSED"
    tolerance: absolute 1e-10
    reference_values: Hand-calculated from α_eq(ice/vapor) at T=233K and T=253K
    blocking: true

  - id: SC8
    phase: testing
    description: Mass conservation across each hook — total heavy-isotope mass (sum over source + destination phases for each isotopologue) is conserved to machine precision after hook application
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_equilibrium_fractionation_hooks.*mass_conservation.*PASSED"
    tolerance: absolute 1e-14
    assertion: For each isotopologue tracer and each hook, |q_iso_src_after + q_iso_dst_after - q_iso_src_before - q_iso_dst_before| < 1e-14
    blocking: true

  - id: SC9
    phase: testing
    description: q_min floor (1e-20 kg/kg) prevents NaN/Inf when source phase is at or near zero; tendency falls through to ratio q_iso/max(q_bulk, q_min) without divide-by-zero
    criterion_type: assertion
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_equilibrium_fractionation_hooks.*qmin_floor.*PASSED"
    assertion: With q_bulk_source set to 0 and 1e-25 in test inputs, all hook outputs remain finite (no NaN, no Inf)
    blocking: true

  - id: SC10
    phase: testing
    description: α=1 limit recovers the existing N=2 passive-copy behavior for the special case q_iso == q_bulk (validation against legacy passive-copy test, where applicable)
    criterion_type: assertion
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_equilibrium_fractionation_hooks.*alpha_one_limit.*PASSED"
    assertion: For tracer where q_iso == q_bulk at hook entry, dq_iso = dq_bulk (1e-14 tolerance)
    blocking: false

  - id: SC11
    phase: integration
    description: Human review of code changes for adherence to scientific-modeling-conventions (units in comments, no spurious abstractions, ratio-preserving formula clearly documented at call site) and of REGISTRY_README.md updates for accuracy
    criterion_type: human-review
    reviewer: rfiorella
    blocking: true

# Out of scope
out_of_scope:
  - Kinetic fractionation (deferred to future spec — will add AlphaKin* formulation)
  - Evaporation hook (kinetic-dominated; deferred)
  - Sublimation hook (kinetic-dominated; deferred)
  - Liquid/ice equilibrium fractionation (melting and freezing use α=1 placeholder in this spec; real α_liq_ice deferred)
  - Tag-specific source-attribution semantics (this spec only enforces ratio-preservation for tags; tag evolution by region/source is a future spec)
  - Performance optimization / Kokkos kernel refactoring of the hook bodies
  - Surface flux (ocean/land) isotope-aware coupling
  - 5-day standalone run / Tier 1 conservation diagnostics (deferred — Tier 0 only for now per compute budget)
  - Modifications to P3, SHOC, or any source outside components/eamxx/src/physics/water_tracers/

# Resolved decisions
resolved_decisions:
  - decision: This spec covers equilibrium-fractionation hooks only — condensation, deposition, melting, freezing. Evaporation and sublimation are kinetic-dominated and deferred.
    rationale: User direction (Q1). Kinetic fractionation requires a separate α formulation and additional state (e.g., effective humidity over water/ice) not yet wired through.
    date: 2026-05-13

  - decision: All four in-scope hooks (condensation, deposition, melting, freezing) share a single unified ratio-preserving formula. α is supplied per-process and per-isotopologue.
    rationale: User direction (Q9). H2O and H216O have α=1 by construction in water_isotopes.hpp, so they pass through the same code path consistently with tags. Avoids splitting passive-copy and fractionation code paths.
    date: 2026-05-13

  - decision: Formula — dq_iso = α(iso, T) * (q_iso_src / max(q_bulk_src, q_min)) * dq_bulk_tend, applied symmetrically to source and destination phases.
    rationale: Closed-system small-f approximation, matching the bulk scheme's own saturation-adjustment assumption (Q7). α = 1 for tags, H2O, H216O. α = AlphaEqLiquidVapor(iso, T) for condensation; AlphaEqIceVapor(iso, T) for deposition; α = 1 for melting and freezing in this spec (Q6).
    date: 2026-05-13

  - decision: No sub-stepping inside the hook. α is evaluated at the temperature passed in, against the instantaneous vapor (or ice) ratio at hook entry. Single tendency application per hook call.
    rationale: User direction (Q7). Matches the closure assumption of the bulk parameterization that called the hook. Sub-stepping would require revisiting bulk-scheme integration.
    date: 2026-05-13

  - decision: Sign-preserving application — the per-isotope tendency carries the same sign as dq_bulk_tend with no guard. If the bulk scheme delivers a negative tendency through (e.g.) the condensation hook, isotopes follow.
    rationale: User direction (Q8). Keeps isotope behavior tightly coupled to bulk behavior; any sign-handling issue surfaces at the bulk level, not in isotope code.
    date: 2026-05-13

  - decision: q_min floor of 1e-20 kg/kg applied to the source-phase bulk mass in the ratio q_iso_src / q_bulk_src. No cap on tendency magnitude (conservation is enforced elsewhere).
    rationale: User direction (Q8). 1e-20 is well below any physically meaningful water mixing ratio; prevents NaN/Inf without distorting realistic regimes.
    date: 2026-05-13

  - decision: Retire the passive_copy_* function names; rename to eq_condensation, eq_deposition, eq_melting, eq_freezing. initialize_water_tracer_hooks() registers the new names. No backwards-compatibility shim.
    rationale: User direction (Q9). The new code path supersedes passive-copy semantically (ratio-preserving rather than mass-tendency-copy); keeping the old name would be misleading.
    date: 2026-05-13

  - decision: Validation tier 0 only — compile + unit tests. No 5-day run, no Tier 1 conservation diagnostics in this spec.
    rationale: User direction (Q5). Compute budget. The unit tests cover conservation invariants directly without a model run.
    date: 2026-05-13

  - decision: Unit-test reference data combines hand-calculated point values (HDO at 273.15K, 253.15K for condensation; 253K, 233K for deposition) with invariant checks (ratio-preservation, mass conservation, α=1 limit, finite-output under q_min stress).
    rationale: User direction (Q12). Point values catch formula errors; invariant checks catch logic errors. Tolerance 1e-10 for point comparisons (consistent with the prior fractionation-function tests).
    date: 2026-05-13

  - decision: Work is confined to components/eamxx/src/physics/water_tracers/. P3 and SHOC call sites are not modified in this spec — the hook table is the seam.
    rationale: User direction (Q11). The hook signatures and call sites were established in the prior water-vars spec; this spec only swaps the implementations.
    date: 2026-05-13

# Ask-before actions (project-specific additions to global policy)
ask_before:
  - Modifying any file outside components/eamxx/src/physics/water_tracers/ (e.g., P3 or SHOC source, components/eamxx/CMakeLists.txt, share/ headers)
  - Adding new dependencies or external libraries
  - Changing the hook function-pointer signatures in water_tracer_hooks.hpp (would break the established interface)
  - Changing the public API of water_isotopes.hpp (AlphaEqLiquidVapor / AlphaEqIceVapor signatures or behavior)
  - Removing or renaming the existing water_mass_passive_copy_n2 test if it would lose coverage

# Parallelization
allow_parallelization: false

# Post-completion review
request_performance_review: false
request_code_review: false
---

# Apply Equilibrium Isotopic Fractionation in Water-Tracer Phase-Change Hooks

## Context

Prior specs in the EAMXX-wiso project established three pieces of infrastructure:

1. **Compile-time tracer registry** (spec `2026-05-06-eamxx-water-tracer-isotopologue-registry`) — provides `WTRC_MAX_CNST`, `tracer_isotopologue(i)`, `tracer_name(i)`, and `tracer_is_tag(i)`. The registry distinguishes bulk water (catalog index 0 = H2O), heavy isotopologues (HDO, H218O, H217O, HTO), and passive tag tracers.
2. **Equilibrium fractionation functions** (spec `2026-05-07-test-equilibrium-fractionation`) — `AlphaEqLiquidVapor(isotope, T)` and `AlphaEqIceVapor(isotope, T)` in `water_isotopes.hpp`, validated to 1e-10 tolerance against hand-calculated values from Horita & Wesolowski (1994) and Merlivat (1978). By construction these return α = 1 for H2O and H216O.
3. **Process-specific hook table** (spec `2026-05-06-eamxx-water-vars-add-tracer-dim`) — `water_tracer_hooks.hpp` defines six hooks (`condensation_hook`, `evaporation_hook`, `deposition_hook`, `sublimation_hook`, `freezing_hook`, `melting_hook`) as function pointers initialized to no-ops in default builds. Under `SCREAM_TRACE_WATER=ON`, `initialize_water_tracer_hooks()` currently wires them to `passive_copy_*` implementations that add the bulk tendency to every additional tracer. The hooks are already called from P3 and SHOC at the established phase-change sites.

This spec replaces the passive-copy hook bodies with **equilibrium fractionation physics** for the four equilibrium-dominated processes (condensation, deposition, melting, freezing). Evaporation and sublimation remain on the passive-copy / no-op path until a future spec adds kinetic fractionation.

The work is the first time real isotope physics actually runs inside the EAMxx atmosphere. Tier 0 validation (compile + unit tests) is sufficient for this stage given compute budget; physical-reasonableness diagnostics in a short standalone run will follow in a separate Tier 1 spec.

## Approach

### Unified ratio-preserving formulation

All four in-scope hooks use the same expression for the per-tracer tendency:

```
dq_iso = α(iso, T) * (q_iso_src / max(q_bulk_src, q_min)) * dq_bulk_tend
```

where:
- `q_iso_src` and `q_bulk_src` are the source-phase mass mixing ratios at hook entry
- `q_min = 1e-20 kg/kg` floors the denominator to prevent divide-by-zero
- `dq_bulk_tend` is the bulk phase-change tendency provided by the calling scheme (P3, SHOC, etc.)
- `α(iso, T)` is the fractionation factor — looked up from the registry and `water_isotopes.hpp` as follows:

| Hook | α source | α value |
|---|---|---|
| `eq_condensation` | `AlphaEqLiquidVapor(iso, T)` | 1.0 for H2O/H216O/tags; physical α for HDO/H218O/H217O/HTO |
| `eq_deposition` | `AlphaEqIceVapor(iso, T)` | as above |
| `eq_melting` | constant 1.0 (placeholder; α_liq_ice deferred) | 1.0 for all |
| `eq_freezing` | constant 1.0 (placeholder; α_liq_ice deferred) | 1.0 for all |

The tendency is applied symmetrically: the source phase loses `dq_iso` and the destination phase gains `dq_iso`, sign-preserved with `dq_bulk_tend`.

### Why this unifies tags and isotopologues

For a tag tracer with α = 1, the formula collapses to a strictly ratio-preserving update:

```
q_tag(after) / q_bulk(after) = q_tag(before) / q_bulk(before)
```

For H2O and H216O (catalog indices 0 and 1), `AlphaEqLiquidVapor` and `AlphaEqIceVapor` already return 1.0, so isotopologue tracers for these species also flow through the same code path consistently with tags. The only tracers that experience non-unit α are the heavy isotopologues (HDO, H218O, H217O, HTO) — and they pass through the same loop, just with a different α value.

This replaces the existing `passive_copy_*` implementations, which used `q_iso += dq_bulk` (mass-tendency-copy, not ratio-preserving). The names `passive_copy_*` are retired in favor of `eq_*`.

### Closed-system formulation — what this matches

The formula above is the closed-system / small-f / saturation-adjustment limit of Rayleigh distillation. For condensation, this matches the closure assumption of P3 and SHOC: a single bulk tendency is delivered, and within that tendency the new condensate is in instantaneous equilibrium with the source vapor at temperature T. No sub-stepping is performed inside the hook — α is evaluated once, against the vapor ratio at hook entry.

For deposition, the same closed-system treatment applies, with α from `AlphaEqIceVapor`. Note that deposition in nature is often kinetic (supersaturation over ice drives Jouzel-Merlivat-style kinetic enhancement); this spec uses the pure equilibrium α as a first-pass implementation and flags the kinetic upgrade as future work.

For melting and freezing, α = 1 is used because we do not yet have a `AlphaEqLiquidIce` function. The hooks are wired through the same code path for architectural uniformity (no `passive_copy` / `fractionation` split), making it trivial to swap in α_liq_ice later.

### Edge cases

- **q_min floor.** When `q_bulk_src` is at or below `q_min = 1e-20 kg/kg`, the denominator is clamped. Because `q_iso_src` is bounded above by `q_bulk_src` for any physical isotopologue, the resulting ratio is well-defined and bounded.
- **Sign of `dq_bulk_tend`.** Applied as-given. If the calling scheme produces a negative-sign tendency through (e.g.) the condensation hook, isotopes follow with the same sign. No guard.
- **Tendency exceeds source-phase mass.** Not capped at the hook level. Downstream conservation handling (in the bulk scheme or in a post-hook clipper, future work) is responsible.

### Decomposition

1. **`eamxx_water_tracers.cpp`** — replace six `passive_copy_*` functions with four `eq_*` functions (condensation, deposition, melting, freezing). Evaporation and sublimation retain a passive-copy implementation (or rename to `passive_copy_*` and keep, pending kinetic spec) — to be confirmed during implementation; default is to leave them as the existing passive-copy bodies, untouched.
2. **`initialize_water_tracer_hooks()`** — register `eq_condensation`, `eq_deposition`, `eq_melting`, `eq_freezing`. Evaporation and sublimation continue to point at the existing passive-copy implementations.
3. **Unit test** — `tests/test_equilibrium_fractionation_hooks.cpp`. Synthetic `(qv, qc, qi, T, p, dq_tend)` inputs with a 4-tracer registry: H2O, H216O (or alternative for variety), HDO, tag. Verify:
   - SC4: ratio q_tag/q_bulk preserved across all four hooks
   - SC5: H2O/H216O tracers identical to α=1 ratio-preserving result
   - SC6: HDO condensation matches hand-calculated reference at 273.15K and 253.15K (tolerance 1e-10)
   - SC7: HDO deposition matches hand-calculated reference at 233K and 253K
   - SC8: mass conservation per isotopologue across each hook
   - SC9: finite output under q_min stress (q_bulk_src ∈ {0, 1e-25})
   - SC10: α=1 special case (q_iso == q_bulk) recovers dq_iso = dq_bulk
4. **CMakeLists.txt** — register the new test with CTest, following the pattern used for `test_equilibrium_fractionation` and `water_tracer_registry_n4_test`.
5. **Docs** — update `REGISTRY_README.md` / `REGISTRY_FOLLOWUP.md` to mark equilibrium fractionation as implemented and to list remaining work (kinetic, evaporation/sublimation, real liquid-ice α, Tier 1 standalone run).

Relevant skills to load during implementation: `e3sm-build-and-test`, `e3sm-platforms`, `scientific-modeling-conventions`.

### Risks

- **Existing N=2 passive-copy test (`water_mass_passive_copy_n2.cpp`)** may need to be updated: its current pass criterion assumes mass-tendency-copy, which is equivalent to ratio-preserving only when `q_iso == q_bulk`. If the test sets up `q_iso = q_bulk` initially (likely), it should still pass; if not, it needs adjustment. Verify during implementation.
- **Closed-system vs. Rayleigh discrepancy.** Some bulk schemes (e.g., precipitation removal in P3) are physically Rayleigh-like, but in this spec we only intercept the in-step phase change, not precipitation. The Rayleigh case will manifest naturally over multiple steps. Document this assumption in the source.
- **Floating-point identity for α=1.** `AlphaEqLiquidVapor(H2O, T)` must return exactly 1.0 (not 1.0 + tiny FP noise) for SC5 to hit 1e-14. The prior fractionation test (SC2 in `2026-05-07-test-equilibrium-fractionation`) verified this. Re-verify if `water_isotopes.hpp` has changed.

## References

- **Prior specs (this project):**
  - `2026-05-06-eamxx-water-tracer-isotopologue-registry-COMPLETE.md` — registry API
  - `2026-05-06-eamxx-water-vars-add-tracer-dim-hook-design.md` — hook interface design
  - `2026-05-07-test-equilibrium-fractionation.md` — α function unit tests
- **Source files:**
  - `components/eamxx/src/physics/water_tracers/water_isotopes.hpp`
  - `components/eamxx/src/physics/water_tracers/water_tracer_hooks.hpp`
  - `components/eamxx/src/physics/water_tracers/water_tracer_registry.hpp`
  - `components/eamxx/src/physics/water_tracers/eamxx_water_tracers.cpp`
- **Scientific references:**
  - Horita, J. & Wesolowski, D.J. (1994). Liquid-vapor fractionation of oxygen and hydrogen isotopes of water from the freezing to the critical temperature. *Geochim. Cosmochim. Acta* 58, 3425–3437.
  - Merlivat, L. (1978). Molecular diffusivities of H2 16O, HD 16O, and H2 18O in gases. *J. Chem. Phys.* 69, 2864–2871.
  - Jouzel, J. & Merlivat, L. (1984). Deuterium and oxygen 18 in precipitation: Modeling of the isotopic effects during snow formation. *JGR Atmospheres* 89, 11749–11757. (For future kinetic-deposition spec.)

## Notes

- This spec is intentionally narrow: equilibrium-only, four hooks, Tier 0 validation. The architectural goal is to land the first real isotope physics in EAMxx and verify correctness through targeted unit tests before committing compute to Tier 1 runs.
- Once merged, the unblocked follow-ups are: (i) kinetic fractionation for evaporation/sublimation; (ii) real liquid-ice α; (iii) Tier 1 standalone run with conservation diagnostics; (iv) tag-source attribution semantics.
- The `eq_*` naming is forward-compatible with future `kin_*` hooks for evaporation/sublimation, making the dispatch table self-documenting about which physics applies where.
