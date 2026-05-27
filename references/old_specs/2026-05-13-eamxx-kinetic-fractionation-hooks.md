---
spec_id: 2026-05-13-eamxx-kinetic-fractionation-hooks
spec_type: model-e3sm
spec_version: 1.0
title: Add kinetic fractionation to hydrometeor evaporation and WBF/supersaturation deposition hooks; complete equilibrium sublimation
created: 2026-05-13
author: rfiorella
project: EAMXX-wiso
work_summary: Kinetic fractionation for hydrometeor evaporation (Stewart 1975, two compile-selectable forms) and WBF/supersaturation deposition (Jouzel-Merlivat 1984); equilibrium sublimation completion.
priority: normal
estimated_effort_hours: 10

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
    description: Contains existing equilibrium α functions (AlphaEqLiquidVapor, AlphaEqIceVapor), kinetic stubs to replace (AlphaKineticEvap, AlphaKineticDepo), constants (difrm, dkfac, fkhum, recrit), and the surface-flux constants aksmc/akrfa/akrfb (out of scope here — documented as surface-flux per Q6 investigation).
    format: C++ header
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_tracer_hooks.hpp
    description: Hook function-pointer table. This spec extends evaporation_hook and deposition_hook signatures to carry RH-over-water (h) and saturation-ratio-over-ice (s_i) fields respectively.
    format: C++ header
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/eamxx_water_tracers.cpp
    description: Implements hook bodies. This spec adds kin_evaporation, eq_sublimation, and kin_deposition; the prior eq_deposition is removed (replaced by kin_deposition, which reduces to eq when S_i=1).
    format: C++ source
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_tracer_registry.hpp
    description: Provides tracer_is_tag(i) and tracer_isotopologue(i) used to dispatch tag (α=1 ratio-preserving) vs isotopologue (fractionating) paths.
    format: C++ header
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/CMakeLists.txt
    description: Build-system file. Adds SCREAM_WISO_EVAP_KINETIC_FORM string option (values: multiplicative, stewart_exact; default: stewart_exact) and propagates as a compile definition.
    format: CMake
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/p3
    description: P3 microphysics call sites for evaporation and deposition hooks. Call sites need to pass h (over water) and s_i (over ice) to the updated hook signatures. Ask-before for modifications.
    format: C++ source directory
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/shoc
    description: SHOC call sites for evaporation hook (turbulent moistening tendencies). Needs to pass h. Ask-before for modifications.
    format: C++ source directory
    required: true

# Deliverables
deliverables:
  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_isotopes.hpp
    description: |
      Replace empty stubs AlphaKineticEvap and AlphaKineticDepo with three new template functions returning the new vapor isotope ratio R_E (uniform return contract — not an α):
        - AlphaKinLiquidVapor_multiplicative(iso, T, h) — Stewart 1975 effective-humidity multiplicative form (α_eff = α_eq * α_kin_only with α_kin_only from dkfac/fkhum/difrm); returns R_E = α_eff * R_L
        - AlphaKinLiquidVapor_stewart_exact(iso, T, h, R_L, R_a) — Stewart 1975 eq. 8 closed form using ambient vapor ratio R_a
        - AlphaKinIceVapor(iso, T, S_i) — Jouzel-Merlivat 1984 WBF; collapses to AlphaEqIceVapor at S_i=1
      Update "not sure! RPF" comments for aksmc/akrfa/akrfb with proper attribution as surface-flux (ocean evaporation) kinetic constants — not used in this spec but kept for the future surface-flux pathway.
    format: C++ header
    validation_method: compiles; functions return correct values vs hand-calculation; old stubs absent

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_tracer_hooks.hpp
    description: |
      Extend EvaporationHookFn signature to include `const view_2d& h` (RH over water at grid cell, dimensionless 0–1). Extend DepositionHookFn signature to include `const view_2d& s_i` (saturation ratio over ice, dimensionless). SublimationHookFn signature unchanged (equilibrium-only). Update no-op default implementations to match new signatures.
    format: C++ header
    validation_method: compiles in both default (OFF) and TRACE_WATER=ON builds

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/eamxx_water_tracers.cpp
    description: |
      - Remove eq_deposition; add kin_deposition using AlphaKinIceVapor (Jouzel-Merlivat).
      - Remove existing passive_copy_evaporation and passive_copy_sublimation; add kin_evaporation (uniform ratio-preserving for tags, Stewart form for isotopologues, dispatch via #ifdef on SCREAM_WISO_EVAP_KINETIC_FORM_STEWART_EXACT) and eq_sublimation (uses AlphaEqIceVapor).
      - Update initialize_water_tracer_hooks() to register kin_evaporation, eq_sublimation, kin_deposition alongside the unchanged eq_condensation, eq_melting, eq_freezing.
    format: C++ source
    validation_method: compiles for both CMake form options; unit tests pass

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/CMakeLists.txt
    description: |
      Add option SCREAM_WISO_EVAP_KINETIC_FORM (string, default "stewart_exact", valid values: "multiplicative" | "stewart_exact"). Configure-time validation rejects unknown values with clear error. Propagate to source as SCREAM_WISO_EVAP_KINETIC_FORM_STEWART_EXACT compile definition when value == "stewart_exact".
    format: CMake
    validation_method: cmake -L shows the option; invalid value rejected at configure time

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/tests/test_kinetic_fractionation_hooks.cpp
    description: |
      Catch2 unit tests covering: (a) AlphaKinLiquidVapor_multiplicative hand-calc verification; (b) AlphaKinLiquidVapor_stewart_exact hand-calc verification; (c) AlphaKinIceVapor hand-calc verification; (d) collapse-to-equilibrium invariants (h=1 for evap, S_i=1 for depo); (e) hook-level ratio preservation for tags; (f) hook-level mass conservation per isotopologue; (g) kin_deposition at S_i=1 matches the prior eq_deposition computation (regression).
    format: C++ source
    validation_method: compiles and all REQUIRE statements pass under both CMake form values

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/tests/CMakeLists.txt
    description: Register test_kinetic_fractionation_hooks with CTest.
    format: CMake
    validation_method: test appears in test-all-eamxx output

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/p3
    description: |
      Update P3 call sites that invoke evaporation_hook and deposition_hook to pass the new h and s_i fields. h is computed from cell qv and saturation specific humidity over water at cell T; s_i is computed from qv and saturation specific humidity over ice at cell T. Use existing P3 saturation functions (Polysvp1 or equivalent) — no new physics calls. ASK BEFORE editing any P3 file.
    format: C++ source
    validation_method: build succeeds with new signatures; unit-test infrastructure still passes

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/shoc
    description: |
      Update SHOC call sites that invoke evaporation_hook (turbulent moistening with sub-saturation tendencies) to pass h. ASK BEFORE editing any SHOC file.
    format: C++ source
    validation_method: build succeeds; turbulent-mixing path compiles

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/REGISTRY_README.md
    description: Update follow-up notes to mark kinetic evaporation (hydrometeor) and WBF/supersaturation deposition as implemented; identify remaining work (surface evaporation via AlphaKMol/coupler, real liquid/ice α, Tier 1 standalone run).
    format: Markdown
    validation_method: manual review

# Success criteria
success_criteria:
  - id: SC1
    phase: implementation
    description: Default build (SCREAM_TRACE_WATER=OFF) compiles unchanged; no-op hook defaults updated for new signatures but produce no runtime cost
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | tee build_default.log && \
      ! grep -E "error:|Error " build_default.log
    expected_output: Build completes; no error lines from water_tracers translation unit
    blocking: true

  - id: SC2
    phase: implementation
    description: Build succeeds under SCREAM_TRACE_WATER=ON with N=4 registry and SCREAM_WISO_EVAP_KINETIC_FORM=stewart_exact (default)
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      cmake -DSCREAM_TRACE_WATER=ON \
            -DSCREAM_WATER_TRACERS_FILE=cmake/water_tracers/registry_n4.cmake \
            -DSCREAM_WISO_EVAP_KINETIC_FORM=stewart_exact \
            -S . -B build_n4_stewart && \
      cmake --build build_n4_stewart --target water_tracers 2>&1 | tee build_n4_stewart.log
    expected_output: libwater_tracers.a builds; kin_evaporation_stewart_exact symbol bound
    blocking: true

  - id: SC3
    phase: implementation
    description: Build succeeds with SCREAM_WISO_EVAP_KINETIC_FORM=multiplicative
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      cmake -DSCREAM_TRACE_WATER=ON \
            -DSCREAM_WATER_TRACERS_FILE=cmake/water_tracers/registry_n4.cmake \
            -DSCREAM_WISO_EVAP_KINETIC_FORM=multiplicative \
            -S . -B build_n4_mult && \
      cmake --build build_n4_mult --target water_tracers 2>&1 | tee build_n4_mult.log
    expected_output: libwater_tracers.a builds; kin_evaporation_multiplicative symbol bound
    blocking: true

  - id: SC4
    phase: implementation
    description: Invalid value for SCREAM_WISO_EVAP_KINETIC_FORM produces a clear configure-time error
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ! cmake -DSCREAM_TRACE_WATER=ON \
              -DSCREAM_WATER_TRACERS_FILE=cmake/water_tracers/registry_n4.cmake \
              -DSCREAM_WISO_EVAP_KINETIC_FORM=invalid_value \
              -S . -B build_bad 2>&1 | tee bad.log && \
      grep -E "SCREAM_WISO_EVAP_KINETIC_FORM.*(multiplicative|stewart_exact)" bad.log
    expected_output: CMake fails; error message lists the two valid values
    blocking: true

  - id: SC5
    phase: testing
    description: AlphaKinLiquidVapor_multiplicative at (T=298.15K, h=0.7) for HDO matches hand-calculated R_E from the multiplicative form (α_eff = α_eq * (1 + (1-h_eff)*(1/difrm - 1)*dkfac)^-1, h_eff = fkhum*h + (1-fkhum), R_E = α_eff * R_L_VSMOW)
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_kinetic_fractionation_hooks.*alpha_kin_lv_multiplicative_hdo.*PASSED"
    tolerance: absolute 1e-10
    reference_values: Hand-calculated from Stewart 1975 multiplicative form using difrm(HDO)=0.9757, dkfac=0.58, fkhum=0.25, AlphaEqLiquidVapor(HDO, 298.15K) from Horita & Wesolowski (1994).
    blocking: true

  - id: SC6
    phase: testing
    description: AlphaKinLiquidVapor_stewart_exact at (T=298.15K, h=0.7, R_L=R_VSMOW(HDO), R_a=0.5*R_VSMOW(HDO)) matches hand-calculated R_E from Stewart 1975 eq. 8
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_kinetic_fractionation_hooks.*alpha_kin_lv_stewart_hdo.*PASSED"
    tolerance: absolute 1e-10
    reference_values: Hand-calculated from Stewart 1975 eq. 8 (closed form with explicit R_a dependence) at given inputs.
    blocking: true

  - id: SC7
    phase: testing
    description: AlphaKinIceVapor at (T=253K, S_i=1.2) for HDO and H218O matches hand-calculated R_E from Jouzel-Merlivat 1984 WBF formulation
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_kinetic_fractionation_hooks.*alpha_kin_iv_jm84.*PASSED"
    tolerance: absolute 1e-10
    reference_values: Hand-calculated from Jouzel & Merlivat (1984) α_eff = α_eq * S_i / (α_eq * (1/difrm) * (S_i - 1) + 1) using difrm and AlphaEqIceVapor at given (T, isotope).
    blocking: true

  - id: SC8
    phase: testing
    description: Collapse-to-equilibrium invariants — AlphaKinLiquidVapor (both forms) at h=1 reproduces AlphaEqLiquidVapor * R_L; AlphaKinIceVapor at S_i=1 reproduces AlphaEqIceVapor * R_v
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_kinetic_fractionation_hooks.*collapse_to_eq.*PASSED"
    tolerance: absolute 1e-12
    assertion: |
      For all isotopologues and a range of temperatures:
        |AlphaKinLiquidVapor_multiplicative(iso, T, h=1.0) - AlphaEqLiquidVapor(iso, T) * R_L| < 1e-12
        |AlphaKinLiquidVapor_stewart_exact(iso, T, h=1.0, R_L, R_a=R_L) - AlphaEqLiquidVapor(iso, T) * R_L| < 1e-12
        |AlphaKinIceVapor(iso, T, S_i=1.0) - AlphaEqIceVapor(iso, T) * R_v| < 1e-12
    blocking: true

  - id: SC9
    phase: testing
    description: kin_evaporation hook ratio-preserves for tags across the four-tracer test configuration
    criterion_type: assertion
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_kinetic_fractionation_hooks.*kin_evap_tag_ratio.*PASSED"
    assertion: For tags through kin_evaporation, |R_tag_after/R_bulk_after - R_tag_before/R_bulk_before| < 1e-12 in both source and destination phases
    blocking: true

  - id: SC10
    phase: testing
    description: Mass conservation per isotopologue across kin_evaporation, eq_sublimation, and kin_deposition hooks
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_kinetic_fractionation_hooks.*mass_conservation.*PASSED"
    tolerance: absolute 1e-14
    assertion: For each isotopologue and each of the three new/modified hooks, |q_iso_src_after + q_iso_dst_after - q_iso_src_before - q_iso_dst_before| < 1e-14
    blocking: true

  - id: SC11
    phase: testing
    description: kin_deposition at S_i=1 produces output bit-for-bit equivalent (or within 1e-14) to the prior spec's eq_deposition at the same inputs — confirms the Jouzel-Merlivat formulation reduces correctly to the equilibrium case
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_kinetic_fractionation_hooks.*depo_regression_si_one.*PASSED"
    tolerance: absolute 1e-14
    reference_values: Computed in-test by also calling AlphaEqIceVapor and the ratio-preserving formula from the prior spec.
    blocking: true

  - id: SC12
    phase: testing
    description: eq_sublimation hook produces correct fractionation for HDO and H218O at T=233K, T=253K — δ-value of newly sublimated vapor is depleted relative to source ice
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_kinetic_fractionation_hooks.*eq_sublimation.*PASSED"
    tolerance: absolute 1e-10
    reference_values: Hand-calculated from AlphaEqIceVapor (ice → vapor direction; vapor ratio = R_ice / α_eq).
    blocking: true

  - id: SC13
    phase: integration
    description: Human review of hook-signature changes (water_tracer_hooks.hpp), P3/SHOC call-site edits, CMake option implementation, and updated REGISTRY_README.md. Confirm the surface-evaporation pathway (AlphaKMol, aksmc/akrfa/akrfb) is unmodified and that the "not sure! RPF" comments have been replaced with proper Q6-attribution text.
    criterion_type: human-review
    reviewer: rfiorella
    blocking: true

# Out of scope
out_of_scope:
  - Surface evaporation (ocean/land flux) kinetic fractionation — handled by AlphaKMol and the wind-speed-based aksmc/akrfa/akrfb constants; this pathway requires coupler-level integration and is a separate spec
  - Modifying AlphaKMol or the surface-flux constants (aksmc/akrfa/akrfb) — comments will be updated for attribution only; numerical values and the AlphaKMol function are not touched
  - Real liquid/ice equilibrium fractionation (eq_melting and eq_freezing remain α=1 from the prior spec; α_liq_ice deferred)
  - Kinetic fractionation for sublimation (this spec treats sublimation as equilibrium-only; eq_sublimation is the first real sublimation hook implementation)
  - Tag-specific source-attribution semantics
  - Lee et al. 2007 raindrop-fall-distance corrections to evaporation kinetic
  - Tier 1 standalone run / 5-day conservation diagnostics
  - Performance optimization or Kokkos-kernel refactoring of hook bodies
  - Reconsidering the Q6 surface-flux constants origin beyond comment updates (a separate spec may revisit aksmc/akrfa/akrfb when wiring the surface pathway)

# Resolved decisions
resolved_decisions:
  - decision: Bundle three pieces into one spec — hydrometeor evaporation (kinetic), sublimation (equilibrium-only completion), WBF/supersaturation deposition (kinetic).
    rationale: User direction (Q1). They share the molecular-diffusivity machinery (difrm) and α-function additions; sublimation completes the equilibrium baseline that the prior spec left open.
    date: 2026-05-13

  - decision: Use the existing constants in water_isotopes.hpp WaterIsotopologues struct — `difrm` (Merlivat 1978 + Schoenemann 2014), `dkfac` (0.58, Stewart 1975), `fkhum` (0.25), `recrit` (1.0). No new constants needed.
    rationale: User direction (Q2). Constants were already placed in the struct by the original developer, clearly for Stewart-style kinetic. Reuse keeps a single source of truth.
    date: 2026-05-13

  - decision: aksmc/akrfa/akrfb are wind-speed-dependent surface-flux kinetic constants for ocean evaporation (smooth regime v10<7 m/s uses aksmc; rough regime uses akrfa*v10+akrfb). Origin traced to the original CAM-wiso Fortran in share/csm_share/shr/water_isotopes.F90 (commit 3d72660c01). These are NOT touched in this spec; their comments are updated to remove "not sure! RPF" and replace with proper attribution.
    rationale: Q6 investigation finding. They belong to the same surface-evaporation pathway as AlphaKMol, which Q10 declared out of scope. Documenting attribution now while we have it.
    date: 2026-05-13

  - decision: Extend evaporation_hook signature with `const view_2d& h` (RH over water at grid cell, 0–1) and deposition_hook signature with `const view_2d& s_i` (saturation ratio over ice). Sublimation_hook signature unchanged.
    rationale: User direction (Q3 option a). Kinetic fractionation needs the saturation state; passing it through the hook keeps consistency with the bulk scheme's own calculation rather than recomputing inside the hook with potentially divergent saturation functions.
    date: 2026-05-13

  - decision: kin_deposition replaces eq_deposition from the prior spec. The Jouzel-Merlivat 1984 formulation reduces analytically to AlphaEqIceVapor at S_i=1, so this is a strict generalization, not a different physics choice.
    rationale: User direction (Q4 option a). Avoids a runtime branch on S_i and a "two-deposition-hooks" architecture. SC11 verifies the reduction empirically.
    date: 2026-05-13

  - decision: Validation tier 0 — compile + unit tests. No 5-day run; defer Tier 1 to a separate "first integrated isotope physics run" spec.
    rationale: User direction (Q5). Compute budget.
    date: 2026-05-13

  - decision: For atmospheric hydrometeor evaporation use Stewart 1975 formulation, NOT Craig-Gordon. Constants dkfac=0.58 and fkhum=0.25 already in WaterIsotopologues struct are Stewart-style and used directly.
    rationale: User direction (Q7). Matches the original developer intent embodied in the existing constants.
    date: 2026-05-13

  - decision: Replace empty stubs AlphaKineticEvap and AlphaKineticDepo with new functions named to match the AlphaEqLiquidVapor/AlphaEqIceVapor convention — AlphaKinLiquidVapor_multiplicative, AlphaKinLiquidVapor_stewart_exact, AlphaKinIceVapor. Stubs are deleted, not preserved.
    rationale: User direction (Q8). Consistent naming aids readability; empty stubs have no value once real implementations land.
    date: 2026-05-13

  - decision: One evaporation hook for hydrometeor evaporation only (cloud water and rain re-evaporation in P3, sub-saturation tendencies in SHOC). Surface evaporation (ocean/land flux) is OUT OF SCOPE and is the future AlphaKMol-pathway spec.
    rationale: User direction (Q9 / Q10). Hydrometeor evaporation uses ambient grid-cell RH because droplet/raindrop temperature is not directly available and Stewart 1975's effective-humidity factor fkhum is meant to absorb that uncertainty. Surface evaporation has known surface T and saturation conditions and uses a different parameterization pathway.
    date: 2026-05-13

  - decision: Implement BOTH Stewart 1975 forms — multiplicative (α_eff = α_eq * α_kin_only) and exact (Stewart 1975 eq. 8 closed form with explicit R_a). CMake string option SCREAM_WISO_EVAP_KINETIC_FORM (values multiplicative | stewart_exact, default stewart_exact) selects which is wired into kin_evaporation_hook via #ifdef. Both functions are always compiled (so both are independently unit-testable).
    rationale: User direction (Q11/Q12). Default to the physically more complete form; let users compare with the convention used by other isotope-enabled GCMs (multiplicative) via build flag. String option (not bool) is self-documenting in cmake -L output.
    date: 2026-05-13

  - decision: Uniform return contract — both AlphaKinLiquidVapor forms and AlphaKinIceVapor return the new vapor isotope ratio R_E directly (not an α). The hook then computes dq_iso = R_E * dq_bulk_evap (or dq_bulk_subl/_depo as appropriate).
    rationale: Multiplicative and Stewart-exact forms differ in whether α depends on R_a; returning R_E uniformly hides that asymmetry from the hook and gives a single dispatch code path. Tags (α=1) still flow through the ratio-preserving formula from the prior spec.
    date: 2026-05-13

  - decision: AlphaKMol is not modified in this spec. Surface evaporation kinetic fractionation (ocean flux) remains under AlphaKMol or the wind-speed-based aksmc/akrfa/akrfb constants; both are unchanged.
    rationale: User direction (Q10). Separate concern from in-atmosphere phase changes.
    date: 2026-05-13

# Ask-before actions (project-specific additions to global policy)
ask_before:
  - Modifying any source under components/eamxx/src/physics/p3/ (call-site updates for new evaporation_hook and deposition_hook signatures)
  - Modifying any source under components/eamxx/src/physics/shoc/ (call-site updates for new evaporation_hook signature)
  - Modifying components/eamxx/CMakeLists.txt (adding SCREAM_WISO_EVAP_KINETIC_FORM option and validation)
  - Modifying or deleting AlphaKMol in water_isotopes.hpp (out of scope)
  - Changing the surface-flux constants aksmc/akrfa/akrfb numerical values (only comment updates allowed)
  - Adding new external dependencies or libraries
  - Renaming or relocating the existing eq_condensation, eq_melting, or eq_freezing hooks from the prior spec
  - Restructuring the hook function-pointer table beyond the documented signature extensions

# Parallelization
allow_parallelization: false

# Post-completion review
request_performance_review: false
request_code_review: false
---

# Add Kinetic Fractionation to Hydrometeor Evaporation and WBF/Supersaturation Deposition Hooks; Complete Equilibrium Sublimation

## Context

The prior spec (`2026-05-13-eamxx-equilibrium-fractionation-hooks`) wired equilibrium fractionation for four phase-change hooks (condensation, deposition, melting, freezing). It explicitly deferred evaporation and sublimation to a "kinetic" spec, and treated deposition as a pure-equilibrium first pass. In nature, three of these processes are dominated by kinetic effects:

1. **Hydrometeor evaporation** (cloud water and rain re-evaporating into sub-saturated air) — kinetic fractionation enhances depletion of the residual liquid and modifies the isotope signature of the re-evaporated vapor. Standard treatment is Stewart (1975), using an effective humidity formulation.
2. **Vapor deposition under supersaturation** (mixed-phase clouds and pure ice clouds where S_i > 1) — the Wegener-Bergeron-Findeisen mechanism preferentially grows ice at the expense of liquid. Standard treatment is Jouzel-Merlivat (1984), which reduces the effective α to less than α_eq because of the ice-side kinetic correction.
3. **Sublimation** (ice → vapor) — *not* this spec. Per user direction (Q1), sublimation is treated as equilibrium-only for now, leaving room for a future spec to add kinetic effects relevant to e.g. snowflake sublimation in dry layers.

This spec also closes a small gap: the prior equilibrium spec did not include a `eq_sublimation` hook because sublimation was kinetic-dominated. With sublimation now treated as equilibrium-only (per user revision), this spec adds `eq_sublimation` using `AlphaEqIceVapor`.

### Provenance investigation (Q6)

During scoping, the existing constants `aksmc`, `akrfa`, `akrfb` in `WaterIsotopologues` were marked "not sure! RPF" in the source. Git archaeology traced them to commit `3d72660c01`, where they appear in the original CAM-wiso Fortran (`share/csm_share/shr/water_isotopes.F90`) as **surface kinetic exchange constants** for the ocean evaporation flux, parameterized by 10m wind speed:

```fortran
if (v10 < 7.0)  kmol = aksmc(isp)                      ! smooth regime
else            kmol = akrfa(isp)*v10 + akrfb(isp)     ! rough regime
alphakn = 1 - kmol
```

This places them in the same physical domain as `AlphaKMol` — both parameterize ocean-surface evaporation kinetic fractionation, the only difference being inputs (ustar/zbot/rbot vs v10). They are **not** for in-atmosphere phase changes. Per Q10, the entire surface-evaporation pathway is out of scope for this spec. The constants are preserved untouched; this spec only updates their comments to replace "not sure! RPF" with proper attribution.

## Approach

### α-function additions to water_isotopes.hpp

Three new template functions, replacing the empty `AlphaKineticEvap` and `AlphaKineticDepo` stubs:

**(1) AlphaKinLiquidVapor_multiplicative(iso, T, h)**

Stewart 1975 multiplicative form. Returns the new vapor isotope ratio:

```
h_eff   = fkhum * h + (1 - fkhum)
α_eq    = AlphaEqLiquidVapor(iso, T)
α_kin   = 1 / (1 + (1 - h_eff) * (1 / difrm(iso) - 1) * dkfac)
α_eff   = α_eq * α_kin
R_E     = α_eff * R_L                          ← returned value, R_L taken from source ratio
```

For tags (and isotopologues with α_eq = 1 = α_kin), this reduces to R_E = R_L — no fractionation, ratio preserved.

**(2) AlphaKinLiquidVapor_stewart_exact(iso, T, h, R_L, R_a)**

Stewart 1975 eq. 8 closed form. Explicitly threads the ambient vapor isotope ratio R_a through:

```
α_eq   = AlphaEqLiquidVapor(iso, T)
R_E    = α_eq^(-1) * (R_L - h * R_a) / ((1/difrm(iso))^dkfac * (1 - h))
```

At h = 1 (saturated, no evaporation): the formula is degenerate (0/0), but physically dq_evap = 0 in that limit, so the hook does not call the function. The unit test (SC8) verifies the limit by approaching h → 1 from below with the convention that R_E → α_eq * R_L (closed-system equilibrium limit, achieved when R_a → α_eq * R_L). At R_a = R_L (ambient and source-phase ratios equal): the formula simplifies and matches the multiplicative-form result.

**(3) AlphaKinIceVapor(iso, T, S_i)**

Jouzel-Merlivat 1984 WBF formulation. Returns the new vapor isotope ratio for deposition under supersaturation:

```
α_eq   = AlphaEqIceVapor(iso, T)
α_kin  = S_i / (α_eq * (1/difrm(iso)) * (S_i - 1) + 1)
α_eff  = α_eq * α_kin
R_v_new = α_eff * R_v                            ← returned, R_v = source vapor ratio
```

Note: for deposition the source is vapor and destination is ice. The "new vapor ratio" returned here refers to the *deposit's* equivalent vapor ratio at deposition time — the destination phase formed from R_v with α_eff. The hook applies this as: `dq_i_iso = α_eff * R_v * dq_dep_bulk`.

At S_i = 1: α_kin = 1, α_eff = α_eq, recovering the prior spec's eq_deposition exactly (SC11 verifies). The Jouzel-Merlivat formula is a strict generalization.

### Hook implementations (eamxx_water_tracers.cpp)

**kin_evaporation (replaces passive_copy_evaporation):**

```cpp
void kin_evaporation(const view_3d& qv, const view_3d& qc,
                     const view_2d& T, const view_2d& p,
                     const view_2d& dqc_tend,
                     const view_2d& h,                   // NEW
                     int ncol, int nlev) {
  if (WTRC_MAX_CNST == 1) return;
  for (int i = 1; i < WTRC_MAX_CNST; ++i) {
    const auto iso_idx = tracer_isotopologue(i);
    const auto is_tag  = tracer_is_tag(i);
    for (int icol = 0; icol < ncol; ++icol) {
      for (int ilev = 0; ilev < nlev; ++ilev) {
        const Real qc_src = qc(icol, 0, ilev);
        const Real R_L    = qc(icol, i, ilev) / std::max(qc_src, q_min);
        Real R_E;
        if (is_tag || iso_idx == 0 || iso_idx == 1) {
          R_E = R_L;                              // ratio-preserving
        } else {
          const auto iso_name = WaterIsotopologueNames[iso_idx];
#ifdef SCREAM_WISO_EVAP_KINETIC_FORM_STEWART_EXACT
          const Real qv_src = qv(icol, 0, ilev);
          const Real R_a    = qv(icol, i, ilev) / std::max(qv_src, q_min);
          R_E = AlphaKinLiquidVapor_stewart_exact(iso_name, T(icol, ilev),
                                                  h(icol, ilev), R_L, R_a);
#else
          R_E = AlphaKinLiquidVapor_multiplicative(iso_name, T(icol, ilev),
                                                   h(icol, ilev)) * R_L; // returns α_eff * R_L
#endif
        }
        const Real dq_iso = R_E * dqc_tend(icol, ilev);   // dqc_tend < 0 for evap
        qc(icol, i, ilev) += dqc_tend(icol, ilev) > 0 ? dq_iso : -dq_iso; // sign-preserving
        qv(icol, i, ilev) -= dqc_tend(icol, ilev) > 0 ? dq_iso : -dq_iso;
      }
    }
  }
}
```

(Sign-handling sketched here; final implementation will track sign conventions cleanly.)

**eq_sublimation (new — equilibrium-only):**

```cpp
void eq_sublimation(const view_3d& qi, const view_3d& qv,
                    const view_2d& T, const view_2d& p,
                    const view_2d& dqi_tend, int ncol, int nlev) {
  // Same ratio-preserving + α structure as eq_condensation, but with:
  //   source = qi, destination = qv
  //   α = AlphaEqIceVapor(iso_name, T)
  //   R_E = α * R_i  (sublimated vapor ratio derived from source ice ratio)
  // Wait — for sublimation, the new vapor ratio is R_v = R_i / α_eq (since ice is heavier).
  // The hook applies: dq_v_iso = (R_i / α_eq) * dq_subl_bulk
}
```

(Detailed convention: α_eq(ice/vap) is defined so R_ice/R_vap = α_eq > 1. For sublimation, new vapor has R_v_new = R_ice_source / α_eq.)

**kin_deposition (replaces eq_deposition):**

```cpp
void kin_deposition(const view_3d& qv, const view_3d& qi,
                    const view_2d& T, const view_2d& p,
                    const view_2d& dqv_tend,
                    const view_2d& s_i,                  // NEW
                    int ncol, int nlev) {
  // Same structure as eq_condensation but with α_eff from AlphaKinIceVapor(iso, T, S_i)
  // which collapses to AlphaEqIceVapor at S_i=1.
}
```

### Hook-signature extensions and call-site impact

`water_tracer_hooks.hpp` signature changes:

```cpp
using EvaporationHookFn = void(*)(
  const view_3d& qv, const view_3d& qc,
  const view_2d& T, const view_2d& p,
  const view_2d& dqc_tend,
  const view_2d& h,                          // ADDED
  int ncol, int nlev
);

using DepositionHookFn = void(*)(
  const view_3d& qv, const view_3d& qi,
  const view_2d& T, const view_2d& p,
  const view_2d& dqv_tend,
  const view_2d& s_i,                        // ADDED
  int ncol, int nlev
);
```

`SublimationHookFn` is unchanged. The default no-op functions and the existing call sites in P3 and SHOC must be updated to match. P3/SHOC modifications are ask-before per project policy.

**At call sites, the new fields are computed from existing P3/SHOC state:**
- `h = qv / qsat_water(T, p)` using P3's existing saturation function (Polysvp1 over water)
- `s_i = qv / qsat_ice(T, p)` using P3's existing saturation function (Polysvp1 over ice)

No new physics calls — only re-using already-computed grid-cell saturation diagnostics. The values are passed as views over (COL, LEV).

### CMake option

```cmake
set(SCREAM_WISO_EVAP_KINETIC_FORM "stewart_exact" CACHE STRING
    "Form of Stewart 1975 kinetic evaporation: multiplicative | stewart_exact")
set_property(CACHE SCREAM_WISO_EVAP_KINETIC_FORM PROPERTY STRINGS
             "multiplicative" "stewart_exact")

if (NOT (SCREAM_WISO_EVAP_KINETIC_FORM STREQUAL "multiplicative" OR
         SCREAM_WISO_EVAP_KINETIC_FORM STREQUAL "stewart_exact"))
  message(FATAL_ERROR
    "SCREAM_WISO_EVAP_KINETIC_FORM must be one of: multiplicative, stewart_exact. "
    "Got: ${SCREAM_WISO_EVAP_KINETIC_FORM}")
endif()

if (SCREAM_WISO_EVAP_KINETIC_FORM STREQUAL "stewart_exact")
  target_compile_definitions(water_tracers PRIVATE SCREAM_WISO_EVAP_KINETIC_FORM_STEWART_EXACT)
endif()
```

### Decomposition

1. `water_isotopes.hpp` — add three new template α functions; delete stubs; update `aksmc/akrfa/akrfb` attribution comments
2. `water_tracer_hooks.hpp` — extend `EvaporationHookFn` and `DepositionHookFn` signatures; update no-op defaults
3. `eamxx_water_tracers.cpp` — add `kin_evaporation`, `eq_sublimation`, `kin_deposition`; remove `eq_deposition`; update `initialize_water_tracer_hooks()`; CMake-conditional dispatch
4. `CMakeLists.txt` — add `SCREAM_WISO_EVAP_KINETIC_FORM` option + validation
5. Unit test `test_kinetic_fractionation_hooks.cpp` — covers all functions, both CMake forms, hook-level invariants
6. **P3 call-site updates** (ASK BEFORE) — pass `h` and `s_i` to the relevant hook calls
7. **SHOC call-site updates** (ASK BEFORE) — pass `h` to the relevant hook calls
8. `REGISTRY_README.md` — mark kinetic evap and WBF depo as implemented; update follow-up list

### Risks

- **Stewart-exact at h → 1.** The formula is 0/0 at h=1. The hook never calls the function in that limit (dq_evap = 0 by construction in the bulk scheme), but the unit test (SC8) must approach the limit carefully — recommended pattern: pick h = 0.9999 and verify the result is within 1e-8 of α_eq * R_L.
- **Sign convention on dqc_tend in evaporation.** Per the prior spec, sign is preserved from the bulk tendency. `dqc_tend < 0` means evaporation (qc decreasing, qv increasing). The hook must apply isotope tendency consistently — both qc decreases and qv increases proportionally, with the new vapor ratio R_E driving the *amount* moving.
- **Saturation-function call-site source.** Using P3's saturation function at the call site (not recomputing in the hook) ensures consistency with P3's own decision about whether the cell is supersaturated. The downside is that SHOC may use a different saturation function — verify during call-site updates and either standardize or document the small discrepancy.
- **kin_deposition at S_i < 1.** Subsaturated-over-ice cases. The Jouzel-Merlivat formula stays well-defined for S_i ∈ (0, 1], but physically deposition would not occur there (the bulk scheme delivers dq_dep_bulk = 0). Document that the hook is only meaningful for S_i ≥ 1 and that other cases produce zero tendency.

Relevant skills to load during implementation: `e3sm-build-and-test`, `e3sm-platforms`, `scientific-modeling-conventions`.

## References

- **Prior specs (this project):**
  - `2026-05-13-eamxx-equilibrium-fractionation-hooks.md` — equilibrium-α baseline; this spec replaces its `eq_deposition` with `kin_deposition` and adds `kin_evaporation`, `eq_sublimation`
  - `2026-05-07-test-equilibrium-fractionation.md` — α function unit tests
  - `2026-05-06-eamxx-water-tracer-isotopologue-registry-COMPLETE.md` — registry API
  - `2026-05-06-eamxx-water-vars-add-tracer-dim-hook-design.md` — hook interface design
- **Source files:**
  - `components/eamxx/src/physics/water_tracers/water_isotopes.hpp`
  - `components/eamxx/src/physics/water_tracers/water_tracer_hooks.hpp`
  - `components/eamxx/src/physics/water_tracers/eamxx_water_tracers.cpp`
  - `components/eamxx/src/physics/water_tracers/water_tracer_registry.hpp`
  - Original Fortran provenance (git history): `share/csm_share/shr/water_isotopes.F90` in commit `3d72660c01`
- **Scientific references:**
  - Stewart, M. K. (1975). Stable isotope fractionation due to evaporation and isotopic exchange of falling waterdrops: Applications to atmospheric processes and evaporation of lakes. *JGR* 80, 1133–1146.
  - Merlivat, L. (1978). Molecular diffusivities of H2 16O, HD 16O, and H2 18O in gases. *J. Chem. Phys.* 69, 2864–2871.
  - Jouzel, J. & Merlivat, L. (1984). Deuterium and oxygen 18 in precipitation: Modeling of the isotopic effects during snow formation. *JGR Atmos.* 89, 11749–11757.
  - Horita, J. & Wesolowski, D.J. (1994). *Geochim. Cosmochim. Acta* 58, 3425–3437. (Equilibrium α from prior spec, used here at h=1 limits.)
  - Schoenemann, S. W. et al. (2014). Triple oxygen isotope analysis of polar firn — context for `difrm` diffusivity ratios.

## Notes

- This is the second "real isotope physics" spec for EAMxx, following the equilibrium hooks. After this lands, the in-atmosphere phase-change isotope physics is feature-complete except for: (i) real liquid/ice α (small enhancement for melting/freezing); (ii) Lee et al. 2007 raindrop-fall corrections; (iii) the surface-evaporation pathway via AlphaKMol / coupler.
- The CMake compile-selectable evaporation form is a deliberate research-code feature, not a production toggle. The default (`stewart_exact`) is the physically more complete formulation; switching to `multiplicative` is intended primarily for inter-model comparisons against CAM-wiso, LMDZ-iso, etc.
- The `aksmc/akrfa/akrfb` provenance was a small finding worth noting: an originally-undocumented set of constants in the C++ port turned out to encode the wind-speed-based surface flux kinetic. Updating the comments here closes that loose end without otherwise touching the surface-flux pathway.
