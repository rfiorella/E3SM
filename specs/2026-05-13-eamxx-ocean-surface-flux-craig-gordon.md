---
spec_id: 2026-05-13-eamxx-ocean-surface-flux-craig-gordon
spec_type: model-e3sm
spec_version: 1.0
title: Ocean surface flux of water isotopologues via Craig-Gordon kinetic (atmospheric side, in-EAMxx)
created: 2026-05-13
author: rfiorella
project: EAMXX-wiso
work_summary: Add a surface_flux_hook to the water-tracer framework and an in-EAMxx Craig-Gordon implementation that computes per-isotopologue ocean surface fluxes locally (using AlphaKMol + AlphaEqLiquidVapor + R_VSMOW), with a passive-VSMOW stopgap for land/sea-ice fractions.
priority: normal
estimated_effort_hours: 8

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
    description: |
      Source of α functions used by the Craig-Gordon implementation:
      - AlphaEqLiquidVapor(iso, T) — equilibrium liquid-vapor fractionation at SST
      - AlphaKMol(iso, rbot, zbot, ustar) — Brutsaert/Merlivat-Jouzel 1979 kinetic for surface layer (already implemented; used as-is)
      - R_VSMOW values from WaterIsotopologues::rnat (= {1.0, 0.9976, 155.76e-6, 2005.2e-6, 402e-6, 77.88e-6} for H2O, H216O, HDO, H218O, H217O, HTO)
      - boce array (= 1.0 everywhere; "Mean ocean surface enrichment relative to VSMOW") — multiplied with R_VSMOW for R_ocean
    format: C++ header
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_tracer_hooks.hpp
    description: Hook function-pointer table. This spec adds a new SurfaceFluxHookFn signature and surface_flux_hook function pointer alongside the existing phase-change hooks. No-op default added.
    format: C++ header
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_tracer_registry.hpp
    description: Provides tracer_is_tag(i) and tracer_isotopologue(i) used to dispatch tag (passive-copy) vs isotopologue (Craig-Gordon-fractionated) paths at the surface.
    format: C++ header
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/eamxx_water_tracers.cpp
    description: |
      Implementation file. This spec adds ocean_craig_gordon_surface_flux() and updates initialize_water_tracer_hooks() to register it. CMake-gated via #ifdef on SCREAM_WISO_SURFACE_FLUX_MODE_OCEAN_CRAIG_GORDON.
    format: C++ source
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    description: |
      Imports surf_evap (line 44, registered as Computed kg/m²/s, scalar2d). After import, this spec calls surface_flux_hook to distribute the bulk flux to per-isotopologue surface flux fields. The importer already has access to surf_evap, sst, wind_speed_10m, ocnfrac, landfrac, icefrac, fv, ram1 — all the inputs we need. ASK BEFORE editing.
    format: C++ source
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
    description: |
      Consumes surf_evap and applies it to the vapor field at the lowest model level (vapor_flux(i) = surf_evap(i) at eamxx_shoc_process_interface.hpp:310). For per-isotopologue tracers, the analogous tendency application happens through the new per-tracer surf_evap fields plumbed via the field manager. May require small modification to apply per-tracer surface flux at lowest level. ASK BEFORE editing.
    format: C++ source
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/CMakeLists.txt
    description: Add SCREAM_WISO_SURFACE_FLUX_MODE string option (values none | ocean_craig_gordon, default ocean_craig_gordon) and propagate via compile definition. ASK BEFORE editing.
    format: CMake
    required: true

# Deliverables
deliverables:
  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_tracer_hooks.hpp
    description: |
      Add SurfaceFluxHookFn typedef and surface_flux_hook function pointer.

      Signature:
        using SurfaceFluxHookFn = void(*)(
          const view_3d& qv,               // (COL, CMP, LEV) — atmospheric water vapor tracers
          const view_2d& surf_evap_iso,    // (COL, CMP) — OUTPUT: per-tracer surface flux to apply
          const view_1d& surf_evap_bulk,   // (COL) — bulk surface evap (kg/m²/s); convention surf_evap > 0 upward
          const view_1d& sst,              // (COL) — ocean surface temperature (K)
          const view_1d& ustar,            // (COL) — friction velocity (m/s); derivable from fv/ram1 if not direct
          const view_1d& z_bot,            // (COL) — height of lowest model level (m)
          const view_1d& rho_bot,          // (COL) — density at lowest level (kg/m³)
          const view_1d& h,                // (COL) — RH over water just above surface (qv_lowest / qsat_water(T_lowest))
          const view_1d& ocnfrac,          // (COL) — ocean fraction (0–1)
          int ncol, int nlev
        );

      Add no-op default implementation. Sublimation_hook, evaporation_hook, deposition_hook, etc. from prior specs are unchanged.
    format: C++ header
    validation_method: compiles in default and SCREAM_TRACE_WATER=ON builds

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/eamxx_water_tracers.cpp
    description: |
      Implement ocean_craig_gordon_surface_flux().

      Per-column loop. For each tracer i (CMP > 0):
        if WTRC_MAX_CNST == 1: return
        ocnfrac_col = ocnfrac(icol)
        bulk_flux   = surf_evap_bulk(icol)
        if tag(i):
          # Passive: tag surface flux = bulk flux (no ratio modification)
          surf_evap_iso(icol, i) = bulk_flux
        else:
          iso_idx  = tracer_isotopologue(i)
          iso_name = WaterIsotopologueNames[iso_idx]
          # OCEAN portion — Craig-Gordon kinetic
          if ocnfrac_col > 0:
            alpha_eq  = AlphaEqLiquidVapor(iso_name, sst(icol))
            alpha_kin = AlphaKMol(iso_name, rho_bot(icol), z_bot(icol), ustar(icol))
            R_ocean   = rnat(iso_idx) * boce(iso_idx)
            qv_bulk_lowest = qv(icol, 0, nlev-1)
            qv_iso_lowest  = qv(icol, i, nlev-1)
            R_a       = qv_iso_lowest / max(qv_bulk_lowest, q_min)
            h_col     = max(min(h(icol), 0.999), 0.0)     # bound h to (0, 1) for the formula
            R_E_ocean = alpha_kin * (R_ocean / alpha_eq - h_col * R_a) / (1.0 - h_col)
            flux_ocean = R_E_ocean * bulk_flux * ocnfrac_col
          else:
            flux_ocean = 0.0
          # NON-OCEAN portion (land + sea ice combined) — passive VSMOW stopgap (Q5 option b)
          R_nonocean = rnat(iso_idx)
          flux_nonocean = R_nonocean * bulk_flux * (1.0 - ocnfrac_col)
          # Combined
          surf_evap_iso(icol, i) = flux_ocean + flux_nonocean

      Update initialize_water_tracer_hooks() to register surface_flux_hook under #ifdef on SCREAM_WISO_SURFACE_FLUX_MODE_OCEAN_CRAIG_GORDON. When that macro is not defined (mode == none), the no-op default stays installed and no per-isotopologue surface flux is applied.
    format: C++ source
    validation_method: compiles for both CMake mode values; unit tests pass

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    description: |
      After importing surf_evap (existing line 44), add per-tracer surf_evap_iso fields to the field manager (one per CMP > 0; sized COL×CMP) and call surface_flux_hook to populate them. The downstream consumer (SHOC, via eamxx_shoc_process_interface.cpp) then applies these tendencies to per-isotopologue qv at the lowest level analogous to its existing handling of surf_evap.

      No changes to the importer's handling of surf_evap itself — the new code is additive (computes per-tracer fluxes alongside).

      ASK BEFORE editing.
    format: C++ source
    validation_method: compiles; new per-tracer fields appear in field-manager listings; SHOC sees them at expected dimensions

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
    description: |
      Extend the existing surf_evap → vapor_flux(i) at lowest level (line 310 of .hpp) to additionally apply per-tracer surf_evap_iso fluxes to qv(icol, i, nlev-1) for i in [1, WTRC_MAX_CNST). Use the same dt × gravit × rpdel conversion. Inherit the existing safety clipping for negativity.

      ASK BEFORE editing.
    format: C++ source
    validation_method: compiles; per-tracer surface tendencies applied correctly at lowest level; mass conservation respected

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/CMakeLists.txt
    description: |
      Add SCREAM_WISO_SURFACE_FLUX_MODE cache string option:
        - Values: none | ocean_craig_gordon (default: ocean_craig_gordon)
        - Configure-time validation rejects invalid values with clear error
        - When value == "ocean_craig_gordon": defines SCREAM_WISO_SURFACE_FLUX_MODE_OCEAN_CRAIG_GORDON compile macro
        - When value == "none": no macro defined; no-op surface_flux_hook stays installed
      ASK BEFORE editing.
    format: CMake
    validation_method: cmake -L lists the option; invalid value rejected at configure

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/tests/test_ocean_surface_flux.cpp
    description: |
      Catch2 unit tests covering:
        - α=1 limit: H2O and H216O tracers (rnat = 1.0 and 0.9976) and tags produce R_E ≈ 1.0 → surf_evap_iso ≈ R_ocean × bulk_flux (within FP)
        - HDO Craig-Gordon hand-calc at (SST=298.15K, ustar=0.3 m/s, z_bot=20m, rho_bot=1.18 kg/m³, ocnfrac=1.0, h=0.7, R_a=R_VSMOW(HDO))
        - H218O Craig-Gordon hand-calc at the same inputs
        - Ocean-fraction scaling: ocnfrac=0 → flux for non-tag = R_VSMOW × bulk × 1.0 (all non-ocean stopgap); ocnfrac=0.5 → half-and-half blend; ocnfrac=1 → all Craig-Gordon
        - Mass conservation: sum of per-tracer surf_evap_iso over CMP > 0 stays bounded (tracer-by-tracer; bulk is separately conserved by SHOC's existing path)
        - h bounding: h = 1.0 input is clamped to 0.999 internally to avoid the 1/(1-h) singularity; result is finite
        - q_min protection: qv_bulk_lowest near 0 doesn't produce NaN in R_a
        - Tags: passive copy independent of any ratio (tag initialized at 0.5×bulk drifts toward 1.0 after a sequence of surface flux applications — documented expected drift)
        - none-mode invariant: when SCREAM_WISO_SURFACE_FLUX_MODE_OCEAN_CRAIG_GORDON is not defined, surface_flux_hook is the no-op default and surf_evap_iso fields are not modified
    format: C++ source
    validation_method: compiles and all REQUIRE statements pass under both CMake mode values

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/tests/CMakeLists.txt
    description: Register test_ocean_surface_flux with CTest.
    format: CMake
    validation_method: test appears in test-all-eamxx output

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/REGISTRY_README.md
    description: |
      Document the new surface_flux_hook framework. Note the per-component (ocean / land / sea-ice) future split, the tag-attribution future spec, and the in-EAMxx-only nature of this spec (CIME-coupled path via eamxx_cpl_indices.F90 is separate future work).
    format: Markdown
    validation_method: manual review

# Success criteria
success_criteria:
  - id: SC1
    phase: planning
    description: Hard prerequisites — verify both prior specs (equilibrium-fractionation-hooks and kinetic-fractionation-hooks) are implemented and merged before starting this spec's implementation phase.
    criterion_type: human-review
    reviewer: rfiorella
    blocking: true

  - id: SC2
    phase: implementation
    description: Default build (SCREAM_TRACE_WATER=OFF) compiles unchanged
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | tee build_default.log && \
      ! grep -E "error:|Error " build_default.log
    expected_output: Build completes; no errors
    blocking: true

  - id: SC3
    phase: implementation
    description: Build under SCREAM_TRACE_WATER=ON with N=4 registry and SCREAM_WISO_SURFACE_FLUX_MODE=ocean_craig_gordon (default)
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      cmake -DSCREAM_TRACE_WATER=ON \
            -DSCREAM_WATER_TRACERS_FILE=cmake/water_tracers/registry_n4.cmake \
            -DSCREAM_WISO_SURFACE_FLUX_MODE=ocean_craig_gordon \
            -S . -B build_sf_cg && \
      cmake --build build_sf_cg --target water_tracers 2>&1 | tee build_sf_cg.log
    expected_output: libwater_tracers.a builds; ocean_craig_gordon_surface_flux symbol bound
    blocking: true

  - id: SC4
    phase: implementation
    description: Build with SCREAM_WISO_SURFACE_FLUX_MODE=none — no-op hook active
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      cmake -DSCREAM_TRACE_WATER=ON \
            -DSCREAM_WATER_TRACERS_FILE=cmake/water_tracers/registry_n4.cmake \
            -DSCREAM_WISO_SURFACE_FLUX_MODE=none \
            -S . -B build_sf_none && \
      cmake --build build_sf_none --target water_tracers 2>&1 | tee build_sf_none.log
    expected_output: libwater_tracers.a builds; no-op default surface_flux_hook installed
    blocking: true

  - id: SC5
    phase: implementation
    description: Invalid SCREAM_WISO_SURFACE_FLUX_MODE value produces clear configure-time error listing valid values
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ! cmake -DSCREAM_TRACE_WATER=ON \
              -DSCREAM_WATER_TRACERS_FILE=cmake/water_tracers/registry_n4.cmake \
              -DSCREAM_WISO_SURFACE_FLUX_MODE=invalid \
              -S . -B build_bad 2>&1 | tee bad.log && \
      grep -E "SCREAM_WISO_SURFACE_FLUX_MODE.*(none|ocean_craig_gordon)" bad.log
    expected_output: CMake fails; error message lists the two valid values
    blocking: true

  - id: SC6
    phase: testing
    description: α=1 limit — H2O and H216O tracers and tags produce surf_evap_iso equal to bulk surf_evap × R_ocean (where R_ocean = rnat(iso)*boce(iso) = rnat(iso) since boce=1 in this spec)
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_ocean_surface_flux.*alpha_one_limit.*PASSED"
    tolerance: absolute 1e-14
    assertion: For H2O (idx 0) and H216O (idx 1) and tags, |surf_evap_iso - rnat(iso) * surf_evap_bulk| < 1e-14 at any (SST, ustar, ocnfrac, h, R_a)
    blocking: true

  - id: SC7
    phase: testing
    description: HDO Craig-Gordon at canonical inputs (SST=298.15K, ustar=0.3 m/s, z_bot=20m, rho_bot=1.18 kg/m³, ocnfrac=1, h=0.7, R_a=R_VSMOW(HDO)) matches hand-calculated R_E from textbook Craig-Gordon
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_ocean_surface_flux.*hdo_cg_canonical.*PASSED"
    tolerance: absolute 1e-10
    reference_values: Hand-calculated R_E = α_kin × (R_VSMOW(HDO) / α_eq − h × R_a) / (1 − h) with α_eq = AlphaEqLiquidVapor(HDO, 298.15K) from Horita & Wesolowski (1994) and α_kin = AlphaKMol(HDO, 1.18, 20, 0.3) from the existing implementation.
    blocking: true

  - id: SC8
    phase: testing
    description: H218O Craig-Gordon at same canonical inputs matches hand-calculation
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_ocean_surface_flux.*h218o_cg_canonical.*PASSED"
    tolerance: absolute 1e-10
    reference_values: Hand-calculated R_E using HOH18O equilibrium α and HOH18O kinetic from the same canonical state.
    blocking: true

  - id: SC9
    phase: testing
    description: Ocean fraction scaling — ocnfrac sweep through {0.0, 0.25, 0.5, 0.75, 1.0} produces correctly blended fluxes; at ocnfrac=0 the result equals R_VSMOW(iso) × bulk_flux (all non-ocean stopgap); at ocnfrac=1 the result equals R_E_CG × bulk_flux (all ocean Craig-Gordon)
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_ocean_surface_flux.*ocnfrac_blend.*PASSED"
    tolerance: absolute 1e-12
    assertion: surf_evap_iso(ocnfrac) = ocnfrac × R_E_CG × bulk + (1-ocnfrac) × R_VSMOW × bulk for each isotopologue and each ocnfrac value
    blocking: true

  - id: SC10
    phase: testing
    description: h-bounding edge case — input h = 1.0 (saturated) is clamped to 0.999 internally; result is finite and within physical bounds
    criterion_type: assertion
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_ocean_surface_flux.*h_clamp.*PASSED"
    assertion: With h = 1.0 input, surf_evap_iso is finite (no NaN, no Inf) and matches the h=0.999 case within tolerance 1e-8
    blocking: true

  - id: SC11
    phase: testing
    description: q_min protection — qv_bulk_lowest near zero (q_bulk = 0 or 1e-25) does not produce NaN in R_a; surf_evap_iso remains finite
    criterion_type: assertion
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_ocean_surface_flux.*qmin_protection.*PASSED"
    assertion: With qv_bulk(icol, 0, nlev-1) in {0.0, 1e-25}, all surf_evap_iso outputs are finite
    blocking: true

  - id: SC12
    phase: testing
    description: Tag passive-copy behavior — tag surf_evap = bulk surf_evap regardless of any ratio; documented drift toward bulk-equivalent over many applications
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_ocean_surface_flux.*tag_passive.*PASSED"
    tolerance: absolute 1e-14
    assertion: For tag tracers, surf_evap_iso = surf_evap_bulk (no ratio multiplication, no Craig-Gordon)
    blocking: true

  - id: SC13
    phase: testing
    description: none-mode invariant — when SCREAM_WISO_SURFACE_FLUX_MODE=none, the surface_flux_hook is the no-op default and surf_evap_iso fields remain unchanged from their initial (zero) state
    criterion_type: assertion
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -E "test_ocean_surface_flux.*none_mode_noop.*PASSED"
    assertion: In the none-mode build, surf_evap_iso is bit-identical to its pre-hook (zero) value after the hook call
    blocking: true

  - id: SC14
    phase: integration
    description: |
      Human review of: (a) the new surface_flux_hook signature; (b) the atmosphere_surface_coupling_importer.cpp edits; (c) the SHOC consumer edits; (d) CMake option implementation; (e) REGISTRY_README.md updates; (f) confirmation that eamxx_cpl_indices.F90 is untouched (CIME-coupled path is future work); (g) confirmation that AlphaKMol is used as-is with no modification.
    criterion_type: human-review
    reviewer: rfiorella
    blocking: true

# Out of scope
out_of_scope:
  - Land-surface isotope fractionation (evapotranspiration, soil water, plant transpiration) — handled by the passive-VSMOW stopgap on (1 − ocnfrac); proper land surface treatment is a future spec
  - Sea-ice surface isotope fractionation (sublimation/melt of sea ice with its own isotope ratio) — also stopgap'd to passive-VSMOW within (1 − ocnfrac); proper sea-ice fractionation is a future spec
  - CIME-coupled F-compset isotope plumbing — extending eamxx_cpl_indices.F90 to declare HDO/H218O/H217O/HTO import/export fields, coupler-side isotope flux provision, ocean model isotope flux computation. This spec computes ocean fluxes locally in EAMxx instead.
  - Modifying AlphaKMol itself — used as-is with no functional changes
  - Modifying the existing surf_evap_h216o coupler plumbing — preserved untouched as the future CIME pathway
  - Tag-source attribution semantics (e.g., tag #1 = "Tropical-Pacific ocean source") — tags get pure passive-copy at the surface in this spec; future spec defines tag region/source schemes
  - Boce factor (mean ocean enrichment relative to VSMOW) tuning — the existing values (1.0 everywhere) are used; future spec may regionalize
  - Restart/checkpoint handling for the new per-tracer surf_evap_iso fields
  - Tier 1 standalone-DP runs exercising this surface-flux machinery (separate spec — Tier 1 DP-SCREAM BOMEX/RICO/MPACE runs)
  - Performance optimization or Kokkos-kernel refactoring of the surface-flux hook body

# Resolved decisions
resolved_decisions:
  - decision: Implement ocean surface flux fractionation locally inside EAMxx (not via coupler-provided isotope fluxes), using the existing AlphaKMol + AlphaEqLiquidVapor functions with R_VSMOW as the ocean source ratio.
    rationale: User direction (Q3 from refinement round). Required for standalone DP-SCREAM Tier 1 runs which do not have a real coupled ocean providing isotope fluxes. The CIME-coupled pathway (extending eamxx_cpl_indices.F90 + ocean-side isotope physics) is a separate future spec.
    date: 2026-05-13

  - decision: New hook in the existing water_tracer_hooks framework, named surface_flux_hook, with function-pointer + no-op-default + #ifdef-gated implementation pattern matching the existing phase-change hooks.
    rationale: Architectural consistency. Keeps all isotope-physics dispatch in one place; surface flux is logically a "where mass enters" hook analogous to the "where mass changes phase" hooks. User direction (Q1 of refinement round).
    date: 2026-05-13

  - decision: Focus on OCEAN fraction for proper Craig-Gordon fractionation; apply a passive-VSMOW stopgap (R_ocean × surf_evap × (1 − ocnfrac)) for the non-ocean fraction. Land and sea-ice get the same VSMOW treatment as a placeholder.
    rationale: User direction (Q5 option b). Avoids a small but real mass-budget artifact at the surface boundary — total per-isotopologue flux summed across all surfaces equals R_VSMOW × bulk_flux to leading order (with Craig-Gordon refinement only on the ocean fraction). Land and sea-ice get their own proper treatment in future specs.
    date: 2026-05-13

  - decision: Textbook Craig-Gordon formulation — R_E = α_kin × (R_ocean / α_eq − h × R_a) / (1 − h) — with α_eq from AlphaEqLiquidVapor(iso, SST), α_kin from AlphaKMol, R_ocean from rnat × boce, R_a from the tracer ratio at the lowest model level.
    rationale: User direction (Q1 of refinement round). Stewart-style effective humidity was for hydrometeor evaporation where droplet T is unknown; for the surface flux we have SST and atmospheric state directly, so textbook Craig-Gordon is the natural fit.
    date: 2026-05-13

  - decision: R_a (ambient atmospheric vapor ratio just above surface) is computed from the tracer field at the lowest model level: R_a = q_iso(icol, i, nlev-1) / max(q_bulk(icol, nlev-1), q_min).
    rationale: User direction (Q2 of refinement round). Self-consistent with the atmospheric state being simulated.
    date: 2026-05-13

  - decision: surf_evap sign convention: positive upward (atmosphere gaining water). The hook applies +R_E × surf_evap to the per-tracer fields, matching the existing bulk surf_evap convention in EAMxx.
    rationale: User direction (Q3 of refinement round). No special sign handling; matches the SHOC ingestion path.
    date: 2026-05-13

  - decision: Tags get pure passive copy at the surface — surf_evap_tag = surf_evap_bulk (no R multiplication, no Craig-Gordon). Document that this causes tag ratio to drift toward 1.0 over a long run as fresh surface flux dilutes tags initialized at other ratios.
    rationale: User direction (Q4 of refinement round). Tag source-attribution semantics are intentionally deferred to a future spec; passive copy at surface is the simplest defensible default in the interim.
    date: 2026-05-13

  - decision: h (RH over water just above the surface) is clamped to [0, 0.999] before entering the Craig-Gordon denominator to avoid the (1 − h) singularity at saturation. Physically, surface evaporation goes to zero at h = 1, so the clamp does not produce a non-physical flux.
    rationale: Numerical safety. The bulk surf_evap is delivered by the coupler and may be small but nonzero even at near-saturation; without the clamp R_E diverges. h = 0.999 cap is a conservative choice; unit test SC10 verifies finiteness.
    date: 2026-05-13

  - decision: CMake option SCREAM_WISO_SURFACE_FLUX_MODE — string, values {none, ocean_craig_gordon}, default ocean_craig_gordon. The `none` mode is kept for development/debug — useful for isolating whether observed surface flux behavior is the hook's contribution.
    rationale: User direction (Q5a of refinement round). String option self-documents in cmake -L; default to physically meaningful mode.
    date: 2026-05-13

  - decision: Call site for surface_flux_hook is atmosphere_surface_coupling_importer.cpp (control layer), after surf_evap is imported and before downstream consumers see per-tracer fluxes. Per-tracer surf_evap_iso fields are added to the field manager so SHOC can consume them analogous to surf_evap.
    rationale: Investigation finding. The importer is the gateway for all surface fluxes from the coupler and already has access to surf_evap, sst, wind_speed_10m, ocnfrac, landfrac, icefrac, fv, ram1 — no field-manager additions for inputs are needed beyond the per-tracer output fields.
    date: 2026-05-13

  - decision: eamxx_cpl_indices.F90 is NOT modified in this spec. The existing surf_evap_h216o pattern is preserved untouched; HDO/H218O/H217O/HTO are not added to the coupler import/export lists. The CIME-coupled pathway is future work.
    rationale: User direction. Standalone DP-SCREAM runs are the target; the coupler pathway is meaningful only for full F-compset runs which require additional coupler-side work in the ocean model.
    date: 2026-05-13

  - decision: AlphaKMol is used as-is. The existing implementation (Brutsaert/Merlivat-Jouzel 1979 surface-layer closure using rbot, zbot, ustar, with Charnock z0 and smooth/rough regime branching at the critical Reynolds number) is the kinetic component of the Craig-Gordon formula.
    rationale: Pre-existing implementation; no modification needed. Q10 of prior kinetic spec already declared AlphaKMol out of scope for modification.
    date: 2026-05-13

  - decision: Both prior specs (equilibrium-fractionation-hooks and kinetic-fractionation-hooks) are hard prerequisites. SC1 is a planning-phase human-review checkpoint that verifies both are merged before this spec's implementation phase begins.
    rationale: User direction (Q0). This spec's surface flux uses AlphaEqLiquidVapor which lands with the equilibrium spec; AlphaKMol predates both but the hooks framework infrastructure (water_tracer_hooks.hpp, registry) is built up incrementally.
    date: 2026-05-13

  - decision: Validation tier 0 — compile + unit tests. No 5-day run, no end-to-end DP-SCREAM run in this spec. Tier 1 DP-SCREAM BOMEX/RICO/MPACE runs are the next spec (3b) and depend on this one.
    rationale: User direction (Q6 of refinement round). Keeps scope focused; end-to-end runs go to spec 3b.
    date: 2026-05-13

# Ask-before actions (project-specific additions to global policy)
ask_before:
  - Modifying components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp (adding the surface_flux_hook call after surf_evap import)
  - Modifying components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp (extending the lowest-level qv-tracer tendency application)
  - Adding new fields to the EAMxx field manager (per-tracer surf_evap_iso)
  - Adding new entries to components/eamxx/src/share/util/io_metadata/io_metadata.csv
  - Modifying components/eamxx/CMakeLists.txt (SCREAM_WISO_SURFACE_FLUX_MODE option)
  - Modifying or removing AlphaKMol in water_isotopes.hpp (must remain unchanged)
  - Modifying components/eamxx/src/mct_coupling/eamxx_cpl_indices.F90 (CIME-coupled path; out of scope)
  - Modifying or removing the existing surf_evap_h216o coupler plumbing
  - Changing the boce or rnat numerical values in WaterIsotopologues
  - Restructuring the existing hook function-pointer table beyond the documented signature additions (surface_flux_hook addition only)
  - Adding new external dependencies or libraries

# Parallelization
allow_parallelization: false

# Post-completion review
request_performance_review: false
request_code_review: false
---

# Ocean Surface Flux of Water Isotopologues via Craig-Gordon Kinetic (Atmospheric Side, In-EAMxx)

## Context

The two prior specs — `2026-05-13-eamxx-equilibrium-fractionation-hooks` and `2026-05-13-eamxx-kinetic-fractionation-hooks` — wired isotope fractionation into the atmospheric phase-change pathways (condensation, evaporation, deposition, sublimation, melting, freezing). Both explicitly deferred the **surface-flux pathway** (water entering the atmosphere from ocean evaporation, land evapotranspiration, sea-ice sublimation) to a future spec on coupler-level integration.

Without a surface flux, multi-hour DP-SCREAM runs are interpretively limited: precipitation removes isotopologue mass from the atmosphere with no resupply, so the simulation depletes toward zero. To make Tier 1 runs meaningful, this spec adds an **in-EAMxx** surface flux for water isotopologues, scoped to the ocean fraction only and using the existing `AlphaKMol` (Brutsaert/Merlivat-Jouzel 1979) and `AlphaEqLiquidVapor` (Horita & Wesolowski 1994) functions.

### Investigation findings

Searching `components/eamxx/src/` revealed two relevant pieces of pre-existing infrastructure:

1. **Bulk surface flux pathway is wired.** `atmosphere_surface_coupling_importer.cpp` imports `surf_evap` (kg/m²/s, positive upward) from the coupler at line 44, and SHOC consumes it via `vapor_flux(i) = surf_evap(i)` at `eamxx_shoc_process_interface.hpp:310`. Mass enters the atmosphere at the lowest model level through this path.

2. **Partial isotope coupling exists, for ONE tracer only.** `eamxx_cpl_indices.F90:101` declares `surf_evap_h216o` (mapped to coupler's `Faxx_evap_h216o`), and lines 194–195 export `Faxa_rainl_h216o` / `Faxa_snowl_h216o`. H216O is the "labeled bulk" tracer (α=1 by construction) — only useful as a sanity check, not real isotope fractionation. The other isotopologues (HDO, H218O, H217O, HTO) are not plumbed through the coupler. Extending the coupler infrastructure to include them would require ocean-model-side work (the ocean must provide the per-isotopologue fluxes), which is out of scope here.

The pragmatic path for DP-SCREAM Tier 1 runs is therefore: **compute the per-isotopologue ocean fluxes locally inside EAMxx**, using the bulk `surf_evap` and the already-imported ocean state (`sst`, `wind_speed_10m`, `ocnfrac`, `fv`, `ram1`).

### Scope

- New `surface_flux_hook` in the existing water-tracer hook framework
- `ocean_craig_gordon_surface_flux` implementation for the **ocean fraction**, using textbook Craig-Gordon (Q1 of refinement round)
- **Passive-VSMOW stopgap** for the non-ocean fraction (land + sea ice combined; per Q5 option b) to avoid a mass-budget artifact while proper land and sea-ice fractionation are deferred
- **Tags** get pure passive copy at the surface — `surf_evap_tag = surf_evap_bulk` — with documented drift toward bulk-equivalent over many surface-flux applications (Q4)
- CMake option `SCREAM_WISO_SURFACE_FLUX_MODE` with values `none` (development) and `ocean_craig_gordon` (default)
- Validation tier 0 — compile + unit tests
- The CIME-coupled pathway (extending `eamxx_cpl_indices.F90` to include HDO/H218O/H217O/HTO + ocean-side isotope physics) is **out of scope** and is a separate future spec

## Approach

### Hook signature

```cpp
// water_tracer_hooks.hpp — new hook
using SurfaceFluxHookFn = void(*)(
  const view_3d& qv,               // (COL, CMP, LEV) — water vapor tracer fields
  const view_2d& surf_evap_iso,    // (COL, CMP) — OUTPUT: per-tracer surface flux
  const view_1d& surf_evap_bulk,   // (COL) — bulk surface evap (kg/m²/s), > 0 upward
  const view_1d& sst,              // (COL) — sea surface temperature (K)
  const view_1d& ustar,            // (COL) — friction velocity (m/s)
  const view_1d& z_bot,            // (COL) — height of lowest model level (m)
  const view_1d& rho_bot,          // (COL) — density at lowest level (kg/m³)
  const view_1d& h,                // (COL) — RH over water near surface
  const view_1d& ocnfrac,          // (COL) — ocean fraction in [0, 1]
  int ncol, int nlev
);

inline SurfaceFluxHookFn surface_flux_hook = default_noop_surface_flux;

inline void default_noop_surface_flux(
  const view_3d&, const view_2d&, const view_1d&, const view_1d&,
  const view_1d&, const view_1d&, const view_1d&, const view_1d&,
  const view_1d&, int, int) {}
```

Tag and isotopologue dispatch happens *inside* the implementation, consistent with the existing phase-change hooks.

### Implementation — ocean_craig_gordon_surface_flux

For each column, for each CMP index `i > 0`:

```
if tag(i):                                          # tags: passive copy, no ratio modification
  surf_evap_iso(icol, i) = surf_evap_bulk(icol)
  continue

iso_idx  = tracer_isotopologue(i)
iso_name = WaterIsotopologueNames[iso_idx]

# --- OCEAN portion: Craig-Gordon kinetic ---
α_eq  = AlphaEqLiquidVapor(iso_name, sst(icol))                       # Horita-Wesolowski 1994
α_kin = AlphaKMol(iso_name, rho_bot(icol), z_bot(icol), ustar(icol))  # Merlivat-Jouzel 1979
R_ocean = rnat[iso_idx] * boce[iso_idx]                               # boce = 1 currently
qv_bulk_low = qv(icol, 0, nlev-1)
qv_iso_low  = qv(icol, i, nlev-1)
R_a = qv_iso_low / max(qv_bulk_low, q_min)
h_clamped = clamp(h(icol), 0.0, 0.999)                                # avoid 1/(1-h) singularity at h=1
R_E_ocean = α_kin * (R_ocean / α_eq - h_clamped * R_a) / (1.0 - h_clamped)
flux_ocean = R_E_ocean * surf_evap_bulk(icol) * ocnfrac(icol)

# --- NON-OCEAN portion: passive-VSMOW stopgap ---
R_nonocean = rnat[iso_idx]                                            # land + sea ice combined
flux_nonocean = R_nonocean * surf_evap_bulk(icol) * (1.0 - ocnfrac(icol))

# --- Combined ---
surf_evap_iso(icol, i) = flux_ocean + flux_nonocean
```

Sign convention follows `surf_evap` (positive upward = atmosphere gaining water).

### Call-site changes

**`atmosphere_surface_coupling_importer.cpp`** — after the existing `surf_evap` import block, register per-tracer `surf_evap_iso` fields with the field manager (dimensions COL × CMP, units kg/m²/s) and call `surface_flux_hook` to populate them. No changes to the existing `surf_evap` handling.

**`eamxx_shoc_process_interface.cpp`** — the existing `vapor_flux(i) = surf_evap(i)` block applies bulk `surf_evap` to `qv` at the lowest model level. Extend to additionally apply each `surf_evap_iso(icol, i)` for `i ∈ [1, WTRC_MAX_CNST)` to `qv(icol, i, nlev-1)` using the same `dt × gravit × rpdel` conversion. Preserve the existing safety clipping for negativity.

Both call-site files are outside `components/eamxx/src/physics/water_tracers/` and are flagged `ask-before` per project policy.

### CMake option

```cmake
set(SCREAM_WISO_SURFACE_FLUX_MODE "ocean_craig_gordon" CACHE STRING
    "Surface flux of water isotopologues: none | ocean_craig_gordon")
set_property(CACHE SCREAM_WISO_SURFACE_FLUX_MODE PROPERTY STRINGS
             "none" "ocean_craig_gordon")

if (NOT (SCREAM_WISO_SURFACE_FLUX_MODE STREQUAL "none" OR
         SCREAM_WISO_SURFACE_FLUX_MODE STREQUAL "ocean_craig_gordon"))
  message(FATAL_ERROR
    "SCREAM_WISO_SURFACE_FLUX_MODE must be one of: none, ocean_craig_gordon. "
    "Got: ${SCREAM_WISO_SURFACE_FLUX_MODE}")
endif()

if (SCREAM_WISO_SURFACE_FLUX_MODE STREQUAL "ocean_craig_gordon")
  target_compile_definitions(water_tracers PRIVATE SCREAM_WISO_SURFACE_FLUX_MODE_OCEAN_CRAIG_GORDON)
endif()
```

`initialize_water_tracer_hooks()` registers `ocean_craig_gordon_surface_flux` only when the macro is defined. In `none` mode the no-op default remains installed.

### Decomposition

1. `water_tracer_hooks.hpp` — add `SurfaceFluxHookFn` typedef, `surface_flux_hook` function pointer, no-op default
2. `eamxx_water_tracers.cpp` — implement `ocean_craig_gordon_surface_flux`; update `initialize_water_tracer_hooks()` with `#ifdef` dispatch
3. `CMakeLists.txt` — add `SCREAM_WISO_SURFACE_FLUX_MODE` option + validation + compile-definition propagation
4. **Importer call-site update (ASK BEFORE)** — `atmosphere_surface_coupling_importer.cpp`
5. **SHOC consumer update (ASK BEFORE)** — `eamxx_shoc_process_interface.cpp`
6. Unit test `test_ocean_surface_flux.cpp` — covers α=1 limit, HDO/H218O canonical points, ocnfrac blending, h-clamping, q_min protection, tag passive-copy, none-mode invariant
7. `REGISTRY_README.md` — document the new hook and the per-component (ocean / land / sea-ice / CIME-coupled) future-work split

### Risks

- **h-clamp at 0.999.** The clamp prevents the Craig-Gordon singularity at saturation but slightly biases R_E in the near-saturated regime. Quantitatively, the bias is small because the bulk `surf_evap` itself approaches zero at h → 1 (the surface model returns small evaporation rates), so the *flux* through the formula stays bounded. SC10 verifies finiteness; physical sensitivity to the clamp value is a future-work investigation.
- **R_a self-consistency in the first timestep.** At t=0, the atmospheric tracer field starts with whatever IC the run prescribes (typically `q_iso = R_VSMOW × q_bulk` per the future Tier 1 spec). R_a = R_VSMOW initially, which makes R_E_ocean = α_kin × R_VSMOW × (1/α_eq - h)/(1-h) — physically reasonable. After the first step, R_a evolves and Craig-Gordon adjusts accordingly. No special-casing needed.
- **AlphaKMol's `recrit` regime branching at the critical Reynolds number.** This produces a slope discontinuity in α_kin at the regime boundary. Unit tests at canonical inputs avoid the boundary, but the discontinuity could surface in long DP-SCREAM runs if wind speed straddles the critical value. Flagged for future-work investigation if observed.
- **Per-tracer `surf_evap_iso` field-manager registration.** The 2D field (COL × CMP) shape is non-standard if the field manager assumes scalar2d. May require adding a new 2D-with-CMP field type or registering one field per CMP index. Investigate during implementation; the choice does not affect physics.
- **SHOC's existing surf_evap clipping** (line 642–658 of `eamxx_shoc_process_interface.cpp` — "Total mass of column vapor should be greater than mass of surf_evap"). This safety check must be replicated per-tracer; if a tracer's column mass is less than its surf_evap_iso, the tracer flux must be clipped consistently with the bulk path. Implementation must mirror the existing pattern.

Relevant skills to load during implementation: `e3sm-build-and-test`, `e3sm-platforms`, `scientific-modeling-conventions`.

## References

- **Prior specs (this project):**
  - `2026-05-13-eamxx-equilibrium-fractionation-hooks.md` — provides AlphaEqLiquidVapor (used here for α_eq at SST)
  - `2026-05-13-eamxx-kinetic-fractionation-hooks.md` — sibling spec for atmospheric phase-change kinetic; AlphaKMol unchanged across both
  - `2026-05-06-eamxx-water-tracer-isotopologue-registry-COMPLETE.md` — registry API for tracer dispatch
- **Source files referenced:**
  - `components/eamxx/src/physics/water_tracers/water_isotopes.hpp` — AlphaEqLiquidVapor, AlphaKMol, rnat, boce
  - `components/eamxx/src/physics/water_tracers/water_tracer_hooks.hpp` — existing hook framework
  - `components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp` — surf_evap import (line 44), all required ocean state
  - `components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp` — surf_evap → vapor_flux at lowest level (line 310 of .hpp)
  - `components/eamxx/src/mct_coupling/eamxx_cpl_indices.F90` — existing surf_evap_h216o pattern (lines 101, 133) — NOT modified in this spec
- **Scientific references:**
  - Craig, H. & Gordon, L. I. (1965). Deuterium and oxygen 18 variations in the ocean and the marine atmosphere. In *Stable Isotopes in Oceanographic Studies and Paleotemperatures*, ed. E. Tongiorgi, 9–130.
  - Merlivat, L. & Jouzel, J. (1979). Global climatic interpretation of the deuterium-oxygen 18 relationship for precipitation. *JGR* 84, 5029–5033.
  - Horita, J. & Wesolowski, D.J. (1994). *Geochim. Cosmochim. Acta* 58, 3425–3437. (For α_eq at SST.)
  - Brutsaert, W. (1975). A theory for local evaporation (or heat transfer) from rough and smooth surfaces at ground level. *Water Resour. Res.* 11, 543–550. (Surface-layer formulation underpinning AlphaKMol.)

## Notes

- This spec is intentionally narrow: ocean only, in-EAMxx, Tier 0 validation. It establishes the surface-flux hook framework and a defensible default implementation; future specs add land (ET partitioning into transpiration vs soil evaporation), sea ice (sublimation with sea-ice isotope composition), and the CIME-coupled pathway with ocean-model-side isotope physics.
- Once this lands, the prerequisite chain for the Tier 1 DP-SCREAM runs spec (3b) is complete: equilibrium hooks + kinetic hooks + surface flux → meaningful multi-hour runs with quasi-steady-state atmospheric isotope distributions.
- The `none` CMake mode is genuinely useful for debugging: running with `none` for a few hours, observing the depletion pattern as precipitation removes isotopologues with no resupply, validates that the atmospheric phase-change hooks are doing what they should and isolates the surface flux's contribution.
- Tag drift toward 1.0 under passive-copy at the surface is a known limitation of this spec. Documenting it here closes a small loop: the tag-attribution spec will replace the passive-copy line with a tag-source-dependent ratio source (e.g., tag #1 = Tropical-Pacific water; surf_evap_tag1 = R_TropPac × surf_evap_bulk × tropical_pacific_mask).
