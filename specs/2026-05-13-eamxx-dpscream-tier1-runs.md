---
spec_id: 2026-05-13-eamxx-dpscream-tier1-runs
spec_type: model-e3sm
spec_version: 1.0
title: Tier 1 DP-SCREAM validation runs (dycomsrf01, arm97, comble) with the integrated water-isotope physics
created: 2026-05-13
author: rfiorella
project: EAMXX-wiso
work_summary: First end-to-end Tier 1 validation of the water-isotope physics by running three DP-SCREAM cases (dycomsrf01, arm97, comble) at 5x5 / 3km with N=4 tracer registry, verifying bulk-water BFB, per-isotopologue mass conservation, δ-value plausibility, and tag drift behavior.
priority: normal
estimated_effort_hours: 16

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
  build_mode: eamxx-cime
  target_compset: F2010-SCREAMv1-DP
  target_resolution: ne4pg2_ne4pg2     # placeholder; actual DP testmods set their own resolution via shell_commands
  platform: docker-local
  container_image: rfiorella/model-containers:e3sm-openmpi-dev-latest
  validation_tier: 1

# Inputs
inputs:
  - path: /code/E3SM/EAMXX-wiso/components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/dpxx/dycomsrf01
    description: |
      In-tree DP-SCREAM testmod for DYCOMS-II RF01 (marine stratocumulus). Exercises warm-phase shallow cloud physics: eq_condensation + kin_evaporation. Standard duration ~6h sim time. Substitutes for BOMEX in the originally-proposed test plan.
    format: shell_commands + namelist overrides
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/dpxx/arm97
    description: |
      In-tree DP-SCREAM testmod for ARM SGP 1997 IOP (continental deep convection). Exercises convection + precipitation pathway: eq_condensation, kin_evaporation (rain re-evap), kin_deposition (deep cloud anvils), eq_sublimation. Standard duration ~24h sim time. Substitutes for RICO.
    format: shell_commands + namelist overrides
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/dpxx/comble
    description: |
      In-tree DP-SCREAM testmod for COMBLE (Arctic cold-air-outbreak mixed-phase). Exercises mixed-phase cloud physics + WBF/supersaturation deposition: kin_deposition (S_i > 1), eq_sublimation, eq_freezing, eq_melting. Standard duration ~12h sim time. Substitutes for MPACE.
    format: shell_commands + namelist overrides
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/cmake/water_tracers/registry_n4.cmake
    description: |
      N=4 tracer registry: [bulk H2O (CMP=0), HDO (CMP=1), H218O (CMP=2), passive tag (CMP=3)]. Drives the SCREAM_WATER_TRACERS_FILE option at build configure.
    format: CMake config file
    required: true

  - path: /code/E3SM/EAMXX-wiso/specs/2026-05-13-eamxx-equilibrium-fractionation-hooks.md
    description: Prior spec — equilibrium fractionation hooks for condensation, deposition (later replaced by kin_deposition), melting, freezing. Hard prerequisite.
    format: spec
    required: true

  - path: /code/E3SM/EAMXX-wiso/specs/2026-05-13-eamxx-kinetic-fractionation-hooks.md
    description: Prior spec — kinetic fractionation for hydrometeor evaporation, WBF/supersaturation deposition (replaces eq_deposition), equilibrium sublimation. Hard prerequisite.
    format: spec
    required: true

  - path: /code/E3SM/EAMXX-wiso/specs/2026-05-13-eamxx-ocean-surface-flux-craig-gordon.md
    description: Prior spec — ocean surface flux of water isotopologues via Craig-Gordon kinetic, with passive-VSMOW stopgap for non-ocean fraction. Hard prerequisite.
    format: spec
    required: true

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/scripts/test-all-eamxx
    description: Existing test driver for standalone EAMxx. NOT used here; spec uses CIME create_test instead. Listed for context.
    format: bash script
    required: false

# Deliverables
deliverables:
  - path: /code/E3SM/EAMXX-wiso/components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/dpxx/dycomsrf01_wiso
    description: |
      New testmod directory derived from dycomsrf01, adding water-tracer namelist settings. shell_commands sets SCREAM_TRACE_WATER=ON and points SCREAM_WATER_TRACERS_FILE at registry_n4.cmake; output handle adds per-tracer qv/qc/qi/precip fields at hourly cadence. Three sibling testmods will be created: dycomsrf01_wiso, arm97_wiso, comble_wiso.
    format: shell_commands + namelist overrides + output yaml
    validation_method: CIME create_test recognizes the testmod; build succeeds; namelists include water-tracer settings

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/dpxx/arm97_wiso
    description: Sibling of dycomsrf01_wiso for the ARM97 case.
    format: shell_commands + namelist overrides + output yaml
    validation_method: As above

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/dpxx/comble_wiso
    description: Sibling of dycomsrf01_wiso for the COMBLE case.
    format: shell_commands + namelist overrides + output yaml
    validation_method: As above

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/scripts/run-wiso-dpscream-tier1
    description: |
      New bash driver script that:
        1. Verifies prerequisite specs are merged (greps for sentinel markers in source — eq_condensation symbol, kin_evaporation symbol, ocean_craig_gordon_surface_flux symbol)
        2. For each case in {dycomsrf01, arm97, comble}, runs three configurations: TRACE_WATER=OFF, TRACE_WATER=ON+N=1, TRACE_WATER=ON+N=4
        3. Uses CIME create_test with the existing F2010-SCREAMv1-DP compset and the per-case wiso testmod
        4. Stores outputs under /data/runs/wiso_dpscream_tier1/<case>/<config>/
        5. Stores baselines (first run) under /data/baselines/EAMXX-wiso/wiso_dpscream_tier1/
        6. Invokes the analysis script after all 9 runs complete (3 cases × 3 configs)
      Configurable via env vars: WISO_DPSCREAM_DURATION_OVERRIDE (override case duration in seconds; default = published-protocol duration); WISO_DPSCREAM_OUTPUT_FREQ_HOURS (default 1); WISO_DPSCREAM_REGEN_BASELINES (default 0).
    format: bash script
    validation_method: script runs all 9 case-config combinations without crash on the docker-local platform within 16-core / 128 GB budget

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/scripts/analyze-wiso-dpscream-tier1
    description: |
      New Python analysis script that loads each case's three netCDF outputs (TRACE_WATER=OFF, ON+N=1, ON+N=4) and computes:
        (1) Bulk-water BFB comparison across configs (CMP=0 in N>=1 configs vs the OFF baseline)
        (2) Column-integrated isotope mass over time, per tracer
        (3) δD and δ18O time series at multiple levels and column-integrated
        (4) Tag drift: column-mean q_tag/q_bulk vs time
        (5) Conservation residual: d(column isotope mass)/dt vs (surface flux - precipitation flux)
      Writes a per-case markdown report under /data/runs/wiso_dpscream_tier1/<case>/report.md and a top-level summary under /data/runs/wiso_dpscream_tier1/SUMMARY.md.
    format: Python script (xarray + numpy)
    validation_method: report files generated; numerical content matches success-criteria expectations

  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/REGISTRY_README.md
    description: |
      Append a "Tier 1 validation" section documenting the three cases, the configurations run, success criteria, and how to reproduce. Pointer to the driver and analysis scripts.
    format: Markdown
    validation_method: manual review

# Success criteria
success_criteria:
  - id: SC1
    phase: planning
    description: Hard prerequisites — verify all three prior specs are implemented and merged before starting this spec's implementation phase. The driver script's prerequisite check (greps for eq_condensation, kin_evaporation, ocean_craig_gordon_surface_flux symbols in the water_tracers source) must pass.
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers && \
      grep -q "eq_condensation" eamxx_water_tracers.cpp && \
      grep -q "kin_evaporation" eamxx_water_tracers.cpp && \
      grep -q "ocean_craig_gordon_surface_flux" eamxx_water_tracers.cpp
    expected_output: All three symbols present; exit code 0
    blocking: true

  - id: SC2
    phase: implementation
    description: Three new wiso testmod directories (dycomsrf01_wiso, arm97_wiso, comble_wiso) are recognized by CIME and contain the expected namelist additions
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/dpxx && \
      for c in dycomsrf01_wiso arm97_wiso comble_wiso; do \
        test -d "$c" && \
        grep -q "SCREAM_TRACE_WATER" "$c/shell_commands" && \
        grep -q "registry_n4" "$c/shell_commands" || exit 1; \
      done
    expected_output: All three testmods exist and contain water-tracer namelist settings
    blocking: true

  - id: SC3
    phase: implementation
    description: Driver script run-wiso-dpscream-tier1 exists, is executable, and the prerequisite-check pre-flight passes (does not actually run cases yet — just self-check)
    criterion_type: shell
    command: |
      test -x /code/E3SM/EAMXX-wiso/components/eamxx/scripts/run-wiso-dpscream-tier1 && \
      /code/E3SM/EAMXX-wiso/components/eamxx/scripts/run-wiso-dpscream-tier1 --dry-run --check-prereqs
    expected_output: Script exists, executable, dry-run prereq check passes
    blocking: true

  - id: SC4
    phase: implementation
    description: Analysis script analyze-wiso-dpscream-tier1 exists and runs on a synthetic test output (small fake netCDF) to verify it doesn't crash on its own framework
    criterion_type: shell
    command: |
      test -x /code/E3SM/EAMXX-wiso/components/eamxx/scripts/analyze-wiso-dpscream-tier1 && \
      /code/E3SM/EAMXX-wiso/components/eamxx/scripts/analyze-wiso-dpscream-tier1 --self-test
    expected_output: Script exists, executable, self-test passes
    blocking: true

  - id: SC5
    phase: testing
    description: |
      DYCOMSrf01 runs to completion for all three configurations (TRACE_WATER=OFF, ON+N=1, ON+N=4) within the 16-core / 128 GB budget. Standard ~6h sim time.
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx/scripts && \
      ./run-wiso-dpscream-tier1 --case dycomsrf01 2>&1 | tee /data/runs/wiso_dpscream_tier1/dycomsrf01.log && \
      grep -q "RUN COMPLETED" /data/runs/wiso_dpscream_tier1/dycomsrf01.log
    expected_output: All three configs complete; netCDF outputs exist under /data/runs/wiso_dpscream_tier1/dycomsrf01/
    blocking: true

  - id: SC6
    phase: testing
    description: ARM97 runs to completion for all three configurations. Standard ~24h sim time.
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx/scripts && \
      ./run-wiso-dpscream-tier1 --case arm97 2>&1 | tee /data/runs/wiso_dpscream_tier1/arm97.log && \
      grep -q "RUN COMPLETED" /data/runs/wiso_dpscream_tier1/arm97.log
    expected_output: All three configs complete; netCDF outputs exist
    blocking: true

  - id: SC7
    phase: testing
    description: COMBLE runs to completion for all three configurations. Standard ~12h sim time.
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx/scripts && \
      ./run-wiso-dpscream-tier1 --case comble 2>&1 | tee /data/runs/wiso_dpscream_tier1/comble.log && \
      grep -q "RUN COMPLETED" /data/runs/wiso_dpscream_tier1/comble.log
    expected_output: All three configs complete; netCDF outputs exist
    blocking: true

  - id: SC8
    phase: testing
    description: |
      Bulk-water BFB across all three configurations for each case. The bulk water (CMP=0 in N>=1) is bit-for-bit identical to the SCREAM_TRACE_WATER=OFF baseline. This must hold at every output time step and every grid cell.
    criterion_type: shell
    command: |
      /code/E3SM/EAMXX-wiso/components/eamxx/scripts/analyze-wiso-dpscream-tier1 --check bfb_bulk 2>&1 | tee bfb.log && \
      grep -q "BFB_BULK: PASS for dycomsrf01 arm97 comble" bfb.log
    expected_output: BFB pass for all three cases
    blocking: true

  - id: SC9
    phase: testing
    description: |
      Per-isotopologue column-integrated mass conservation. For each isotopologue tracer i and each case, the residual |d(M_iso)/dt - (F_surface_iso - F_precip_iso)| relative to |M_iso| stays below 1e-10 throughout the run.
    criterion_type: tolerance
    command: |
      /code/E3SM/EAMXX-wiso/components/eamxx/scripts/analyze-wiso-dpscream-tier1 --check mass_conservation 2>&1 | tee mc.log && \
      grep -q "MASS_CONSERVATION: PASS" mc.log
    tolerance: relative 1e-10
    reference_values: Computed in-script from netCDF output time series of column-integrated isotope mass, surface flux, and precipitation flux.
    blocking: true

  - id: SC10
    phase: testing
    description: |
      δ-value plausibility bounds across all three cases. δD ∈ [-300, +50] ‰ and δ18O ∈ [-50, +10] ‰ in vapor and condensate at all output times and all model levels.
    criterion_type: tolerance
    command: |
      /code/E3SM/EAMXX-wiso/components/eamxx/scripts/analyze-wiso-dpscream-tier1 --check delta_bounds 2>&1 | tee delta.log && \
      grep -q "DELTA_BOUNDS: PASS" delta.log
    tolerance: |
      δD ∈ [-300, +50] ‰
      δ18O ∈ [-50, +10] ‰
    reference_values: |
      Plausibility bounds chosen as conservative envelope from published isotope-enabled GCM studies. Values outside these bounds indicate either a model bug or a regime not covered by typical published cases (in which case bounds may need revision).
    blocking: true

  - id: SC11
    phase: testing
    description: |
      Tag drift documentation. The tag tracer initialized at ratio 1.0 (q_tag = q_bulk at t=0 everywhere) is tracked over the run. Column-mean q_tag/q_bulk is expected to drift by no more than 1e-4 across the simulation, since tag fractionation in the hooks is α=1 (ratio-preserving) and surface flux is passive-copy. Larger drift indicates a code bug in the tag-handling dispatch.
    criterion_type: tolerance
    command: |
      /code/E3SM/EAMXX-wiso/components/eamxx/scripts/analyze-wiso-dpscream-tier1 --check tag_drift 2>&1 | tee tag.log && \
      grep -q "TAG_DRIFT: PASS" tag.log
    tolerance: absolute 1e-4
    reference_values: |
      Column-mean q_tag/q_bulk = 1.0 at t=0. Drift across run < 1e-4 expected because: (a) atmospheric phase-change hooks use ratio-preserving formula with α=1 for tags; (b) surface flux uses passive-copy (dq_tag = dq_bulk_flux); (c) all numerical operations should preserve this to FP order-of-operations noise.
    blocking: true

  - id: SC12
    phase: testing
    description: |
      Each case produces a per-case markdown report and the top-level summary report. Reports contain the diagnostic time series, BFB result, mass-conservation residual, δ-value ranges, tag drift, and pass/fail per success criterion.
    criterion_type: shell
    command: |
      for c in dycomsrf01 arm97 comble; do \
        test -f /data/runs/wiso_dpscream_tier1/$c/report.md || exit 1; \
      done && \
      test -f /data/runs/wiso_dpscream_tier1/SUMMARY.md
    expected_output: All four report files exist
    blocking: true

  - id: SC13
    phase: testing
    description: |
      Wall-time and memory budget — total wall time for all 9 runs (3 cases × 3 configs) fits within an 8-hour interactive session on the 16-core / 128 GB docker-local platform. Peak memory usage stays below 100 GB.
    criterion_type: tolerance
    command: |
      /code/E3SM/EAMXX-wiso/components/eamxx/scripts/analyze-wiso-dpscream-tier1 --check resources 2>&1 | tee res.log && \
      grep -q "WALL_TIME_OK: PASS" res.log && grep -q "MEMORY_OK: PASS" res.log
    tolerance: |
      wall_time < 28800 seconds (8h)
      peak_memory < 100 GB
    blocking: false

  - id: SC14
    phase: integration
    description: |
      Human review of the SUMMARY.md report and per-case reports for: (a) physical reasonableness of δ-value patterns (e.g., COMBLE should show δD decreasing with altitude through the mixed-phase region; ARM97 convective anvils should be substantially depleted; DYCOMSrf01 stratocumulus layer should show characteristic δ structure); (b) confirmation that the BFB-bulk check actually validates the bulk pathway is untouched; (c) confirmation that mass-conservation residual is not just below tolerance due to compensating errors. Determine whether observed behavior warrants any new specs (e.g., addressing systematic biases observed).
    criterion_type: human-review
    reviewer: rfiorella
    blocking: true

# Out of scope
out_of_scope:
  - BOMEX, RICO, MPACE per se — substituted with the in-tree dycomsrf01, arm97, comble testmods which exercise the same physics pathways. The published-case BOMEX/RICO/MPACE configurations would require either external scmlib scripts or new in-tree testmods, both substantial additional work.
  - Land surface isotope physics — the passive-VSMOW stopgap from spec 3a covers (1 - ocnfrac); proper land ET partitioning, soil-water δ profile, plant transpiration are future specs
  - Sea-ice surface isotope physics — same stopgap as land; proper sea-ice sublimation with sea-ice δ is future
  - CIME-coupled F-compset runs — these DP-SCREAM runs use prescribed surface conditions from the case testmod, not a coupled ocean. CIME-coupled isotope plumbing (extending eamxx_cpl_indices.F90) is a future spec.
  - Tag source-attribution semantics — tag is initialized at ratio 1.0 (bulk-equivalent) and stays there as a passive-copy check; the tag-attribution spec defines region-source-specific tag semantics
  - H217O and HTO — registry_n4.cmake uses bulk + HDO + H218O + tag. These two isotopologues are functional in the code but not exercised by this validation; a future spec may add a wider-N validation
  - Comparison against published δ-value observations — this spec checks plausibility bounds only; quantitative comparison to e.g. NEEM, GISP2, or ground-based GNIP data is a separate analysis spec
  - Performance profiling and Kokkos kernel optimization
  - Restart/checkpoint testing — runs go to completion in one batch; restart-correctness testing is a future spec
  - Long-duration runs (multi-day) — case durations match published protocol; longer integrations require evaluating drift behaviors not addressed here
  - Resolution sweep — 5x5 / 3km is the only configuration tested; sensitivity to dx or column count is future work
  - Re-running with SCREAM_WISO_EVAP_KINETIC_FORM=multiplicative (Stewart-multiplicative form rather than the default stewart_exact) — the default form is used; comparison is a future spec
  - Re-running with SCREAM_WISO_SURFACE_FLUX_MODE=none — only the ocean_craig_gordon mode is exercised; the none mode is for development/debug

# Resolved decisions
resolved_decisions:
  - decision: Three hard prerequisites — equilibrium-fractionation-hooks, kinetic-fractionation-hooks, and ocean-surface-flux-craig-gordon specs must all be implemented and merged before this spec's implementation phase begins. The driver script enforces this via symbol-presence checks at startup.
    rationale: User direction (Q0 of refinement round). Without all three, the integrated runs are meaningless — the prior atmospheric phase-change physics alone would deplete the atmosphere over time without a surface source.
    date: 2026-05-13

  - decision: In-tree DP-SCREAM testmods (dycomsrf01, arm97, comble) are used instead of the originally-named BOMEX, RICO, MPACE. New sibling testmods (dycomsrf01_wiso, arm97_wiso, comble_wiso) are created to add water-tracer namelist settings.
    rationale: User direction. Investigation showed BOMEX/RICO/MPACE are not in-tree; they're managed via the external scmlib resource. The in-tree testmods exercise the same physical regimes: marine stratocumulus + warm-phase shallow cloud (DYCOMSrf01 ≈ BOMEX), continental deep convection + precipitation (ARM97 ≈ RICO), Arctic mixed-phase + WBF (COMBLE ≈ MPACE). For validation purposes the physics coverage is equivalent.
    date: 2026-05-13

  - decision: Each case is run in three configurations — TRACE_WATER=OFF (bulk-only baseline), TRACE_WATER=ON with N=1 (bulk-only via the registry path), TRACE_WATER=ON with N=4 (bulk + HDO + H218O + tag). Total 9 runs.
    rationale: User direction (Q3 of original round = N=4 + Q6 = strict BFB-bulk). The three configurations together validate: (i) tracer infrastructure doesn't perturb bulk physics in OFF→ON+N=1; (ii) adding isotope tracers doesn't perturb bulk physics in N=1→N=4. Both are SC8.
    date: 2026-05-13

  - decision: Standard published case durations — dycomsrf01 ~6h, arm97 ~24h, comble ~12h. Run-duration override via WISO_DPSCREAM_DURATION_OVERRIDE env var (in seconds) is provided for future flexibility but defaults to published protocol.
    rationale: User direction. Published durations are chosen to capture the established quasi-steady-state regimes for each case. Standard durations also fit within the compute budget without truncation.
    date: 2026-05-13

  - decision: Build pipeline is CIME create_test (option (i) of Q1 in the original round; all three cases use CIME, since this is the natural path for the in-tree dpxx testmods).
    rationale: User direction (Q1). The in-tree dpxx testmods are wired into CIME testing; standalone DP-SCREAM is not separately wired in this repo. Sticking with CIME for all three cases avoids a build-pipeline split.
    date: 2026-05-13

  - decision: Initial conditions for isotope tracers — q_iso(t=0) = R_VSMOW(iso) × q_bulk(t=0) at every grid cell. Tag initialized at q_tag = q_bulk (ratio 1.0).
    rationale: Simplest physically-reasonable initial condition. Spec 3a's ocean Craig-Gordon surface flux + atmospheric phase-change hooks then drive the atmosphere toward a quasi-steady-state composition over the first few hours.
    date: 2026-05-13

  - decision: Output frequency — hourly netCDF per tracer (qv, qc, qi at all levels; column-integrated values; surface flux per tracer; precipitation flux per tracer).
    rationale: Hourly is standard for DP-SCREAM analysis and balances diagnostic resolution against storage. For dycomsrf01 (6h) → 6 output frames; for arm97 (24h) → 24; for comble (12h) → 12. Output-frequency override via WISO_DPSCREAM_OUTPUT_FREQ_HOURS env var.
    date: 2026-05-13

  - decision: Strict BFB-bulk requirement (SC8) — bulk water (CMP=0) is bit-for-bit identical across all three configurations for each case. Any deviation is a code bug, not FP noise.
    rationale: User direction (Q6 of original round). The hooks redistribute tendencies among tracers without touching the bulk; correct implementation preserves bit-identical bulk-water behavior.
    date: 2026-05-13

  - decision: Mass conservation tolerance — column-integrated isotope mass change rate matches surface flux − precipitation flux to relative tolerance 1e-10 (SC9). This is the only conservation check in this spec.
    rationale: User direction (Q7 of original round). δ-value plausibility (SC10) and tag drift (SC11) are separate plausibility checks, not strict conservation; documented as such in their criteria.
    date: 2026-05-13

  - decision: δ-value plausibility bounds — δD ∈ [-300, +50] ‰ and δ18O ∈ [-50, +10] ‰ in both vapor and condensate at all levels and times.
    rationale: Conservative envelope of expected values from published isotope-GCM studies (NEEM, CAM-wiso, LMDZ-iso). Tight enough to catch implementation bugs (unphysical δ values), loose enough to accommodate range of physical conditions across the three cases. May need revision if observed during runs that a particular case routinely exits these bounds for legitimate physical reasons.
    date: 2026-05-13

  - decision: Tag drift tolerance — column-mean q_tag/q_bulk drift across run < 1e-4 (SC11). Tighter than the mass-conservation tolerance because tag drift accumulates more slowly (ratio-preserving formula has near-zero contribution per timestep).
    rationale: Tags are exercised by atmospheric phase-change hooks (α=1 ratio-preserving formula from prior specs) and by surface flux (passive-copy from spec 3a). Both pathways should preserve tag ratio to FP order-of-operations noise. 1e-4 allows for accumulated FP noise over thousands of timesteps but rejects systematic bias.
    date: 2026-05-13

  - decision: Compute budget — 16 cores, 128 GB, 8 hours wall time for all 9 runs combined. Memory ceiling 100 GB (leaves headroom for analysis script + system).
    rationale: User direction. 5x5 / 3km DP-SCREAM with N=4 tracer dim is small (25 columns × ~128 levels × 4 tracers = 12800 tracer cells per column-level pair). Memory and wall time are not tight constraints; SC13 verifies but is non-blocking.
    date: 2026-05-13

  - decision: Baselines — first-run results are stored as Tier 1 baselines under /data/baselines/EAMXX-wiso/wiso_dpscream_tier1/. Future runs of this spec compare against these baselines (BFB or near-BFB depending on configuration). Regeneration via WISO_DPSCREAM_REGEN_BASELINES=1 env var.
    rationale: First Tier 1 establishes the reference; subsequent re-runs (after any future isotope-physics changes) compare to detect regressions.
    date: 2026-05-13

  - decision: Analysis script in Python (xarray + numpy) reading netCDF outputs. Separate from the runner so it can be re-invoked on existing outputs without re-running models.
    rationale: Standard scientific Python ecosystem for netCDF time-series analysis. xarray's labelled-array model fits per-tracer / per-level diagnostic computations cleanly.
    date: 2026-05-13

  - decision: Validation tier 1 (DP-SCREAM variant) — per-case-duration runs + BFB-bulk + mass conservation + δ-value plausibility + tag drift + resource bounds + human-review of summary. Codifies the project-specific Tier 1 definition from Q8 of the earlier round.
    rationale: Adapts the global E3SM Tier 1 ("5-day run + conservation") to DP-SCREAM context where case durations are case-protocol-specific rather than a uniform 5 days. Documented as a project Tier 1 variant.
    date: 2026-05-13

# Ask-before actions (project-specific additions to global policy)
ask_before:
  - Modifying any existing testmod under components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/dpxx/ (dycomsrf01, arm97, comble) — only adding *_wiso siblings is in scope
  - Modifying components/eamxx/cime_config/config_compsets.xml or compset-defining files
  - Modifying components/eamxx/CMakeLists.txt
  - Modifying components/eamxx/scripts/test-all-eamxx (this spec uses CIME create_test instead)
  - Adding fields to the field manager or IO metadata (output yaml additions are in scope; broader field-manager changes are not)
  - Modifying components/eamxx/src/control or components/eamxx/src/share — these spec runs use the integrated implementation from prior specs; the runner does not modify physics code
  - Modifying any source under components/eamxx/src/physics/water_tracers/ — physics is locked from prior specs at this point
  - Re-tuning the boce, rnat, dkfac, fkhum, or recrit constants
  - Storing run outputs or baselines outside /data/runs/wiso_dpscream_tier1/ or /data/baselines/EAMXX-wiso/wiso_dpscream_tier1/
  - Running with non-default SCREAM_WISO_EVAP_KINETIC_FORM or SCREAM_WISO_SURFACE_FLUX_MODE values (the spec uses both defaults)

# Parallelization
allow_parallelization: false

# Post-completion review
request_performance_review: false
request_code_review: false
---

# Tier 1 DP-SCREAM Validation Runs (DYCOMSrf01, ARM97, COMBLE) with the Integrated Water-Isotope Physics

## Context

This is the first end-to-end Tier 1 validation of the EAMXX-wiso water-isotope physics. The three prior specs landed on the same day (2026-05-13) cover, respectively:

1. **Equilibrium fractionation hooks** (`2026-05-13-eamxx-equilibrium-fractionation-hooks`) — `eq_condensation`, `eq_melting`, `eq_freezing` (and originally `eq_deposition`, since superseded)
2. **Kinetic fractionation hooks** (`2026-05-13-eamxx-kinetic-fractionation-hooks`) — `kin_evaporation` (Stewart 1975 hydrometeor evap, CMake-selectable form), `kin_deposition` (Jouzel-Merlivat 1984 WBF/supersaturation; replaces `eq_deposition`), `eq_sublimation`
3. **Ocean surface flux via Craig-Gordon** (`2026-05-13-eamxx-ocean-surface-flux-craig-gordon`) — `surface_flux_hook` with `ocean_craig_gordon_surface_flux` implementation using `AlphaKMol` + `AlphaEqLiquidVapor(SST)` + R_VSMOW for the ocean fraction; passive-VSMOW stopgap for non-ocean fraction

Together these establish atmospheric phase-change fractionation (six processes) and a defensible surface-flux source. Without this Tier 1 validation, the physics is unverified end-to-end: unit tests confirm correctness in isolation, but the integrated runtime behavior — conservation across many timesteps, plausibility under realistic forcings, BFB invariance of the bulk path — has not been exercised.

### Why these three cases

The originally-proposed BOMEX, RICO, and MPACE are not in-tree DP-SCREAM testmods in this repo; they are managed via the external scmlib resource (`https://github.com/E3SM-Project/scmlib/wiki/Doubly-Periodic-SCREAM-Home`). The in-tree dpxx testmods are dycomsrf01, arm97, and comble — which substitute on physics grounds:

- **dycomsrf01** (marine stratocumulus, ~6h) replaces BOMEX → exercises eq_condensation in shallow warm cloud, kin_evaporation in drizzle drift
- **arm97** (continental deep convection, ~24h) replaces RICO → exercises strongly kin_evaporation in re-evaporating precipitation, kin_deposition in deep-convective anvil ice
- **comble** (Arctic cold-air-outbreak mixed-phase, ~12h) replaces MPACE → the closest analog, exercising kin_deposition under S_i > 1 (WBF), eq_sublimation, and eq_freezing/eq_melting

The case names changed; the physics coverage did not.

### Compute fit

5x5 columns at 3km horizontal, full vertical levels, N=4 tracer dimension → 25 columns × ~128 levels × 4 tracers = a small problem by EAMxx standards. On 16 cores / 128 GB the limiting factor is wall time of the 24h ARM97 case; with 3 configurations per case and 3 cases, total wall time of ~6h is expected.

## Approach

### Run matrix

| Case | Configuration | TRACE_WATER | N tracers | Purpose |
|---|---|---|---|---|
| dycomsrf01 | bulk baseline | OFF | n/a | Reference for BFB check |
| dycomsrf01 | N=1 path | ON | 1 | Validates registry/hook plumbing doesn't perturb bulk |
| dycomsrf01 | N=4 isotopes | ON | 4 | Real isotope physics; primary validation |
| arm97 | (same three) | — | — | — |
| comble | (same three) | — | — | — |

Total: 9 runs. The three configurations per case are necessary to enforce SC8 (BFB-bulk across all three).

### What we expect to see (physically)

**dycomsrf01 (marine stratocumulus, 6h):**
- Surface flux from Craig-Gordon over ocean keeps atmospheric δD near a quasi-steady value (-70 to -110 ‰ in BL vapor depending on h and ustar)
- δD in stratocumulus condensate is heavier than vapor by α_eq − 1 ≈ +80‰
- Drizzle re-evaporation through kin_evaporation depletes the residual liquid; vapor below cloud reflects this
- δ18O follows δD with d-excess ≈ 10‰ in BL vapor

**arm97 (continental deep convection, 24h):**
- Surface flux limited to ocnfrac fraction of grid (ARM SGP is continental → ocnfrac ≈ 0; mostly passive-VSMOW from non-ocean stopgap)
- Deep convective updraft produces strong vertical δD gradient (BL ~-100‰, upper troposphere ~-300‰ from Rayleigh distillation)
- Rain re-evaporation (kin_evaporation) modifies near-surface vapor toward heavier values during precipitation
- Anvil ice condensate δD even more depleted than upper-troposphere vapor

**comble (Arctic mixed-phase, 12h):**
- Surface flux from Craig-Gordon over open ocean; very cold SST → α_eq large for HDO and H218O
- Mixed-phase clouds with S_i > 1 → kin_deposition reduces effective fractionation vs equilibrium
- δD in vapor depleted strongly due to cold-region effect (-200 to -300‰)
- Ice condensate at cloud top heavily depleted

### Diagnostics computed by the analysis script

For each case and configuration:

1. **BFB-bulk** — compare bulk q_vap, q_cloud, q_ice, q_rain, q_snow time series across configurations. Strict bit-equality required.
2. **Per-isotopologue conservation** — compute column-integrated mass M_iso(t) = ∫ ρ × q_iso × dz over each column at each output time. Verify d(M_iso)/dt ≈ F_surf_iso − F_precip_iso to relative tolerance 1e-10.
3. **δ-value time series** — for each output time and each level, compute δD = (R/R_VSMOW − 1) × 1000 for vapor and condensate. Verify within plausibility bounds.
4. **Tag drift** — column-mean q_tag/q_bulk vs time. Verify drift < 1e-4 across run.
5. **Resource usage** — wall time per run, peak memory.

Outputs per case under `/data/runs/wiso_dpscream_tier1/<case>/` include:
- Raw netCDF from each of the 3 configurations
- `report.md` with all diagnostics
- A top-level `SUMMARY.md` aggregates pass/fail across cases.

### Driver and analysis scripts

`run-wiso-dpscream-tier1` (bash):
- Verifies prerequisite symbols are present in source
- For each case in {dycomsrf01, arm97, comble}:
  - For each config in {OFF, ON_N1, ON_N4}:
    - `create_test SMS_Ld<dur>.<res>.F2010-SCREAMv1-DP.docker_gnu --testmod-dirs <case>_wiso --user-mods-dirs <wiso config>`
    - Runs the case, stores output under standard CIME location
    - Copies netCDF to `/data/runs/wiso_dpscream_tier1/<case>/<config>/`
- After all 9 runs, invokes `analyze-wiso-dpscream-tier1`

`analyze-wiso-dpscream-tier1` (Python):
- Loads outputs via xarray
- Runs BFB check, mass conservation, δ bounds, tag drift, resource diagnostics
- Writes per-case and summary reports

### Decomposition

1. Three new testmod directories (`dycomsrf01_wiso`, `arm97_wiso`, `comble_wiso`) under `cime_config/testdefs/testmods_dirs/eamxx/dpxx/` — each adds water-tracer namelist settings on top of the existing case namelist
2. `run-wiso-dpscream-tier1` bash driver
3. `analyze-wiso-dpscream-tier1` Python analysis
4. Output yaml additions per testmod (per-tracer field output at hourly cadence)
5. `REGISTRY_README.md` Tier-1-validation section

### Risks

- **Per-isotopologue output field registration.** The N=4 configuration must produce netCDF output that distinguishes the per-tracer dimension (CMP=0..3). Verify the field manager and IO yaml support this. If only bulk fields are output by default, the analysis script will be unable to compute δ-values. Investigate during implementation.
- **CIME testmod compatibility with SCREAM_TRACE_WATER.** The dpxx testmods set EAMxx namelist defaults; SCREAM_TRACE_WATER and SCREAM_WATER_TRACERS_FILE are CMake-configure options, not runtime namelist. The wiso siblings need to set them via shell_commands at case-build time before `case.build`.
- **Bit-for-bit bulk water under FP-ordering changes.** If the hook loop introduces FP-ordering changes (e.g., new sum-reduction order in tendency aggregation) that don't touch the bulk semantics, the bulk could still drift at 1e-15. Strict BFB requires that no such reordering occurs. The hooks were designed to leave bulk untouched, but verify during the first run; if non-BFB observed, investigate before relaxing the criterion.
- **Tag drift at FP noise.** 1e-4 is a generous tolerance for tag drift, but if it's systematically violated, the tag-handling dispatch in the hooks (passive copy at α=1) likely has a subtle bug. The unit tests from prior specs should have caught most issues; this is the first end-to-end verification.
- **arm97 wall-time exceeding the 8-hour budget.** 24h sim time on 25 columns is small, but CIME overhead, initial IC processing, and three configurations multiplied could compound. Monitor during first run.
- **δ-value plausibility bounds revision.** The proposed bounds [-300, +50] for δD and [-50, +10] for δ18O may be too tight for COMBLE upper-cloud levels or arm97 anvils where Rayleigh distillation produces very depleted values. If observed, the bounds need adjustment with documented justification — not the implementation.

Relevant skills to load during implementation: `e3sm-build-and-test`, `e3sm-platforms`, `python-analysis-conventions`.

## References

- **Prior specs (this project):**
  - `2026-05-13-eamxx-equilibrium-fractionation-hooks.md`
  - `2026-05-13-eamxx-kinetic-fractionation-hooks.md`
  - `2026-05-13-eamxx-ocean-surface-flux-craig-gordon.md`
- **EAMxx DP documentation:**
  - `components/eamxx/docs/user/dp_eamxx.md` — pointer to scmlib resource
  - External: https://github.com/E3SM-Project/scmlib/wiki/Doubly-Periodic-SCREAM-Home
- **In-tree DP-SCREAM testmods (substituted for BOMEX/RICO/MPACE):**
  - `components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/dpxx/dycomsrf01/`
  - `components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/dpxx/arm97/`
  - `components/eamxx/cime_config/testdefs/testmods_dirs/eamxx/dpxx/comble/`
- **Case-specific scientific references:**
  - **DYCOMSrf01:** Stevens, B. et al. (2003). On entrainment rates in nocturnal marine stratocumulus. *Q.J.R. Meteorol. Soc.* 129, 3469–3493. Plus Ackerman, A. S. et al. (2009). Large-eddy simulations of a drizzling, stratocumulus-topped marine boundary layer. *Mon. Weather Rev.* 137, 1083–1110.
  - **ARM97:** Xie, S. et al. (2005). Simulations of midlatitude frontal clouds by single-column and cloud-resolving models during the ARM March 2000 cloud IOP. *J. Geophys. Res.* 110, D15S03. (ARM97 IOP context.)
  - **COMBLE:** Geerts, B. et al. (2022). The COMBLE campaign. *Bull. Amer. Meteorol. Soc.* 103, E1371–E1396.
- **Isotope-validation literature for comparison context:**
  - Risi, C. et al. (2010). LMDZ-iso evaluation against GNIP. *J. Geophys. Res.* 115, D24123.
  - Werner, M. et al. (2011). ECHAM5-wiso evaluation. *J. Geophys. Res.* 116, D15109.

## Notes

- This spec closes the loop on a sequence of four specs landed today (the three prior implementation specs plus this validation spec). After this lands, the project has its first quantitatively-validated water-isotope-enabled EAMxx and is positioned for follow-on work: land/sea-ice surface flux, full CIME-coupled F-compset runs, tag-source attribution, restart correctness, and comparison against observational δ-value datasets.
- The case-name substitution (DYCOMSrf01/ARM97/COMBLE for BOMEX/RICO/MPACE) is a pragmatic adaptation to the in-tree testmod inventory. If a future spec wants to add the originally-named cases either via scmlib integration or new in-tree testmods, the driver and analysis scripts produced here generalize: they take a `--case <name>` argument and assume only that the case has a sibling `<name>_wiso` testmod and produces netCDF output with the standard tracer dimensions.
- The "Tier 1 (DP-SCREAM variant)" designation is documented in the resolved_decisions to distinguish from the global E3SM Tier 1 (5-day ne4 + conservation). The criteria are equivalent in spirit (run completes + conservation + plausibility) but adapted to DP-SCREAM context.
- Output stored under `/data/runs/wiso_dpscream_tier1/` and baselines under `/data/baselines/EAMXX-wiso/wiso_dpscream_tier1/`. These directories are inside the docker container and persist via the standard docker-local volume mounts.
