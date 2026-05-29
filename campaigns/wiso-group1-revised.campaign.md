---
campaign_id: 2026-05-29-wiso-group1-revised
title: "EAMxx Water-Tracer Campaign: Group 1 — Structural Array Extension (Revised)"
created: 2026-05-29
author: rfiorella
project: EAMxx-wiso

repo_path: /vast/home/rfiorella/E3SM
base_branch: eamxx-wiso-dev
remote: origin

branch_prefix: wiso

pr:
  draft: true
  labels:
    - wiso
    - water-tracers
    - group1
  size_soft_limit:
    files: 40
    lines: 1500

specs:
  # PR 1: Metadata, types, and BFB feasibility gate (Tier 0)
  - path: specs/2026-05-28-water-tracer-metadata-and-gate.md
  
  # PR 2a: Field infrastructure (FieldTag, grid layouts, accessor patterns)
  - path: specs/2026-05-28-extend-qv-tracer-2a-infrastructure.md
  
  # PR 2b: SHOC process (PBL turbulence)
  - path: specs/2026-05-28-extend-qv-tracer-2b-shoc.md
  
  # PR 2c: P3 process (microphysics)
  - path: specs/2026-05-28-extend-qv-tracer-2c-p3.md
  
  # PR 2d: Remaining processes (RRTMGP, ZM, HOMME, surface coupling)
  - path: specs/2026-05-28-extend-qv-tracer-2d-remaining.md
  
  # PR 3: Extend qc, qi, qm to tracer dimension (Tier 2)
  - path: specs/2026-05-28-extend-cloud-tracer.md
  
  # PR 4: Extend qr, qs, rain, snow to tracer dimension (Tier 2)
  - path: specs/2026-05-28-extend-precip-tracer.md
  
  # PR 5: Tracer ratio utility and validation test (Tier 1)
  - path: specs/2026-05-28-tracer-validation.md

execution:
  mode: stacked
  on_spec_halt: stop-chain
  pause_between_specs: true
---

# EAMxx Water-Tracer Campaign: Group 1 — Structural Array Extension (Revised)

## Context

Group 1 extends all EAMxx water reservoir fields to support multiple tracer
constituents along a new leading dimension. This is the foundational campaign
for water tracers and isotope ratios in EAMxx.

At completion:
- All water mass fields (`qv`, `qc`, `qi`, `qr`, `qs`, `qm`, `rain`, `snow`) 
  have shape `(tracer, col, lev)` instead of `(col, lev)`
- Users register water tracers via `add_water_tracer(NAME <name> LONGNAME <desc>)` in CMake
- Slot 0 is canonical bulk water and remains bit-for-bit identical to pre-campaign baseline
- Tracer ratio validation test proves passive transport correctness

**Critical success criterion:** Existing EAMxx regression tests pass BFB when 
`SCREAM_NUM_TRACERS=1` (bulk water only in slot 0).

## Revision Notes (2026-05-29)

**Original campaign** had 5 specs (PRs 1-5). During execution, spec 2 ("Extend qv to tracer dimension") was determined to be too large in scope and was **split into 4 sub-specs (2a-2d)**:

- **Spec 2a**: Field infrastructure only (FieldTag, grid layouts, accessor patterns)
- **Spec 2b**: SHOC process (first physics process, establishes pattern)
- **Spec 2c**: P3 microphysics (validates pattern for complex multi-kernel process)
- **Spec 2d**: Remaining processes (RRTMGP, ZM, HOMME, surface coupling)

**Rationale for split:**
1. Original spec 2 scope was ~40+ files vs 14 listed deliverables
2. Architectural decisions (accessor pattern selection) required validation in infrastructure before applying to all processes
3. Incremental approach enables faster BFB debugging and easier code review
4. Each sub-spec is now ~10-15 files and <1500 lines (within PR size guidelines)

**Campaign now has 8 PRs** instead of 5. Total scope unchanged, just better structured.

## Ordering Rationale

- **PR 1 (Tier 0)** — Establishes types, metadata, CMake function, and runs 
  BFB feasibility gate. Must complete before PRs 2a-5. Gate determines whether 
  exact BFB is achievable (result: YES, passed).

- **PR 2a (Tier 1)** — Field infrastructure: `TRACER` FieldTag, `get_3d_tracer_layout()`, 
  and both accessor patterns (explicit indexing and Kokkos subview) with BFB validation 
  for each. Critical gate determines which accessor pattern to use in PRs 2b-2d.

- **PR 2b (Tier 2)** — Extends SHOC (PBL turbulence) to tracer-aware qv. First 
  physics process modification. Establishes pattern for PRs 2c-2d. SHOC chosen 
  because it has standalone test suite for isolated BFB validation.

- **PR 2c (Tier 2)** — Extends P3 (microphysics) to tracer-aware qv. Most complex 
  qv consumer (all phase transitions). Validates pattern for complex multi-kernel 
  processes. Critical BFB test.

- **PR 2d (Tier 2)** — Extends remaining processes (RRTMGP, ZM, HOMME, surface 
  coupling) to tracer-aware qv. Each is smaller scope than SHOC/P3. Full-atmosphere 
  integration test validates all processes work together. After this PR, qv tracer 
  extension is complete.

- **PR 3 (Tier 2)** — Extends `qc`, `qi`, `qm` to tracer dimension. Follows 
  pattern from PRs 2b-2d. More processes touch cloud water than vapor (radiation, 
  aerosols, cloud fraction).

- **PR 4 (Tier 2)** — Extends `qr`, `qs`, `rain`, `snow` to tracer dimension. 
  Follows pattern from PRs 2b-2d. Primarily P3 microphysics and sedimentation.

- **PR 5 (Tier 1)** — Adds tracer ratio utility and validation test. Test 
  scales surface evaporation flux by 0.5 for tracer 1, verifies precipitation 
  ratios converge to 0.5 within numerical precision (rtol=1e-12). This proves 
  passive tracer transport is correct.

## Test Strategy

- **Tier 0 (PR 1)**: Header/types compilation, CMake config, prototype 
  performance and BFB tests
- **Tier 1 (PRs 2a, 5)**: Unit tests and doubly-periodic tests with tracer validation
- **Tier 2 (PRs 2b-2d, 3-4)**: Full EAMxx test suite with baseline comparison, 
  BFB requirement enforced

## Architectural Decisions (Documented in PR 2a)

These decisions apply to PRs 2a-5:

**Decision 1: Compile-time tracer count**  
Use `SCREAM_NUM_TRACERS` as compile-time constant. Rationale: PR 1 BFB gate passed with this approach. Simpler implementation and easier BFB maintenance.

**Decision 2: New grid layout method**  
Add `get_3d_tracer_layout(ntracers)` instead of repurposing vector layout. Rationale: Clearer semantics, better documentation, minimal cost.

**Decision 3: Hybrid accessor pattern validation**  
Implement both explicit indexing (`qv(0, icol, ilev)`) and Kokkos subview (`subview(qv, 0, ALL, ALL)`) in PR 2a. Test both for BFB. Use the pattern that passes for PRs 2b-2d. Rationale: Subview offers better performance/readability but BFB is non-negotiable. Testing both minimizes risk.

**Accessor pattern result (from PR 2a)**: TBD - will be documented in PR 2a's resolved_decisions and carried forward to PRs 2b-2d.

## Group 1 Boundary

After PR 5 merges:
1. Regenerate baseline with `SCREAM_NUM_TRACERS=1`
2. Compare against pre-campaign baseline (must be BFB)
3. Run tracer scaling validation test
4. Verify performance overhead < 2% for `SCREAM_NUM_TRACERS=1`
5. Verify performance overhead < 5% for `SCREAM_NUM_TRACERS=2`

## References

- `wiso_group1_campaign_revision.md` — detailed campaign plan (this repo root)
- `.claude/project-config.yml` — project configuration
- `claude-workflows/skills/spec-schema/SKILL.md` — spec contract
- `claude-workflows/skills/tracer-conservation/SKILL.md` — invariants
- `claude-workflows/skills/regression-baseline/SKILL.md` — BFB proof strategy
- `claude-workflows/skills/eamxx-cpp-conventions/SKILL.md` — C++/Kokkos patterns
- Original campaign: `campaigns/wiso-group1.campaign.md` (superseded by this revision)
