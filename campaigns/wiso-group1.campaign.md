---
campaign_id: 2026-05-28-wiso-group1
title: "EAMxx Water-Tracer Campaign: Group 1 — Structural Array Extension"
created: 2026-05-28
author: rfiorella
project: EAMxx-wiso

repo_path: /code/E3SM/EAMXX-wiso
base_branch: wiso-dev
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
  
  # PR 2: Extend qv to tracer dimension (Tier 2)
  - path: specs/2026-05-28-extend-qv-tracer.md
  
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

# EAMxx Water-Tracer Campaign: Group 1 — Structural Array Extension

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

## Ordering rationale

- **PR 1 (Tier 0)** — Establishes types, metadata, CMake function, and runs 
  BFB feasibility gate. Must complete before PRs 2-5. Gate determines whether 
  exact BFB is achievable or requires architectural fallback.

- **PR 2 (Tier 2)** — Extends `qv` to tracer dimension. Pattern established 
  here applies to all other water species. Includes field layout extension, 
  field manager integration, kernel updates across all physics that touch `qv`.

- **PR 3 (Tier 2)** — Extends `qc`, `qi`, `qm` to tracer dimension. Follows 
  PR 2 pattern. More processes touch cloud water than vapor (radiation, 
  aerosols, cloud fraction).

- **PR 4 (Tier 2)** — Extends `qr`, `qs`, `rain`, `snow` to tracer dimension. 
  Follows PR 2 pattern. Primarily P3 microphysics and sedimentation.

- **PR 5 (Tier 1)** — Adds tracer ratio utility and validation test. Test 
  scales surface evaporation flux by 0.5 for tracer 1, verifies precipitation 
  ratios converge to 0.5 within numerical precision (rtol=1e-12). This proves 
  passive tracer transport is correct.

## Test strategy

- **Tier 0 (PRs 1)**: Header/types compilation, CMake config, prototype 
  performance and BFB tests
- **Tier 2 (PRs 2-4)**: Full EAMxx test suite with baseline comparison, 
  BFB requirement enforced
- **Tier 1 (PR 5)**: Doubly-periodic test with tracer ratio validation

## Group 1 boundary

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
