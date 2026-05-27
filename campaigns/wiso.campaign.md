---
campaign_id: 2026-05-26-wiso
title: "EAMxx water tracers and isotope ratios"
created: 2026-05-26
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
  size_soft_limit:
    files: 40
    lines: 1500                              # raised: array-dim PRs in group 1 will exceed default 600

specs:
  # Group 1: structural array extension (PRs 1-4)
  - path: specs/2026-05-26-water-species-concept.md
  - path: specs/2026-05-26-extend-qv-multi-tracer.md
  - path: specs/2026-05-26-extend-cloud-multi-tracer.md
  - path: specs/2026-05-26-extend-precip-multi-tracer.md

  # Group 2: fractionation primitives (PRs 6-9). PR 5 lives off-chain.
  - path: specs/2026-05-26-equilibrium-fractionation.md
  - path: specs/2026-05-26-alpha-diff-functions.md
  - path: specs/2026-05-26-net-fractionation.md
  - path: specs/2026-05-26-wtrc-ratio-utility.md

  # DP harness brought forward so group 3 specs can use it.
  - path: specs/2026-05-26-dp-harness.md

  # Group 3: parameterization hooks (PRs 10-14)
  - path: specs/2026-05-26-ocean-evap-hook.md
  - path: specs/2026-05-26-shoc-hook.md
  - path: specs/2026-05-26-p3-hook.md
  - path: specs/2026-05-26-zm-hook.md
  - path: specs/2026-05-26-rrtmgp-radiation-audit.md

  # Group 4: auxiliary + tagged tracers (PRs 15-20)
  - path: specs/2026-05-26-ch4-oxidation-hdo.md
  - path: specs/2026-05-26-region-tagged-evap.md
  - path: specs/2026-05-26-sh-decomp-prototype.md
  - path: specs/2026-05-26-sh-decomp-impl.md
  - path: specs/2026-05-26-evap-tracer-multiplicative.md
  - path: specs/2026-05-26-condensation-tagged.md
  - path: specs/2026-05-26-parcel-integrated.md

  # Group 5: comprehensive guide
  - path: specs/2026-05-26-wiso-user-guide.md

execution:
  mode: stacked
  on_spec_halt: stop-chain
  pause_between_specs: true                  # true: explicit checkpoint at group boundaries; downgrade per-spec in resolved_decisions if not desired
---

# EAMxx water tracers and isotope ratios

## Context

This campaign implements water tracers and isotope ratios in EAMxx,
adapting the architecture used in CAM (see `references/` for the CAM
implementation) to EAMxx's C++ / Kokkos / EKAT idioms. The work is
delivered as a stacked chain of pull requests because:

- The array-dimension extension PRs (group 1) touch the entire
  microphysics + transport surface and are large enough that a single
  monolithic PR would be unreviewable.
- The fractionation primitives (group 2) are independently testable and
  benefit from review separation from the array extension work.
- The parameterization hook PRs (group 3) are per-component, naturally
  reviewable one parameterization at a time.
- The tagged-tracer extensions (group 4) are research-grade and may
  evolve based on experience with the earlier groups.

The campaign plan, conservation invariants, group boundaries, and
baseline policy live in `wiso_campaign_plan.md` in the repo root.

## Ordering rationale

- **Group 1 first (PRs 1-4)** — establishes the data structures. Every
  later PR depends on the species-dim array layout. Each PR in this
  group is a pure refactor that must be slice-0 BFB against the
  pre-campaign baseline. Order: enum + add_tracer → qv → cloud → precip,
  because the cloud/precip array shapes follow the qv pattern.

- **Group 2 (PRs 6-9)** — pure functions; no group-1 dependency on
  fractionation arithmetic. Could run in parallel with group 1, but
  stacked-PR mode runs them after. PR 9 (`wtrc_ratio` utility) is a
  blocking dependency for group 3.

- **DP harness PR before group 3** — group-3 hooks validate against the
  DP harness; harness must merge first.

- **Group 3 (PRs 10-14)** — depends on group 1 (array shapes) and PR 9
  (ratio utility). Order matches frequency of phase-change calls:
  ocean evap → SHOC → P3 → ZM → RRTMGP audit.

- **Group 4 (PRs 15-20)** — depends on the full stack. PR 17 split into
  17a (analysis) and 17b (implementation) because spherical-harmonic
  decomposition is research-grade and benefits from a prototype phase.

- **Group 5 (PR 21b)** — documentation pass at the end.

## References

- `wiso_campaign_plan.md` — the campaign source of truth
- `wiso_implementation_order.md` — original prose roadmap (retained for
  history; do not edit)
- `references/` — CAM wiso implementation, prior specs, publications
- `claude-workflows/skills/spec-schema/SKILL.md` — spec contract
- `claude-workflows/skills/tracer-conservation/SKILL.md` — invariants
- `claude-workflows/skills/regression-baseline/SKILL.md` — BFB proof
