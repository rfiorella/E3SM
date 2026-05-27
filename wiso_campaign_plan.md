# EAMxx water-tracer / isotope campaign plan

This document is the source of truth for the EAMxx water-tracer + isotope
campaign. Each section maps onto a concrete artifact in this repo: the
project config, the campaign manifest, the per-step specs, and the
baselines.

This file replaces the prose-only `wiso_implementation_order.md`. The old
file is retained for history; do not edit it.

## Pre-campaign bootstrap

Before any spec runs, the following must exist in this repo:

1. **`.claude/project-config.yml`** — captures defaults inherited by
   every spec. See the template at the bottom of this doc and the
   committed copy in this repo. Required fields:
   - `host_source_path: /code/E3SM/EAMXX-wiso`
   - `container_image: rfiorella/model-containers:e3sm-openmpi-dev-latest`
   - `target_compset` and `target_resolution` (full-grid Tier-2 baseline target)
   - `dp_compset` and `dp_resolution` (DP test harness — see
     `claude-workflows/skills/e3sm-platforms/references/eamxx-doubly-periodic.md`)
   - `baseline_root: /data/baselines`
   - `spec_mode: mutable`
   - `branch_prefix: wiso`
   - `base_branch: wiso-dev`

2. **`references/` folder** — populate with:
   - CAM water-tracer / wiso implementation extracted from a tagged CAM SHA
     (note the SHA in `references/CAM_SHA.txt`)
   - Prior wiso draft specs (any earlier attempts in this repo)
   - PDFs of the publications cited in steps 6-8 and 15-20
   - A short `references/README.md` indexing the above

3. **`wiso-dev` long-lived branch** — created off `master` at the campaign
   start. All PRs land on `wiso-dev`. Periodically rebased onto `master`
   via `/restack-campaign`.

4. **Pre-campaign Tier-2 baseline** — generate at the SHA where `wiso-dev`
   branches from `master`. This is the regression target for every PR that
   extends an array dimension. Procedure: see
   `claude-workflows/skills/regression-baseline/SKILL.md`.

## Design principles (binding)

- Water-tracer / isotope-ratio code lives in `scream::water_tracers`
  namespace. New files under `components/eamxx/src/physics/water_tracers/`.
- Multi-tracer arrays use a leading tracer-species dimension. Species
  index 0 = bulk H2O. **Slice-0 values must remain bit-for-bit identical
  to the pre-extension scalar values for any single-tracer build.** This
  is the `slice-1-unchanged` invariant defined in
  `claude-workflows/skills/tracer-conservation/SKILL.md`.
- Existing tests are never modified. Tests for new behavior are added,
  not edited into existing ones.
- CMake-selectable schemes are declared once and documented in each spec's
  `resolved_decisions` block.

## Skills referenced by every spec in this campaign

- `e3sm-build-and-test` — validation tiers, lint, cprnc, conservation
- `e3sm-platforms` (with `references/docker-local.md` and
  `references/eamxx-doubly-periodic.md`)
- `eamxx-cpp-conventions` — C++17, Kokkos, EKAT, field-manager tracers
- `scientific-modeling-conventions` — provenance, units, comments
- `tracer-conservation` — invariants and tolerance defaults
- `regression-baseline` — proving array-dim extensions are non-functional
- `spec-schema` — spec format the validator enforces

Each spec body's "Approach" section names the subset it relies on.

## Conservation invariants (canonical)

These are referenced by spec success criteria via the `tracer-conservation`
skill. Repeated here for visibility:

| Invariant | rtol | atol | Verified by |
|---|---|---|---|
| slice-0 BFB vs pre-extension baseline | 0 | 0 | `cprnc` |
| water-mass closure (column, timestep) | 1e-12 | 1e-15 | `scripts/check_water_balance.py` |
| isotope-mass closure (per species) | 1e-10 | 1e-15 | `scripts/check_isotope_balance.py` (deliverable of step 9) |
| sum-of-species ≡ bulk H2O | 1e-8 | 1e-15 | same script |
| alpha(H2O) = alpha(H216O) = 1 exactly | 0 | 0 | unit tests |
| fractionation vs published values | 1e-6 | 1e-9 | unit tests over [233 K, 313 K] |

## Branch and baseline model

- **Stacked chain.** Only PR 1 bases on `wiso-dev`. Every subsequent PR
  bases on the previous PR's branch, because each spec depends on the
  state produced by the spec before it (the array-dim machinery added in
  PR 2 is required to build PR 3; the ratio utility from PR 9 is required
  by every group-3 hook; etc.). Branch naming `wiso-NN-<spec-slug>` with
  NN = 1-indexed position. The campaign-orchestrator creates each branch
  off the preceding one automatically when `execution.mode: stacked`.
- Merge order matches chain order. When PR N merges into `wiso-dev`, PR
  N+1's base is rewritten from `wiso-NN-<spec>` to `wiso-dev` and
  rebased; downstream PRs cascade the same way. The campaign-orchestrator's
  restack mode handles this.
- `wiso-dev` itself rebased onto `master` no more than once per month,
  immediately after a clean group boundary (see groups below). Restack
  via `/restack-campaign campaigns/wiso.campaign.md`.
- Baselines named `<branch>-<compset>-<resolution>` under `baseline_root`.
  Each baseline directory carries a `BASELINE.txt` recording parent SHA,
  date, compset, resolution, stop options.
- Baseline refresh decision tree lives in `regression-baseline`. Bump
  with `-rN` suffix; never overwrite.
- PR size policy: default soft limit is 25 files / 600 lines. Steps 2-4
  (array-dim extensions) will exceed. Manifest sets
  `pr.size_soft_limit.lines: 1500`. Specs that exceed even that split into
  `Na`/`Nb` sub-specs.

## Test layering

| Layer | Run cost | What it proves |
|---|---|---|
| Unit (Catch2) | seconds | Pure-function fractionation factors match published values |
| Component | minutes | A single parameterization's tracer hook is internally consistent |
| DP integration | ~10 min | End-to-end run with hooks active, conservation holds |
| Full-grid Tier-2 | hours | BFB / non-BFB-with-tolerance vs baseline |

Steps 6-9 stop at Unit. Steps 9a-9b (namelist + mass-fixer) add
Component. Steps 10-14 add Component + DP. Group boundaries (after
step 4, after step 9b, after step 14b, after step 20, after step 20c)
trigger Tier-2.

## Spec groups

Specs grouped so that group boundaries are natural Tier-2 checkpoints.
`pause_between_specs: false` within a group; `true` at group boundaries.

### Group 1 — structural array extension (PRs 1-4)

Highest risk. Every PR ≥ 2 must include slice-0 BFB regression vs the
pre-campaign baseline.

| PR | Spec id | Title | Tier | Notes |
|---|---|---|---|---|
| 1 | `2026-MM-DD-water-species-concept` | Water-species enum + `add_tracer` in `scream::water_tracers` | 0 | Header-only + CMake. No run. |
| 2 | `2026-MM-DD-extend-qv-multi-tracer` | Extend qv to species dim | 2 | Slice-0 BFB required. |
| 3 | `2026-MM-DD-extend-cloud-multi-tracer` | Extend qc, qi to species dim | 2 | Same pattern as PR 2. |
| 4 | `2026-MM-DD-extend-precip-multi-tracer` | Extend qr, qs to species dim | 2 | Same pattern. |

Group 1 boundary: Tier-2 full-grid baseline regenerated with the array
machinery in place but `tracer_count = 1`. This becomes the baseline for
groups 2-3.

### Group 2 — fractionation primitives + utilities (PRs 6-9b)

PR 5 (phase-change enumeration) is moved out of the chain — see
"Out-of-chain work" below.

| PR | Spec id | Title | Tier | Notes |
|---|---|---|---|---|
| 6 | `2026-MM-DD-equilibrium-fractionation` | Equilibrium fractionation functions | 0 | Unit tests vs published tabulated values. CMake-selectable schemes. **Deliverables include** `components/eamxx/tests/water_tracers/data/horita-wesolowski-1994.csv`, `data/majoube-1971.csv`, `data/merlivat-nief-1967.csv` per `tracer-conservation` skill convention. |
| 7 | `2026-MM-DD-alpha-diff-functions` | Molecular-diffusivity fractionation | 0 | Unit tests with exact expected values. **Deliverables include** `data/merlivat-1978.csv`, `data/schoenemann-2014.csv` (molecular-diffusivity coefficients per isotopologue). |
| 8 | `2026-MM-DD-net-fractionation` | Net fractionation (Brutsaert + Craig-Gordon, Stewart 1975, Ciais-Jouzel 1994) | 0 | Unit tests. CMake-selectable. **In scope: partial-equilibration kinetics for falling hydrometeors**, mirroring CAM `wtrc_equil_time` (`references/CAM/CAM/src/physics/cam/water_tracers.F90:5636`). Separate from the macrostep evaporation closure. |
| 9 | `2026-MM-DD-wtrc-ratio-utility` | `wtrc_ratio` utility + `check_isotope_balance.py` deliverable | 0 | Used by everything in group 3. Blocking. |
| 9a | `2026-MM-DD-wiso-namelist` | YAML / runtime config for the ~10 wiso flags (`trace_water`, `wisotope`, `wtrc_lh2oadj`, `wtrc_lzmlin`, `wtrc_warn_only`, `wtrc_add_cvprecip`, `wtrc_add_stprecip`, `wtrc_alpha_kinetic`, `wtrc_check_total_h2o`, `wtrc_detrain_in_macrop`, `wtrc_use_ice_supsat`) currently sitting as `inline bool` defaults in `water_tracers.hpp:20-31` | 1 | Mirrors CAM `wtrc_readnl`. Blocking for group-3 hooks that read these flags. |
| 9b | `2026-MM-DD-wtrc-mass-fixer` | Wiso-aware negative-mass / mass-conservation fixer | 1 | Mirrors CAM `wtrc_mass_fixer` (`water_tracers.F90:4233`) + `qneg_module` wiso entries. Required before production science runs — small concentrations + roundoff produce negative isotope masses that break sum-of-species closure. Hooks into EAMxx's existing `qneg`-style mechanism along the species dim. |

### Group 3 — parameterization hooks (PRs 10-14)

Each PR uses DP harness for Component validation; isotope-mass closure
required. Group boundary triggers Tier-2.

| PR | Spec id | Title | Tier |
|---|---|---|---|
| 10 | `2026-MM-DD-ocean-evap-hook` | Ocean evaporation isotope partitioning | 1 |
| 10b | `2026-MM-DD-surface-flux-inputs` | Coupler-side delivery of wiso surface fluxes (land transpiration from ELM, sea-ice sublimation, ocean from data file in F-compsets). Mirrors CAM `camsrfexch.F90` + `atm_import_export.F90` wiso paths. | 1 |
| 11 | `2026-MM-DD-shoc-hook` | SHOC phase-change hooks | 1 |
| 12 | `2026-MM-DD-p3-hook` | P3 microphysics phase-change hooks | 1 |
| 13 | `2026-MM-DD-zm-hook` | Zhang-McFarlane phase-change hooks | 1 |
| 14a | `2026-MM-DD-rrtmgp-radiation-audit` | Audit RRTMGP / aerosol activation for unexpected phase changes | 0 |
| 14b | (created if 14a finds gaps) | Hooks for parameterizations identified by 14a | 1 |

Group 3 boundary: full-grid Tier-2 run with `tracer_count = 4` (H2O,
HDO, H218O, H217O). Compare against group-2 boundary baseline (slice-0)
and verify isotope-mass closure and sum-of-species closure.

### Group 4 — auxiliary tracers and tagged tracers (PRs 14c, 15-20)

| PR | Spec id | Title | Tier |
|---|---|---|---|
| 14c | `2026-MM-DD-tritium-decay` | HTO radioactive decay using `WaterIsotopologues<Scalar>::hlhto` half-life constant. Mirrors CAM `wtrc_rad_decay` (`water_tracers.F90:6738`). Small standalone PR. | 0 |
| 15 | `2026-MM-DD-ch4-oxidation-hdo` | CH4 + OH → HDO impact, CMake flag | 1 |
| 16 | `2026-MM-DD-region-tagged-evap` | Region-tagged evaporation tracer (lat/lon box + shapefile) | 1 |
| 17a | `2026-MM-DD-sh-decomp-prototype` | Spherical-harmonic decomp prototype (spec_type: analysis) | n/a |
| 17b | `2026-MM-DD-sh-decomp-impl` | Production implementation of SH / needlet tagged tracers | 1 |
| 18 | `2026-MM-DD-evap-tracer-multiplicative` | Multiplicative-factor evap tracers (time, T, etc.) | 1 |
| 19 | `2026-MM-DD-condensation-tagged` | Condensation-tagged tracers (T, p at phase change) | 1 |
| 20 | `2026-MM-DD-parcel-integrated` | Parcel-integrated quantities (Q·dt etc.) | 1 |

Step 17 split into 17a (prototype / analysis) + 17b (implementation) to
avoid mixing research and implementation in a single PR.

### Group 4.5 — production support (PRs 20a-20c)

Required before production science runs but not blocking for the
campaign's scientific build-out. Group 4.5 runs after group 4 closes
and before the final user-guide pass. Group boundary triggers Tier-2.

| PR | Spec id | Title | Tier | Notes |
|---|---|---|---|---|
| 20a | `2026-MM-DD-wiso-initial-conditions` | IC dataset + reader for wiso fields | 1 | Mirrors CAM `const_init.F90` wiso paths. Default IC = VSMOW ratios applied to bulk; optional file input for non-trivial spin-up state. Without this, every wiso run starts from arbitrary state → wasted spin-up. |
| 20b | `2026-MM-DD-wiso-diagnostics` | δ-value derived diagnostics, R/Rvsmow conversions, column-integrated isotope quantities, history field registration | 0 | Mirrors CAM `cam_diagnostics.F90` wiso paths. Output convention: store mass mixing ratios in history; provide derived-field utilities and a Python helper for δ conversion in `python-analysis-conventions` style. |
| 20c | `2026-MM-DD-wiso-restart` | Restart support for wiso fields | 1 | Required for any multi-segment run. Group 1 specs explicitly defer restart; this PR closes that gap. |

### Group 5 — comprehensive harness (PR 21, split)

| PR | Spec id | Title | Tier |
|---|---|---|---|
| 21a | `2026-MM-DD-dp-harness` | DP test cases for the wiso campaign | 1 |
| 21b | `2026-MM-DD-wiso-user-guide` | Comprehensive wiso user guide + test recipes | 0 |

21a is moved earlier in scheduling — it is created in parallel with
group 2 so that group 3 hooks can use it. 21b is the documentation pass
at the end.

## Out-of-chain work

These items live outside the campaign chain because they don't ship code
on a PR cadence:

- **Step 5 (phase-change enumeration)** — produced as
  `docs/wiso/phase_changes.md` via a single `spec_type: analysis` spec
  before group 3. Output is a markdown table grouped by phase-change type
  and parameterization. Used as the input checklist for PRs 10-14. The
  doc must include a CAM → EAMxx mapping section since EAMxx's SHOC
  unifies what CAM splits across `clubb_intr.F90`, `macrop_driver.F90`,
  `convect_shallow.F90`, `uwshcu.F90`, and `vertical_diffusion.F90`;
  PR 11 (shoc-hook) must reference that mapping to confirm coverage.
- **References folder population** — one-time setup; not a PR.
- **Tier-2 baseline regeneration** — happens at group boundaries; the
  orchestrator runs it but it is not its own spec.

## Restack policy

Run `/restack-campaign campaigns/wiso.campaign.md` when:
- a PR in the chain merges (every downstream PR's base must move up one
  link; the orchestrator handles the cascade),
- `master` advances and contains a change that materially affects EAMxx
  physics or build system (rebase `wiso-dev` onto `master`, then cascade
  every PR in the chain onto the new `wiso-dev`),
- a group boundary has just been cleared (natural restack point — clean
  point to absorb upstream changes).

The campaign-orchestrator's restack mode handles the mechanical rebase
across the whole stack. Do not rebase individual PR branches by hand;
the chain integrity check (C10/C11) verifies each branch's base matches
the manifest order.

## Post-completion review

Specs in groups 1 and 4 enable `post_completion_review` with
`performance-specialist` (adding array dims has perf risk in Kokkos
kernels; tagged-tracer accounting can blow up memory). Specs in groups
2, 3, and 4.5 enable `code-reviewer` only.

## Validator expectations

Each spec must satisfy validator hard-fails (H1-H12 in `spec-schema`):
- ≥3 runnable success criteria
- Every shell/comparison/tolerance criterion has `cmd`
- Every tolerance criterion has rtol or atol > 0
- Non-empty `out_of_scope`
- `ask_before` includes destructive ops (CMake API changes touching
  shared field-manager headers count)
- No `TODO(spec):` left in the body

Soft-warns to watch (S8, S10, S11): every criterion's `cmd` should
reference a path that lives in `deliverables`. No lint-only specs (S10).
Tests cited under `pytest <path>` or `<path>_tests` must exist on disk
or be listed in deliverables (S11).

## Project-config template

The file at `.claude/project-config.yml` follows this shape (see the
committed copy in this repo for the actual values):

```yaml
project_id: EAMxx-wiso
host_source_path: /code/E3SM/EAMXX-wiso
container_image: rfiorella/model-containers:e3sm-openmpi-dev-latest
default_platform: docker-local
target_compset: <full-grid Tier-2 compset>
target_resolution: <full-grid Tier-2 resolution>
dp_compset: <DP compset>
dp_resolution: <DP resolution>
baseline_root: /data/baselines
spec_mode: mutable
branch_prefix: wiso
base_branch: wiso-dev
remote: origin
```

The `<...>` placeholders are filled in by `/init-project`.
