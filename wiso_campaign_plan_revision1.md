# EAMxx water-tracer / isotope campaign plan

This document is the source of truth for the EAMxx water-tracer + isotope
campaign. Each section maps onto a concrete artifact in this repo: the
project config, the campaign manifest, the per-step specs, and the
baselines.

This file replaces the prose-only `wiso_implementation_order.md`. The old
file is retained for history; do not edit it.

## Revision note

This version addresses the two highest-risk issues identified in review:

1. **Tracer semantics are now explicit.** Bulk water, isotope tracers, and
   tagged tracers are distinct concepts. Bulk H2O is not reconstructed by
   summing all isotope tracer slots. The `sum-of-species == bulk H2O`
   invariant applies only to explicitly declared closed tagged-tracer
   partition groups.

2. **Bit-for-bit preservation is now a gated compatibility target.** Exact
   BFB is required for the single-tracer compatibility path only if PR 1
   proves the selected field-layout strategy can support it without
   unacceptable architectural or performance cost. If not, the campaign
   must pause and adopt the documented non-BFB fallback before PRs 2-4
   proceed.

## Pre-campaign bootstrap

Before any spec runs, the following must exist in this repo:

1. **`.claude/project-config.yml`** — captures defaults inherited by
   every spec. See the template at the bottom of this doc and the
   committed copy in this repo. Required fields:

   - `host_source_path: /code/E3SM/EAMXX-wiso`
   - `container_image: rfiorella/model-containers:e3sm-openmpi-dev-latest`
   - `target_compset` and `target_resolution` meaning the full-grid Tier-2
     baseline target
   - `dp_compset` and `dp_resolution` meaning the DP test harness; see
     `claude-workflows/skills/e3sm-platforms/references/eamxx-doubly-periodic.md`
   - `baseline_root: /data/baselines`
   - `spec_mode: mutable`
   - `branch_prefix: wiso`
   - `base_branch: wiso-dev`

2. **`references/` folder** — populate with:

   - CAM water-tracer / wiso implementation extracted from a tagged CAM SHA
     with the SHA recorded in `references/CAM_SHA.txt`
   - prior wiso draft specs, if any
   - PDFs of the publications cited in steps 6-8 and 15-20
   - a short `references/README.md` indexing the above

3. **`docs/wiso/tracer_data_model.md`** — required before PR 2 starts.
   PR 1 creates this document and the campaign treats it as binding. It
   must define:

   - the difference between bulk water, isotope tracers, and tagged tracers
   - whether each tracer slot is prognostic, diagnostic, isotope, tagged, or
     auxiliary
   - storage units for isotope tracers
   - how isotope ratios are computed
   - whether H216O is explicit or implicit
   - which tracer groups, if any, form a closed partition of bulk water
   - the BFB compatibility result from the PR 1 feasibility gate

4. **`wiso-dev` long-lived branch** — created off `master` at the campaign
   start. All PRs land on `wiso-dev`. Periodically rebased onto `master`
   via `/restack-campaign`.

5. **Pre-campaign Tier-2 baseline** — generate at the SHA where `wiso-dev`
   branches from `master`. This is the regression target for every PR that
   extends an array dimension. Procedure: see
   `claude-workflows/skills/regression-baseline/SKILL.md`.

## Design principles, binding

### Namespaces and file placement

- Water-tracer / isotope-ratio code lives in the
  `scream::water_tracers` namespace.
- New files live under:

  ```text
  components/eamxx/src/physics/water_tracers/
  ```

- CMake-selectable schemes are declared once and documented in each spec's
  `resolved_decisions` block.

### Binding tracer data model

The campaign uses a slot-based water-tracer axis, but the slots do not all
have the same physical meaning.

#### Slot 0: canonical bulk water

- Slot 0 is the canonical bulk water amount for each water reservoir:
  `qv`, `qc`, `qi`, `qr`, and `qs`.
- Slot 0 is the only slot that represents the total physical water amount
  used by legacy EAMxx physics.
- Existing EAMxx water physics must continue to read and update the slot-0
  bulk field.
- Bulk water is not computed by summing all tracer slots.

#### Isotope tracer slots

For the initial isotope configuration:

```text
tracer_count = 4
slot 0 = bulk H2O
slot 1 = HDO water-equivalent tracer mass
slot 2 = H2^18O water-equivalent tracer mass
slot 3 = H2^17O water-equivalent tracer mass
```

The notation above names isotope tracer identities. It does **not** mean
that bulk water is the sum of slots 0-3.

Rules:

- Isotope tracer slots store water-equivalent mass mixing ratios unless
  `docs/wiso/tracer_data_model.md` explicitly chooses a different unit.
- H216O is implicit by default, not a separately prognosed slot.
- Isotope ratios are derived from isotope tracer mass and bulk water using
  `scream::water_tracers::wtrc_ratio`.
- Isotope tracers are independently transported and fractionated but do not
  form a complete partition of bulk water.
- Therefore, this invariant is invalid for isotope mode:

  ```text
  sum(all slots) == bulk H2O
  ```

#### Tagged tracer slots

Tagged tracers are different from isotope tracers.

- A tagged-tracer group may be declared as a closed partition of bulk water
  only when its metadata says so.
- For a closed tagged partition, the valid invariant is:

  ```text
  sum(slots in that tagged group) == slot 0 bulk H2O
  ```

- For non-partition tags, such as multiplicative diagnostic tags or
  parcel-integrated quantities, no sum-to-bulk invariant is implied.

#### Required tracer metadata

Every nonzero tracer slot must carry metadata sufficient to determine which
invariants apply:

```text
name
kind: isotope | tagged_partition | tagged_diagnostic | auxiliary
reservoirs: qv | qc | qi | qr | qs | subset
units
ratio_standard, if isotope
partition_group_id, if tagged_partition
conserved_independently: true | false
```

PR 1 must introduce the enum/types necessary to represent this metadata.

### BFB compatibility policy

The campaign distinguishes three cases.

#### 1. Single-tracer compatibility mode

```text
tracer_count = 1
slot 0 only
```

In this mode, the campaign target is exact bit-for-bit identity against the
pre-extension scalar baseline.

#### 2. Multi-tracer isotope mode

```text
tracer_count > 1
isotope slots active
```

In this mode, exact BFB of slot 0 is desirable but not required. Slot-0 bulk
water must instead satisfy the non-BFB scientific regression tolerances
defined in the relevant spec and must pass water-mass closure.

#### 3. Multi-tracer tagged mode

```text
tracer_count > 1
tagged tracers active
```

BFB is not expected. The relevant tagged-partition or tagged-diagnostic
invariants apply.

### PR 1 BFB feasibility gate

PR 1 is a blocking design/architecture gate.

Before PRs 2-4 extend `qv`, `qc`, `qi`, `qr`, or `qs`, PR 1 must prototype
the field-layout strategy enough to answer:

- Can `tracer_count = 1` preserve exact BFB against the scalar baseline?
- Does this require retaining a scalar compatibility path?
- Does a leading tracer dimension change Kokkos pack/vectorization order?
- Does field registration or restart/history metadata change output order?
- What is the expected performance and memory impact?

If PR 1 proves that exact BFB for `tracer_count = 1` is feasible, then PRs
2-4 must enforce it as a hard gate.

If PR 1 proves that exact BFB is not feasible without unacceptable
architectural or performance cost, then PR 1 must update these files before
PR 2 starts:

```text
docs/wiso/tracer_data_model.md
claude-workflows/skills/tracer-conservation/SKILL.md
campaigns/wiso.campaign.md
```

Those updates must replace exact BFB criteria in PRs 2-4 with documented
non-BFB regression criteria. This is a campaign pause point and requires
human approval.

### Existing-test policy

Existing test intent must not be weakened.

- Prefer adding new tests rather than editing existing tests.
- Existing tests may be updated only when required by intentional interface
  changes such as field-rank or metadata changes.
- Any updated existing test must preserve the original scientific or
  software-contract assertion in an equivalent form.

## Skills referenced by every spec in this campaign

- `e3sm-build-and-test` — validation tiers, lint, cprnc, conservation
- `e3sm-platforms` with `references/docker-local.md` and
  `references/eamxx-doubly-periodic.md`
- `eamxx-cpp-conventions` — C++17, Kokkos, EKAT, field-manager tracers
- `scientific-modeling-conventions` — provenance, units, comments
- `tracer-conservation` — invariants and tolerance defaults
- `regression-baseline` — proving array-dim extensions are non-functional
  in single-tracer compatibility mode
- `spec-schema` — spec format enforced by the validator

Each spec body's "Approach" section names the subset it relies on.

## Conservation and regression invariants, canonical

These are referenced by spec success criteria through the
`tracer-conservation` skill.

| Invariant | Applies to | rtol | atol | Verified by |
|---|---|---:|---:|---|
| single-tracer slot-0 BFB vs pre-extension baseline | `tracer_count = 1`, if PR 1 BFB gate passes | 0 | 0 | `cprnc` |
| single-tracer slot-0 non-BFB fallback | only if PR 1 BFB gate explicitly fails and human approval is recorded | spec-defined | spec-defined | `cprnc` plus field norms |
| water-mass closure, column/timestep | all modes | 1e-12 | 1e-15 | `scripts/check_water_balance.py` |
| isotope-mass closure, per isotope tracer | isotope mode | 1e-10 | 1e-15 | `scripts/check_isotope_balance.py` |
| isotope ratio finite and physically bounded | isotope mode | spec-defined | spec-defined | `scripts/check_isotope_balance.py` |
| closed tagged-partition sum equals bulk H2O | only tracer groups with `kind = tagged_partition` and shared `partition_group_id` | 1e-8 | 1e-15 | `scripts/check_tagged_partition_balance.py` |
| alpha(H2O) = alpha(H216O) = 1 exactly | fractionation unit tests | 0 | 0 | unit tests |
| fractionation vs published values | fractionation unit tests | 1e-6 | 1e-9 | unit tests over `[233 K, 313 K]` |

Important exclusions:

- Isotope slots do **not** satisfy `sum(all species) == bulk H2O`.
- Bulk water is not reconstructed by summing isotope tracer slots.
- H216O is implicit unless the data-model document explicitly changes this.

## Branch and baseline model

- **Stacked chain.** Only PR 1 bases on `wiso-dev`. Every subsequent PR
  bases on the previous PR's branch, because each spec depends on the state
  produced by the spec before it.
- Branch naming:

  ```text
  wiso-NN-<spec-slug>
  ```

  where `NN` is the 1-indexed campaign position.

- Merge order matches chain order.
- When PR N merges into `wiso-dev`, PR N+1's base is rewritten from
  `wiso-NN-<spec>` to `wiso-dev` and rebased. Downstream PRs cascade the
  same way.
- `wiso-dev` itself is rebased onto `master` no more than once per month,
  immediately after a clean group boundary.
- Restack via:

  ```text
  /restack-campaign campaigns/wiso.campaign.md
  ```

- Baselines are named:

  ```text
  <branch>-<compset>-<resolution>
  ```

  under `baseline_root`.

- Each baseline directory carries a `BASELINE.txt` recording:

  - parent SHA
  - date
  - compset
  - resolution
  - stop options

- Baseline refresh decision tree lives in `regression-baseline`.
- Bump refreshed baselines with a `-rN` suffix; never overwrite.
- PR size policy:

  - default soft limit: 25 files / 600 lines
  - steps 2-4 may exceed this
  - manifest sets `pr.size_soft_limit.lines: 1500`
  - specs exceeding even that split into `Na` / `Nb` sub-specs

## Test layering

| Layer | Run cost | What it proves |
|---|---|---|
| Unit, Catch2 | seconds | Pure-function fractionation factors match published values |
| Component | minutes | A single parameterization's tracer hook is internally consistent |
| DP integration | about 10 min | End-to-end run with hooks active and conservation checked |
| Full-grid Tier-2 | hours | single-tracer BFB or approved non-BFB regression; multi-tracer conservation and scientific sanity |

Steps 6-9 stop at Unit.

Steps 9a-9b add Component validation.

Steps 10-14 add Component + DP validation.

Group boundaries trigger Tier-2 after:

```text
step 4
step 9b
step 14b
step 20
step 20c
```

## Spec groups

Specs are grouped so that group boundaries are natural Tier-2 checkpoints.
`pause_between_specs: false` within a group; `true` at group boundaries.

### Group 1 — tracer data model and structural array extension, PRs 1-4

Highest risk. PR 1 is a blocking feasibility gate. PRs 2-4 may not start
until PR 1 has landed `docs/wiso/tracer_data_model.md` and resolved the
single-tracer BFB strategy.

| PR | Spec id | Title | Tier | Notes |
|---|---|---|---:|---|
| 1 | `2026-MM-DD-water-tracer-data-model-and-bfb-gate` | Water-tracer metadata, `add_tracer`, data-model doc, and BFB feasibility gate | 0 | Header/types plus minimal field-layout prototype. Must decide exact-BFB path or approved fallback before PR 2. |
| 2 | `2026-MM-DD-extend-qv-multi-tracer` | Extend `qv` to tracer slot dimension | 2 | Enforce PR 1 BFB decision for `tracer_count = 1`. No isotope sum-to-bulk check. |
| 3 | `2026-MM-DD-extend-cloud-multi-tracer` | Extend `qc`, `qi` to tracer slot dimension | 2 | Same pattern as PR 2. |
| 4 | `2026-MM-DD-extend-precip-multi-tracer` | Extend `qr`, `qs` to tracer slot dimension | 2 | Same pattern as PR 2. |

Group 1 boundary:

- Tier-2 full-grid baseline regenerated with array machinery in place and
  `tracer_count = 1`.
- If exact BFB was approved by PR 1, this baseline must be BFB-equivalent
  to the pre-campaign baseline.
- If exact BFB fallback was approved by PR 1, this baseline records the
  accepted non-BFB field norms and becomes the campaign structural baseline
  for groups 2-3.

### Group 2 — fractionation primitives and utilities, PRs 6-9b

PR 5, phase-change enumeration, is moved out of the chain. See
"Out-of-chain work" below.

| PR | Spec id | Title | Tier | Notes |
|---|---|---|---:|---|
| 6 | `2026-MM-DD-equilibrium-fractionation` | Equilibrium fractionation functions | 0 | Unit tests vs published tabulated values. CMake-selectable schemes. Deliverables include `components/eamxx/tests/water_tracers/data/horita-wesolowski-1994.csv`, `data/majoube-1971.csv`, and `data/merlivat-nief-1967.csv`. |
| 7 | `2026-MM-DD-alpha-diff-functions` | Molecular-diffusivity fractionation | 0 | Unit tests with exact expected values. Deliverables include `data/merlivat-1978.csv` and `data/schoenemann-2014.csv`. |
| 8 | `2026-MM-DD-net-fractionation` | Net fractionation: Brutsaert + Craig-Gordon, Stewart 1975, Ciais-Jouzel 1994 | 0 | Unit tests. CMake-selectable. Includes partial-equilibration kinetics for falling hydrometeors, mirroring CAM `wtrc_equil_time`. |
| 9 | `2026-MM-DD-wtrc-ratio-utility` | `wtrc_ratio` utility plus `check_isotope_balance.py` | 0 | Blocking. Must implement isotope ratios using bulk slot 0 plus isotope slots, not species summation. |
| 9a | `2026-MM-DD-wiso-namelist` | YAML / runtime config for wiso flags | 1 | Mirrors CAM `wtrc_readnl`. Blocking for group-3 hooks that read these flags. |
| 9b | `2026-MM-DD-wtrc-mass-fixer` | Wiso-aware negative-mass / mass-conservation fixer | 1 | Must document whether corrections conserve isotope mass globally, columnwise, or ratio-preservingly. Must emit diagnostics for corrections. |

Required flags for PR 9a include:

```text
trace_water
wisotope
wtrc_lh2oadj
wtrc_lzmlin
wtrc_warn_only
wtrc_add_cvprecip
wtrc_add_stprecip
wtrc_alpha_kinetic
wtrc_check_total_h2o
wtrc_detrain_in_macrop
wtrc_use_ice_supsat
```

Mass-fixer policy for PR 9b:

- Do not silently enforce isotope sum-to-bulk closure.
- Do not independently clip isotope tracers in a way that creates
  unbounded isotope-ratio excursions.
- Record fixer tendencies or correction totals in diagnostics.
- Prefer ratio-preserving local fixes when physically justified.
- Clearly document any nonlocal conservation correction.

### Group 3 — parameterization hooks, PRs 10-14

Each PR uses the DP harness for Component validation. Isotope-mass closure
is required. Group boundary triggers Tier-2.

| PR | Spec id | Title | Tier |
|---|---|---|---:|
| 10 | `2026-MM-DD-ocean-evap-hook` | Ocean evaporation isotope partitioning | 1 |
| 10b | `2026-MM-DD-surface-flux-inputs` | Coupler-side delivery of wiso surface fluxes | 1 |
| 11 | `2026-MM-DD-shoc-hook` | SHOC phase-change hooks | 1 |
| 12 | `2026-MM-DD-p3-hook` | P3 microphysics phase-change hooks | 1 |
| 13 | `2026-MM-DD-zm-hook` | Zhang-McFarlane phase-change hooks, if present in target EAMxx configuration | 1 |
| 14a | `2026-MM-DD-rrtmgp-radiation-audit` | Audit RRTMGP / aerosol activation for unexpected phase changes | 0 |
| 14b | created if 14a finds gaps | Hooks for parameterizations identified by 14a | 1 |

Surface-flux PR 10b is a major integration milestone, not a small adjunct.
It must verify:

- coupler field names
- units
- sign conventions
- land/ocean/ice fractional areas
- behavior for F-compsets with prescribed ocean data
- fallback isotope composition when a surface component cannot provide
  isotope-aware fluxes

Group 3 boundary:

- Full-grid Tier-2 run with:

  ```text
  tracer_count = 4
  slots = bulk H2O, HDO, H2^18O, H2^17O
  ```

- Compare slot-0 bulk water against the group-2 boundary baseline using
  the PR 1-approved policy:

  - exact BFB only if the campaign explicitly retained BFB for this mode
  - otherwise non-BFB field norms and water closure

- Verify isotope-mass closure per isotope tracer.
- Verify isotope ratios are finite and within spec-defined physical bounds.
- Do **not** verify `sum(all slots) == bulk H2O` for isotope mode.

### Group 4 — auxiliary tracers and tagged tracers, PRs 14c, 15-20

| PR | Spec id | Title | Tier |
|---|---|---|---:|
| 14c | `2026-MM-DD-tritium-decay` | HTO radioactive decay using `WaterIsotopologues<Scalar>::hlhto` half-life constant | 0 |
| 15 | `2026-MM-DD-ch4-oxidation-hdo` | CH4 + OH to HDO impact, CMake flag | 1 |
| 16 | `2026-MM-DD-region-tagged-evap` | Region-tagged evaporation tracer, lat/lon box plus shapefile | 1 |
| 17a | `2026-MM-DD-sh-decomp-prototype` | Spherical-harmonic decomp prototype, `spec_type: analysis` | n/a |
| 17b | `2026-MM-DD-sh-decomp-impl` | Production implementation of SH / needlet tagged tracers | 1 |
| 18 | `2026-MM-DD-evap-tracer-multiplicative` | Multiplicative-factor evap tracers, time, T, etc. | 1 |
| 19 | `2026-MM-DD-condensation-tagged` | Condensation-tagged tracers, T and p at phase change | 1 |
| 20 | `2026-MM-DD-parcel-integrated` | Parcel-integrated quantities, such as `Q*dt` | 1 |

Step 17 is split into 17a, prototype / analysis, and 17b,
implementation, to avoid mixing research and implementation in one PR.

For every tagged-tracer PR:

- declare whether the tracer is `tagged_partition` or `tagged_diagnostic`
- if `tagged_partition`, declare `partition_group_id`
- only closed partition groups use the sum-to-bulk invariant
- diagnostic tags must define their own conservation or boundedness checks

### Group 4.5 — production support, PRs 20a-20c

Required before production science runs. This group runs after group 4 and
before the final user-guide pass. Group boundary triggers Tier-2.

| PR | Spec id | Title | Tier | Notes |
|---|---|---|---:|---|
| 20a | `2026-MM-DD-wiso-initial-conditions` | IC dataset plus reader for wiso fields | 1 | Mirrors CAM `const_init.F90` wiso paths. Default IC = VSMOW ratios applied to bulk; optional file input for spin-up state. |
| 20b | `2026-MM-DD-wiso-diagnostics` | delta-value derived diagnostics, R/RVSMOW conversions, column-integrated isotope quantities, history field registration | 0 | Store mass mixing ratios in history; provide derived-field utilities and Python helper for delta conversion. |
| 20c | `2026-MM-DD-wiso-restart` | Restart support for wiso fields | 1 | Required for any multi-segment run. |

Although listed here, minimal diagnostic hooks needed for debugging isotope
mass, isotope ratios, and mass-fixer corrections may be introduced earlier
by PR 9 or PR 9b.

### Group 5 — comprehensive harness, PR 21 split

| PR | Spec id | Title | Tier |
|---|---|---|---:|
| 21a | `2026-MM-DD-dp-harness` | DP test cases for the wiso campaign | 1 |
| 21b | `2026-MM-DD-wiso-user-guide` | Comprehensive wiso user guide plus test recipes | 0 |

PR 21a is scheduled in parallel with group 2 so that group 3 hooks can use
it. PR 21b is the final documentation pass.

## Out-of-chain work

These items live outside the campaign chain because they do not ship code
on a PR cadence:

- **Step 5, phase-change enumeration** — produced as:

  ```text
  docs/wiso/phase_changes.md
  ```

  through a single `spec_type: analysis` spec before group 3.

  Output is a markdown table grouped by phase-change type and
  parameterization. Used as the input checklist for PRs 10-14.

  The doc must include a CAM-to-EAMxx mapping section because EAMxx's SHOC
  unifies what CAM splits across:

  ```text
  clubb_intr.F90
  macrop_driver.F90
  convect_shallow.F90
  uwshcu.F90
  vertical_diffusion.F90
  ```

  PR 11, `shoc-hook`, must reference that mapping to confirm coverage.

- **References folder population** — one-time setup; not a PR.
- **Tier-2 baseline regeneration** — happens at group boundaries; the
  orchestrator runs it but it is not its own spec.

## Restack policy

Run:

```text
/restack-campaign campaigns/wiso.campaign.md
```

when:

- a PR in the chain merges
- `master` advances and contains a change that materially affects EAMxx
  physics or the build system
- a group boundary has just been cleared

The campaign-orchestrator's restack mode handles the mechanical rebase
across the whole stack.

Do not rebase individual PR branches by hand. The chain integrity check
verifies that each branch's base matches the manifest order.

## Post-completion review

Specs in group 1 enable `post_completion_review` with:

```text
performance-specialist
code-reviewer
```

Reason: adding tracer dimensions can affect Kokkos memory layout,
vectorization, field registration, restart/history size, and GPU memory
pressure.

Specs in group 4 also enable `performance-specialist`.

Specs in groups 2, 3, and 4.5 enable `code-reviewer` by default, and may
enable `performance-specialist` when a PR touches hot kernels.

## Validator expectations

Each spec must satisfy validator hard-fails H1-H12 in `spec-schema`:

- at least 3 runnable success criteria
- every shell/comparison/tolerance criterion has `cmd`
- every tolerance criterion has `rtol` or `atol` greater than 0
- BFB comparisons with `rtol = 0` and `atol = 0` must be represented as
  explicit `cprnc_exact` or equivalent non-tolerance criteria
- non-empty `out_of_scope`
- `ask_before` includes destructive operations
- CMake API changes touching shared field-manager headers count as
  `ask_before`
- no `TODO(spec):` left in the body

Soft-warns to watch:

- every criterion's `cmd` should reference a path that lives in
  `deliverables`
- no lint-only specs
- tests cited under `pytest <path>` or `<path>_tests` must exist on disk
  or be listed in deliverables

## Project-config template

The file at `.claude/project-config.yml` follows this shape. See the
committed copy in this repo for actual values.

```yaml
project_id: EAMxx-wiso
host_source_path: /code/E3SM/EAMXX-wiso
container_image: rfiorella/model-containers:e3sm-openmpi-dev-latest
default_platform: docker-local
target_compset: 
target_resolution: 
dp_compset: 
dp_resolution: 
baseline_root: /data/baselines
spec_mode: mutable
branch_prefix: wiso
base_branch: wiso-dev
remote: origin
```

The `<...>` placeholders are filled in by `/init-project`.