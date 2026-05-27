---
spec_id: 2026-05-26-phase-change-enumeration
spec_type: analysis
spec_version: 1
title: "EAMxx phase-change enumeration: catalog of water phase-change sites per parameterization"
created: 2026-05-26T00:00:00-06:00
author: rfiorella
project: EAMxx-wiso

inputs:
  source_files:
    - components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
    - components/eamxx/src/physics/p3/impl/p3_main_impl.hpp
    - components/eamxx/src/physics/p3/impl/p3_main_impl_part2.hpp
    - components/eamxx/src/physics/p3/impl/p3_main_impl_part3.hpp
    - components/eamxx/src/physics/shoc/eamxx_shoc_process_interface.cpp
    - components/eamxx/src/physics/shoc/impl/shoc_main_impl.hpp
    - components/eamxx/src/physics/zm/eamxx_zm_process_interface.cpp
    - components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.cpp
    - components/eamxx/src/physics/rrtmgp/eamxx_rrtmgp_process_interface.cpp
    - components/eamxx/src/physics/surface_coupling/surface_coupling_utils.hpp
    - wiso_campaign_plan.md
    - campaigns/wiso.campaign.md
  data: []
  baseline: null

deliverables:
  - "docs/wiso/phase_changes.md — markdown table enumerating every water phase-change emission site in EAMxx atmosphere physics, grouped by parameterization (P3, SHOC, ZM, ocean evap, RRTMGP audit) and phase-change type (vapor↔liquid, vapor↔ice, liquid↔ice, surface-flux vapor source). Required columns: parameterization, phase-change type, source species, destination species, source file, line number, function/kernel name, callsite context (in macrostep vs in substep), notes on whether mid-step T/p trajectory is preserved at the callsite (load-bearing for fractionation correctness). Header section cites the wiso_campaign_plan.md group-3 PRs that consume this checklist (PRs 10-14)."
  - "scripts/wiso/audit_phase_changes.py — Python script that greps the EAMxx source tree for candidate phase-change tendency emissions (regex on tendency-variable naming conventions: qX2qY, qX_to_qY, evap_, sublim_, cond_, dep_, melt_, freeze_, riming_), writes a candidate list to docs/wiso/phase_changes.candidates.csv, and prints a coverage report comparing candidates to the curated docs/wiso/phase_changes.md. Exit nonzero if a candidate is not categorized in the curated doc (either as in-scope or explicitly excluded). The doc-vs-candidates diff is the regression check that catches new phase-change paths added upstream."
  - "scripts/wiso/test_audit_phase_changes.py — pytest tests for the audit script: (a) a fixture EAMxx source snippet contains qX2qY tendency assignments, audit script lists them; (b) a curated doc covering all candidates exits zero; (c) a curated doc missing one candidate exits nonzero and names the missing site."
  - "manifest.lock — pip-tools style frozen environment file pinning the analysis environment (python, pytest, ruff). Lives under scripts/wiso/ alongside the script."

success_criteria:
  - id: doc-exists
    type: shell
    cmd: "test -f docs/wiso/phase_changes.md"
    expect: exit_zero
    phase: implementation
    verifies:
      - deliverable: docs/wiso/phase_changes.md

  - id: doc-has-required-sections
    type: assertion
    cmd: "grep -c '^## ' docs/wiso/phase_changes.md"
    expect: "5"
    phase: implementation
    verifies:
      - claim: "doc has one second-level section per group-3 parameterization target (P3, SHOC, ZM, ocean-evap, RRTMGP-audit)"

  - id: doc-cites-each-group3-pr
    type: shell
    cmd: "grep -qE 'PR ?10' docs/wiso/phase_changes.md && grep -qE 'PR ?11' docs/wiso/phase_changes.md && grep -qE 'PR ?12' docs/wiso/phase_changes.md && grep -qE 'PR ?13' docs/wiso/phase_changes.md && grep -qE 'PR ?14' docs/wiso/phase_changes.md"
    expect: exit_zero
    phase: implementation

  - id: audit-script-exists
    type: shell
    cmd: "test -x scripts/wiso/audit_phase_changes.py"
    expect: exit_zero
    phase: implementation
    verifies:
      - deliverable: scripts/wiso/audit_phase_changes.py

  - id: audit-script-runs-clean
    type: shell
    cmd: "python scripts/wiso/audit_phase_changes.py --repo-root . --doc docs/wiso/phase_changes.md"
    expect: exit_zero
    phase: testing
    verifies:
      - claim: "every grepped candidate appears in the curated doc (either categorized or explicitly excluded)"

  - id: pytest-audit-script
    type: shell
    cmd: "pytest scripts/wiso/test_audit_phase_changes.py -v"
    expect: exit_zero
    phase: testing
    verifies:
      - deliverable: scripts/wiso/test_audit_phase_changes.py

  - id: ruff-lint
    type: shell
    cmd: "ruff check scripts/wiso/"
    expect: exit_zero
    gate: advisory
    on_fail: skip
    phase: implementation

  - id: p3-substep-callsites-flagged
    type: shell
    cmd: "grep -qE 'update_prognostic_(ice|liquid)' docs/wiso/phase_changes.md && grep -qE 'mid-step T|temperature trajectory|substep' docs/wiso/phase_changes.md"
    expect: exit_zero
    phase: implementation
    verifies:
      - claim: "doc explicitly flags P3 update_prognostic_ice / update_prognostic_liquid as mid-substep T-mutating callsites (the load-bearing reason group 3 hooks intercept inside substeps rather than reading macrostep tendency fields)"

  - id: reproducibility-manifest-present
    type: shell
    cmd: "test -f scripts/wiso/manifest.lock"
    expect: exit_zero
    phase: implementation

  - id: doc-completeness-checkpoint
    type: human-review
    description: "Read docs/wiso/phase_changes.md end-to-end. Confirm: (a) every parameterization slated for a group-3 hook PR has at least one phase-change entry; (b) the in-substep vs macrostep distinction is correctly annotated; (c) ocean evaporation (surface flux source for vapor) is included even though it is not a microphysical phase change in the traditional sense; (d) RRTMGP audit row is honest about what was checked and what was concluded; (e) any phase-change path the audit script flagged but the doc excluded includes a one-line reason."
    phase: integration

out_of_scope:
  - "Any implementation of fractionation hooks. This spec produces the input checklist; PRs 10-14 consume it."
  - "Changes to P3 / SHOC / ZM / RRTMGP source code (read-only audit)."
  - "Quantitative phase-change rate estimates (the doc lists where phase changes happen, not how often)."
  - "EAM (legacy Fortran) phase-change cataloging — this spec is EAMxx-only."
  - "Coupling phase changes (e.g., land surface evaporation routed through the coupler) — flagged in the doc but cataloged in the coupled-component documentation, not here."
  - "Performance or vectorization characterization of phase-change kernels."

resolved_decisions:
  - "spec_type is analysis per the campaign plan (wiso_campaign_plan.md, 'Out-of-chain work' section). Deliverable is markdown documentation backed by an auditable extraction script — the script makes the doc verifiable rather than aspirational."
  - "Audit script is a candidate-finder, not a categorizer. Categorization (which phase-change paths are in-scope for group-3 hooks vs. excluded) is a judgment call documented in the curated doc. The script's job is to ensure the human-curated doc covers every grepped candidate; that lets the doc stay current as upstream evolves."
  - "P3 update_prognostic_ice / update_prognostic_liquid sites get explicit flagging in the doc. Per references/old_specs/2026-05-06-eamxx-water-condensates-add-tracer-dim.md resolved_decisions, these calls update th_atm mid-substep, so macrostep tendency fields have lost the temperature trajectory needed for fractionation. The doc's annotation is what group-3 hook specs cite to justify intercepting in-substep rather than post-hoc."
  - "Surface evaporation (ocean → atmosphere vapor flux) is included even though it is not a microphysical phase change. Reason: from the isotope perspective it is the dominant vapor source globally and the kinetic fractionation that occurs at the ocean-air interface is what PR 10 (ocean-evap-hook) implements."
  - "Per python-analysis-conventions skill (loaded implicitly via spec_type=analysis): pytest test discipline, ruff lint, frozen environment manifest, no global state in the script."

ask_before:
  - Modifying any EAMxx source file (this is a read-only audit; if a phase-change site is genuinely missing or misnamed, file a separate issue rather than patching here).
  - Adding any new analysis runtime to scripts/wiso/ beyond Python + pytest + ruff (the environment manifest is small by design).
  - Removing or relocating docs/wiso/ outputs after they exist (downstream group-3 specs cite the path).
  - Pushing the audit script to any shared CI workflow beyond the spec's own success criteria.

execution:
  mode: checkpoint
  checkpoints:
    - after: planning
      requires: user-confirmation
    - after: implementation
      requires: tests-pass
    - after: integration
      requires: user-confirmation
  max_iterations_per_phase: 5
  parallelization:
    allowed: false
    max_subagents: 1
    plan: null

post_completion_review:
  enabled: false
  reviewers: []
  blocking: false

analysis_specific:
  language: python
  environment: e3sm-py311
  test_strategy: test-along
  manifest_required: true
  performance_budget: null
  output_type: script
---

# EAMxx phase-change enumeration

## Context

This spec is the out-of-chain "Step 5" from
`wiso_campaign_plan.md`. It does not produce a code change to
EAMxx; it produces the input checklist that the four group-3 hook
PRs (PRs 10-14: ocean-evap, SHOC, P3, ZM, RRTMGP audit) consume.
Without this enumeration each group-3 spec has to re-derive which
phase-change sites need an isotope hook and which do not, and the
enumeration is likely to drift between specs.

The work is read-only: grep the EAMxx physics tree for water-
phase-change tendency emissions, categorize each site by phase-
change type (vapor↔liquid, vapor↔ice, liquid↔ice, surface-flux
vapor source) and parameterization, and capture the callsite
context that matters for isotope physics — most importantly
whether the callsite is inside a substep that mutates temperature.

The substep distinction is load-bearing. P3's level loop calls
`update_prognostic_ice` and `update_prognostic_liquid` and these
update `th_atm` — hence `T_atm` — mid-substep
(`p3_main_impl_part2.hpp` around lines 410-450). The registered
macrostep tendency fields (`qr2qv_evap`, `qi2qv_sublim`,
`qc2qr_accret`, etc., at `eamxx_p3_process_interface.cpp:129-143`)
are aggregates that have already lost the temperature trajectory
across substeps. Equilibrium and kinetic fractionation factors
depend on `T`, so a hook that reads macrostep tendencies and
applies fractionation at macrostep `T` gives the wrong isotopic
delta. Group-3 hook PRs need to know which paths are macrostep-
safe and which require in-substep interception; the doc
delivered here is what tells them.

Why an audit script and not just a hand-curated markdown table?
The doc rots the moment a new phase-change tendency lands
upstream. A grep-based candidate finder keeps the doc honest: if
the script flags a candidate the doc does not cover, CI (or the
campaign-orchestrator) fails and the doc has to be updated. The
script does not categorize — categorization is judgment — but it
ensures coverage.

## Approach

1. **Grep pass over EAMxx physics.** `audit_phase_changes.py`
   walks `components/eamxx/src/physics/` looking for tendency
   variable names matching the convention regex:
   `qX2qY` (e.g., `qr2qv_evap`), `qX_to_qY`, `evap_`, `sublim_`,
   `cond_`, `dep_`, `melt_`, `freeze_`, `riming_`. For each match
   it records source file, line number, surrounding function /
   kernel name, and the immediate parameterization directory.
   Output is `docs/wiso/phase_changes.candidates.csv`.

2. **Categorize each candidate.** Hand-curate
   `docs/wiso/phase_changes.md` with one second-level section per
   group-3 parameterization target (P3, SHOC, ZM, ocean-evap,
   RRTMGP-audit), one table row per phase-change site, columns:
   parameterization, phase-change type, source species,
   destination species, source file, line, function name,
   callsite context (macrostep vs in-substep), notes on
   temperature-trajectory preservation, target group-3 PR.

3. **Annotate substep callsites.** Explicitly flag
   `update_prognostic_ice` / `update_prognostic_liquid` and any
   other site that mutates `T` mid-substep. Cite the line in
   `p3_main_impl_part2.hpp`. This is the load-bearing annotation
   for PR 12 (p3-hook).

4. **Ocean-evap row.** Surface vapor flux from the ocean is not a
   microphysical phase change in the traditional sense, but it is
   the dominant kinetic-fractionation site at the air-sea
   interface. Include it as its own section with a pointer to PR
   10.

5. **RRTMGP audit row.** PR 14a is itself an audit ("any
   unexpected phase changes in RRTMGP / aerosol activation").
   This doc's RRTMGP section is the structured handoff: list
   every place RRTMGP touches a water species and annotate what
   was concluded (in-scope, out-of-scope, needs-hook).

6. **Coverage check.** `audit_phase_changes.py` re-runs after the
   doc is curated, comparing candidates to doc entries. Exits
   nonzero if a candidate is uncovered (neither categorized nor
   explicitly excluded with a reason). Exits zero when coverage
   is complete.

7. **Pytest discipline.** `test_audit_phase_changes.py` covers
   the script's edge cases: empty candidate file, candidate not
   in doc, candidate in doc but flagged excluded, doc with
   bogus excluded reason.

8. **Frozen environment.** `manifest.lock` pins Python, pytest,
   ruff versions for reproducibility per the
   `python-analysis-conventions` skill.

Risks:

- **False negatives in the regex.** A phase-change site that uses
  unconventional naming (no `qX2qY` pattern, no
  `evap/sublim/cond/dep/melt/freeze/riming` token) will not show
  up as a candidate, will not be flagged, and may quietly miss
  isotope hooks in group 3. Mitigation: the regex is documented
  in the script, the curated doc has a "manually identified, not
  caught by regex" section for known additions, and the human-
  review checkpoint at the end is the catch-all.
- **Doc drift.** Once group 3 starts, every hook PR adds to or
  changes the doc. The coverage script keeps the doc honest
  against the source tree, but it does not prevent group-3 PRs
  from contradicting each other. Mitigation: each group-3 spec
  has the doc in its `inputs.source_files` so changes surface in
  diffs.
- **Substep distinction misread.** If the doc misclassifies a
  callsite as macrostep-safe when it is actually in-substep, a
  group-3 hook PR will silently produce wrong isotope physics.
  The human-review checkpoint covers this; the cited line
  numbers in `p3_main_impl_part2.hpp` are the auditable
  evidence.

Skills relied on: `python-analysis-conventions` (test discipline,
ruff, manifest), `scientific-modeling-conventions` (provenance,
no magic numbers, named constants), `spec-schema` (this spec's
structure).

## References

- `wiso_campaign_plan.md` — "Out-of-chain work" section
  defining this spec.
- `campaigns/wiso.campaign.md` — group 3 PR list (PRs 10-14)
  that consume the doc.
- `components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp:129-143`
  — P3 phase-change tendency registration list (`qr2qv_evap`,
  `qi2qv_sublim`, `qc2qr_accret`, `qc2qr_autoconv`,
  `qv2qi_vapdep`, `qc2qi_berg`, `qc2qr_ice_shed`,
  `qc2qi_collect`, `qr2qi_collect`, `qc2qi_hetero_freeze`,
  `qr2qi_immers_freeze`, `qi2qr_melt`, plus sedimentation).
- `components/eamxx/src/physics/p3/impl/p3_main_impl_part2.hpp`
  (around 410-450) — `update_prognostic_ice` /
  `update_prognostic_liquid`; cited as the load-bearing
  in-substep T-mutating callsites.
- `components/eamxx/src/physics/p3/impl/p3_main_impl.hpp:295, 303, 310`
  — adaptive sedimentation substepping comments.
- `references/old_specs/2026-05-06-eamxx-water-condensates-add-tracer-dim.md`
  resolved_decisions — prior articulation of the macrostep-vs-
  substep distinction inherited here.
- `claude-workflows/skills/python-analysis-conventions/SKILL.md`.
