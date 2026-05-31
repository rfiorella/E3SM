# Campaign progress: 2026-05-29-wiso-group1-revised

## Environment
- Repo: /vast/home/rfiorella/E3SM
- Remote: origin (rfiorella/E3SM)
- Base branch: eamxx-wiso-dev
- Repo branch at start: wiso/02-extend-qv-tracer
- Started: 2026-05-29T00:00:00Z

## Specs

### 1/8 - 2026-05-28-water-tracer-metadata-and-gate (COMPLETE)
- Spec: specs/2026-05-28-water-tracer-metadata-and-gate.md
- Branch: wiso/01-water-tracer-metadata-and-gate (based on eamxx-wiso-dev)
- PR: (previously completed, not part of this revised campaign execution)
- Status: COMPLETE
- Notes: Completed before campaign revision. Gates passed: BFB feasible, performance acceptable.

### 2/8 - 2026-05-28-extend-qv-tracer-2a-infrastructure (COMPLETE)
- Spec: specs/2026-05-28-extend-qv-tracer-2a-infrastructure.md
- Branch: wiso/02a-extend-qv-tracer-infrastructure (based on wiso/01-water-tracer-metadata-and-gate)
- PR: (completed prior to campaign orchestrator execution)
- Status: COMPLETE
- Commit: 5bb940daeb Add tracer field infrastructure for EAMxx water tracers

### 3/8 - 2026-05-28-extend-qv-tracer-2b-shoc (COMPLETE)
- Spec: specs/2026-05-28-extend-qv-tracer-2b-shoc.md
- Branch: wiso/02b-extend-qv-tracer-shoc (based on wiso/02a-extend-qv-tracer-infrastructure)
- PR: (completed prior to campaign orchestrator execution)
- Status: COMPLETE
- Commit: b91dbc6653 Extend qv to tracer dimension - Part 2b: SHOC Process

### 4/8 - 2026-05-28-extend-qv-tracer-2c-p3 (COMPLETE)
- Spec: specs/2026-05-28-extend-qv-tracer-2c-p3.md
- Branch: wiso/02c-extend-qv-tracer-p3 (based on wiso/02b-extend-qv-tracer-shoc)
- PR: (completed prior to campaign orchestrator execution)
- Status: COMPLETE
- Commit: 797dac004e Extend qv to tracer dimension - Part 2c: P3 Microphysics

### 5/8 - 2026-05-28-extend-qv-tracer-2d-remaining (COMPLETE)
- Spec: specs/2026-05-28-extend-qv-tracer-2d-remaining.md
- Branch: wiso/02d-extend-qv-tracer-remaining (based on wiso/02c-extend-qv-tracer-p3)
- PR: (completed prior to campaign orchestrator execution)
- Status: COMPLETE
- Commit: a5ea3cd1d6 Extend qv to tracer dimension - Part 2d: Remaining Processes

### 6/8 - 2026-05-28-extend-cloud-tracer (COMPLETE)
- Spec: specs/2026-05-28-extend-cloud-tracer.md
- Branch: wiso/03-extend-cloud-tracer (based on wiso/02d-extend-qv-tracer-remaining)
- PR: (completed prior to campaign orchestrator execution)
- Status: COMPLETE
- Commit: a12c3206a0 Extend qc, qi, qm to tracer dimension

### 7/8 - 2026-05-28-extend-precip-tracer (COMPLETE)
- Spec: specs/2026-05-28-extend-precip-tracer.md
- Branch: wiso/04-extend-precip-tracer (based on wiso/03-extend-cloud-tracer)
- PR: https://github.com/rfiorella/E3SM/pull/9
- Status: COMPLETE
- Commit: 05da046ec4 Extend qr and surface precipitation to tracer dimension
- Notes: Spec corrected 2026-05-30 - scope reduced to qr (3D), precip_liq_surf_mass (2D), precip_ice_surf_mass (2D). P3 does not use separate qs/rain/snow fields.

### 8/8 - 2026-05-28-tracer-validation (PENDING)
- Spec: specs/2026-05-28-tracer-validation.md
- Branch: (not yet created)
- PR: (not yet created)
- Status: PENDING
