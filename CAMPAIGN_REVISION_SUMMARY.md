# Campaign Revision Summary: wiso-group1

**Date**: 2026-05-29  
**Reason**: Spec 2 scope mismatch discovered during investigation phase

## What Happened

Campaign `wiso-group1` was launched with 5 specs (PRs 1-5). During execution:

1. **PR 1 COMPLETED** ✅
   - Water tracer metadata and BFB feasibility gate
   - Both gates PASSED (BFB exact match, no performance overhead)
   - Branch: `wiso/01-water-tracer-metadata-and-gate`
   - GitHub PR: https://github.com/rfiorella/E3SM/pull/3

2. **PR 2 HALTED** 🛑
   - "Extend qv to tracer dimension" 
   - Investigation phase revealed scope mismatch:
     - Spec listed 14 deliverables
     - Investigation found ~40+ files need modification
     - Architectural decisions needed before implementation
   - Campaign halted per orchestrator protocol

## Resolution: Split Spec 2 into Sub-Specs

You selected **Option C** (split into sub-specs) with **hybrid accessor pattern validation**.

### New Structure: Spec 2 → Specs 2a, 2b, 2c, 2d

| Spec | Title | Scope | Purpose |
|------|-------|-------|---------|
| **2a** | Field Infrastructure | FieldTag, grid layouts, accessor patterns | Establish infrastructure, validate both EXPLICIT and SUBVIEW patterns for BFB, choose one for 2b-2d |
| **2b** | SHOC Process | PBL turbulence | First physics process, establishes pattern |
| **2c** | P3 Process | Microphysics | Validates pattern for complex multi-kernel processes |
| **2d** | Remaining Processes | RRTMGP, ZM, HOMME, surface | Completes qv extension, full-atmosphere integration test |

### Campaign Structure: 5 PRs → 8 PRs

**Original:**
1. Metadata + gate
2. Extend qv (MONOLITHIC)
3. Extend cloud water
4. Extend precip
5. Validation test

**Revised:**
1. Metadata + gate ✅ DONE
2a. Field infrastructure
2b. SHOC process
2c. P3 process
2d. Remaining processes
3. Extend cloud water (unchanged)
4. Extend precip (unchanged)
5. Validation test (unchanged)

## Architectural Decisions Documented

These apply to PRs 2a-5:

### Decision 1: Compile-time tracer count
**Choice**: `SCREAM_NUM_TRACERS` as compile-time constant  
**Rationale**: PR 1 BFB gate passed with this. Simpler implementation.

### Decision 2: Grid layout method
**Choice**: New `get_3d_tracer_layout()` method  
**Rationale**: Clearer semantics than repurposing vector layout.

### Decision 3: Accessor pattern (HYBRID)
**Choice**: Implement BOTH patterns in spec 2a, test both for BFB, use winner in 2b-2d  
**Patterns**:
- **EXPLICIT**: `qv(0, icol, ilev)` everywhere (safer for BFB)
- **SUBVIEW**: `Kokkos::subview(qv, 0, ALL, ALL)` (cleaner, potentially faster)

**Selection criteria**:
- If both pass BFB → use SUBVIEW (performance + readability)
- If only EXPLICIT passes → use EXPLICIT (safety)
- If neither passes → HALT (architectural issue)

**Result**: TBD - will be determined in spec 2a execution

## Files Created

### New Specs
- `specs/2026-05-28-extend-qv-tracer-2a-infrastructure.md`
- `specs/2026-05-28-extend-qv-tracer-2b-shoc.md`
- `specs/2026-05-28-extend-qv-tracer-2c-p3.md`
- `specs/2026-05-28-extend-qv-tracer-2d-remaining.md`

### Updated Campaign
- `campaigns/wiso-group1-revised.campaign.md` (new manifest)
- Original: `campaigns/wiso-group1.campaign.md` (superseded)

### Investigation Artifacts
- `specs/2026-05-28-extend-qv-tracer.progress.md` (investigation findings)
- `campaigns/2026-05-28-wiso-group1.progress.md` (campaign ledger)

## Next Steps

### Option A: Resume Campaign with Revised Manifest
```bash
/run-campaign campaigns/wiso-group1-revised.campaign.md
```

This will:
1. Skip PR 1 (already complete)
2. Execute specs 2a, 2b, 2c, 2d in sequence
3. Continue with PRs 3, 4, 5 (unchanged)

### Option B: Manual Execution of Sub-Specs
Execute each sub-spec individually:
```bash
/run-spec specs/2026-05-28-extend-qv-tracer-2a-infrastructure.md
# Review results, then:
/run-spec specs/2026-05-28-extend-qv-tracer-2b-shoc.md
# etc.
```

**Recommendation**: **Option A** - let the campaign orchestrator handle the stacked PR chain.

## Benefits of This Revision

1. **Better code review**: Each PR is now ~10-15 files, <1500 lines (within guidelines)
2. **Faster BFB debugging**: If BFB breaks, smaller scope makes bisection easier
3. **Incremental validation**: Infrastructure validated before applying to all processes
4. **Pattern establishment**: SHOC proves pattern before applying to 6 other processes
5. **Risk management**: P3 (most complex) gets its own PR with focused attention

## Current Repository State

- **Branch**: `wiso/02-extend-qv-tracer` (created but empty - can be deleted)
- **Working tree**: Clean
- **Ready to resume**: Yes, with revised campaign manifest
- **PR 1 status**: Merged? Draft? (check GitHub)

## Questions for You

Before resuming:

1. **Should we resume with the revised campaign now?** (Option A)
2. **Do you want to review the new specs first?**
3. **Should PR 1 be merged before proceeding?** (It's currently draft)
4. **Any concerns about the 5→8 PR expansion?**

Let me know how you'd like to proceed!
