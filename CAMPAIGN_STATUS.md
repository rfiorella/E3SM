# Campaign Status: wiso-group1-revised

**Last Updated**: 2026-05-29  
**Campaign**: wiso-group1-revised (8 specs, revised from original 5)

## Progress Summary

| Spec | Status | Branch | PR | Notes |
|------|--------|--------|---|-------|
| 1 | ✅ COMPLETE | wiso/01-water-tracer-metadata-and-gate | #3 | Metadata + BFB gates PASSED |
| 2a | 🔄 CODE COMPLETE | wiso/02a-extend-qv-tracer-infrastructure | #4 | Needs validation on E3SM machine |
| 2b | ⏸️ BLOCKED | - | - | Waiting for 2a validation results |
| 2c | ⏸️ PENDING | - | - | Depends on 2b |
| 2d | ⏸️ PENDING | - | - | Depends on 2c |
| 3 | ⏸️ PENDING | - | - | Depends on 2d |
| 4 | ⏸️ PENDING | - | - | Depends on 3 |
| 5 | ⏸️ PENDING | - | - | Depends on 4 |

## Current Blocker: Spec 2a Validation

### What's Needed

Spec 2a code is complete (PR #4) but requires validation on an E3SM-configured machine to determine the accessor pattern for subsequent specs.

**Validation script**: `/vast/home/rfiorella/E3SM/components/eamxx/tests/water_tracers/validate_implementation.sh`

**Requirements**:
- E3SM machine with proper configuration (CMake, Kokkos, MPI, build environment)
- Machine examples: Perlmutter, Chrysalis, Compy, or docker-local setup
- Current machine (darwin-fe3.lanl.gov) lacks E3SM build environment

### What the Validation Tests

The script runs BFB tests for both accessor patterns:
1. **EXPLICIT**: `qv(0, icol, ilev)` indexing
2. **SUBVIEW**: `Kokkos::subview(qv, 0, ALL, ALL)` accessor

**Result determines implementation strategy for specs 2b-2d.**

## Options to Continue

### Option 1: Run Validation on E3SM Machine (Recommended)

On a configured E3SM machine (e.g., Perlmutter):
```bash
cd /vast/home/rfiorella/E3SM
git checkout wiso/02a-extend-qv-tracer-infrastructure
cd components/eamxx
./tests/water_tracers/validate_implementation.sh
```

Then return here with results and continue campaign.

### Option 2: Default to EXPLICIT Pattern

Proceed with EXPLICIT accessor pattern (conservative choice) without validation:
- Lower risk for BFB preservation
- More verbose code
- Can be changed later if SUBVIEW proves better

### Option 3: Parallel Development

While waiting for validation machine access:
- Continue implementing specs 2b-2d with placeholder pattern
- Update once validation complete
- Risk: may need to refactor if wrong pattern chosen

## What Happens After Validation

Once validation completes and accessor pattern is determined:

1. Update spec 2b-2d with chosen pattern (replace `{{PATTERN}}` placeholder)
2. Continue campaign execution:
   - Spec 2b: SHOC process
   - Spec 2c: P3 microphysics
   - Spec 2d: Remaining processes
   - Specs 3-4: Extend other water species
   - Spec 5: Validation test

## Resuming Campaign

After validation, resume with:
```bash
cd /vast/home/rfiorella/E3SM
# Update spec 2b with accessor pattern result
# Then continue:
/run-spec specs/2026-05-28-extend-qv-tracer-2b-shoc.md
# Or resume full campaign:
/run-campaign campaigns/wiso-group1-revised.campaign.md
```

## Files & Artifacts

- **Campaign manifest**: `campaigns/wiso-group1-revised.campaign.md`
- **Spec files**: `specs/2026-05-28-extend-qv-tracer-2*.md` (2a-2d)
- **Progress ledgers**: `specs/*.progress.md`, `campaigns/*.progress.md`
- **Summary documents**: `CAMPAIGN_REVISION_SUMMARY.md`, this file

## Contact

If you have access to an E3SM-configured machine and can run the validation, please do so and report results. The campaign can then continue immediately.
