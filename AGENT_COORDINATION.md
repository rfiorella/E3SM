# Agent Coordination Framework

**Project**: P3 Water Isotope Implementation  
**Strategy**: 3-Agent Parallel Development  
**Expected Timeline**: 5-6 months  
**Created**: 2026-04-23

---

## 🎯 Agent Assignments

### Agent 1: "Stratiform" (Critical Path Lead)
**GitHub**: `@agent-stratiform`  
**Branch**: `agent1-stratiform`  
**Phases**: 2 → 3 → 4 → 7  
**Duration**: 75-95 days (15-19 weeks)  
**Priority**: 🔴 HIGHEST

**Responsibilities**:
- Phase 2: Rain evaporation (Stewart model)
- Phase 3: Vapor-liquid phase changes
- Phase 4: Vapor-ice phase changes
- Phase 7: Sedimentation

**Primary Files**:
- `components/eam/src/physics/p3/eam/micro_p3.F90`
- `components/eam/src/physics/p3/eam/micro_p3_interface.F90`
- `components/eam/src/physics/p3/eam/water_tracers_p3.F90` (create)
- `components/eam/src/physics/cam/water_tracers.F90` (extend)

**Key Deliverable**: Rain evaporation fractionation working by Week 4 (enables validation)

---

### Agent 2: "Convection" (Independent Track)
**GitHub**: `@agent-convection`  
**Branch**: `agent2-convection`  
**Phase**: 8  
**Duration**: 12-16 days (2-3 weeks)  
**Priority**: 🟠 HIGH

**Responsibilities**:
- Phase 8: Deep convection (ZM scheme)
- Phase 8: Shallow convection (UW scheme)

**Primary Files** (NO OVERLAP with Agent 1):
- `components/eam/src/physics/cam/zm_conv_intr.F90`
- `components/eam/src/physics/cam/uwshcu.F90` or `convect_shallow.F90`
- `components/eam/src/physics/cam/water_tracers.F90` (different functions)

**Key Deliverable**: Convection isotopes ported by Week 6 (then waits for integration)

---

### Agent 3: "Collection" (Delayed Start)
**GitHub**: `@agent-collection`  
**Branch**: `agent3-collection`  
**Phases**: 5 → 6  
**Duration**: 30-40 days (6-8 weeks)  
**Start**: After Agent 1 completes Phase 4 (~Week 12)  
**Priority**: 🟡 MEDIUM

**Responsibilities**:
- Phase 5: Freezing and melting
- Phase 6: Collection and riming

**Primary Files** (PARTIAL OVERLAP with Agent 1):
- `components/eam/src/physics/p3/eam/micro_p3.F90` (different functions)
- `components/eam/src/physics/p3/eam/water_tracers_p3.F90` (extend Agent 1's work)

**Key Deliverable**: Freezing/collection processes working by Week 18

---

## 📅 Timeline & Milestones

### Week 1-3: Phase 1 (ALL AGENTS)
**Status**: BLOCKING - All agents must complete Phase 1 together

**Daily Standup Topics**:
- Module porting progress
- Compilation issues
- Unit test results

**Week 3 Deliverable**: 
- [ ] All Phase 1 modules ported
- [ ] Model compiles with `wisotope=.true.`
- [ ] Unit tests passing
- [ ] **READY TO SPLIT INTO PARALLEL TRACKS**

---

### Week 4-7: Parallel Work Begins

#### Agent 1: Phase 2 (Rain Evaporation)
**Critical Path Item** - Highest priority

**Week 4**: 
- Modify `evaporate_rain()` function
- Create `wtrc_p3_inter_phase2()` interface

**Week 5-6**: 
- Port Stewart model functions
- Implement equilibration calculation

**Week 7**: 
- Testing and validation
- **Milestone**: Rain evaporation depletes δD correctly

#### Agent 2: Phase 8 (Convection)
**Independent Track** - Can start immediately after Phase 1

**Week 4**: 
- Port `wtrc_q1q2_pjr()` for ZM deep convection

**Week 5**: 
- Port `wtrc_shallow()` for shallow convection

**Week 6**: 
- Unit testing
- **Milestone**: Convection isotopes ported
- **Status**: WAITING FOR INTEGRATION (Phase 9)

**Week 7+**: Agent 2 is done, can help with documentation or Phase 9 prep

---

### Week 8-11: Agent 1 Continues

#### Agent 1: Phase 3 (Vapor-Liquid)
**Week 8-9**: 
- Modify autoconversion and accretion
- Port Rayleigh distillation
- **Milestone**: Condensation fractionation working

#### Agent 1: Phase 4 (Vapor-Ice)
**Week 10-11**: 
- Modify ice nucleation and deposition
- Implement ice-vapor fractionation
- **Milestone**: Ice processes working
- **TRIGGER**: Agent 3 can now start

---

### Week 12-17: Three Agents Working in Parallel

#### Agent 1: Phase 7 (Sedimentation)
**Week 12-16**: 
- Modify rain/ice sedimentation functions
- Port `wtrc_sediment()` function (complex!)
- Layer-by-layer fractionation
- **Milestone**: Precipitation isotopes match observations qualitatively

#### Agent 3: Phases 5+6 (Starts Week 12)
**Week 12-14**: Phase 5 - Freezing and melting
**Week 15-17**: Phase 6 - Collection and riming
**Milestone**: Phase transitions and collection working

**Coordination Required**: Agents 1 & 3 both modifying `micro_p3.F90`
- See "Merge Conflict Resolution" section below

---

### Week 18-20: Integration (Phase 9)

**ALL AGENTS** reconvene for final integration

**Week 18**: 
- Merge all agent branches
- Resolve conflicts
- Integration testing

**Week 19-20**: 
- Port conservation checks
- Port diagnostic output
- Global validation simulation

**Final Milestone**: 
- [ ] Mass conservation < 0.1%
- [ ] Precipitation isotopes match GNIP observations
- [ ] δD-δ18O meteoric water line slope ~8

---

## 🔄 Daily Coordination

### Daily Standup (Async - 5 min each agent)
**Time**: 9:00 AM local time  
**Format**: Post to shared Slack channel or GitHub Discussions

**Template**:
```
Agent: [1/2/3]
Phase: [Current phase number]
Yesterday: [What I completed]
Today: [What I'm working on]
Blockers: [Any issues or questions]
Files Modified: [List of files touched yesterday]
```

**Example**:
```
Agent: 1
Phase: 2 (Rain Evaporation)
Yesterday: 
  - Ported stewart_isoevap() function
  - Created unit tests for equilibration calculation
Today: 
  - Integrate Stewart model into wtrc_apply_rates()
  - Test with single-column case
Blockers: 
  - Need clarification on parameter φ (tuning factor) - see ISOTOPE_PHYSICS.md
Files Modified: 
  - components/eam/src/physics/cam/water_tracers.F90 (+150 lines)
```

---

## 📞 Weekly Sync Meeting

**Time**: Friday 2:00 PM  
**Duration**: 30 minutes  
**Format**: Video call (required) + written notes

### Agenda Template

**1. Progress Review** (10 min)
- Agent 1: Status on critical path
- Agent 2: Status on convection
- Agent 3: Status on collection (if started)

**2. Integration Check** (10 min)
- Merge conflicts this week?
- Integration test results
- Any interface changes that affect other agents?

**3. Next Week Planning** (5 min)
- What will each agent start/finish?
- Any coordination needed?
- Upcoming milestones

**4. Blockers & Help Needed** (5 min)
- Science questions (defer to ISOTOPE_PHYSICS.md)
- Technical issues
- Resource needs

---

## 🌿 Branch Strategy

### Main Branches
```
main (or master)
  └── feature/p3-isotopes
        ├── feature/p3-isotopes-phase1 (Phase 1 base - Week 1-3)
        │
        ├── agent1-stratiform (Agent 1's work)
        │     ├── agent1-phase2
        │     ├── agent1-phase3
        │     ├── agent1-phase4
        │     └── agent1-phase7
        │
        ├── agent2-convection (Agent 2's work)
        │
        └── agent3-collection (Agent 3's work)
              ├── agent3-phase5
              └── agent3-phase6
```

### Branch Naming Convention
```
agent[1-3]-[phase-name]      # Main agent branch
agent[1-3]-phase[N]          # Specific phase branch (optional)
agent[1-3]-bugfix-[issue]    # Bug fixes
agent[1-3]-test-[feature]    # Experimental work
```

---

## 🔀 Merge Protocol

### Daily: Pull from Base
**Every morning before starting work**:
```bash
git checkout agent1-stratiform  # Your branch
git fetch origin
git merge origin/feature/p3-isotopes-phase1  # Get Phase 1 updates
git push origin agent1-stratiform
```

### Every 2-3 Days: Share Your Work
```bash
# Commit your work
git add -A
git commit -m "Phase 2: Implement Stewart evaporation model"
git push origin agent1-stratiform

# Post in Slack/Discussions:
# "Agent 1: Pushed updates to agent1-stratiform. Modified water_tracers.F90 (+150 lines). No interface changes."
```

### Weekly: Integration Test (Friday)
**Integration Manager** (rotating role):

```bash
# Create integration branch
git checkout -b integration-week-[N]
git merge origin/agent1-stratiform
git merge origin/agent2-convection
git merge origin/agent3-collection  # if started

# Resolve conflicts (see below)

# Build and test
cd ~/cases/integration_test_week_N
./case.build
./case.submit --test XS_2x5_ndays

# If successful:
# - Post results to team
# - Tag integration-week-N
# - Agents continue on their branches

# If failures:
# - Identify which agent's code caused issue
# - That agent fixes on their branch
# - Retry next day
```

---

## ⚔️ Merge Conflict Resolution

### High-Risk Files (Require Coordination)

#### File: `micro_p3.F90`
**Agents**: 1 and 3 both modify this file

**Conflict Prevention**:
- **Agent 1**: Modifies functions in Phase 2-4,7
  - `evaporate_rain()` (line ~2800)
  - `ice_nucleation()` (line ~1200)
  - `ice_deposition_sublimation()` (line ~2200)
  - `rain_sedimentation()` (line ~3200)

- **Agent 3**: Modifies functions in Phase 5-6
  - `cldliq_immersion_freezing()` (line ~1800)
  - `rain_immersion_freezing()` (line ~1900)
  - `ice_melting()` (line ~2500)
  - `ice_cldliq_collection()` (line ~2000)

**Coordination Rule**: 
- Post in daily standup which functions you're modifying
- If two agents need same function, coordinate who goes first

**Conflict Resolution**:
```bash
# Agent 3 merges Agent 1's work regularly
git checkout agent3-collection
git merge origin/agent1-stratiform

# If conflict in micro_p3.F90:
git diff  # Check conflicting lines
# Usually: Both agents added different functions
# Resolution: Keep both additions

# Mark as resolved
git add micro_p3.F90
git commit -m "Merge agent1: Keep both function modifications"
```

---

#### File: `water_tracers_p3.F90`
**Agents**: 1 creates, 3 extends

**Conflict Prevention**:
- **Agent 1** (Phase 2): Creates file with basic structure
  ```fortran
  subroutine wtrc_p3_inter_phase2(...)
    ! Only handles Phases 2-4,7 processes
  end subroutine
  ```

- **Agent 3** (Phase 5-6): Extends with more processes
  ```fortran
  subroutine wtrc_p3_inter_full(...)
    ! Adds Phases 5-6 processes
    call wtrc_p3_inter_phase2(...)  ! Call Agent 1's work
    ! Add Phase 5-6 rates
  end subroutine
  ```

**Coordination Rule**: 
- Agent 1 defines clear interfaces in Phase 2
- Agent 3 extends without modifying Agent 1's functions
- Use function overloading or separate entry points

---

#### File: `water_tracers.F90`
**Agents**: All agents add different functions

**Low Risk**: 
- Each agent adds new functions to end of file
- No modification of existing functions

**Conflict Resolution**: 
- Usually auto-mergeable
- Manual: Keep all additions from all agents

---

### Conflict Resolution Workflow

**Step 1: Identify Conflict Type**
```bash
git merge origin/agent2-convection
# CONFLICT in micro_p3.F90
```

**Step 2: Examine Conflict**
```bash
git diff micro_p3.F90
# <<<<<<< HEAD (your changes)
#   call isotope_evap_rain(qr_evap_rate)  # Agent 1
# =======
#   call isotope_freeze_rain(qr_frz_rate)  # Agent 3
# >>>>>>> origin/agent2-convection
```

**Step 3: Decide Resolution**
- **Both different functions**: Keep both (most common)
- **Same function modified**: Coordinate with other agent
- **Interface change**: Highest priority - escalate immediately

**Step 4: Resolve and Test**
```bash
# Edit file to keep both changes
vim micro_p3.F90

# Mark resolved
git add micro_p3.F90
git commit -m "Merge agent3: Keep both evap and freeze calls"

# TEST IMMEDIATELY
./case.build
make test_phase2_phase5
```

**Step 5: Notify Team**
```
Slack: "Resolved merge conflict in micro_p3.F90. Kept both evap_rain and freeze_rain modifications. Tests passing."
```

---

## 🧪 Testing Coordination

### Unit Tests (Independent)
Each agent runs their own unit tests without coordination:

```bash
# Agent 1
cd components/eam/test
make test_stewart_model
make test_rayleigh_distillation
make test_ice_fractionation

# Agent 2
make test_zm_convection
make test_shallow_convection

# Agent 3
make test_freezing_processes
make test_collection_processes
```

**Rule**: All unit tests must pass before pushing to your branch

---

### Integration Tests (Weekly)

**Friday Integration Test**:
```bash
# Integration Manager creates test case
cd $CIME_ROOT/scripts
./create_newcase --case ~/cases/integration_week_N \
                 --compset WCYCL1850 \
                 --res ne30pg2_r05_IcoswISC30E3r5 \
                 --machine pm-cpu

cd ~/cases/integration_week_N

# Enable isotopes
echo "wisotope = .true." >> user_nl_eam

# Short test
./case.setup
./case.build
./xmlchange STOP_N=2
./xmlchange STOP_OPTION=ndays
./case.submit --no-batch

# Check results
if [ $? -eq 0 ]; then
  echo "✅ Week N integration test PASSED"
  # Post results to team
else
  echo "❌ Week N integration test FAILED"
  # Debug and identify which agent's code caused failure
fi
```

**Results Posted to Team**:
```
Week 6 Integration Test Results
================================
✅ Compilation: SUCCESS
✅ Model run:   SUCCESS (completed 2 days)
✅ Mass conservation: 0.05% error (within tolerance)
✅ Phase 2 (Agent 1): Rain evaporation working
✅ Phase 8 (Agent 2): Convection ported (waiting for full integration)
⚠️  Warning: δD values seem high in stratosphere (expected for Phase 6 issue)

Next Steps:
- Agent 1: Continue with Phase 3
- Agent 2: Documentation and Phase 9 prep
- Integration test Week 7: Same protocol
```

---

## 📊 Progress Tracking

### Shared Project Board (GitHub Projects)

**Columns**:
- **Backlog**: Issues not yet started
- **In Progress - Agent 1**
- **In Progress - Agent 2**
- **In Progress - Agent 3**
- **Review**: Waiting for code review
- **Testing**: In integration testing
- **Done**: Completed and merged

**Issue Labels**:
- `agent-1`, `agent-2`, `agent-3`: Assigned agent
- `phase-2`, `phase-3`, etc.: Phase number
- `blocked`: Waiting for another agent or decision
- `merge-conflict`: Needs conflict resolution
- `high-priority`: Critical path item
- `question`: Needs science/design decision

---

### Weekly Progress Report Template

**Posted every Friday after integration test**:

```markdown
# Week N Progress Report (Date Range)

## Agent 1 (Stratiform)
**Phase**: [N]
**Completed This Week**:
- [ ] Task 1
- [ ] Task 2

**Next Week Goals**:
- [ ] Task 3
- [ ] Task 4

**Blockers**: None / [Description]

## Agent 2 (Convection)
**Phase**: [N]
**Status**: [In Progress / Waiting for Integration]
**Completed This Week**:
- [ ] Task 1

**Blockers**: None

## Agent 3 (Collection)
**Phase**: [N or Not Started]
**Status**: [Not Started / In Progress]
**Start Date**: [TBD or Date]

## Integration Status
✅ **Weekly Integration Test**: PASSED / FAILED
⚠️  **Warnings**: [Any issues]
📊 **Test Results**: [Link to output]

## Overall Timeline
On Track / 1 Week Behind / 2 Weeks Ahead

## Decisions Needed
1. [Decision 1] - See ISOTOPE_PHYSICS.md Question #N
2. [Decision 2]
```

---

## 🚨 Escalation Protocol

### When to Escalate

**Immediate Escalation** (Post in #urgent Slack channel):
- Build broken for >2 hours
- Integration test fails 2 weeks in a row
- Blocking dependency discovered between agents
- Science question needs expert input (can't proceed)

**Next-Day Escalation** (Raise in Friday meeting):
- Falling behind schedule by >1 week
- Repeated merge conflicts in same file
- Design decision needed
- Need help from another agent

**Eventually Escalate** (Document in GitHub issue):
- Minor questions about implementation
- Performance optimization ideas
- Documentation gaps

---

### Escalation Response Time

| Severity | Response Time | Responsible |
|----------|---------------|-------------|
| 🔴 Build Broken | 2 hours | Agent who broke it |
| 🟠 Integration Failure | 1 day | Integration Manager |
| 🟡 Behind Schedule | 1 week | Team Lead |
| 🟢 Minor Question | Next sync | Any team member |

---

## 📝 Communication Channels

### Primary: GitHub
- **Issues**: Task tracking and bug reports
- **Pull Requests**: Code reviews
- **Discussions**: Questions and design decisions
- **Projects**: Progress tracking board

### Secondary: Slack (or equivalent)
- **#p3-isotopes-general**: General discussion
- **#p3-isotopes-daily**: Daily standups
- **#p3-isotopes-urgent**: Build breaks, blockers
- **#p3-isotopes-science**: Science questions (link to ISOTOPE_PHYSICS.md)

### Tertiary: Email
- Weekly progress reports (archive)
- Formal milestone approvals

---

## ✅ Definition of Done

### Phase Complete Checklist

Before marking a phase as "done" and moving to next:

**Code**:
- [ ] All functions implemented
- [ ] Code compiles without warnings
- [ ] Unit tests written and passing
- [ ] Code reviewed by at least one other person
- [ ] No commented-out code or TODO markers

**Documentation**:
- [ ] Functions documented (purpose, inputs, outputs)
- [ ] README updated if needed
- [ ] Known issues documented

**Testing**:
- [ ] Unit tests pass
- [ ] Integration test passes
- [ ] Physics validation (if applicable to phase)

**Coordination**:
- [ ] Other agents notified of completion
- [ ] Merge conflicts resolved
- [ ] Branch pushed to remote

---

## 🎓 Lessons Learned (Update After Each Phase)

### Phase 1 Lessons
*To be filled after Phase 1 completion*

### Phase 2 Lessons
*To be filled after Phase 2 completion*

### Phase 3 Lessons
*To be filled after Phase 3 completion*

---

## 🏁 Final Integration (Phase 9)

### Integration Week Procedure

**Monday**: 
- All agents create final PRs to integration branch
- Code freeze on all agent branches

**Tuesday**: 
- Integration Manager merges all PRs
- Resolve all conflicts
- Full build and unit test suite

**Wednesday**: 
- Short global test (5 days)
- Check for runtime errors

**Thursday**: 
- Start multi-month validation simulation
- Begin conservation checks

**Friday**: 
- Review preliminary results
- Plan Phase 9 conservation and diagnostics work

---

## 📞 Contact Information

| Role | Agent | Contact | Availability |
|------|-------|---------|--------------|
| Agent 1 (Stratiform) | TBD | TBD | TBD |
| Agent 2 (Convection) | TBD | TBD | TBD |
| Agent 3 (Collection) | TBD | TBD | TBD |
| Integration Manager | Rotating | TBD | TBD |
| Science Lead | TBD | TBD | TBD |
| Project Manager | TBD | TBD | TBD |

---

## 🔄 This Document

**Version**: 1.0  
**Last Updated**: 2026-04-23  
**Next Review**: After Phase 1 completion  

**Update Frequency**: 
- After each phase completion
- When processes change
- When new team members join

**Owner**: Project Manager

---

**Ready to coordinate? Let's start Phase 1!** 🚀
