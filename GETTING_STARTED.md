# Getting Started with P3 Water Isotope Implementation

**Quick-start guide for beginning Phase 1 implementation**

---

## 📋 Prerequisites Checklist

Before starting Phase 1, ensure you have:

- [ ] Access to both codebases:
  - [ ] EAMv2-wiso: `/Users/rfiorella/code/E3SM/EAMv2-wiso` (source)
  - [ ] EAMv3-wiso: `/Users/rfiorella/code/E3SM/EAMv3-wiso` (target)
- [ ] Read the documentation:
  - [ ] `P3_ISOTOPE_README.md` - Project overview
  - [ ] `P3_ISOTOPE_IMPLEMENTATION_PLAN.md` - Detailed plan
  - [ ] `CLAUDE.md` - Repository guide
- [ ] Development environment ready:
  - [ ] HPC account with E3SM access
  - [ ] Git configured
  - [ ] Fortran compiler working
  - [ ] CIME build system functional
- [ ] Team alignment:
  - [ ] Phase 0 review complete
  - [ ] Design decisions finalized
  - [ ] Personnel assigned

---

## 🚀 Phase 1 Quick Start (Next 3 Weeks)

### Week 1: Module Porting

**Monday - Wednesday** (Days 1-3): Issue #1.1

```bash
# Set up paths
export EAMv2_ROOT=/Users/rfiorella/code/E3SM/EAMv2-wiso
export EAMv3_ROOT=/Users/rfiorella/code/E3SM/EAMv3-wiso

# Create feature branch
cd $EAMv3_ROOT
git checkout -b feature/p3-isotopes-phase1

# Copy water_isotopes.F90
cp $EAMv2_ROOT/share/util/water_isotopes.F90 \
   $EAMv3_ROOT/share/util/

# Test compilation
cd $EAMv3_ROOT/share/util
# Add to build system (see CMakeLists.txt in that directory)

# Create unit test
cd $EAMv3_ROOT/components/eam/test
# Create test_water_isotopes.F90
```

**Thursday** (Day 4): Issue #1.2

```bash
# Copy water_types.F90
cp $EAMv2_ROOT/share/util/water_types.F90 \
   $EAMv3_ROOT/share/util/

# Test compilation
```

**Friday** (Day 5): Issue #1.3

```bash
# Copy water_tracer_vars.F90
cp $EAMv2_ROOT/components/eam/src/physics/cam/water_tracer_vars.F90 \
   $EAMv3_ROOT/components/eam/src/physics/cam/

# Update build system
```

---

### Week 2: Core Framework

**Monday - Friday** (Days 8-12): Issue #1.4

```bash
# Create skeleton water_tracers.F90
cd $EAMv3_ROOT/components/eam/src/physics/cam

# Port functions (in order):
# 1. wtrc_readnl() - Read namelist
# 2. wtrc_init() - Initialize module
# 3. wtrc_register() - Register constituents
# 4. wtrc_init_rates() - Initialize rate arrays
# 5. wtrc_add_rates() - Add process rates
# 6. wtrc_get_alpha() - Fractionation wrapper
# 7. wtrc_ratio() - Calculate isotope ratios

# IMPORTANT: Only port Phase 1 subset!
# Do NOT port:
#   - wtrc_apply_rates() (Phase 2)
#   - stewart_isoevap() (Phase 2)
#   - wtrc_sediment() (Phase 7)
#   - Convection functions (Phase 8)
```

---

### Week 3: Integration

**Monday - Wednesday** (Days 15-17): Issue #1.5

```fortran
! File: components/eam/src/physics/cam/physics_types.F90

! Add to physics_state type:
type physics_state
   ! ... existing fields ...
   
   ! Water isotope tracers (if wisotope=.true.)
   real(r8), pointer :: qtrc(:,:,:)    ! (pcols, pver, num_tracers)
   ! qtrc indices will be set by wtrc_register()
end type physics_state

! Add to physics_ptend type:
type physics_ptend
   ! ... existing fields ...
   
   ! Water isotope tendencies (if wisotope=.true.)
   real(r8), pointer :: qtrc_tend(:,:,:)
end type physics_ptend
```

**Thursday** (Day 18): Issue #1.6

```fortran
! File: components/eam/src/physics/cam/physpkg.F90

! In phys_init() subroutine, add:
use water_tracers, only: wtrc_init, wtrc_register

! After constituent registration:
if (wisotope) then
   call wtrc_register()
   call wtrc_init()
endif
```

**Friday** (Day 19): Issue #1.7

```xml
<!-- File: components/eam/bld/namelist_files/namelist_definition.xml -->

<!-- Add water_tracer_nl namelist group -->
<namelist_definition>
  <entry id="wisotope" type="logical" 
         category="water_tracers" group="water_tracer_nl">
    Enable water isotope tracers (.false. = disabled, .true. = enabled)
    Default: .false.
  </entry>
  
  <entry id="wtrc_alpha_kinetic" type="logical"
         category="water_tracers" group="water_tracer_nl">
    Enable kinetic fractionation effects
    Default: .true.
  </entry>
</namelist_definition>
```

---

### Week 3 End: Phase 1 Testing

**Day 20**: Full Phase 1 test

```bash
# Create test case
cd $CIME_ROOT/scripts
./create_newcase --case ~/cases/test_phase1_isotopes_off \
                 --compset WCYCL1850 \
                 --res ne30pg2_r05_IcoswISC30E3r5 \
                 --machine pm-cpu

cd ~/cases/test_phase1_isotopes_off
./case.setup
./case.build    # Should succeed

# Test with isotopes enabled
cd $CIME_ROOT/scripts
./create_newcase --case ~/cases/test_phase1_isotopes_on \
                 --compset WCYCL1850 \
                 --res ne30pg2_r05_IcoswISC30E3r5 \
                 --machine pm-cpu

cd ~/cases/test_phase1_isotopes_on

# Add to user_nl_eam:
echo "wisotope = .true." >> user_nl_eam

./case.setup
./case.build    # Should succeed

# Run short test (2 days)
./xmlchange STOP_N=2
./xmlchange STOP_OPTION=ndays
./case.submit

# Check results
# - Model should complete successfully
# - Isotope tracers should be initialized
# - No impact on standard physics yet
```

---

## 📁 File Organization

```
EAMv3-wiso/
├── CLAUDE.md                           # Repository guide
├── P3_ISOTOPE_README.md                # Project overview
├── P3_ISOTOPE_IMPLEMENTATION_SPEC.md   # Technical specification
├── P3_ISOTOPE_FUNCTION_MAPPING.md      # MG2→P3 mapping
├── P3_ISOTOPE_IMPLEMENTATION_PLAN.md   # This detailed plan
├── P3_ISOTOPE_ISSUES.md                # 50 detailed issues
├── P3_ISOTOPE_DEPENDENCIES.md          # Dependency diagrams
├── GETTING_STARTED.md                  # Quick-start guide (this file)
│
├── share/util/
│   ├── water_isotopes.F90              # ← Port in Phase 1 (Day 1-3)
│   ├── water_types.F90                 # ← Port in Phase 1 (Day 4-5)
│   └── CMakeLists.txt                  # ← Update for new modules
│
└── components/eam/src/physics/
    ├── cam/
    │   ├── water_tracer_vars.F90       # ← Port in Phase 1 (Day 6-7)
    │   ├── water_tracers.F90           # ← Port Phase 1 subset (Day 8-14)
    │   ├── physics_types.F90           # ← Modify in Phase 1 (Day 15-17)
    │   ├── physpkg.F90                 # ← Modify in Phase 1 (Day 18-19)
    │   └── ...
    └── p3/eam/
        ├── micro_p3.F90                # ← Modify in Phase 2 (Day 21-23)
        ├── micro_p3_interface.F90      # ← Modify in Phase 2 (Day 41-43)
        └── water_tracers_p3.F90        # ← Create in Phase 2 (Day 24-27)
```

---

## 🧪 Testing Workflow

### Unit Tests (Each Function)

```fortran
! Example: test_wiso_alpl.F90
program test_wiso_alpl
  use water_isotopes, only: wiso_alpl, isphdo, isph218o
  implicit none
  
  real(r8) :: T, alpha_hdo, alpha_h218o
  real(r8), parameter :: tol = 1.0e-3
  
  ! Test at 0°C (273.15 K)
  T = 273.15_r8
  
  alpha_hdo = wiso_alpl(isphdo, T)
  print *, "HDO α at 0°C:", alpha_hdo
  ! Expected: ~1.084
  if (abs(alpha_hdo - 1.084_r8) > tol) then
    print *, "FAIL: HDO fractionation incorrect"
    stop 1
  endif
  
  alpha_h218o = wiso_alpl(isph218o, T)
  print *, "H218O α at 0°C:", alpha_h218o
  ! Expected: ~1.0098
  if (abs(alpha_h218o - 1.0098_r8) > tol) then
    print *, "FAIL: H218O fractionation incorrect"
    stop 1
  endif
  
  print *, "PASS: Liquid-vapor fractionation correct"
end program
```

### Integration Test (Phase Complete)

```bash
#!/bin/bash
# test_phase1_complete.sh

# Test 1: Compile without isotopes
./case.build --clean-all
./case.build
if [ $? -ne 0 ]; then
  echo "FAIL: Build without isotopes failed"
  exit 1
fi

# Test 2: Compile with isotopes
echo "wisotope = .true." > user_nl_eam
./case.build --clean-all
./case.build
if [ $? -ne 0 ]; then
  echo "FAIL: Build with isotopes failed"
  exit 1
fi

# Test 3: Run short simulation
./xmlchange STOP_N=1
./xmlchange STOP_OPTION=ndays
./case.submit --no-batch
if [ $? -ne 0 ]; then
  echo "FAIL: Run with isotopes failed"
  exit 1
fi

echo "PASS: Phase 1 integration tests complete"
```

---

## 🐛 Common Issues & Solutions

### Issue 1: Module compilation errors

**Error**: `module water_isotopes not found`

**Solution**:
```bash
# Check CMakeLists.txt in share/util/
# Ensure water_isotopes.F90 is listed:
set(SOURCES
    water_isotopes.F90
    water_types.F90
    # ... other files
)

# Rebuild from scratch
./case.build --clean-all
./case.build
```

---

### Issue 2: Namelist not recognized

**Error**: `unknown namelist variable: wisotope`

**Solution**:
```bash
# Check namelist_definition.xml has water_tracer_nl group defined
# Rebuild configure:
cd components/eam/bld
./configure
```

---

### Issue 3: Constituent registration fails

**Error**: `ERROR: too many constituents`

**Solution**:
```fortran
! In physics_types.F90, increase pcnst parameter
! Each isotope adds 6 species × 3-5 phases = 18-30 constituents
! Default pcnst is usually ~40, may need to increase to 70
```

---

## 📚 Key References (Keep Handy)

### Scientific Papers
1. **Stewart (1975)**: "Stable isotope fractionation due to evaporation"
   - Core physics for rain evaporation (Phase 2)

2. **Horita & Wesolowski (1994)**: "Liquid-vapor fractionation of oxygen and hydrogen isotopes"
   - Equilibrium fractionation factors for liquid-vapor

3. **Merlivat & Nief (1967)**: "Fractionnement isotopique lors des changements d'état solide-vapeur"
   - Ice-vapor equilibrium fractionation

4. **Morrison & Milbrandt (2015)**: "Parameterization of cloud microphysics based on the prediction of bulk ice particle properties"
   - P3 microphysics scheme (J. Atmos. Sci.)

### Code References
- **EAMv2-wiso repository**: `/Users/rfiorella/code/E3SM/EAMv2-wiso`
  - Reference implementation with MG2 microphysics
  - Key file: `components/eam/src/physics/cam/water_tracers.F90` (7,813 lines)

- **E3SM Documentation**: https://docs.e3sm.org/E3SM/
- **CIME Documentation**: https://esmci.github.io/cime/

---

## ✅ Phase 1 Completion Checklist

Before moving to Phase 2, verify:

### Code
- [ ] All Phase 1 modules compile without warnings
- [ ] Unit tests pass for all fractionation functions
- [ ] Model builds successfully with `wisotope=.true.`
- [ ] Model builds successfully with `wisotope=.false.` (backward compatible)
- [ ] No memory leaks detected (valgrind if available)

### Testing
- [ ] Short test run (1-2 days) completes successfully
- [ ] Isotope tracers initialized correctly
- [ ] Initial isotope values reasonable (not NaN or Inf)
- [ ] History output includes isotope fields (if requested)
- [ ] Standard physics unchanged when isotopes disabled

### Documentation
- [ ] Code comments added to modified files
- [ ] Function documentation updated
- [ ] Build instructions documented
- [ ] Known issues documented

### Team
- [ ] Code review completed
- [ ] Phase 1 review meeting held
- [ ] Approval to proceed to Phase 2
- [ ] Phase 2 personnel assigned

---

## 🎯 Success Criteria Summary

**Phase 1 Success** = Code compiles, model runs, isotopes initialized (not yet active in physics)

**Phase 2 Success** = Rain evaporation depletes heavy isotopes correctly (δD becomes more negative)

**Phase 9 Success** = Precipitation isotopes match GNIP observations within ±20‰

**Final Success** = Publication-ready isotope-enabled E3SM v3 with P3 microphysics

---

## 💬 Getting Help

### Internal Resources
- **Documentation**: All `P3_ISOTOPE_*.md` files in repository
- **Reference code**: EAMv2-wiso implementation
- **Issues tracker**: See `P3_ISOTOPE_ISSUES.md` for detailed issue templates

### External Resources
- **E3SM Forums**: https://github.com/E3SM-Project/E3SM/discussions
- **CIME Support**: https://esmci.github.io/cime/
- **E3SM Documentation**: https://docs.e3sm.org/

### Questions to Address in Phase 0 Review
1. Is the phased approach acceptable?
2. Should Phase 8 (convection) be moved earlier?
3. What validation datasets are priorities?
4. What is the target completion date?
5. Who will be involved in coding, testing, validation?

---

## 🏁 Ready to Start?

If you've completed all items in the Prerequisites Checklist, you're ready to begin Phase 1!

**First task**: Issue #1.1 - Port water_isotopes.F90 (Days 1-3)

```bash
# Let's go!
cd /Users/rfiorella/code/E3SM/EAMv3-wiso
git checkout -b feature/p3-isotopes-phase1
git push -u origin feature/p3-isotopes-phase1

# Start Issue #1.1
cp /Users/rfiorella/code/E3SM/EAMv2-wiso/share/util/water_isotopes.F90 \
   share/util/
```

Good luck! 🚀

---

**Document Version**: 1.0  
**Created**: 2026-04-23  
**Next Review**: After Phase 1 completion (Week 3)
