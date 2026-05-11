# Water Tracers Implementation Status

**Date:** February 17, 2026  
**Author:** Rich Fiorella  
**Source:** water_tracers.F90 (7,813 lines Fortran) → C++/Kokkos

## Implementation Progress Summary

### ✅ COMPLETED FUNCTIONS (22 functions - ~40% by count)

#### Core Calculation Functions (Highest Usage)
1. ✅ **wtrc_ratio** (lines 6966-6995) - **CRITICAL - 80+ calls**
   - Calculate isotope ratios with numerical precision handling
   - Returns standard ratio when denominator < wtrc_qmin
   
2. ✅ **wtrc_get_alpha** (lines 7042-7140) - **30+ calls**
   - Get fractionation factor for phase transitions
   - Delegates to water_types::wtype_get_alpha
   
3. ✅ **wtrc_liqvap_equil** (lines 6289-6341) - **15+ calls**
   - Equilibrate liquid and vapor isotopes
   - Uses implicit solution to avoid instabilities
   
4. ✅ **wtrc_add_rate** (lines 901-943)
   - Add single process rate to tracking matrix
   - Handles forward and reverse rates
   
5. ✅ **wtrc_add_rates** (lines 947-991) - **15+ calls** - **JUST COMPLETED**
   - Add array of process rates
   - Handles special freezing rain to cloud ice case (MG2)

#### Fractionation & Physics
6. ✅ **wtrc_vap_distil** (lines 6390-6422)
   - Rayleigh distillation for vapor phase
   - Uses classic f^alpha formula

7. ✅ **wtrc_efac** (lines 6721-6748)
   - Equilibrium implicit factor calculation
   
8. ✅ **wtrc_dqequil** (lines 6344-6387)
   - Calculate equilibrium change for distillation

9. ✅ **wtrc_equil_time** (lines 6548-6718)
   - Calculate rain equilibration fraction
   - Based on Stewart (1975) raindrop exchange model

10. ✅ **wtrc_init_rates** (lines 995-1016)
    - Initialize process rate matrices to zero

11. ✅ **wtrc_get_rstd** (wrapper to water_isotopes)
    - Get standard isotopic ratios

12. ✅ **wtrc_get_rao2** (lines 7242-7267)
    - Get O2 fractionation for chemistry

#### Query Functions (All 10 implemented)
13-22. ✅ **wtrc_is_wtrc, wtrc_is_vap, wtrc_is_liq, wtrc_is_ice, wtrc_is_cvrain, wtrc_is_cvsnow, wtrc_is_strain, wtrc_is_stsnow, wtrc_is_tagged, wtrc_get_icnst**
    - Simple boolean type checking functions
    - All inline for performance

---

## ⚠️ HIGH PRIORITY REMAINING WORK (8 functions - ~2,750 lines)

### 1. ⚠️ wtrc_apply_rates_mg1 (Lines 1020-1809, ~789 lines) 
**Status:** Interface declared, implementation needed  
**Complexity:** VERY HIGH  
**Usage:** Called by MG1 microphysics scheme

**Key Components:**
- Iterative equilibration loops (wtrc_niter iterations)
- Pre-sedimentation process application
- Sedimentation handling
- Post-sedimentation process application  
- Precipitation phase changes (rain/snow freezing/melting)
- Rayleigh distillation for ice/snow formation
- Bergeron process handling
- Vapor-liquid equilibration
- Rain equilibration below cloud base
- Numerical precision corrections (two-pass error fixing)
- Mass conservation calculations

**Conversion Challenge:**
- Deeply nested loops (iter, k, i, iwset, isrctype, idsttype)
- State-dependent fractionation (temperature, humidity)
- Iterative convergence
- Multiple optional arguments
- Integration with physics buffer for precipitation

**Lines of Code:** ~789 lines

---

### 2. ⚠️ wtrc_apply_rates (Lines 2068-2697, ~629 lines)
**Status:** Interface declared, implementation needed  
**Complexity:** VERY HIGH  
**Usage:** Called by MG2 microphysics scheme

**Key Components:**
- Similar to wtrc_apply_rates_mg1 but for MG2
- Additional handling for stratiform rain/snow
- More sophisticated precipitation tracking
- Integration with wtrc_sediment (not wtrc_sediment_mg1)
- Phase change tracking (frzpre, mltpre)
- Rain->cloud ice conversion logic

**Differences from MG1:**
- Uses wtrc_sediment instead of wtrc_sediment_mg1
- More explicit precip phase change handling
- Different precipitation variable names (precr, preci outputs)
- Uses pre_rates_in instead of pre_rates

**Lines of Code:** ~629 lines

---

### 3. ⚠️ wtrc_sediment_mg1 (Lines 2850-3136, ~286 lines)
**Status:** Interface declared, implementation needed  
**Complexity:** HIGH  
**Usage:** Called by wtrc_apply_rates_mg1

**Key Components:**
- CFL-stable sub-stepping for sedimentation
- Fall velocity calculations (fc, fi)
- Cloud fraction handling (liqcldf, icecldf)
- Flux divergence calculations
- Instant evaporation/sublimation of flux entering clear-sky
- Isotopic equilibration during sedimentation
- Temperature updates from latent heating
- Precipitation accumulation

**Algorithm:**
1. Calculate sub-steps for CFL stability: `nstep = max(int(rgvm*dtime/pdel+1), nstep)`
2. Loop over sub-steps
3. Calculate sedimentation fluxes: `falout = fall_velocity * q`
4. Apply flux divergence with cloud fraction weighting
5. Equilibrate vapor/liquid for isotopes
6. Update temperature from phase changes
7. Accumulate surface precipitation

**Lines of Code:** ~286 lines

---

### 4. ⚠️ wtrc_sediment (Lines 3138-3623, ~485 lines)
**Status:** Interface declared, implementation needed  
**Complexity:** HIGH  
**Usage:** Called by wtrc_apply_rates (MG2)

**Key Components:**
- Similar to wtrc_sediment_mg1 but handles 4 hydrometeor types:
  - Cloud liquid (fc)
  - Cloud ice (fi)
  - Stratiform rain (fr)
  - Stratiform snow (fs)
- More complex CFL calculation across all 4 types
- Separate precipitation tracking for rain vs snow
- Rain equilibration below cloud base (Stewart isotope exchange)
- Different output structure (precr, preci arrays)

**Additional Complexity vs MG1:**
- 4 fall velocity arrays instead of 2
- Separate rain/snow mass tracking
- More sophisticated equilibration logic
- Stewart isotope evaporation model integration

**Lines of Code:** ~485 lines

---

### 5. ⚠️ wtrc_mg_inter (Lines 1812-2064, ~252 lines)
**Status:** Interface declared, implementation needed  
**Complexity:** MEDIUM-HIGH  
**Usage:** Interface between MG2 microphysics and wtrc_apply_rates

**Key Components:**
- Setup pre-sedimentation rates from MG2 process rates:
  - Deposition/sublimation (pcmei/ncmei)
  - Rain/snow evaporation (preo, prdso)
  - Liquid freezing (mnuccco, mnuccto, msacwio)
  - Autoconversion/accretion (prao, prco, psacwso, pracso)
  - Bergeron processes (bergo, bergso)
  - Ice processes (praio, prcio, mnuccro)
- Setup post-sedimentation rates:
  - Condensation/deposition (qcreso, qireso)
  - Freezing (homoo)
  - Melting (melto)
  - Precipitation phase changes (frzrpst, meltspst)
- Call wtrc_apply_rates with prepared rate matrices
- Return precipitation output (precr, preci)

**Algorithm:**
1. Initialize pre_rates and post_rates arrays
2. Split MG2 rates into positive/negative components
3. Map MG2 process names to wtrc rate matrix indices
4. Call wtrc_add_rates for each process
5. Call wtrc_apply_rates
6. Return precipitation

**Lines of Code:** ~252 lines

---

### 6. ⚠️ wtrc_q1q2_pjr (Lines 6752-7063, ~311 lines)
**Status:** Interface declared, implementation needed  
**Complexity:** VERY HIGH  
**Usage:** Called by ZM deep convection scheme

**Key Components:**
- Calculate water tracer tendencies for ZM convection
- Updraft calculations:
  - Interface values (qhat) using log interpolation
  - Updraft humidity with condensation
  - Isotopic fractionation in updraft
  - Cloud liquid production (ql)
  - Precipitation production (rprd)
- Downdraft calculations:
  - Downdraft humidity
  - Rain evaporation
  - Isotopic equilibration (optional, currently disabled)
- Detrainment tendencies (wtdlf)

**Algorithm Sections:**
1. Exit early if no convection (lengath==0)
2. Initialize variables for all wsets
3. Calculate interface values (qhat)
4. Calculate updraft with condensation and fractionation
5. Calculate precipitation production
6. Calculate downdraft with evaporation
7. Apply precipitation evaporation
8. Calculate final tendencies
9. Handle detrainment
10. Handle subsidence

**Special Features:**
- Handles compressed array indexing (ideep)
- Mass fixing for numerical errors (uqdiff, dqdiff)
- Optional rain equilibration (currently commented out)
- Ratio tracking throughout vertical column

**Lines of Code:** ~311 lines

---

### 7. ⚠️ wtrc_mass_fixer (Lines 4928-5195, ~267 lines)
**Status:** Interface declared, implementation needed  
**Complexity:** MEDIUM  
**Usage:** Called by physics driver

**Key Components:**
- Reset standard (H2O) water tracer to match bulk Q
- Adjust isotope tracers proportionally
- Handle all 7 water types
- Fix negative values
- Conservation enforcement

**Algorithm:**
1. Loop over all water types
2. For H2O tracer, set equal to bulk water
3. Calculate correction needed
4. Apply proportional correction to isotope tracers
5. Ensure non-negative values

**Lines of Code:** ~267 lines

---

### 8. ⚠️ wtrc_check_h2o (Lines 4372-4703, ~331 lines)
**Status:** Interface declared, implementation needed  
**Complexity:** MEDIUM  
**Usage:** ~10+ calls for mass conservation checking

**Key Components:**
- Check total H2O conservation
- Compare H2O tracer vs bulk water
- Report errors with detailed diagnostics
- Optional tendency checking
- Configurable warning vs abort

**Algorithm:**
1. Sum all H2O tracer components
2. Sum all bulk water components  
3. Compare totals
4. If mismatch exceeds threshold:
   - Print diagnostic information
   - Optionally abort or just warn
5. Check individual water types

**Lines of Code:** ~331 lines

---

## 📊 Implementation Statistics

### By Function Count
- **Completed:** 22 / 54 = **41%**
- **High Priority Remaining:** 8 functions
- **Medium/Low Priority Remaining:** 24 functions

### By Lines of Code
- **Completed:** ~1,600 lines implemented
- **High Priority Remaining:** ~2,750 lines
- **Total Source:** ~7,813 lines
- **Progress:** ~20% by line count

### By Complexity
- **Simple (Query functions):** 10/10 = 100% complete
- **Medium (Utilities):** 6/12 = 50% complete
- **Complex (Physics integration):** 6/32 = 19% complete

---

## 🔧 Implementation Strategy for Remaining Functions

### Phase 1: Core Microphysics Integration (Highest Priority)
**Estimated Effort:** 2-3 weeks

1. **wtrc_mg_inter** (252 lines)
   - Start here - it's the interface, so understanding it helps with others
   - Map MG2 processes to rate matrices
   
2. **wtrc_sediment_mg1** (286 lines)
   - Simpler than wtrc_sediment (only 2 hydrometeor types)
   - Core algorithm for CFL sub-stepping
   
3. **wtrc_apply_rates_mg1** (789 lines)
   - Most complex but most important
   - Integrates everything together

### Phase 2: MG2 Support (High Priority)
**Estimated Effort:** 2 weeks

4. **wtrc_sediment** (485 lines)
   - Extends wtrc_sediment_mg1 to 4 hydrometeor types
   
5. **wtrc_apply_rates** (629 lines)
   - Similar to wtrc_apply_rates_mg1 but for MG2

### Phase 3: Convection & Mass Fixing (High Priority)
**Estimated Effort:** 1-2 weeks

6. **wtrc_q1q2_pjr** (311 lines)
   - ZM convection scheme interface
   - Complex but well-documented
   
7. **wtrc_mass_fixer** (267 lines)
   - Enforces conservation
   
8. **wtrc_check_h2o** (331 lines)
   - Diagnostic checking

### Phase 4: Remaining Functions (Medium/Low Priority)
**Estimated Effort:** 1-2 weeks

- Initialization functions (wtrc_readnl, wtrc_init, wtrc_register)
- Precipitation functions
- Diagnostics and checking functions
- Chemistry functions
- Unused functions (marked as such)

---

## 🎯 Next Immediate Steps

1. **Implement wtrc_mg_inter** (~252 lines)
   - Provides template for how MG2 integrates
   - Relatively straightforward mapping of process rates
   
2. **Implement wtrc_sediment_mg1** (~286 lines)
   - Core sedimentation algorithm
   - Needed by wtrc_apply_rates_mg1
   
3. **Implement wtrc_apply_rates_mg1** (~789 lines)
   - Main integration routine
   - Most complex but unlocks full MG1 functionality

---

## 📝 Notes on Implementation

### Kokkos Patterns Needed
- **Nested parallel_for:** For iter loops over columns/levels
- **parallel_reduce:** For sum/max operations (CFL calculations)
- **MDRangePolicy:** For multi-dimensional loops
- **Views:** Replace Fortran arrays
- **TeamPolicy:** Potentially for column-level parallelism

### EAMxx Infrastructure Needed
- Physics state structures
- Physics tendency structures (ptend)
- Field manager for precipitation output
- Saturation vapor pressure functions (qsat_water, qsat_ice)
- Physics buffer replacement

### Testing Strategy
- Unit test each function with simple inputs
- Compare against Fortran for identical inputs
- Single-column tests with full microphysics
- Check mass conservation
- GPU compatibility verification

---

## 📚 References

See WATER_TRACERS_CONVERSION_NOTES.md for:
- Complete function catalog
- Scientific references
- Usage analysis
- Design patterns

---

**Document Status:** Current as of February 17, 2026  
**Next Update:** After completing Phase 1 (wtrc_mg_inter, wtrc_sediment_mg1, wtrc_apply_rates_mg1)
