# EAMv2 Water Isotope Testing - Future Validation Phases

**Purpose:** Document future validation and testing phases for the EAMv2 water isotope implementation.

**Status:** Phase 1 (Basic Smoke Tests) completed. Phases 2-5 planned for future implementation.

---

## Phase 1: Basic Smoke Tests (✓ COMPLETED)

**Status:** Implemented and operational

**What it does:**
- Builds EAM with 4 different isotope tracer packages
- Runs 5-day simulations with Fi2000climo compset at ne30pg2 resolution
- Validates model completion and presence of isotope tracer fields

**Success Criteria:**
- ✓ Model builds successfully with isotope tracers
- ✓ Simulations complete without crashes
- ✓ Isotope tracer variables present in output
- ✓ Tracer fields contain non-zero values

**Tools:**
- `test_isotopes.sh` - Main test suite
- `validate_isotope_output.py` - Basic output validation
- `run_isotope_tests_docker.sh` - Docker wrapper

---

## Phase 2: Physical Realism Validation (FUTURE)

**Priority:** HIGH  
**Estimated Effort:** 2-3 weeks

### Objectives

Verify that isotope tracers produce physically realistic spatial and temporal patterns consistent with known isotope behavior in the atmosphere.

### Validation Tests

#### 2.1 Isotopic Ratio Calculations

**Implementation:**
- Calculate δD and δ18O from tracer ratios
- Formula: δ = [(R_sample / R_standard) - 1] × 1000 (per mil)
- Standard ratios (SMOW): HDO/H2O = 155.76×10⁻⁶, H218O/H2O = 2005.20×10⁻⁶

**Expected Results:**
- Global mean δD: approximately -25‰ to +5‰
- Global mean δ18O: approximately -5‰ to 0‰
- Ratios should be finite and within physical bounds

**Success Criteria:**
- No NaN or infinite values in calculated deltas
- Ratios fall within observed atmospheric range
- Spatial patterns show expected gradients

#### 2.2 Meteoric Water Line (MWL)

**Implementation:**
- Plot δD vs δ18O from precipitation output
- Fit linear relationship: δD = slope × δ18O + intercept
- Compare to Global Meteoric Water Line (GMWL): δD = 8 × δ18O + 10

**Expected Results:**
- Slope: 7-9 (global mean ~8)
- Intercept: 0-15 (global mean ~10)
- R² > 0.95 for global fit

**Success Criteria:**
- Simulated MWL slope within ±1 of observations
- Strong correlation between δD and δ18O
- Regional variations consistent with climatology

**Tools Needed:**
- Python script to calculate deltas and fit MWL
- Comparison with GNIP (Global Network of Isotopes in Precipitation) data
- Visualization of global MWL and regional variations

#### 2.3 Polar Depletion Gradient

**Implementation:**
- Calculate mean δD and δ18O as function of latitude
- Compare tropical (±30°) vs polar (>60°) values

**Expected Results:**
- Progressive depletion from equator to poles
- Tropical δ18O: -5‰ to +2‰
- Polar δ18O: -30‰ to -15‰
- Temperature effect: ~0.5-0.6‰/°C for δ18O

**Success Criteria:**
- Clear latitudinal gradient in isotope values
- Polar regions show 20-30‰ depletion in δ18O relative to tropics
- Gradient consistent with observations

**Tools Needed:**
- Zonal mean analysis scripts
- Comparison with observed isotope climatologies
- Temperature-isotope correlation analysis

#### 2.4 Deuterium Excess (d-excess)

**Implementation:**
- Calculate d-excess: d = δD - 8 × δ18O
- Analyze spatial and seasonal patterns

**Expected Results:**
- Global mean: ~10‰
- Higher values in ocean evaporation regions (10-20‰)
- Lower values in continental precipitation (5-10‰)
- Positive values throughout (kinetic effects during evaporation)

**Success Criteria:**
- Global mean d-excess: 8-12‰
- Spatial patterns consistent with observations
- No widespread negative d-excess values (unphysical)

**Tools Needed:**
- d-excess calculation script
- Global and regional analysis
- Seasonal cycle validation

#### 2.5 Temporal Patterns

**Implementation:**
- Run longer simulation (1 year or more)
- Analyze seasonal cycles of δD and δ18O
- Compare to observed seasonal variations

**Expected Results:**
- Seasonal amplitude: 5-15‰ for δ18O (varies by region)
- Summer enrichment, winter depletion at mid-latitudes
- Tropical regions show smaller seasonal variations

**Success Criteria:**
- Realistic seasonal amplitude
- Correct phase of seasonal cycle
- Regional patterns match observations

### Deliverables

1. **Physical Validation Script** (`validate_physical_realism.py`)
   - Delta calculations
   - MWL fitting and plotting
   - Latitudinal gradient analysis
   - d-excess calculations
   - Statistical comparisons

2. **Reference Data**
   - GNIP station data for comparison
   - Climatological isotope fields
   - Published isotope patterns from literature

3. **Validation Report Template**
   - Automated generation of validation figures
   - Statistical metrics
   - Pass/fail criteria for each test

4. **Documentation**
   - Update README with Phase 2 instructions
   - Document expected ranges and tolerances
   - Reference key isotope studies

### Success Metrics

- 90% of grid cells have δD and δ18O within observed ranges
- MWL slope within 7.5-8.5
- Latitudinal gradient present and realistic
- d-excess predominantly positive (>95% of grid cells)
- Seasonal cycle amplitude within factor of 2 of observations

---

## Phase 3: Mass Conservation Analysis (FUTURE)

**Priority:** HIGH  
**Estimated Effort:** 2-3 weeks

### Objectives

Verify that water mass is conserved when tracking isotope tracers and that total water equals the sum of individual isotope species.

### Validation Tests

#### 3.1 Total Water Consistency

**Implementation:**
- Sum all isotope tracer mixing ratios
- Compare to standard total water (Q, CLDLIQ, CLDICE)
- Calculate relative errors: |(Q_total - Q_isotopes)/Q_total|

**Expected Results:**
- Relative error < 1e-10 (machine precision)
- No systematic drift over time
- Conservation maintained at each timestep

**Success Criteria:**
- Mean relative error < 1e-8
- Maximum relative error < 1e-6
- No increasing trend in error

#### 3.2 Global Water Budget

**Implementation:**
- Calculate global integrals of water tracers
- Track evaporation, precipitation, and storage
- Verify closure: dW/dt = E - P (where W=storage, E=evap, P=precip)

**Expected Results:**
- Budget closure within 1%
- No unphysical sources or sinks
- Isotope tracers follow same budget as total water

**Success Criteria:**
- Global water budget closes within 0.1%
- No drift in global mean water content
- Isotope budget consistent with total water budget

#### 3.3 Physical Constraints

**Implementation:**
- Check for negative mixing ratios
- Verify ratios R = Q_isotope/Q_total stay within bounds
- Check for NaN, Inf, or other unphysical values

**Expected Results:**
- Zero negative values
- Ratios: 0 < R_HDO/R_H2O < 1.2 (with margin for variations)
- Ratios: 0 < R_H218O/R_H2O < 1.2
- No NaN or Inf values

**Success Criteria:**
- 100% of values are physical (no negatives, no NaN/Inf)
- Isotope ratios stay within theoretical bounds
- No conservation violations flagged by model

#### 3.4 Component Coupling Conservation

**Implementation:**
- Track isotope fluxes across component boundaries
- Verify atmosphere-ocean, atmosphere-land isotope exchanges
- Check that coupler conserves isotope mass

**Expected Results:**
- Fluxes into atmosphere = sum of component fluxes
- No mass loss or gain at component boundaries
- Coupler isotope fields (FLDS_WISO) consistent

**Success Criteria:**
- Component flux balance within 0.1%
- No systematic errors at boundaries
- Land and ocean isotope ratios influence atmosphere correctly

### Deliverables

1. **Conservation Analysis Script** (`validate_mass_conservation.py`)
   - Total water vs isotope sum checks
   - Global budget calculations
   - Physical constraint violations
   - Time series analysis of conservation

2. **Monitoring Tools**
   - Real-time conservation checking during runs
   - Automated alerts for violations
   - Diagnostic output fields

3. **Debugging Tools**
   - Identification of conservation violation sources
   - Component-by-component analysis
   - Process-level tracking

### Success Metrics

- Total water matches isotope sum to machine precision
- Global water budget closes within 0.1%
- Zero negative mixing ratios
- Zero NaN or Inf values
- Component coupling conserves mass within 0.1%

---

## Phase 4: Performance and Scalability (FUTURE)

**Priority:** MEDIUM  
**Estimated Effort:** 1-2 weeks

### Objectives

Quantify the computational cost of isotope tracking and ensure the implementation scales efficiently.

### Validation Tests

#### 4.1 Runtime Overhead

**Implementation:**
- Run identical simulations with and without isotopes
- Measure wall-clock time for each tracer package
- Calculate overhead percentage

**Expected Results:**
- Overhead roughly proportional to number of tracers
- h2o_hdo (14 tracers): ~15-25% overhead
- all_stable_wiso (35 tracers): ~30-50% overhead

**Success Criteria:**
- Overhead less than 2× number of tracers
- No nonlinear scaling with tracer count
- Overhead consistent across different cases

#### 4.2 Memory Usage

**Implementation:**
- Profile memory usage with different tracer packages
- Measure peak memory and average memory
- Compare to non-isotope runs

**Expected Results:**
- Memory increase proportional to tracers
- ~2-5% increase per tracer
- No memory leaks during long runs

**Success Criteria:**
- Memory overhead < 3× tracer percentage
- Stable memory usage (no leaks)
- Sufficient memory available for production runs

#### 4.3 Scalability Testing

**Implementation:**
- Run isotope cases with different processor counts
- Measure strong scaling efficiency
- Test weak scaling with different resolutions

**Expected Results:**
- Isotope tracers scale similarly to base model
- Strong scaling efficiency > 80% up to typical core counts
- No load imbalance introduced by isotopes

**Success Criteria:**
- Scaling similar to non-isotope runs (within 5%)
- Parallel efficiency maintained
- No bottlenecks introduced

#### 4.4 I/O Performance

**Implementation:**
- Measure output file sizes
- Time history file writing
- Test with different I/O configurations (PIO settings)

**Expected Results:**
- Output file size scales with number of tracers
- h2o_hdo: ~2× larger files than standard
- all_stable_wiso: ~5× larger files
- I/O time proportional to data volume

**Success Criteria:**
- I/O time increase < 2× data size increase
- No I/O bottlenecks
- Output frequency not limited by isotopes

### Deliverables

1. **Performance Benchmarking Suite**
   - Automated timing scripts
   - Memory profiling tools
   - Scalability test configurations

2. **Optimization Recommendations**
   - Identify performance bottlenecks
   - Suggest code improvements if needed
   - Document best practices for isotope runs

3. **Performance Documentation**
   - Expected overhead tables
   - Scaling plots
   - Resource requirements guide

### Success Metrics

- Runtime overhead < 50% for all_stable_wiso
- Linear scaling of overhead with tracer count
- Memory usage predictable and reasonable
- Parallel scaling not degraded by isotopes

---

## Phase 5: Scientific Validation (FUTURE)

**Priority:** LOW-MEDIUM  
**Estimated Effort:** 1-2 months

### Objectives

Validate isotope simulations against observational datasets and published studies to ensure scientific accuracy.

### Validation Tests

#### 5.1 GNIP Station Comparison

**Implementation:**
- Extract precipitation isotopes at GNIP station locations
- Compare monthly means to observations
- Calculate bias, RMSE, correlation for each station

**Expected Results:**
- Global mean bias < 2‰ for δ18O
- RMSE < 3-4‰ for δ18O
- Correlation > 0.7 for most stations

**Success Criteria:**
- Bias within range of other isotope-enabled models
- Correlations show skill in reproducing observations
- Regional patterns captured correctly

#### 5.2 Published Model Intercomparison

**Implementation:**
- Run standard intercomparison cases (e.g., SWING2)
- Compare to other isotope-enabled models (ECHAM, CAM5)
- Use same forcing and boundary conditions

**Expected Results:**
- Results within envelope of model spread
- Similar patterns to published E3SM isotope studies
- No obvious biases relative to other models

**Success Criteria:**
- Within model ensemble range
- Reasonable agreement with multi-model mean
- No outlier behavior

#### 5.3 Long-term Climate Simulation

**Implementation:**
- Run multi-year (5-10 year) simulation
- Analyze climatological means
- Check for drift or unphysical trends

**Expected Results:**
- Stable climatology after spinup
- No drift in global mean isotopes
- Realistic interannual variability

**Success Criteria:**
- Climatological mean stable after 2-3 years
- No significant trends in global means
- Variability consistent with observations

#### 5.4 Process-level Validation

**Implementation:**
- Validate specific processes (evaporation, condensation, etc.)
- Check fractionation factors during phase changes
- Compare to theoretical expectations

**Expected Results:**
- Equilibrium fractionation factors match literature values
- Kinetic fractionation active in appropriate conditions
- Process-level isotope behavior consistent with theory

**Success Criteria:**
- Fractionation factors within 5% of expected values
- Physical processes correctly implemented
- Process contributions make sense

### Deliverables

1. **Observational Comparison Tools**
   - GNIP data extraction and comparison scripts
   - Station-by-station validation reports
   - Regional synthesis

2. **Model Intercomparison Analysis**
   - Comparison with published studies
   - Quantification of model spread
   - Performance metrics relative to other models

3. **Long-term Simulation Protocol**
   - Recommended spinup procedures
   - Drift detection and analysis
   - Production run configuration guidance

4. **Scientific Validation Report**
   - Comprehensive comparison with observations
   - Model performance assessment
   - Known biases and limitations
   - Recommendations for scientific applications

### Success Metrics

- Agreement with GNIP observations comparable to other models
- Results consistent with published E3SM isotope studies
- Long-term simulations stable and physically reasonable
- Process-level behavior matches theory
- Ready for scientific applications (paleoclimate, water cycle studies)

---

## Implementation Priority

### Immediate (Next 1-3 months)
1. **Phase 2: Physical Realism** - Essential for scientific credibility
2. **Phase 3: Mass Conservation** - Critical for model integrity

### Near-term (3-6 months)
3. **Phase 4: Performance** - Needed for production runs
4. **Phase 5: GNIP Comparison** - Begin scientific validation

### Long-term (6-12 months)
5. **Phase 5: Full Scientific Validation** - Complete validation for publication
6. **Additional Tests** - Domain-specific validation as needed

---

## Resource Requirements

### Computing
- **Phase 2-3:** Minimal (use existing test output)
- **Phase 4:** Moderate (multiple benchmark runs)
- **Phase 5:** Significant (multi-year simulations)

### Data
- **GNIP database** - Free, ~100 MB
- **Satellite observations** - Varies, ~1-10 GB
- **Model intercomparison data** - Request from colleagues

### Personnel
- **Phase 2-3:** 1 person, 2-3 weeks each
- **Phase 4:** 1 person, 1-2 weeks
- **Phase 5:** 1-2 people, 1-2 months

### Software
- Python analysis tools (netCDF4, xarray, matplotlib)
- Statistical packages (scipy, sklearn)
- Visualization tools (cartopy, basemap)

---

## Documentation Updates Needed

As phases are implemented:

1. Update `README.md` with new validation procedures
2. Create phase-specific documentation
3. Document known issues and limitations
4. Add example usage for each validation tool
5. Maintain changelog of test suite improvements

---

## Contact & References

**Primary Contact:** rfiorella@lanl.gov

**Key References:**
- Noone, D. (2003) - Original isotope tracer implementation
- Nusbaumer, J. (2011) - CAM5 isotope adaptation
- GNIP Database: https://www.iaea.org/services/networks/gnip
- SWING2 Intercomparison: https://www.swing2.org/

**Last Updated:** 2026-04-07  
**Document Version:** 1.0
