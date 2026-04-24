# Isotope Physics Reference & Implementation Questions

**Purpose**: Scientific reference for developers + open questions requiring decisions  
**Audience**: Agent developers, science team, reviewers  
**Status**: Living document - update as questions are resolved

---

## Table of Contents

1. [Core Isotope Physics](#core-isotope-physics)
2. [Fractionation Mechanisms](#fractionation-mechanisms)
3. [Process-Specific Physics](#process-specific-physics)
4. [Implementation Questions](#implementation-questions)
5. [Parameter Values & Tuning](#parameter-values--tuning)
6. [Validation Targets](#validation-targets)

---

## Core Isotope Physics

### Isotope Ratio Notation

**Absolute Ratio**: R = N_heavy / N_light (number of molecules)

**Delta Notation** (per mil, ‰):
```
δ = (R_sample / R_standard - 1) × 1000
```

**Standard Values** (VSMOW - Vienna Standard Mean Ocean Water):
- R_standard(HDO) = 155.76 × 10⁻⁶
- R_standard(H₂¹⁸O) = 2005.20 × 10⁻⁶
- R_standard(H₂¹⁷O) = 379.9 × 10⁻⁶ (approximately)

**Example**: 
- If δD = -100‰, the sample is depleted in deuterium by 10% relative to VSMOW
- Ocean surface: δD ≈ 0‰, δ¹⁸O ≈ 0‰ (by definition)
- Rain in mid-latitudes: δD ≈ -50‰, δ¹⁸O ≈ -8‰
- Antarctic ice: δD ≈ -400‰, δ¹⁸O ≈ -55‰

---

## Fractionation Mechanisms

### 1. Equilibrium Fractionation

**Definition**: Fractionation at thermodynamic equilibrium between phases

**Liquid-Vapor** (Horita & Wesolowski 1994):
```fortran
! Temperature-dependent fractionation factor
α_eq = exp(A/T² + B/T + C)

! For HDO:
A_liq_vap = 1158.8    ! K²
B_liq_vap = -1620.1   ! K
C_liq_vap = 0.5223    ! dimensionless
```

**Ice-Vapor** (Merlivat & Nief 1967, Majoube 1971):
```fortran
! For HDO:
α_ice_vap(T) = exp(16288/T² - 0.0934)
! Note: Ice-vapor fractionation is STRONGER than liquid-vapor
! At -20°C: α_ice_vap ≈ 1.18 vs α_liq_vap ≈ 1.09
```

**Physical Meaning**: Heavy isotopes preferentially partition into condensed phase (liquid or ice) because of stronger bonding.

**❓ QUESTION 1: Temperature Range Validity**
- Horita & Wesolowski formulas valid for T > 273.15 K
- Ice formulas valid for T < 273.15 K
- What to do in supercooled liquid regime (250-273 K)?
- **Current approach in EAMv2-wiso**: Uses liquid formula if liquid present, ice formula if only ice
- **Decision needed**: Should we interpolate? Use liquid formula exclusively for supercooled liquid?

---

### 2. Kinetic Fractionation

**Definition**: Fractionation due to diffusion rate differences (molecular mass dependent)

**Diffusivity Ratios** (difrm):
```fortran
D_H2O / D_HDO   = 1.0251    ! HDO diffuses ~2.5% slower
D_H2O / D_H218O = 1.0285    ! H218O diffuses ~2.8% slower
D_H2O / D_H217O = 1.0145    ! H217O diffuses ~1.4% slower
```

**Kinetic Fractionation Factor**:
```fortran
α_kinetic = (D_H2O / D_iso)^n

! where n = exponent depending on turbulence
!   n = 0.5  : Partially turbulent (typical for atmospheric processes)
!   n = 1.0  : Fully molecular diffusion
!   n = 0.0  : Fully turbulent (no kinetic effect)
```

**❓ QUESTION 2: Kinetic Exponent (n)**
- EAMv2-wiso uses `dkfac = 0.58` (from Stewart 1975)
- Bony et al. (2008) uses n = 0.5
- Jouzel & Merlivat (1984) uses n = 0.5
- **Should we use 0.58 or 0.5?** Or make it a tunable parameter?
- **Does n vary by process** (evaporation vs ice deposition)?

---

### 3. Combined Fractionation

**General Rule**: 
```fortran
α_total = α_equilibrium × α_kinetic
```

**Exception**: Sublimation (ice → vapor)
- Kinetic effects are reversed (vapor is source, not sink)
- α_kinetic = 1 / α_kinetic_deposition (approximately)
- Or neglect kinetic effects for sublimation?

**❓ QUESTION 3: Kinetic Effects in Sublimation**
- EAMv2-wiso: Only uses equilibrium fractionation for sublimation
- Physical argument: Sublimating surface is saturated, no kinetic barrier
- **Should P3 implementation do the same?**
- Or include kinetic effects but reversed?

---

## Process-Specific Physics

### A. Rain Evaporation (Phase 2 - HIGHEST PRIORITY)

**Stewart (1975) Model**: Below-cloud rain evaporation with partial equilibration

**Physics**: As rain falls through unsaturated air, it evaporates. The evaporating water exchanges isotopes with the remaining droplet according to:

```fortran
! Stewart equation:
R_rain_final = ((R_rain_initial - γ·R_vapor)·f^β) + γ·R_vapor

where:
  f = fraction of rain remaining after evaporation
  γ = (α_e·h)/(1 - α_e·α_k·(1-h))  ! Effective fractionation
  β = (1 - α_e·α_k·(1-h))/(α_e·α_k·(1-h))  ! Kinetic exponent
  h = relative humidity (0-1)
  α_e = equilibrium fractionation factor
  α_k = kinetic fractionation factor
```

**Key Physics**:
1. **At high RH (h → 1)**: Approaches equilibrium (γ → α_e, β → 0)
2. **At low RH (h → 0)**: Strong kinetic effects dominate
3. **Evaporation enriches rain** in heavy isotope (makes δD less negative)
4. **Small drops equilibrate more** than large drops

**Partial Equilibration**: Only a fraction of the rain actually equilibrates:
```fortran
f_equil = 1 - exp(-t_fall / t_diffusion)

where:
  t_fall = Δz / v_terminal         ! Time to fall through layer
  t_diffusion = r² / (3·D·(1-h))  ! Diffusion time scale
  r = raindrop radius
  D = water vapor diffusivity
```

**❓ QUESTION 4: Raindrop Size Distribution**
- Stewart uses **Marshall-Palmer** distribution: N(D) = N₀·exp(-Λ·D)
- P3 uses predicted particle properties (more sophisticated)
- **Should we**:
  - A) Use Marshall-Palmer as Stewart did (simple, tested)?
  - B) Use P3's actual raindrop size distribution (consistent, complex)?
  - C) Use mean raindrop size from P3 (compromise)?
- **EAMv2-wiso choice**: Marshall-Palmer (option A)

**❓ QUESTION 5: Stewart Parameter φ (Tuning Factor)**
- Stewart equation includes tuning parameter φ ∈ [0,1]
- φ = 0: No equilibration (kinetic only)
- φ = 1: Full equilibration
- **EAMv2-wiso uses φ = 0.9** (from Bony et al. 2008)
- Bony calibrated this to match observations
- **Should we use same value or re-tune for P3?**
- Expect may need adjustment because P3 evaporation rates differ from MG2

**❓ QUESTION 6: Layer-by-Layer Application**
- Stewart model is for **single layer** evaporation
- Rain falls through **multiple layers** with different T, RH
- **Implementation options**:
  - A) Apply Stewart model independently in each layer (EAMv2-wiso approach)
  - B) Track cumulative fractionation as rain falls
  - C) Use implicit vertical integral
- **EAMv2-wiso choice**: Option A (simplest, most common)

---

### B. Condensation (Phase 3)

**Rayleigh Distillation**: Progressive condensation depletes vapor in heavy isotope

**Classic Rayleigh Equation**:
```fortran
R_vapor_final / R_vapor_initial = f^(α - 1)

where:
  f = fraction of vapor remaining after condensation
  α = equilibrium fractionation factor
```

**Physical Meaning**: 
- First condensate is enriched (δD less negative)
- Remaining vapor becomes progressively depleted (δD more negative)
- Cloud liquid has δD ≈ (α - 1)×1000 ‰ higher than source vapor

**❓ QUESTION 7: Kinetic Effects During Condensation**
- Equilibrium fractionation only? (Standard Rayleigh)
- Or include kinetic effects?
- **Argument for kinetic**: Supersaturation exists, diffusion barrier
- **Argument against**: Cloud droplets are small, equilibrate quickly
- **EAMv2-wiso choice**: Equilibrium only (standard Rayleigh)
- **Should P3 do same?**

**❓ QUESTION 8: P3 Supersaturation Adjustment**
- P3 removes supersaturation **implicitly** (no explicit condensation rate)
- How to extract condensation rate for isotopes?
- **Options**:
  - A) Calculate from Q change before/after saturation adjustment
  - B) Modify P3 to output explicit condensation rate
  - C) Estimate from supersaturation magnitude
- **Recommendation**: Option B (cleaner, explicit)
- **Issue**: Requires modifying P3 core (may affect bit-for-bit reproducibility)

---

### C. Ice Nucleation & Deposition (Phase 4)

**Unique Aspects**:
1. **Stronger fractionation** than liquid: α_ice_vap > α_liq_vap
2. **Supersaturation-dependent kinetics**: Ice can grow in sub-saturated air (w.r.t. liquid)
3. **Different ice formation modes**: Heterogeneous vs homogeneous

**Ice-Vapor Kinetic Fractionation** (Jouzel & Merlivat 1984):
```fortran
α_kinetic_ice = α_diff^n · (S_ice / S_sat)^m

where:
  α_diff = diffusivity ratio
  S_ice = ice supersaturation
  S_sat = saturation vapor pressure over ice
  n, m = empirical exponents
```

**❓ QUESTION 9: Supersaturation Dependence**
- How strongly does supersaturation affect kinetic fractionation?
- Jouzel & Merlivat: Strong effect, m ≈ 0.5-1.0
- Ciais & Jouzel (1994): Weaker effect at low temperatures
- **EAMv2-wiso**: Uses `wiso_akci()` with moderate S_i dependence
- **Should P3 use same formulation?**
- **Issue**: P3 tracks supersaturation explicitly, MG2 did not

**❓ QUESTION 10: Homogeneous vs Heterogeneous Freezing**
- **Heterogeneous nucleation** (slower, warmer): Kinetic effects important
- **Homogeneous nucleation** (fast, colder): No kinetic effects (instantaneous)
- P3 distinguishes these modes
- **Should isotopes treat differently?**
- **EAMv2-wiso**: No distinction (always applied fractionation)
- **Physics suggests**: Homogeneous should be α_kinetic = 1 (no kinetic effect)

---

### D. Bergeron Process (Phase 6)

**Definition**: Vapor transfer from liquid droplets to ice crystals in mixed-phase clouds

**Physics**: 
- Air is supersaturated w.r.t. ice but subsaturated w.r.t. liquid
- Ice grows by deposition while liquid evaporates
- Net result: Liquid → Vapor → Ice

**Two-Step Fractionation**:
```fortran
! Step 1: Liquid evaporation
α_liq_to_vap = 1 / α_eq_liq_vap  ! Inverse (evaporation)

! Step 2: Ice deposition
α_vap_to_ice = α_eq_ice_vap × α_kinetic_ice

! Total:
α_Bergeron = α_vap_to_ice / α_liq_to_vap
           = (α_eq_ice_vap × α_kinetic_ice) / α_eq_liq_vap
```

**❓ QUESTION 11: Bergeron Implementation**
- **MG2 approach**: Explicit rates `bergo` (liquid loss) and `bergso` (ice gain)
- **P3 approach**: Bergeron is implicit in supersaturation adjustment
- **How to implement isotopes?**
  - A) Extract effective Bergeron rate from Q changes
  - B) Treat as evaporation + deposition separately
  - C) Ignore (treat as simple phase transfer)
- **Recommendation**: Option B (most physics-correct)
- **Issue**: Requires identifying when Bergeron is happening

---

### E. Freezing & Melting (Phase 5)

**Key Principle**: Fractionation depends on process speed

**Slow Processes** (equilibrium):
- Immersion freezing at warm temperatures: Small fractionation
- Melting at 0°C: No fractionation (equilibrium)

**Fast Processes** (kinetic):
- Homogeneous freezing at T < -37°C: No fractionation (instantaneous)
- Rapid melting: No fractionation

**❓ QUESTION 12: Freezing Fractionation**
- **Physical expectation**: Ice formation from liquid has small fractionation
- **Literature**: δ(ice) ≈ δ(liquid) + 1 to 3‰ (small enrichment in ice)
- **EAMv2-wiso**: Applied small fractionation for immersion freezing
- **Should P3 do same?**
- **Alternative**: Treat as conservative (α = 1) since effect is small
- **Recommendation**: Include small effect for slow freezing, neglect for fast

**❓ QUESTION 13: Melting Fractionation**
- **Physical expectation**: At 0°C, liquid and ice are in equilibrium, α = 1
- **EAMv2-wiso**: No fractionation during melting
- **Should P3 do same?** (Almost certainly yes)

---

### F. Riming (Phase 6)

**Definition**: Supercooled liquid droplets freeze on contact with ice

**Physics**:
- **Instant freezing**: Droplet freezes immediately on contact
- **Kinetic fractionation**: During impact and freezing
- **Rime mass tracking**: P3 tracks rime mass separately (important!)

**❓ QUESTION 14: Riming Fractionation Factor**
- Literature: Variable, depends on impact velocity and droplet size
- Typical: α_rime ≈ 1.00 to 1.02 (very small fractionation)
- **EAMv2-wiso**: Applied equilibrium fractionation
- **Physical argument**: Fast process → kinetic effects
- **Should P3 use equilibrium, kinetic, or intermediate?**

**❓ QUESTION 15: Rime Mass in P3**
- P3 predicts rime mass fraction (unique feature!)
- **Should isotopes track rime mass separately from ice?**
- **Option A**: Yes - separate isotope tracer for rime mass
- **Option B**: No - combine rime into total ice isotope ratio
- **Recommendation**: Option B (simpler, no observations to validate separate tracking)

---

### G. Sedimentation (Phase 7)

**Critical Process**: Rain/snow falls through multiple layers, evaporating/sublimating

**Complexity**: 
- CFL condition requires sub-stepping
- Layer-by-layer fractionation accumulates
- Droplet/crystal size evolves during fall

**❓ QUESTION 16: CFL Sub-Stepping**
- P3 uses CFL-limited time stepping for sedimentation
- Each sub-step may have evaporation
- **Should isotopes apply fractionation at each sub-step?**
- **EAMv2-wiso**: Yes, sub-step integration
- **Alternative**: Assume steady state and integrate analytically?
- **Recommendation**: Follow EAMv2-wiso (more accurate)

**❓ QUESTION 17: Clear-Sky vs In-Cloud**
- Fractionation is **stronger in clear sky** (RH < 100%)
- Should track separately or use cloud fraction?
- **EAMv2-wiso**: Uses cloud fraction weighting
- **P3 implementation**: Same approach?

---

### H. Convection (Phase 8)

**Physical Differences from Stratiform**:
- **Fast updrafts**: Less time for equilibration
- **Entrainment/detrainment**: Mixing affects isotope ratios
- **Convective precipitation**: Distinct from stratiform

**❓ QUESTION 18: Convective Fractionation**
- **Updraft condensation**: Use Rayleigh model?
- **Downdraft evaporation**: Use Stewart model?
- **EAMv2-wiso**: Rayleigh in updraft, partial evap in downdraft
- **Should P3 convection use same?** (Almost certainly yes)

**❓ QUESTION 19: Convective vs Stratiform Tracking**
- **Why separate?**: Observations show different isotope signatures
- **Convective rain**: Often less depleted (rapid process)
- **Stratiform rain**: More depleted (slow equilibration)
- **Implementation**: Requires separate tracers `iwtcvrain` and `iwtstrain`
- **Validation**: Can compare convective/stratiform δD separately

---

## Implementation Questions Summary

### 🔴 HIGH PRIORITY (Must Decide for Phase 2)

**Q4**: Raindrop size distribution for Stewart model?
- **Options**: Marshall-Palmer / P3 DSD / Mean size
- **Recommendation**: Marshall-Palmer (tested, simple)

**Q5**: Stewart parameter φ value?
- **Options**: 0.9 (EAMv2) / Re-tune / Make configurable
- **Recommendation**: Start with 0.9, re-tune if needed

**Q6**: Layer-by-layer Stewart application?
- **Options**: Independent per layer / Cumulative / Implicit
- **Recommendation**: Independent per layer (EAMv2 approach)

**Q2**: Kinetic exponent n?
- **Options**: 0.5 / 0.58 / Configurable
- **Recommendation**: 0.58 (Stewart value), make namelist parameter

---

### 🟠 MEDIUM PRIORITY (Decide for Phase 3-4)

**Q7**: Kinetic effects during condensation?
- **Options**: Equilibrium only / Include kinetic
- **Recommendation**: Equilibrium only (standard)

**Q8**: P3 supersaturation adjustment for isotopes?
- **Options**: Calculate from ΔQ / Modify P3 / Estimate
- **Recommendation**: Modify P3 to output rate

**Q9**: Supersaturation-dependent kinetic fractionation?
- **Options**: Strong / Weak / Match EAMv2
- **Recommendation**: Match EAMv2 initially

**Q10**: Homogeneous vs heterogeneous nucleation?
- **Options**: Treat differently / Treat same
- **Recommendation**: Different (homogeneous → α_k = 1)

---

### 🟡 LOW PRIORITY (Decide for Phase 5-6)

**Q11**: Bergeron process implementation?
- **Options**: Separate steps / Combined / Ignore
- **Recommendation**: Separate evap + deposition

**Q12**: Freezing fractionation?
- **Options**: Include small effect / Neglect
- **Recommendation**: Include for slow, neglect for fast

**Q13**: Melting fractionation?
- **Options**: Include / Neglect
- **Recommendation**: Neglect (α = 1 at 0°C)

**Q14**: Riming fractionation?
- **Options**: Equilibrium / Kinetic / Intermediate
- **Recommendation**: Small kinetic effect

**Q15**: Separate rime mass isotope tracking?
- **Options**: Yes / No
- **Recommendation**: No (combine into ice)

---

### 🟢 DEFERRED (Can decide later or make configurable)

**Q1**: Temperature range validity for fractionation formulas?
**Q3**: Kinetic effects in sublimation?
**Q16**: CFL sub-stepping for isotope fractionation?
**Q17**: Clear-sky vs in-cloud fractionation?
**Q18**: Convective fractionation details?
**Q19**: Convective vs stratiform tracking approach?

---

## Parameter Values & Tuning

### Fixed Parameters (From Literature)

```fortran
! Equilibrium fractionation (Horita & Wesolowski 1994)
! HDO liquid-vapor: α_eq(T) = exp(1158.8/T² - 1620.1/T + 0.5223)
real(r8), parameter :: alp_hdo_a = 1158.8_r8
real(r8), parameter :: alp_hdo_b = -1620.1_r8
real(r8), parameter :: alp_hdo_c = 0.5223_r8

! Diffusivity ratios (Craig & Gordon 1965, Merlivat 1978)
real(r8), parameter :: difr_hdo   = 0.9755_r8  ! D_H2O / D_HDO
real(r8), parameter :: difr_h218o = 0.9723_r8  ! D_H2O / D_H218O
real(r8), parameter :: difr_h217o = 0.9847_r8  ! D_H2O / D_H217O

! Standard ratios (VSMOW)
real(r8), parameter :: rstd_hdo   = 155.76e-6_r8
real(r8), parameter :: rstd_h218o = 2005.20e-6_r8
real(r8), parameter :: rstd_h217o = 379.9e-6_r8
```

### Tunable Parameters (May Need Adjustment)

```fortran
! Kinetic fractionation exponent
real(r8) :: dkfac = 0.58_r8  ! Default: Stewart (1975)
! Range: 0.5 (Jouzel) to 1.0 (molecular)
! Namelist: wtrc_dkfac

! Stewart model equilibration parameter
real(r8) :: phi = 0.9_r8  ! Default: Bony et al. (2008)
! Range: 0.0 (no equilibration) to 1.0 (full)
! Namelist: wtrc_phi

! Temperature for kinetic effects (ice appears)
real(r8) :: tkini = 253.15_r8  ! -20°C
! Below this T: kinetic effects turn on for ice
! Range: 253-258 K
! Namelist: wtrc_tkini

! Minimum tracer threshold (prevent negative values)
real(r8) :: qmin_tracer = 1.0e-12_r8
! Namelist: wtrc_qmin
```

### Sensitivity Tests Needed

After Phase 7, run sensitivity tests:
1. **dkfac**: 0.5, 0.58, 0.7 → Impact on δD-RH relationship
2. **phi**: 0.7, 0.9, 1.0 → Impact on precipitation δD
3. **tkini**: 250, 253, 258 K → Impact on high-altitude/polar clouds

---

## Validation Targets

### Precipitation (Surface)

**GNIP Network** (Global Network of Isotopes in Precipitation):
- 1000+ stations worldwide
- Monthly δD and δ¹⁸O measurements
- Validation metrics:
  - **Global mean**: δD ≈ -20‰, δ¹⁸O ≈ -3‰
  - **Temperature gradient**: δ¹⁸O ≈ 0.5-0.7‰ per °C
  - **Latitude gradient**: More depleted toward poles
  - **Altitude effect**: δ¹⁸O ≈ -0.2‰ per 100m
  - **Seasonal cycle**: ±50‰ δD at high latitudes

**Meteoric Water Line**:
```
δD = 8 × δ¹⁸O + 10  (global average)

Slope should be ~8 (Rayleigh fractionation)
Intercept ~10 (deuterium excess)
```

### Vapor (Atmospheric Column)

**TES Satellite** (Tropospheric Emission Spectrometer):
- Upper troposphere δD measurements
- Vertical profiles
- Validation: δD_vapor ≈ -200 to -400‰ at 200-300 hPa

**AIRS Satellite**:
- Similar to TES, broader coverage
- δD in middle/upper troposphere

### Aircraft Campaigns

**HIPPO, ATom**: Vertical profiles from aircraft
- δD and δ¹⁸O in vapor
- Multiple latitudes and seasons

---

## Key References

### Fractionation Theory
1. **Horita & Wesolowski (1994)**: "Liquid-vapor fractionation of oxygen and hydrogen isotopes of water from the freezing to the critical temperature"
   - Definitive liquid-vapor equilibrium fractionation

2. **Merlivat & Nief (1967)**: "Fractionnement isotopique lors des changements d'état solide-vapeur et liquide-vapeur de l'eau"
   - Ice-vapor equilibrium fractionation

3. **Majoube (1971)**: "Fractionnement en oxygène 18 et en deutérium entre l'eau et sa vapeur"
   - Temperature dependence of fractionation

### Rain Evaporation
4. **Stewart (1975)**: "Stable isotope fractionation due to evaporation and isotopic exchange of falling waterdrops"
   - THE reference for rain evaporation model

5. **Bony et al. (2008)**: "An assessment of the primary sources of spread of global warming estimates from coupled atmosphere-ocean models"
   - Tuning parameter φ = 0.9

6. **Jouzel & Merlivat (1984)**: "Deuterium and oxygen 18 in precipitation: modeling of the isotopic effects during snow formation"
   - Ice cloud fractionation, kinetic exponents

### Atmospheric Applications
7. **Galewsky et al. (2016)**: "Stable isotopes in atmospheric water vapor and applications to the hydrologic cycle"
   - Review of atmospheric isotope processes

8. **Noone (2012)**: "Pairing measurements of the water vapor isotope ratio with humidity to deduce atmospheric moistening and dehydration in the tropical midtroposphere"
   - Interpretation of isotope observations

---

## Decision Log

*To be updated as questions are resolved*

| Question | Decision | Date | Rationale |
|----------|----------|------|-----------|
| Q5: φ parameter | Start with 0.9 | TBD | Match EAMv2-wiso, re-tune if needed |
| Q2: dkfac | Use 0.58, make namelist parameter | TBD | Stewart value, allow tuning |
| ... | ... | ... | ... |

---

## Notes for Developers

### When You Encounter Ambiguity:

1. **Check this document first**: See if question is already listed
2. **Check EAMv2-wiso code**: What did MG2 implementation do?
3. **Check literature**: Is there a standard approach?
4. **Document your choice**: Add to Decision Log above
5. **Make it configurable**: If uncertain, create namelist parameter

### Adding New Questions:

```markdown
**❓ QUESTION N: [Brief title]**
- **Issue**: [What's the ambiguity?]
- **Options**: [List options A, B, C]
- **EAMv2-wiso approach**: [What did MG2 do?]
- **Literature**: [What do papers say?]
- **Recommendation**: [Suggested approach]
- **Priority**: [High/Medium/Low]
```

### Physics Checks:

After implementing any process, check:
1. **Mass conservation**: Total water unchanged
2. **Direction of fractionation**: Heavy isotope goes to condensed phase
3. **Magnitude**: δD changes by ~10-50‰, not 1000‰
4. **Temperature dependence**: Stronger fractionation at lower T
5. **RH dependence** (evaporation): Stronger effect at lower RH

---

## Quick Reference: Expected Values

### Typical Isotope Ratios

| Location/Process | δD (‰) | δ¹⁸O (‰) |
|-----------------|--------|----------|
| Ocean surface | 0 | 0 |
| Tropical rain | -20 to -80 | -3 to -12 |
| Mid-latitude rain | -50 to -100 | -7 to -15 |
| Polar precipitation | -200 to -400 | -30 to -55 |
| Stratosphere vapor | -400 to -700 | -60 to -100 |

### Fractionation Factors (Approximate)

| Process | Temperature | α (HDO) | α (H₂¹⁸O) |
|---------|-------------|---------|-----------|
| Liquid-vapor | 20°C | 1.079 | 1.0092 |
| Liquid-vapor | 0°C | 1.084 | 1.0098 |
| Ice-vapor | -20°C | 1.176 | 1.0135 |
| Ice-vapor | -40°C | 1.223 | 1.0156 |

**Rule of Thumb**: Ice fractionation is about twice as strong as liquid

---

**Document Status**: Living - update frequently  
**Last Updated**: 2026-04-23  
**Next Review**: After Phase 2 completion (when first fractionation is implemented)

**Questions?** Add them to this document and discuss in weekly sync!
