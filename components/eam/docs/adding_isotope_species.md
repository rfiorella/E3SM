# Adding Isotope Species to EAMv3 Water Isotope Implementation

This document describes how to extend the EAMv3 water isotope implementation from
the default 2-species configuration (HDO + H2-18O) to include all 6 supported species.

## Supported Species

The water isotope infrastructure (defined in `share/util/water_isotopes.F90`) supports
6 species:

| Index | Species | Name     | Description                        |
|-------|---------|----------|------------------------------------|
| 1     | H2O     | `isph2o` | Standard water (bulk tracer)       |
| 2     | H2-16O  | `isph216o`| Light isotopologue                |
| 3     | HDO     | `isphdo` | Deuterated water                   |
| 4     | H2-18O  | `isph218o`| Oxygen-18 water                   |
| 5     | H2-17O  | `isph217o`| Oxygen-17 water                   |
| 6     | HTO     | `isphto` | Tritiated water                    |

The default configuration uses species 1 (H2O as bulk tracer), 3 (HDO), and 4 (H2-18O).

## Configuration Steps

### 1. Build System Configuration

The number of water tracers is set at configure time via `WTRC_MAX_CNST`. Each species
requires 7 tracers (one per water type: vapor, liquid, ice, stratiform rain, stratiform
snow, convective rain, convective snow).

For the default 2-species configuration (H2O + HDO + H2-18O = 3 species x 7 types):
```
WTRC_MAX_CNST = 21
```

For all 6 species (6 x 7 = 42 tracers):
```
WTRC_MAX_CNST = 42
```

This is set via the `configure` script. Look for `-wtrc_max_cnst` or the
`SET_WTRC_MAX_CNST` preprocessor variable in `src/physics/cam/water_tracer_vars.F90`.

### 2. Namelist Configuration

In the `water_tracer_nl` namelist group, configure the tracer names, species, and types.
The following example adds H2-17O (species 5) to the default HDO + H2-18O configuration:

```fortran
&water_tracer_nl
  trace_water    = .true.
  wisotope       = .true.
  water_tracer_model = 'isotope'

  ! Each entry defines one tracer: name, species, type, surface vapor name, surface precip name
  ! H2O tracers (bulk, species 1)
  wtrc_names(1)  = 'H2OV'
  wtrc_names(2)  = 'H2OL'
  wtrc_names(3)  = 'H2OI'
  wtrc_names(4)  = 'H2ORS'
  wtrc_names(5)  = 'H2OSS'
  wtrc_names(6)  = 'H2ORC'
  wtrc_names(7)  = 'H2OSC'

  wtrc_species_names(1:7)  = 7*'H2O'
  wtrc_type_names(1)  = 'vapor'
  wtrc_type_names(2)  = 'liquid'
  wtrc_type_names(3)  = 'ice'
  wtrc_type_names(4)  = 'strat_rain'
  wtrc_type_names(5)  = 'strat_snow'
  wtrc_type_names(6)  = 'conv_rain'
  wtrc_type_names(7)  = 'conv_snow'

  ! HDO tracers (species 3)
  wtrc_names(8)  = 'HDOV'
  wtrc_names(9)  = 'HDOL'
  wtrc_names(10) = 'HDOI'
  wtrc_names(11) = 'HDORS'
  wtrc_names(12) = 'HDOSS'
  wtrc_names(13) = 'HDORC'
  wtrc_names(14) = 'HDOSC'
  wtrc_species_names(8:14) = 7*'HDO'
  ! ... type_names repeat the same pattern ...

  ! H2-18O tracers (species 4)
  wtrc_names(15) = 'H218V'
  ! ... etc ...

  ! H2-17O tracers (species 5) -- NEW
  wtrc_names(22) = 'H217V'
  wtrc_names(23) = 'H217L'
  wtrc_names(24) = 'H217I'
  wtrc_names(25) = 'H217R'  ! strat rain
  wtrc_names(26) = 'H217S'  ! strat snow
  wtrc_names(27) = 'H21RC'  ! conv rain
  wtrc_names(28) = 'H21SC'  ! conv snow
  wtrc_species_names(22:28) = 7*'H217O'
  ! ... type_names repeat ...
/
```

Note: Tracer names are limited to 8 characters for history file compatibility.
Choose names carefully.

### 3. Species-Specific Fractionation

Each species has unique fractionation behavior defined in `share/util/water_isotopes.F90`:

#### H2-17O (species 5)
- Uses the mass-dependent scaling: alpha_17O = alpha_18O^0.529
- This is the canonical mass-dependent fractionation exponent
- The 17O-excess (Delta-17O = delta-17O - 0.528 * delta-18O) is a tracer of
  kinetic vs. equilibrium fractionation processes
- No code changes needed — the fractionation routines (`wiso_alpl`, `wiso_alpi`)
  already handle H2-17O correctly

#### HTO (Tritiated Water, species 6)
- Uses enhanced fractionation: alpha_HTO = alpha_HDO^2.0 (approximately)
- HTO undergoes radioactive decay with half-life of 12.32 years
- The decay is handled by `wtrc_rad_decay()` in `water_tracers.F90`
- HTO is primarily useful as a tracer of stratospheric water vapor residence time
  and nuclear fallout

### 4. Initial Conditions

Each new species requires initial condition fields on the IC file. The fields should be
named to match the `wtrc_names` entries in the namelist. If fields are not found on the
IC file, they will be initialized to zero.

For a proper initialization at natural abundance:
- H2-17O: Use `rstd` ratio from `water_isotopes.F90` (currently set to 1.0 in model
  coordinates; actual SMOW ratio is defined in `rnat`)
- HTO: Initialize to zero (no natural background; only relevant for bomb-test or
  cosmogenic production studies)

### 5. Performance Considerations

| Configuration | Tracers | Approx. Cost Increase |
|--------------|---------|----------------------|
| HDO only     | 7       | ~5-10%               |
| HDO + H2-18O | 14      | ~10-20%              |
| HDO + H2-18O + H2-17O | 21 | ~15-30%         |
| All 6 species | 42     | ~30-50%              |

The cost increase is primarily from:
1. **Advection**: Each tracer is advected by the dynamical core. This scales linearly
   with tracer count and dominates the cost.
2. **Microphysics**: The `wtrc_niter` iterations in `wtrc_p3_inter()` loop over all
   water sets. Cost scales linearly with number of species.
3. **Convection**: Similar scaling to microphysics.
4. **I/O**: More fields to read/write. Can be mitigated by selective history output.

Recommendations:
- For most climate studies, HDO + H2-18O is sufficient
- Add H2-17O only if studying kinetic fractionation processes (e.g., leaf water,
  stratospheric chemistry)
- HTO is only relevant for specialized tracer studies
- H2-16O is primarily a diagnostic tracer for mass conservation checking

### 6. Validation

When adding a new species, verify:
1. **Mass conservation**: The sum of all isotope species for each water type should
   equal the bulk water (within machine precision for H2O + H2-16O)
2. **Delta values**: Check that delta values are in physically reasonable ranges:
   - delta-17O: approximately 0.529 * delta-18O for equilibrium processes
   - delta-T (HTO): no standard expectation, but should be positive for bomb-test era
3. **BFB**: Adding species should not affect non-isotope simulations when
   `wisotope=.false.`

### 7. Code Changes Required

Adding species within the existing 6 requires **no code changes** — only namelist and
build configuration. The fractionation factors, molecular weights, diffusivity ratios,
and natural abundance ratios are already defined in `water_isotopes.F90` for all 6
species.

To add an entirely new species beyond the existing 6 would require:
1. Incrementing `pwtspec` in `water_isotopes.F90`
2. Adding fractionation factor parameterizations (`wiso_alpl`, `wiso_alpi`)
3. Adding molecular weight, diffusivity ratio, and natural abundance ratio
4. Recompiling with updated `WTRC_MAX_CNST`
