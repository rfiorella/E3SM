# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the **E3SM (Energy Exascale Earth System Model)** repository, a state-of-the-art fully coupled Earth system model for climate research. This branch (`rfiorella/eam-v3/wiso`) focuses on implementing **water isotope tracers in the P3 microphysics scheme** (EAMv3-wiso).

**Water Isotope Project**: Porting water isotope infrastructure from EAMv2-wiso (MG2 microphysics) to EAMv3-wiso (P3 microphysics). See `P3_ISOTOPE_README.md` for comprehensive project documentation, including a 10-phase implementation plan (~190-250 person-days effort).

## Build System (CIME)

E3SM uses **CIME** (Common Infrastructure for Modeling the Earth) for case management, building, and testing.

### Creating and Building a Case

```bash
# 1. Create a new case
cime/scripts/create_newcase --case <CASE_DIR> \
                            --compset <COMPSET> \
                            --res <RESOLUTION> \
                            --machine <MACHINE> \
                            --project <PROJECT>

# Example for standard pre-industrial coupled run:
cime/scripts/create_newcase --case ~/cases/my_test \
                            --compset WCYCL1850 \
                            --res ne30pg2_r05_IcoswISC30E3r5 \
                            --machine pm-cpu

# 2. Navigate to case directory and configure
cd <CASE_DIR>
./case.setup

# 3. Build the model
./case.build

# 4. Submit the run
./case.submit
```

**Common compsets**:
- `WCYCL1850`: Pre-industrial climatological forcing
- `WCYCL20TR`: Historical transient forcing (1850-2014)
- `WCYCLSSP585`: SSP-585 future scenario

**Common resolutions**:
- `ne30pg2_r05_IcoswISC30E3r5`: ~100km atmosphere, 0.5° land, 30km ocean
- `northamericax4v1pg2_r025_IcoswISC30E3r5`: North America regionally refined

### Testing

```bash
# Create a single test
cime/scripts/create_test <TEST_NAME> --machine <MACHINE>

# Example: ERS (exact restart) test
cime/scripts/create_test ERS.f19_g16.IELM --machine pm-cpu

# Create test suite
cime/scripts/create_test --test-suite e3sm_developer --machine <MACHINE>
```

Test suites are defined in `cime_config/tests.py`. Test format: `<test_type>.<grid>.<compset>[.<testmod>]`

### Case Management Commands (run from case directory)

```bash
./case.setup          # Set up the case
./case.build          # Build the executable
./case.build --clean  # Clean build
./case.submit         # Submit the run
./xmlquery <VAR>      # Query case variable
./xmlchange <VAR>=<VALUE>  # Change case variable
```

## Repository Structure

### Major Components

```
components/
├── eam/          # E3SM Atmosphere Model (EAM) - formerly CAM
├── eamxx/        # EAM C++ refactor (next-gen atmosphere)
├── elm/          # E3SM Land Model
├── mpas-ocean/   # MPAS Ocean
├── mpas-seaice/  # MPAS Sea Ice
├── mosart/       # River routing model
├── cice/         # Sea ice model (alternative)
└── *_comps/      # Data/stub component models
```

### EAM Atmosphere Model Structure

```
components/eam/src/
├── physics/
│   ├── cam/            # Core physics (formerly CAM)
│   │   ├── physpkg.F90           # Main physics driver
│   │   ├── physics_types.F90     # Physics state variables
│   │   ├── water_tracers.F90     # Water isotope tracer code (7,813 lines)
│   │   ├── water_tracer_vars.F90 # Water tracer variables/namelist
│   │   ├── zm_conv_intr.F90      # Zhang-McFarlane convection
│   │   └── micro_mg_cam.F90      # MG2 microphysics (old)
│   ├── p3/             # P3 microphysics (current)
│   │   └── eam/
│   │       ├── micro_p3.F90           # Main P3 code (4,619 lines)
│   │       ├── micro_p3_interface.F90 # EAM-P3 interface (1,724 lines)
│   │       └── micro_p3_utils.F90     # P3 utilities
│   ├── clubb/          # CLUBB turbulence/cloud scheme
│   ├── rrtmgp/         # Radiation (RRTMGP)
│   └── cosp2/          # COSP satellite simulator
├── control/           # Model control and configuration
└── cpl/              # Coupler interface
```

### Water Isotope Infrastructure

**Shared utilities** (`share/util/`):
- `water_isotopes.F90`: Fractionation factors, isotope calculations (777 lines)
- `water_types.F90`: Water phase type definitions (171 lines)

**EAM isotope code** (`components/eam/src/physics/cam/`):
- `water_tracer_vars.F90`: Tracer variables, namelist settings (122 lines)
- `water_tracers.F90`: MG2-based isotope implementation (7,813 lines)

**Target for new P3 isotope code** (to be created):
- `components/eam/src/physics/p3/eam/water_tracers_p3.F90` (Phase 2)
- `components/eam/src/physics/p3/eam/micro_p3_isotopes.F90` (Phase 3+)

## Water Isotope Implementation Project

### Current Status
- **Phase 0**: Complete (planning and documentation)
- **Next**: Phase 1 (core infrastructure porting, 15-20 days)

### Key Documentation Files (in repo root)
- `P3_ISOTOPE_README.md`: Overview and quick start
- `P3_ISOTOPE_IMPLEMENTATION_SPEC.md`: Main technical specification (27KB)
- `P3_ISOTOPE_FUNCTION_MAPPING.md`: MG2→P3 function mapping (30KB)
- `P3_ISOTOPE_ISSUES.md`: 50 detailed issue templates (55KB)
- `P3_ISOTOPE_DEPENDENCIES.md`: Dependency diagrams and timeline
- `p3_isotope_kanban.json`: Machine-readable project management data

### Isotope Species Tracked
1. **H2O** - Bulk water (standard)
2. **H216O** - Standard oxygen-16
3. **HDO** - Deuterium (δD)
4. **H218O** - Oxygen-18 (δ18O)
5. **H217O** - Oxygen-17 (δ17O)
6. **HTO** - Tritium (radioactive)

### Water Phase Types
Isotopes tracked in: vapor, cloud liquid, cloud ice, stratiform rain/snow, convective rain/snow

### Phased Implementation (10 phases, ~190-250 person-days)
1. **Phase 1**: Core infrastructure port (15-20 days) - CRITICAL PATH
2. **Phase 2**: Rain evaporation (20-25 days) - HIGHEST PRIORITY
3. **Phase 3**: Vapor-liquid phase changes (15-20 days)
4. **Phase 4**: Vapor-ice phase changes (20-25 days)
5. **Phase 5**: Freezing and melting (15-20 days)
6. **Phase 6**: Collection and riming (15-20 days)
7. **Phase 7**: Sedimentation (20-25 days) - CRITICAL PATH
8. **Phase 8**: Convection interface (12-16 days)
9. **Phase 9**: Conservation & diagnostics (14-19 days) - CRITICAL PATH
10. **Phase 10**: Additional physics (15-22 days) - OPTIONAL

**Critical Path Total**: ~120-160 days (Phases 0→1→2→3→4→7→9)

### Source Code Locations

**Reference code (EAMv2-wiso)**: `/Users/rfiorella/code/E3SM/EAMv2-wiso`
- Water isotope infrastructure already ported from here to `share/util/`
- MG2-based `water_tracers.F90` serves as reference for P3 implementation

**Target code (EAMv3-wiso)**: `/Users/rfiorella/code/E3SM/EAMv3-wiso` (this repo)

## Development Workflow

### Git Workflow
- **Main branch**: `master`
- **Current branch**: `rfiorella/eam-v3/wiso`
- Create feature branches for specific phases/issues
- Follow E3SM contribution guidelines (see `CONTRIBUTING.md`)

### Code Modification Guidelines
1. **Read before modifying**: Always read existing code before suggesting changes
2. **Water isotope namelist**: Check `water_tracer_vars.F90` for configuration flags
3. **P3 modifications**: Changes to P3 should maintain bit-for-bit reproducibility when isotopes are OFF
4. **Conservation**: All water mass must be conserved (bulk + all isotopes)
5. **Fortran conventions**: 
   - Use `r8 => shr_kind_r8` for real precision
   - Follow existing naming conventions (e.g., `q` for water vapor, `qc` for cloud liquid)
   - Use explicit `intent(in/out/inout)` for all subroutine arguments

### Testing Strategy
1. **Unit tests**: Test individual functions for correctness, conservation, numerical stability
2. **Single-column tests**: Idealized cases (tropical, Arctic, stratosphere)
3. **Integration tests**: Multi-month global runs with conservation checks
4. **Validation**: Compare to observations (GNIP, TES/AIRS, aircraft, ice cores)

### File Naming Conventions
- Physics modules: `<scheme>_<component>.F90` (e.g., `micro_p3.F90`)
- Interface modules: `<scheme>_interface.F90`
- Variable modules: `<name>_vars.F90`
- Type definitions: `<name>_types.F90`

## Common Development Patterns

### Physics Parameterization Loop Structure
```fortran
! Typical pattern in EAM physics
do i = 1, ncol  ! Loop over columns
  do k = 1, pver  ! Loop over vertical levels
    ! Physics calculations here
  end do
end do
```

### Water Isotope Tracer Loop Pattern (from water_tracers.F90)
```fortran
! Loop over isotope species
do m = ixwti, ixwtx
  ispec = iwspec(m)
  if (ispec > 0) then
    ! Apply fractionation for this isotope species
    ! ispec values: isph2o, isph216o, isphdo, isph218o, isph217o, isphto
  end if
end do
```

### Fractionation Factor Pattern
```fortran
! Equilibrium fractionation (liquid-vapor)
alpl = wiso_alpl(ispec, temp)

! Kinetic fractionation  
akel = wiso_akel(ispec, relhum, temp)

! Combined fractionation during evaporation
alpha_total = alpl * akel
```

## Key Physics Concepts

### P3 Microphysics
- **P3** = Predicted Particle Properties scheme (Morrison & Milbrandt 2015)
- Predicts ice particle properties (rimed fraction, density, size) rather than fixed categories
- More physically realistic than older schemes (MG2)

### Water Isotope Fractionation
- **Equilibrium fractionation**: Temperature-dependent (Horita & Wesolowski 1994 for liquid, Merlivat & Nief 1967 for ice)
- **Kinetic fractionation**: Depends on diffusivity and relative humidity
- **Stewart (1975) model**: Below-cloud evaporation (dominant control on δD in precipitation)
- **Rayleigh distillation**: Progressive depletion during condensation

### Why Rain Evaporation is Priority
Below-cloud evaporation is the dominant control on precipitation isotope ratios reaching the surface. Getting this right first enables validation against surface observations (GNIP network).

## Important Caveats

1. **Computational cost**: E3SM requires HPC resources; most development is on DOE supercomputers
2. **Input data**: Large input datasets required (managed via CIME input data system)
3. **Compilation time**: Full builds can take 30+ minutes on typical HPC systems
4. **Testing time**: Even short tests (XS_2x5_ndays) require batch job submission
5. **Submodules**: Several components are git submodules (cime, YAKL, scorpio, etc.) - handle with care

## Useful References

### Documentation
- E3SM main site: https://e3sm.org
- E3SM docs: https://docs.e3sm.org/E3SM/
- CIME docs: https://esmci.github.io/cime/
- Running E3SM guide: https://e3sm.org/model/running-e3sm/
- GitHub discussions: https://github.com/E3SM-Project/E3SM/discussions

### Scientific References
- **P3 scheme**: Morrison & Milbrandt (2015, 2016) J. Atmos. Sci.
- **Stewart model**: Stewart (1975) - rain evaporation isotope model
- **Equilibrium fractionation**: Horita & Wesolowski (1994), Merlivat & Nief (1967)
- **Kinetic fractionation**: Merlivat & Jouzel (1979)

### Template Run Script
See `run_e3sm.template.sh` for example case creation and configuration.
