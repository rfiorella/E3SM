# P3 Water Isotope Implementation - Project Documentation

**Project**: Port water isotope infrastructure from EAMv2-wiso (MG2 microphysics) to EAMv3-wiso (P3 microphysics)  
**Version**: 1.0  
**Date**: 2026-04-21  
**Status**: Phase 0 Complete - Ready for Review

---

## Documentation Overview

This project contains comprehensive documentation for implementing water isotope tracers in the P3 microphysics scheme within E3SM Atmosphere Model version 3 (EAMv3-wiso).

### Core Documents

1. **[P3_ISOTOPE_IMPLEMENTATION_SPEC.md](./P3_ISOTOPE_IMPLEMENTATION_SPEC.md)** - Main specification document
   - Project overview and scope
   - Isotope species and water phase types
   - Phased implementation plan (10 phases)
   - Success criteria and dependencies
   - Questions for team review

2. **[P3_ISOTOPE_FUNCTION_MAPPING.md](./P3_ISOTOPE_FUNCTION_MAPPING.md)** - Detailed function mapping
   - MG2 to P3 process mapping tables
   - P3 function requirements (37 functions)
   - Water type transition matrix
   - Fractionation factor requirements
   - EAMv2-wiso function dependencies
   - Implementation checklist

3. **[P3_ISOTOPE_ISSUES.md](./P3_ISOTOPE_ISSUES.md)** - Issue templates
   - 50 detailed issue descriptions
   - Tasks and acceptance criteria
   - Effort estimates (190-250 person-days total)
   - Dependency chains
   - Priority labels

4. **[p3_isotope_kanban.json](./p3_isotope_kanban.json)** - Machine-readable kanban data
   - Importable into project management tools
   - Issue dependencies defined
   - Labels and columns configured
   - Ready for Cline, GitHub Projects, or similar

---

## Quick Start

### For Team Review

1. **Read the main specification**: [P3_ISOTOPE_IMPLEMENTATION_SPEC.md](./P3_ISOTOPE_IMPLEMENTATION_SPEC.md)
2. **Review the phased approach**: Pay special attention to Phases 1-2 (infrastructure and rain evaporation)
3. **Check the questions section**: Provide feedback on design decisions
4. **Review function mapping**: [P3_ISOTOPE_FUNCTION_MAPPING.md](./P3_ISOTOPE_FUNCTION_MAPPING.md) for technical details

### For Implementation

1. **Import kanban board**: Use `p3_isotope_kanban.json` in your project management tool
2. **Start with Phase 1**: Infrastructure porting (Issues #1.1 through #1.7)
3. **Follow dependency chain**: Each issue lists its dependencies
4. **Reference function mapping**: Detailed technical requirements for each function

### For Project Management

1. **Total Effort**: 190-250 person-days (155-205 days for core Phases 1-9)
2. **Critical Path**: Phase 0 → 1 → 2 → 3 → 4 → 7 → 9
3. **Phases**: 10 phases (Phase 10 optional)
4. **Issues**: 50 total issues with detailed specifications

---

## Project Structure

### Phase Breakdown

| Phase | Name | Issues | Est. Days | Priority | Status |
|-------|------|--------|-----------|----------|--------|
| 0 | Infrastructure Setup | 1 | 5 | Critical | ✓ Complete |
| 1 | Core Infrastructure Port | 7 | 15-20 | Critical | To Do |
| 2 | Rain Evaporation | 8 | 20-25 | Critical | To Do |
| 3 | Vapor-Liquid Phase Changes | 6 | 15-20 | High | To Do |
| 4 | Vapor-Ice Phase Changes | 6 | 20-25 | High | To Do |
| 5 | Freezing and Melting | 6 | 15-20 | Medium | To Do |
| 6 | Collection and Riming | 5 | 15-20 | Medium | To Do |
| 7 | Sedimentation | 5 | 20-25 | Critical | To Do |
| 8 | Convection Interface | 3 | 12-16 | High | To Do |
| 9 | Conservation & Diagnostics | 3 | 14-19 | Critical | To Do |
| 10 | Additional Physics | 4 | 15-22 | Low (Optional) | To Do |

### Critical Path

The fastest route to a working implementation:

```
Phase 0 (Planning) ✓
  ↓
Phase 1 (Infrastructure) → 15-20 days
  ↓
Phase 2 (Rain Evaporation) → 20-25 days [HIGHEST PRIORITY]
  ↓
Phase 3 (Vapor-Liquid) → 15-20 days
  ↓
Phase 4 (Vapor-Ice) → 20-25 days
  ↓
Phase 7 (Sedimentation) → 20-25 days
  ↓
Phase 9 (Conservation) → 14-19 days
```

**Total Critical Path**: ~120-160 days

---

## Isotope Species

All isotopes from EAMv2-wiso will be tracked:

1. **H2O** - Bulk water (standard)
2. **H216O** - Standard oxygen-16
3. **HDO** - Deuterium (δD)
4. **H218O** - Oxygen-18 (δ18O)
5. **H217O** - Oxygen-17 (δ17O)
6. **HTO** - Tritium (radioactive)

### Water Phase Types

Isotopes tracked in all phases:

- Water vapor
- Cloud liquid
- Cloud ice
- Stratiform rain
- Stratiform snow
- Convective rain (Phase 8)
- Convective snow (Phase 8)

---

## Key Scientific Processes

### Phase 2: Rain Evaporation (HIGHEST PRIORITY)

**Why First**: Below-cloud evaporation is the dominant control on precipitation isotope ratios reaching the surface.

**Physics**:
- Stewart (1975) rain evaporation model
- Equilibrium fractionation (Horita & Wesolowski 1994)
- Kinetic fractionation (diffusivity-based)
- Partial equilibration (droplet size dependent)

**Expected Behavior**:
- Rain evaporation makes δD more negative
- Effect stronger at lower relative humidity
- Larger drops equilibrate less

### Phases 3-6: Phase Changes and Collection

Progressive addition of:
- Condensation (Rayleigh distillation)
- Ice nucleation and growth (supersaturation effects)
- Freezing and melting
- Riming and collection

### Phase 7: Sedimentation

**Critical**: Layer-by-layer fractionation as precipitation falls

**Physics**:
- Stewart model applied per layer
- CFL-limited sub-stepping
- Raindrop size evolution
- Below-cloud fractionation

### Phase 8: Convection

Separate tracking of convective vs stratiform precipitation:
- ZM deep convection (updraft/downdraft fractionation)
- Shallow convection
- Detrainment to stratiform cloud

---

## Dependencies

### Source Code (EAMv2-wiso)

Located at: `/Users/rfiorella/code/E3SM/EAMv2-wiso`

**Key Files to Port**:
- `share/util/water_isotopes.F90` (777 lines)
- `share/util/water_types.F90` (171 lines)
- `components/eam/src/physics/cam/water_tracer_vars.F90` (122 lines)
- `components/eam/src/physics/cam/water_tracers.F90` (7,813 lines)

### Target Code (EAMv3-wiso)

Located at: `/Users/rfiorella/code/E3SM/EAMv3-wiso`

**Key Files to Modify**:
- `components/eam/src/physics/p3/eam/micro_p3.F90` (4,619 lines)
- `components/eam/src/physics/p3/eam/micro_p3_interface.F90` (1,724 lines)
- `components/eam/src/physics/cam/physpkg.F90`
- `components/eam/src/physics/cam/physics_types.F90`

**New Files to Create**:
- `components/eam/src/physics/p3/eam/water_tracers_p3.F90` (Phase 2)
- `components/eam/src/physics/p3/eam/micro_p3_isotopes.F90` (Phase 3+)

---

## Testing Strategy

### Unit Tests (Per Phase)

Each phase includes unit tests for:
- Individual function correctness
- Mass conservation
- Fractionation factor values
- Numerical stability

### Integration Tests

- Single-column model tests
- Idealized cases (tropical, Arctic, stratosphere)
- Comparison to EAMv2-wiso MG2 results

### Validation (Phase 9)

- Global simulation (multi-month)
- Comparison to observations:
  - GNIP precipitation isotope network
  - TES/AIRS satellite vapor isotopes
  - Aircraft campaigns
  - Ice core records
- Physics checks:
  - δD-δ18O meteoric water line (slope ~8)
  - Temperature correlation
  - Altitude effects
  - Seasonal cycle

---

## Questions for Review

Before proceeding to Phase 1 implementation:

### 1. Scope and Phasing

- ✓ Is the phased approach acceptable?
- ✓ Is rain evaporation first (Phase 2) the right priority?
- ? Should convection (Phase 8) be moved earlier?

### 2. Technical Approach

- ? Create new files (`water_tracers_p3.F90`) or extend existing?
- ? P3 condensation rate extraction approach (implicit in P3)?
- ? Level of documentation required for each phase?

### 3. Resources

- ? Who will be involved in coding, testing, validation?
- ? What is target completion date for each phase?
- ? Computational resources for testing/validation?

### 4. Validation

- ? Which observational datasets are priorities?
- ? Acceptable performance overhead (currently < 20% expected)?
- ? Minimum validation requirements before moving to next phase?

---

## Getting Started

### Next Steps After Review

1. **Phase 0 Review Meeting**: Discuss this documentation with science team
2. **Finalize Design Decisions**: Answer questions listed above
3. **Phase 1 Kickoff**: Begin infrastructure porting
4. **Set Up Development Environment**:
   - Create feature branch
   - Set up testing framework
   - Configure CI/CD if available

### Recommended Team Structure

- **Lead Developer**: Overall coordination, complex modules
- **P3 Specialist**: Modifications to P3 microphysics
- **Isotope Physics Expert**: Fractionation implementations
- **Testing/Validation Lead**: Test design and validation
- **Code Reviewer**: Review all changes

---

## References

### Scientific References

- **Stewart (1975)**: Rain evaporation isotope model
- **Horita & Wesolowski (1994)**: Liquid-vapor equilibrium fractionation
- **Merlivat & Nief (1967)**: Ice-vapor equilibrium fractionation
- **Majoube (1971)**: Temperature dependence of fractionation
- **Bony et al. (2008)**: Below-cloud relative humidity effects
- **Morrison & Milbrandt (2015)**: P3 microphysics scheme (J. Atmos. Sci.)
- **Milbrandt & Morrison (2016)**: P3 predicted particle properties (J. Atmos. Sci.)

### Code References

- **EAMv2-wiso Repository**: `/Users/rfiorella/code/E3SM/EAMv2-wiso`
- **EAMv3-wiso Repository**: `/Users/rfiorella/code/E3SM/EAMv3-wiso`
- **E3SM Documentation**: https://e3sm.org/
- **P3 Scheme**: Morrison & Milbrandt (2015, 2016)

---

## Support and Contact

For questions about this implementation plan:

1. Review the detailed documents first
2. Check the function mapping for technical details
3. Consult issue templates for specific tasks
4. Contact project leads (TBD during Phase 0 review)

---

## Document History

| Version | Date | Changes | Author |
|---------|------|---------|--------|
| 1.0 | 2026-04-21 | Initial documentation | OpenCode Analysis |

---

## License and Attribution

This implementation plan is based on:
- Water isotope infrastructure from EAMv2-wiso
- P3 microphysics scheme in EAMv3-wiso
- E3SM modeling framework

Please cite the original P3 papers (Morrison & Milbrandt 2015, 2016) and acknowledge the EAMv2-wiso isotope implementation when using this work.

---

**Status**: Phase 0 Complete - Awaiting Team Review ✓

**Next Milestone**: Phase 1 Kickoff - Core Infrastructure Port

---
