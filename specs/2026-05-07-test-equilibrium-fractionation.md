---
spec_type: model-e3sm
spec_version: 1.0
execution_mode: checkpoint
created_date: 2026-05-07
work_summary: Test equilibrium fractionation functions for correctness
priority: normal
estimated_effort_hours: 2

# Contact and ownership
primary_contact: rfiorella
reviewers: []

# Model-specific fields
model_specific:
  subsystem: eamxx
  build_mode: eamxx-standalone-cmake
  target_compset: null
  target_resolution: ne4
  platform: docker-local
  container_image: rfiorella/model-containers:e3sm-openmpi-dev-latest
  validation_tier: unit-test

# Inputs
inputs:
  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/water_isotopes.hpp
    description: Header file containing AlphaEqLiquidVapor() and AlphaEqIceVapor() functions
    format: C++ header
    required: true

# Deliverables
deliverables:
  - path: /code/E3SM/EAMXX-wiso/components/eamxx/src/physics/water_tracers/tests/test_equilibrium_fractionation.cpp
    description: Catch2-based unit test file for equilibrium fractionation functions
    format: C++ source
    validation_method: compiles and passes all tests

# Success criteria
success_criteria:
  - id: SC1
    phase: implementation
    description: Test file compiles successfully in standalone EAMxx build
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | tee build.log && \
      grep -q "test_equilibrium_fractionation" build.log
    expected_output: Build completes without errors and test binary is created
    blocking: true

  - id: SC2
    phase: implementation
    description: H2O and H216O always return 1.0 (no self-fractionation) for both functions across all temperatures
    criterion_type: assertion
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -A 5 "test_equilibrium_fractionation"
    assertion: All test cases for H2O and H216O self-fractionation pass
    blocking: true

  - id: SC3
    phase: implementation
    description: AlphaEqLiquidVapor() produces correct fractionation factors for HDO, H218O, H217O, HTO across temperature range 233-313K
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -A 10 "AlphaEqLiquidVapor"
    tolerance: absolute 1e-10
    reference_values: Hand-calculated from Horita & Wesolowski (1994) polynomial formulas
    blocking: true

  - id: SC4
    phase: implementation
    description: AlphaEqIceVapor() produces correct fractionation factors for HDO, H218O, H217O, HTO across temperature range 233-313K
    criterion_type: tolerance
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | \
      grep -A 10 "AlphaEqIceVapor"
    tolerance: absolute 1e-10
    reference_values: Hand-calculated from Merlivat (1978) and related papers
    blocking: true

  - id: SC5
    phase: implementation
    description: All tests pass in standalone EAMxx test suite
    criterion_type: shell
    command: |
      cd /code/E3SM/EAMXX-wiso/components/eamxx && \
      ./scripts/test-all-eamxx --preserve-env -m docker -t dbg 2>&1 | tee test.log && \
      grep -q "test_equilibrium_fractionation.*PASSED" test.log
    expected_output: Test execution completes with all test cases passing
    blocking: true

# Out of scope
out_of_scope:
  - Kinetic fractionation testing (deferred to separate effort)
  - Integration of fractionation functions into physics parameterizations
  - Performance optimization of fractionation calculations
  - Testing beyond the Earth-relevant temperature range (233-313K)
  - Full-model simulation testing

# Resolved decisions
resolved_decisions:
  - decision: Use hand-calculated reference values from polynomial formulas in published papers
    rationale: Provides independent validation without circular dependency on implementation
    date: 2026-05-07
    
  - decision: Test approximately 10 representative temperature points across the 233-313K range
    rationale: Balances thorough coverage with practical test runtime
    date: 2026-05-07
    
  - decision: Use absolute tolerance of 1e-10 for fractionation factor comparisons
    rationale: Appropriate for double-precision floating-point calculations of physical constants
    date: 2026-05-07
    
  - decision: Place test file in components/eamxx/src/physics/water_tracers/tests/
    rationale: Follows existing EAMxx test organization pattern
    date: 2026-05-07
    
  - decision: Use Catch2 test framework with simple REQUIRE statements
    rationale: Consistent with existing EAMxx unit test patterns
    date: 2026-05-07

# Ask-before actions (project-specific additions to global policy)
ask_before:
  - Modifying any source files outside the tests/ directory
  - Adding new dependencies or external libraries

# Parallelization
allow_parallelization: false

# Post-completion review
request_performance_review: false
request_code_review: false
---

# Test Equilibrium Fractionation Functions for Correctness

## Context

The EAMXX-wiso project adds water isotope tracking to the EAMxx atmosphere component. Two key functions calculate equilibrium fractionation factors between water phases:

- `AlphaEqLiquidVapor(isotope, temperature)` - liquid/vapor equilibrium fractionation
- `AlphaEqIceVapor(isotope, temperature)` - ice/vapor equilibrium fractionation

These functions are defined in `water_isotopes.hpp` and implement published polynomial formulas from Horita & Wesolowski (1994), Merlivat (1978), and related papers.

Before integrating these functions into physics parameterizations, we need unit tests validating that:

1. The functions produce physically correct results (no self-fractionation for H2O and H216O)
2. The polynomial implementations match expected values from the source papers
3. The functions behave correctly across the full Earth-relevant temperature range

This work creates a focused unit test suite for these two functions, independent of the broader physics integration.

## Approach

### Test Structure

Create a new Catch2-based test file at `components/eamxx/src/physics/water_tracers/tests/test_equilibrium_fractionation.cpp`.

The test file will include the following test cases:

1. **Self-fractionation invariant test**: Verify that H2O and H216O always return exactly 1.0 (no fractionation) for both `AlphaEqLiquidVapor()` and `AlphaEqIceVapor()` across all test temperatures.

2. **Liquid-vapor fractionation accuracy test**: For each of HDO, H218O, H217O, and HTO:
   - Test at approximately 10 representative temperatures spanning 233K to 313K
   - Compare computed fractionation factors to hand-calculated reference values
   - Use absolute tolerance of 1e-10

3. **Ice-vapor fractionation accuracy test**: Same structure as liquid-vapor test but calling `AlphaEqIceVapor()`.

### Temperature Sampling Strategy

Select approximately 10 temperature points distributed across the 233-313K range:
- Lower bound: 233K (relevant for cold atmospheric conditions)
- Upper bound: 313K (relevant for warm surface/atmospheric conditions)
- Include key transition points (e.g., 273.15K freezing point)
- Distribute remaining points to achieve good coverage

### Reference Value Generation

For each test temperature and isotopologue:
1. Manually calculate expected fractionation factors using the polynomial formulas from the source papers
2. Document the calculation method and source reference in test comments
3. Hard-code these reference values as constants in the test file

This approach provides independent validation - the test values come from direct implementation of the published formulas, not from running the code being tested.

### Integration with EAMxx Build System

Update `components/eamxx/src/physics/water_tracers/tests/CMakeLists.txt` to:
1. Add the new test file to the compilation list
2. Register the test with CTest so it runs as part of `./scripts/test-all-eamxx`

Follow existing patterns from other EAMxx unit tests in the directory.

### Verification Plan

After implementation:
1. Build the standalone EAMxx test suite in debug mode
2. Run the test using `./scripts/test-all-eamxx --preserve-env -m docker -t dbg`
3. Verify all test cases pass
4. Review test output to confirm expected coverage (all isotopologues, all temperatures)

## Notes

- This work validates the mathematical correctness of the fractionation functions in isolation
- Kinetic fractionation testing is explicitly deferred to a separate effort
- The test focuses on the computational kernel functions, not their integration into broader physics
- Success depends on accurate hand-calculation of reference values from source papers
- The test will become part of the ongoing regression suite for the water isotope implementation
