include(ScreamUtils)

# tracer_ratio and bulk_ocean_flux, incl. FPE safety (run in the fpe build too).
CreateUnitTest(water_tracers_functions
  SOURCES water_tracers_functions_tests.cpp
  LIBS water_tracers
  LABELS water_tracers physics
)

# Generic tracer ocean fluxes, and the invalid-ocean rule.
CreateUnitTest(water_tracers_ocean_flux
  SOURCES water_tracers_ocean_flux_tests.cpp
  LIBS water_tracers
  LABELS water_tracers physics
)

# water_tracers_file parsing and validation. The test writes its YAML fixtures
# into the working directory.
CreateUnitTest(water_tracers_registry
  SOURCES water_tracers_registry_tests.cpp
  LIBS water_tracers
  LABELS water_tracers physics
)

# Isotope-specific unit tests. Isotope math is always compiled into
# water_tracers, but these tests only make sense (and are only registered)
# when the build also accepts isotope entries.
if (EAMXX_ENABLE_WATER_ISOTOPES)
  # Equilibrium water-isotope fractionation factor unit tests.
  CreateUnitTest(water_isotopes_fractionation
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../../water_isotopes/tests/eamxx_water_isotopes_fractionation_tests.cpp
    LIBS water_tracers
    LABELS water_tracers water_isotopes physics
  )

  # Runtime formulation selection unit tests.
  CreateUnitTest(water_isotopes_runtime_options
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../../water_isotopes/tests/water_isotopes_runtime_options_tests.cpp
    LIBS water_tracers
    LABELS water_tracers water_isotopes physics
  )

  # Isotope ocean evaporation flux: equilibrium invariants, composed from
  # alpha_liquid_vapor and water_tracers::ocean_flux.
  CreateUnitTest(water_isotopes_ocean_flux
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../../water_isotopes/tests/water_isotopes_ocean_flux_tests.cpp
    LIBS water_tracers
    LABELS water_tracers water_isotopes physics
  )

  # Kinetic fractionation of ocean evaporation (alpha_kinetic), its use in
  # water_tracers::ocean_flux, and the bottom-level height it reads.
  CreateUnitTest(water_isotopes_kinetic
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../../water_isotopes/tests/water_isotopes_kinetic_tests.cpp
    LIBS water_tracers
    LABELS water_tracers water_isotopes physics
  )

  # Prints alpha_kinetic for a CSV of cases, for comparison against the
  # Fortran wiso_kmol. Not a ctest test.
  add_executable(water_isotopes_kinetic_driver
    ${CMAKE_CURRENT_SOURCE_DIR}/../../water_isotopes/tests/water_isotopes_kinetic_driver.cpp)
  target_link_libraries(water_isotopes_kinetic_driver water_tracers)
endif()
