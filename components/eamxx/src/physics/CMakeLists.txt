add_library(eamxx_physics INTERFACE)

# Water tracer options. They are defined here, before p3 and shoc, because the
# two of them export their bulk water transfers (the input of the
# water_tracer_transfers process) only in builds with water tracers enabled:
# p3 and shoc check EAMXX_ENABLE_WATER_TRACERS when they are configured.
#
#   EAMXX_ENABLE_WATER_TRACERS  - build the water tracer processes (and have
#                                 SHOC and P3 export their water transfers)
#   EAMXX_ENABLE_WATER_ISOTOPES - allow `type: isotope` entries in the
#                                 tracers file (isotope math is header-only
#                                 and always compiled into water_tracers;
#                                 this flag only gates acceptance and the
#                                 isotope-specific unit tests)
#
# Isotopes are implemented on top of the water tracer infrastructure, so they
# can ONLY be enabled when water tracers are also enabled. If isotopes are
# requested without tracers we auto-enable tracers and warn, rather than build
# an inconsistent configuration.
option(EAMXX_ENABLE_WATER_TRACERS  "Whether to build the EAMxx water tracer process"  OFF)
option(EAMXX_ENABLE_WATER_ISOTOPES "Whether to accept isotope entries in the EAMxx water tracer process" OFF)

# Enforce the dependency: isotopes require tracers.
if (EAMXX_ENABLE_WATER_ISOTOPES AND NOT EAMXX_ENABLE_WATER_TRACERS)
  message(STATUS "WARNING: EAMXX_ENABLE_WATER_ISOTOPES=ON requires EAMXX_ENABLE_WATER_TRACERS; "
                 "auto-enabling water tracers.")
  set(EAMXX_ENABLE_WATER_TRACERS ON CACHE BOOL "Whether to build the EAMxx water tracer process" FORCE)
endif()

# Add all physics parametrizations subfolders
add_subdirectory(p3)
if (SCREAM_DOUBLE_PRECISION)
  add_subdirectory(rrtmgp)
  add_subdirectory(cosp)
  add_subdirectory(tms)
  add_subdirectory(iop_forcing)
  add_subdirectory(zm)
else()
  message(STATUS "WARNING: RRTMGP, COSP, TMS, ZM, and IOPForcing only supported for double precision builds; skipping")
endif()
add_subdirectory(shoc)
add_subdirectory(cld_fraction)
add_subdirectory(spa)
add_subdirectory(spc)
add_subdirectory(nudging)
if (SCREAM_ENABLE_MAM)
  add_subdirectory(mam)
endif()
add_subdirectory(gw)
add_subdirectory(specialized_tracers)
