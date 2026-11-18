# Auxiliary tracers for EAMxx.
#
# This folder groups the "auxiliary tracer" physics processes: water_tracers
# and water_tracer_transfers, which accept generic and isotope tracer entries
# in a single water_tracers_file. The build is governed by the CMake options
# EAMXX_ENABLE_WATER_TRACERS and EAMXX_ENABLE_WATER_ISOTOPES, which are defined
# (and documented) in ../CMakeLists.txt, ahead of p3 and shoc.

if (EAMXX_ENABLE_WATER_TRACERS)
  add_subdirectory(water_tracers)
endif()
