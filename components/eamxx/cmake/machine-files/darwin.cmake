include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

# Darwin-specific settings for LANL HPC
set (EKAT_TEST_LAUNCHER_MANAGE_RESOURCES True CACHE BOOL "")

# Compiler flags for gfortran >= 10
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
    AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
  if (CMAKE_Fortran_FLAGS)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch"
        CACHE STRING "Fortran compiler flags" FORCE)
  else()
    set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"
        CACHE STRING "Fortran compiler flags" FORCE)
  endif()
endif()

# MPI launcher settings for darwin
set(EKAT_MPIRUN_EXE "mpirun" CACHE STRING "")
set(EKAT_MPI_NP_FLAG "-n" CACHE STRING "")
set(EKAT_MPI_EXTRA_ARGS "" CACHE STRING "Extra args for mpirun")

# Disable use of deprecated Kokkos 4 APIs
option(Kokkos_ENABLE_DEPRECATED_CODE_4 "" OFF)
