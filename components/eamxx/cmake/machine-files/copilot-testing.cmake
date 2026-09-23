include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

# No resource manager in CI/container environments
set (EKAT_TEST_LAUNCHER_MANAGE_RESOURCES True CACHE BOOL "")

# -fallow-argument-mismatch is needed for gfortran >= 10 to compile legacy
# Fortran code (e.g. HOMME's bndry_mod.F90, which calls MPI_Isend/Irecv with
# inconsistent argument types across call sites). CMAKE_Fortran_COMPILER_ID
# can't be used to guard this: this file is preloaded via ctest's `-C`
# option, which runs before project()/enable_language(), so that variable is
# still unset here. Detect gfortran directly via `mpifort --version` instead.
execute_process(
  COMMAND mpifort --version
  OUTPUT_VARIABLE _mpifort_version_output
  ERROR_VARIABLE _mpifort_version_output
)
if (_mpifort_version_output MATCHES "GNU Fortran")
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"
      CACHE STRING "Fortran compiler flags" FORCE)
endif()

# Input data directory (set by setup-copilot-env.sh or agent)
if (DEFINED ENV{SCREAM_INPUT_ROOT})
  set(SCREAM_INPUT_ROOT "$ENV{SCREAM_INPUT_ROOT}" CACHE PATH "")
endif()

# BLAS/LAPACK - prefer env vars, fall back to system paths
if (DEFINED ENV{BLAS_ROOT})
  # Try lib first, then lib64 (common on HPC/module installs)
  set(_blas_lib_dir "lib")
  if (NOT EXISTS "$ENV{BLAS_ROOT}/${_blas_lib_dir}/libblas.so"
      AND EXISTS "$ENV{BLAS_ROOT}/lib64/libblas.so")
    set(_blas_lib_dir "lib64")
  endif()
  set(BLAS_LIBRARIES "$ENV{BLAS_ROOT}/${_blas_lib_dir}/libblas.so" CACHE STRING "")
  set(LAPACK_LIBRARIES "$ENV{BLAS_ROOT}/${_blas_lib_dir}/liblapack.so" CACHE STRING "")
elseif (EXISTS "/usr/lib/x86_64-linux-gnu/libblas.so")
  set(BLAS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libblas.so" CACHE STRING "")
  set(LAPACK_LIBRARIES "/usr/lib/x86_64-linux-gnu/liblapack.so" CACHE STRING "")
endif()

# Help Scorpio find PnetCDF when pnetcdf-config is not available (e.g., apt install)
if (NOT DEFINED PnetCDF_C_PATH AND EXISTS "/usr/include/pnetcdf.h")
  set(PnetCDF_C_PATH "/usr" CACHE PATH "")
endif()

# MPI launcher settings
set(EKAT_MPIRUN_EXE "mpirun" CACHE STRING "")
set(EKAT_MPI_NP_FLAG "-n" CACHE STRING "")

# Open MPI refuses to run as root and warns about oversubscription unless
# told otherwise; MPICH's mpiexec has neither restriction and errors out on
# these unrecognized flags. Detect which one is present instead of assuming.
execute_process(
  COMMAND mpirun --version
  OUTPUT_VARIABLE _mpirun_version_output
  ERROR_VARIABLE _mpirun_version_output
)
if (_mpirun_version_output MATCHES "Open MPI")
  set(EKAT_MPI_EXTRA_ARGS "--allow-run-as-root --oversubscribe" CACHE STRING "Extra args for mpirun")
else()
  set(EKAT_MPI_EXTRA_ARGS "" CACHE STRING "Extra args for mpirun")
endif()

# Disable use of deprecated Kokkos 4 APIs
option(Kokkos_ENABLE_DEPRECATED_CODE_4 "" OFF)
