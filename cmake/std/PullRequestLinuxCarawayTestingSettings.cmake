# This file contains the options needed to both run the pull request testing
# for Trilinos for the CUDA 10.2.2 pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of modules must be loaded and path must be augmented.
# (See the sems/PullRequestCuda10.2.2TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxCarawayTestingSettings.cmake

set (CMAKE_CXX_STANDARD "14" CACHE STRING "Set C++ standard to C++14")

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

# Options necessary for CUDA build
set (TPL_ENABLE_MPI ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CUDA ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CUDA_UVM ON CACHE BOOL "Set by default for CUDA PR testing")
set (KOKKOS_ARCH Power9,Volta70 CACHE STRING "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Sacado_ENABLE_HIERARCHICAL_DFAD ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CXX11_DISPATCH_LAMBDA ON CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Set by default for CUDA PR testing")
set (Phalanx_KOKKOS_DEVICE_TYPE CUDA CACHE STRING "Set by default for CUDA PR testing")
set (MPI_EXEC_POST_NUMPROCS_FLAGS "-map-by;socket:PE=4" CACHE STRING "Set by default for CUDA PR testing")

set (Trilinos_ENABLE_DEBUG OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ENABLE_DEBUG_SYMBOLS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (BUILD_SHARED_LIBS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Tpetra_INST_SERIAL ON CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ENABLE_SECONDARY_TESTED_CODE OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_Debug_Bounds_Check ON CACHE BOOL "Set by default for CUDA PR testing")
set (KOKKOS_ENABLE_DEBUG ON CACHE BOOL "Set by default for CUDA PR testing")

# TPL settings specific to CUDA build
#set (TPL_ENABLE_CUDA ON CACHE BOOL "Set by default for CUDA PR testing")
#set (TPL_ENABLE_CUSPARSE ON CACHE BOOL "Set by default for CUDA PR testing")
#set (TPL_ENABLE_Pthread OFF CACHE BOOL "Set by default for CUDA PR testing")

set (TPL_ENABLE_BinUtils OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_BLAS ON CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_BLAS_LIBRARIES "-L$ENV{BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp;-lm" CACHE STRING "Set by default for CUDA PR testing")
set (TPL_ENABLE_Boost ON CACHE BOOL "Set by default for CUDA PR testing")
set (Boost_INCLUDE_DIRS "$ENV{BOOST_ROOT}/include" CACHE FILEPATH "Set by default for CUDA PR testing")
set (TPL_Boost_LIBRARIES "$ENV{BOOST_ROOT}/lib/libboost_program_options.a;$ENV{BOOST_ROOT}/lib/libboost_system.a" CACHE FILEPATH  "Set by default for CUDA PR testing")
set (TPL_ENABLE_BoostLib ON CACHE BOOL "Set by default for CUDA PR testing")
set (BoostLib_INCLUDE_DIRS "$ENV{BOOST_ROOT}/include" CACHE FILEPATH "Set by default for CUDA PR testing")
set (TPL_BoostLib_LIBRARIES "$ENV{BOOST_ROOT}/lib/libboost_program_options.a;$ENV{BOOST_ROOT}/lib/libboost_system.a" CACHE FILEPATH "Set by default for CUDA PR testing")
set (TPL_ENABLE_METIS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_ParMETIS OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_Scotch OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_HWLOC OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_LAPACK ON CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_LAPACK_LIBRARIES "-L$ENV{LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp" CACHE STRING "Set by default for CUDA PR testing")
set (TPL_ENABLE_Zlib OFF CACHE BOOL "Set by default for CUDA PR testing")
set (TPL_ENABLE_CGNS OFF CACHE BOOL "Set by default for CUDA PR testing")



#set (TPL_ENABLE_HDF5 ON CACHE BOOL "Set by default for CUDA PR testing")
#set (HDF5_INCLUDE_DIRS "$ENV{HDF5_ROOT}/include" CACHE FILEPATH "Set by default for CUDA PR testing")
#set (TPL_HDF5_LIBRARIES "-L$ENV{HDF5_ROOT}/lib;$ENV{HDF5_ROOT}/lib/libhdf5_hl.a;$ENV{HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl" CACHE FILEPATH "Set by default for CUDA PR testing")
#set (TPL_ENABLE_Netcdf ON CACHE BOOL "Set by default for CUDA PR testing")
#set (Netcdf_INCLUDE_DIRS "$ENV{NETCDF_ROOT}/include" CACHE FILEPATH "Set by default for CUDA PR testing")
#IF ("$ENV{PNETCDF_ROOT}" STREQUAL "")
#  SET(PNETCDF_ROOT "$ENV{NETCDF_ROOT}")
#ELSE()
#  SET(PNETCDF_ROOT "$ENV{PNETCDF_ROOT}")
#ENDIF()
#set (TPL_Netcdf_LIBRARIES "-L$ENV{NETCDF_ROOT}/lib64;$ENV{NETCDF_ROOT}/lib64/libnetcdf.so;${PNETCDF_ROOT}/lib/libpnetcdf.so;${TPL_HDF5_LIBRARIES}" CACHE STRING "Set by default for CUDA PR testing")
# SuperLU and SuperLUDist is available on ride and could be enabled for the CUDA PR build
#set (TPL_ENABLE_SuperLU OFF CACHE BOOL "Set by default for CUDA PR testing")
#set (TPL_ENABLE_SuperLUDist OFF CACHE BOOL "Set by default for CUDA PR testing")
#set (TPL_ENABLE_BoostLib OFF CACHE BOOL "Set by default for CUDA PR testing")
#set (TPL_ENABLE_Matio OFF CACHE BOOL "Set by default for CUDA PR testing")
#set (TPL_DLlib_LIBRARIES "-ldl" CACHE FILEPATH "Set by default for CUDA PR testing")

# The compile times for two Panzer files went up to over 6 hours. This
# turns off one feature that allows these in about 24 minutes. Please remove
# when issue #7532 is resolved.
set (Sacado_NEW_FAD_DESIGN_IS_DEFAULT OFF CACHE BOOL "Temporary fix for issue #7532" )

# Disable some packages that can't be tested with this PR build
set (Trilinos_ENABLE_ShyLU_NodeTacho OFF CACHE BOOL
  "Can't test Tacho with CUDA without RDC" FORCE)

include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")
