# This file contains the options needed to both run the pull request testing
# for Trilinos for the caraway pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of modules must be loaded and path must be augmented.

# Usage: cmake -C PullRequestLinuxCarawayTestingSettings.cmake

set (CMAKE_CXX_STANDARD "14" CACHE STRING "Set C++ standard to C++14")
set (CMAKE_C_COMPILER "mpicc" CACHE STRING "Set by default for caraway PR testing")
set (CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "Set by default for caraway PR testing")
set (CMAKE_FORTRAN_COMPILER "/home/projects/x86-64/gcc/8.2.0/bin/gfortran" CACHE STRING "Set by default for caraway PR testing")
set (CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Set by default for caraway PR testing")

#set (CMAKE_CXX_FLAGS "--gcc-toolchain=/home/projects/x86-64/gcc/8.2.0" CACHE STRING "Set by default for caraway PR testing")

set (BUILD_SHARED_LIBS "OFF" CACHE STRING "Set by default for caraway PR testing")
set (Trilinos_EXTRA_LINK_FLAGS "-L/home/projects/x86-64-rocm/netlib-lapack/3.8.0/lib64 -lgfortran -lm" CACHE STRING "Set by default for caraway PR testing")
set (Trilinos_ENABLE_TESTS ON CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_ALL_PACKAGES OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ENABLE_EXPLICIT_INSTANTIATION ON CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ASSERT_MISSING_PACKAGES OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ALLOW_NO_PACKAGES OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ENABLE_OpenMP OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ENABLE_Amesos2 ON CACHE BOOL "Set by default for caraway PR testing")
set (Amesos2_ENABLE_SuperLU OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Amesos2_ENABLE_KLU2 ON CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Belos ON CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Ifpack2 ON CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Kokkos ON CACHE BOOL "Set by default for caraway PR testing")
set (Kokkos_ARCH_VEGA908 ON CACHE BOOL "Set by default for caraway PR testing")
set (Kokkos_ENABLE_CUDA OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_HIP ON CACHE BOOL "Set by default for caraway PR testing")
set (Kokkos_ENABLE_OPENMP OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Kokkos_ENABLE_DEPRECATED_CODE OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Trilinos_ENABLE_KokkosKernels ON CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_MueLu ON CACHE BOOL "Set by default for caraway PR testing")
set (MueLu_ENABLE_Kokkos_Refactor ON CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Tpetra ON CACHE BOOL "Set by default for caraway PR testing")
set (Tpetra_ENABLE_CUDA OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Tpetra_INST_HIP ON CACHE BOOL "Set by default for caraway PR testing")
set (Tpetra_INST_SERIAL OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Tpetra_INST_OPENMP OFF CACHE BOOL "Set by default for CUDA PR testing")
set (Tpetra_INST_DOUBLE ON CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Gtest ON CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Teuchos ON CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Xpetra ON CACHE BOOL "Set by default for caraway PR testing")
set (Xpetra_ENABLE_Kokkos_Refactor ON CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Zoltan2 ON CACHE BOOL "Set by default for caraway PR testing")
set (TPL_ENABLE_BLAS ON CACHE BOOL "Set by default for caraway PR testing")
set (TPL_BLAS_LIBRARIES "/home/projects/x86-64-rocm/netlib-lapack/3.8.0/lib64/libblas.a" CACHE STRING "Set by default for caraway PR testing")
set (TPL_ENABLE_LAPACK ON CACHE BOOL "Set by default for caraway PR testing")
set (TPL_LAPACK_LIBRARIES "/home/projects/x86-64-rocm/netlib-lapack/3.8.0/lib64/liblapack.a;gfortran;m" CACHE STRING "Set by default for caraway PR testing")
set (TPL_ENABLE_MPI ON CACHE BOOL "Set by default for caraway PR testing")
set (MPI_USE_COMPILER_WRAPPERS ON CACHE BOOL "Set by default for caraway PR testing")
set (MPI_BASE_DIR "/home/projects/x86-64-rocm/openmpi/4.0.5/gcc/8.2.0" CACHE STRING "Set by default for caraway PR testing")
set (MPI_EXEC "mpirun" CACHE STRING "Set by default for caraway PR testing")
set (MPI_EXEC_NUMPROCS_FLAG "-np" CACHE STRING "Set by default for caraway PR testing")
set (KokkosCore_UnitTest_CudaTimingBased_MPI_1_DISABLE ON CACHE BOOL "Set by default for caraway PR testing")

#New stuff forthe full build
set (Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES ON CACHE BOOL "Set by default for PR testing")
set (Trilinos_ENABLE_SEACAS OFF CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_STK OFF CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_ShyLU OFF CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_ShyLU_DD OFF CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_ShyLU_Node OFF CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_MiniTensor OFF CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Stokhos OFF CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Intrepid2 OFF CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Phalanx OFF CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_ROL OFF CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Panzer OFF CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_Moertel OFF CACHE BOOL "Set by default for caraway PR testing")
set (Trilinos_ENABLE_TrilinosCouplings OFF CACHE BOOL "Set by default for caraway PR testing")


set (TPL_ENABLE_LAPACK ON CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_Boost ON CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_BoostLib ON CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_ParMETIS ON CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_Zlib ON CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_HDF5 ON CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_Netcdf ON CACHE BOOL "Set by default for PR testing")
set (Trilinos_TRACE_ADD_TEST ON CACHE BOOL "Set by default for PR testing")


#trouble building these 2 on caraway
set (TPL_ENABLE_SuperLU OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_Scotch OFF CACHE BOOL "Set by default for PR testing")
#set (TPL_ENABLE_Matio OFF CACHE BOOL "Set by default for PR testing")
#set (TPL_ENABLE_X11 OFF CACHE BOOL "Set by default for PR testing")


SET(Trilinos_EXTRA_LINK_FLAGS "-lgomp -lgfortran -lldl -ldl" CACHE STRING "Set by default for PR testing")

SET(Trilinos_ENABLE_PyTrilinos OFF CACHE BOOL "Set by default for PR testing")

SET(Boost_INCLUDE_DIRS "$ENV{BOOST_INC_DIR}" CACHE PATH "Set by default for PR testing")
SET(Boost_LIBRARY_DIRS "$ENV{BOOST_LIB_DIR}" CACHE PATH "Set by default for PR testing")

SET(BoostLib_INCLUDE_DIRS "$ENV{BOOST_INC_DIR}" CACHE PATH "Set by default for PR testing")
SET(BoostLib_LIBRARY_DIRS "$ENV{BOOST_LIB_DIR}" CACHE PATH "Set by default for PR testing")

SET(ParMETIS_INCLUDE_DIRS "$ENV{PARMETIS_INC_DIR}" CACHE PATH "Set by default for PR testing")
SET(ParMETIS_LIBRARY_DIRS "$ENV{PARMETIS_LIB_DIR}" CACHE PATH "Set by default for PR testing")

SET(Zlib_INCLUDE_DIRS "$ENV{ZLIB_INC_DIR}" CACHE PATH "Set by default for PR testing")
SET(Zlib_LIBRARY_DIRS "$ENV{ZLIB_LIB_DIR}" CACHE PATH "Set by default for PR testing")

SET(HDF5_INCLUDE_DIRS "$ENV{HDF5_INC_DIR}" CACHE PATH "Set by default for PR testing")
SET(HDF5_LIBRARY_DIRS "$ENV{HDF5_LIB_DIR}" CACHE PATH "Set by default for PR testing")

SET(Netcdf_INCLUDE_DIRS "$ENV{NETCDF_INC_DIR}" CACHE PATH "Set by default for PR testing")
SET(Netcdf_LIBRARY_DIRS "$ENV{NETCDF_LIB_DIR}" CACHE PATH "Set by default for PR testing")

#SET(SuperLU_INCLUDE_DIRS "$ENV{SUPERLU_INC_DIR}" CACHE PATH "Set by default for PR testing")
#SET(SuperLU_LIBRARY_DIRS "$ENV{SUPERLU_LIB_DIR}" CACHE PATH "Set by default for PR testing")

#set (TPL_Scotch_INCLUDE_DIRS "$ENV{SCOTCH_INC_DIR}" CACHE PATH "Set by default for PR testing")
#set (Scotch_LIBRARY_DIRS "$ENV{SCOTCH_LIB_DIR}" CACHE PATH "Set by default for PR testing")

#include("${CMAKE_CURRENT_LIST_DIR}/PullRequestLinuxCommonTestingSettings.cmake")
