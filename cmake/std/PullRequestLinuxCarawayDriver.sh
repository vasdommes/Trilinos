#!/bin/bash -l

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

set -x
module load cmake/3.19.3
module load rocm/4.3.0

srun -N1 -p MI100 ${WORKSPACE}/Trilinos/cmake/std/PullRequestLinuxDriver.sh
