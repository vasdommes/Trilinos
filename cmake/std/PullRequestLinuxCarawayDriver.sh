#!/bin/bash -l

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

set -x

sbatch -N1 ${WORKSPACE}/Trilinos/cmake/std/PullRequestLinuxDriver.sh
