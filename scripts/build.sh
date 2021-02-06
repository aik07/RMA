#!/bin/bash

#RMA_DIR=$(pwd)
RMA_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo ${RMA_DIR}


# process the command line argument
while getopts d:p: flag
do
    case "${flag}" in
        d) debug_mode=${OPTARG};;
        p) pebbl_build=${OPTARG};;
    esac
done


# set CMAKE_BUILD_TYPE
if [ "${debug_mode}" = "true" ]; then
  export DEBUG_OPTION="-DCMAKE_BUILD_TYPE=Debug"
else
  export DEBUG_OPTION="-DCMAKE_BUILD_TYPE=Release"
fi

echo ${DEBUG_OPTION}


# set CMAKE_BUILD_TYPE
if [ "${pebbl_build}" = "false" ]; then
  echo "NOT BUILDING PEBBEL"
else
  # Build PEBBL
  export PEBBL_BUILD_DIR=${RMA_DIR}"/external/pebbl/build"
  # create a file if it does not exist
  mkdir -p ${PEBBL_BUILD_DIR}
  cd ${PEBBL_BUILD_DIR}
  cmake -Denable_mpi=ON -Denable_examples=OFF ${DEBUG_OPTION} ..
  make -j8
fi


# Build RMA
export RMA_BUILD_DIR=${RMA_DIR}"/build"
# create a file if it does not exist
mkdir -p  ${RMA_BUILD_DIR}
cd ${RMA_BUILD_DIR}
cmake ${DEBUG_OPTION} ..
make -j8
