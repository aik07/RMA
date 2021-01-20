#!/bin/bash

RMA_DIR=$(pwd)
#RMA_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo ${RMA_DIR}


# process the command line argument
while getopts b: flag
do
    case "${flag}" in
        b) build_type=${OPTARG};;
    esac
done

# set CMAKE_BUILD_TYPE
if [ "${build_type}" = "debug" ]; then
  export DEBUG_OPTION="-DCMAKE_BUILD_TYPE=Debug"
else
  export DEBUG_OPTION="-DCMAKE_BUILD_TYPE=Release"
fi


# Build PEBBL
export PEBBL_BUILD_DIR=${RMA_DIR}"/external/pebbl/build"
# create a file if it does not exist
mkdir -p ${PEBBL_BUILD_DIR}
cd ${PEBBL_BUILD_DIR}
cmake -Denable_mpi=ON -Denable_examples=OFF ${DEBUG_OPTION} ..
make


# Build RMA
export RMA_BUILD_DIR=${RMA_DIR}"/build"
# create a file if it does not exist
mkdir -p  ${RMA_BUILD_DIR}
cd ${RMA_BUILD_DIR}
cmake ${DEBUG_OPTION} ..
make
