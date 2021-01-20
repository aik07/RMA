#!/bin/bash

RMA_DIR=$(pwd)

${RMA_DIR}/build/rma ${RMA_DIR}/data/cleveland.dat

mpirun -np 4 ${RMA_DIR}/build/rma ${RMA_DIR}/data/cleveland.dat
