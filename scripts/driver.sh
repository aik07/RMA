#!/bin/bash

RMA_DIR=/home/kagawa/Projects/thesis/rma/RMA

mpirun -np 8 $(RMA_DIR)/rma $(RMA_DIR)/data/servo.data
