# RMA

[![Build Status](https://travis-ci.com/aik7/RMA.svg?branch=devel)](https://travis-ci.com/aik7/RMA)

RMA is a solver to find an axis-parallel box containing the maximum net
weight of positivly minus negativly covered observations or vice versa
when each observation has a positive or negative weight.

## Software Requirement:

### You need to install
* CMake (version >= 3.0)
* C++ compiler (g++)
* OpenMPI 2.1.1 (openmpi-bin, libopenmpi-dev)

### The following packages are installed by running scripts/build.sh as described below
* [PEBBL](https://github.com/PEBBL/pebbl)

### Note
* The build was tested on Ubuntu 18.04 (Bionic) as shown in [our TravisCI file](https://github.com/aik7/RMA/blob/devel/.travis.yml)


## The description and user guide of RMA
* [Presentation](https://github.com/aik7/RMA/blob/master/RMA_slides.pdf)
* [User Guide](https://github.com/aik7/RMA/blob/master/RMA_user_guide.pdf)

## How to clone and build PEBBL and RMA

* Recursive clone the RMA repository (PEBBL is a submodule in this repository)
```
git clone --recursive https://github.com/aik7/RMA.git
```

* Build RMA with PEBBL (Note: assuming that all commands are executed in the RMA main directory)
```
sh scripts/build.sh
```

* If you want to build RMA and PEBBL with the debug mode

```
sh scripts/build.sh -d true
```

## Example commands to run RMA:

### Serial implementation
```
./build/rma <data_filename>
```

### Parallel implementation
```
mpirun -np 4 ./build/rma <data_filename>
```

### Example script
```
sh scripts/driver.sh  # at the RMA root directory
```

## Parameters

| parameters      |      description                               | data type | range         | default value  |
|-----------------|:-----------------------------------------------|:---------:|--------------:|---------------:|
| isPebblRMA   | whether or not to solve RMA using PEBBL            | bool   | true or false | true              |
| isGreedyRMA   | whether or not to solve RMA using the greedy method          | bool   | true or false | true           |
| isInitGuess   | whether or not to set the PEBBL's initial solution as the greedy RMA solution          | bool   | true or false | true           |
| fracCachedCutPts   | the threshold of cached cutpoints caching which is the fraction of the total cutpoints   | double   | [0, infinity) | 0           |
| isMinMaxGreedy   | whether or not to solve the greedy RMA by min-then-max order; if false, it will solve the greedy RMA by max-then-min order   | bool   | true or false | true           |
| isPrintBBdetails   | whether or not to print out the PEBBL's branch and bound information  | bool   | true or false | false           |
| nonUniformWt | if you set a non-uniform weight file for a given dataset, the code solves the RMA problem using the weight | string   |  | ""        |
| delta | it is a factor of the epsilon integerization method  | double   |  | 0.005       |



## Class Diagram

<p align="center">

<img src="https://github.com/aik7/RMA/blob/devel/figures/RMA_class_org.png" width="300">

* A solid arrow indicates an inheritance relationship
* A dashed arrow indicates a composition relationship

## Source files at src directory
```
├── argRMA.cpp       : a file contains RMA argument class
├── argRMA.h
├── baseRMA.cpp      : a file contains a base class for a SolveRMA class
├── baseRMA.h
├── dataRMA.cpp      : a file contains RMA data classes
├── dataRMA.h
├── driver.cpp       : a driver file
├── greedyRMA.cpp    : a file contains greedy RMA class
├── greedyRMA.h
├── parRMA.cpp       : a file contains parallel RMA class
├── parRMA.h
├── serRMA.cpp       : a file contains serial RMA class
├── serRMA.h
├── Time.h           : a header file for time
├── solveRMA.cpp     : a file contains SolveRMA class
├── solveRMA.h
├── utility.cpp      : a file contains RMA utility functions
└── utility.h
```

## Reference

```
@phdthesis{AiThesis,
  author       = {Ai Kagawa},
  title        = {The Rectangular Maximum Agreement Problem: Applications and Parallel Solution},
  school       = {Rutgers University},
  year         = 2018
}
```
