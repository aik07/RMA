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

* Run the build script to build PEBBL and RMA
```
cd RMA
sh scripts/build.sh
```

<!--
* Build PEBBL
```
mkdir RMA/external/pebbl/build
cd RMA/external/pebbl/build
cmake -Denable_mpi=ON -Denable_examples=OFF ..
make
```

* Build RMA
```
cd ../../../  # go back to the RMA root directory
mkdir build
cd build
cmake ..
make
```
-->


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

Please read the user guide about how to use parameters for the RMA solver.


## Class Diagram

<p align="center">

<img src="https://github.com/aik7/RMA/blob/devel/figures/RMA_class_org.png" width="300">

## Source files at src directory
```
├── argRMA.cpp       : a file contains RMA argument class
├── argRMA.h
├── baseRMA.cpp      : a file contains a base class for a RMA driver class
├── baseRMA.h
├── dataRMA.cpp      : a file contains RMA data classes
├── dataRMA.h        
├── driver.cpp       : a driver file
├── driverRMA.cpp    : a file contains RMA driver class
├── driverRMA.h      
├── greedyRMA.cpp    : a file contains greedy RMA class
├── greedyRMA.h      
├── parRMA.cpp       : a file contains parallel RMA class
├── parRMA.h   
├── serRMA.cpp       : a file contains serial RMA class
├── serRMA.h
└── Time.h           : a header file for time
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
