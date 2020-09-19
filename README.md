# RMA

[![Build Status](https://travis-ci.com/aik7/RMA.svg?branch=travis-ci)](https://travis-ci.com/aik7/RMA)

RMA is a solver to find an axis-parallel box containing the maximum net
weight of positivly minus negativly covered observations or vice versa
when each observation has a positive or negative weight.

## Software Requirement:
* [PEBBL](https://github.com/PEBBL/pebbl)
* C++ compiler
* MPI

## The description and user guide of RMA
* [Presentation](https://github.com/aik7/RMA/blob/master/RMA_slides.pdf)
* [User Guide](https://github.com/aik7/RMA/blob/master/RMA_user_guide.pdf)

## How to download and build RMA

* Clone or download this RMA repository
```
git clone https://github.com/aik7/RMA.git
```

* Run the following command for compiling and building applications
```
cd RMA
make
```

## Example commands to run RMA:

### Serial implementation
```
./rma <data_filename>
```

### Parallel implementation
```
mpirun -np 4 ./rma <data_filename>
```

### Example script
```
sh script/driver.sh
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
