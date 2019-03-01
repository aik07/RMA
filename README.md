# RMA

RMA is a solver to find an axis-parallel box containing the maximum net
weight of positivly minus negativly covered observations or vice versa
when each observation has a positive or negative weight. 

## Software Requirement:
* [PEBBL](https://software.sandia.gov/trac/acro/wiki/Example/Building/acro-pebbl)
* C++ compiler
* [Open MPI](https://www.open-mpi.org/)

## The description and user guide of RMA
* [Presentation](https://github.com/aik07/RMA/blob/master/RMA_slides.pdf)
* [User Guide](https://github.com/aik07/RMA/blob/master/RMA_user_guide.pdf)

## Reference

```
@phdthesis{AiThesis,
  author       = {Ai Kagawa}, 
  title        = {The Rectangular Maximum Agreement Problem: Applications and 
                  Parallel Solution},
  school       = {Rutgers University},
  year         = 2018
}
```

## How to download and build RMA

* Clone or download this RMA repository
```
git clone https://github.com/aik07/RMA.git
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

Please read the user guide about how to use parameters for the RMA solver.
