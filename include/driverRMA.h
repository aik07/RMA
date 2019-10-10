/**********************************************************
* File name:   driverRMA.h
* Author:      Ai Kagawa
* Description: a header file for RMA driver class
***********************************************************/

#ifndef RMA_h
#define RMA_h

#include <iostream>

#include <pebbl_config.h>
#include <pebbl/utilib/CommonIO.h>

#include "argRMA.h"
#include "dataRMA.h"
#include "serRMA.h"
#include "greedyRMA.h"

#ifdef ACRO_HAVE_MPI
  #include <pebbl/pbb/parBranching.h>
  #include "parRMA.h"
  #define outstream ucout
  #define IO(action) if (uMPI::iDoIO) { CommonIO: end_tagging(); action; }
#else
  typedef void parRMA;
  #define outstream cout;
  #define IO(action) action;
#endif


namespace rma {

//static double inf = numeric_limits<double>::infinity();
//static int intInf = numeric_limits<int>::max();

using namespace utilib;
using namespace arg;
using namespace data;
using namespace pebblRMA;
using namespace greedyRMA;


class DriverRMA {

public:

  DriverRMA(ArgRMA* args, Data* data);

  ~DriverRMA() {
    #ifdef ACRO_HAVE_MPI
      if (parallel) { CommonIO::end(); uMPI::done(); }
    #endif // ACRO_HAVE_MPI
  }

  void solveRMA();

//private:
  bool          parallel;

  ArgRMA*       args;
  Data*         data;
  GreedyRMA*    grma;
  RMA*          rma ;
  parRMA*       prma;

	Time          tc;
	double        wallTime;
  double        cpuTime;

};

} // end namespace rma

#endif // RMA_h
