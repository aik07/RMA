#ifndef RMA_h
#define RMA_h

#include <iostream>

#include <pebbl_config.h>
#include <ParameterList.h>

#include "data.h"
#include "greeduRMA.h"

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


namespace pebblRMA {

using namepsace utilib;
using namespace greedyRMA;

class DriverRMA {
  
public:
  DriverRMA(): rma(NULL), pram(NULL), parallel(false) {}
  ~DriverRMA() {
    #ifdef ACRO_HAVE_MPI
      if (parallel) { CommonIO::end(); uMPI::done(); }
    #endif // ACRO_HAVE_MPI
  }
  
private:
  Data*      data;
  RMA*       rma;
  parRMA*    prma;
  greedyRMA* grma;

};

} // end namespace pebblRMA

#endif // RMA_h
