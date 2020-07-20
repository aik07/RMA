/**********************************************************
* File name:   driverRMA.h
* Author:      Ai Kagawa
* Description: a header file for RMA driver class
***********************************************************/

#ifndef RMA_h
#define RMA_h

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>

#include <pebbl_config.h>
#include <pebbl/utilib/CommonIO.h>

#include "baseRMA.h"
#include "dataRMA.h"
#include "serRMA.h"
//#include "greedyRMA.h"

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

  using namespace utilib;
  using namespace base;
  using namespace data;
  using namespace pebblRMA;
  //using namespace greedyRMA;


  class DriverRMA : public BaseRMA {

  public:

    DriverRMA(int& argc, char**& argv);

    ~DriverRMA() {
#ifdef ACRO_HAVE_MPI
      if (parallel) { CommonIO::end(); uMPI::done(); }
#endif // ACRO_HAVE_MPI
    }

    void setData(int& argc, char**& argv) {
      data = new DataRMA(argc, argv, (ArgRMA *) this);
    }

    void setupRMA(int& argc, char**& argv);
    void solveRMA();

  private:

    bool          parallel;

    DataRMA*      data;
    //GreedyRMA*    grma;

    RMA*          rma ;
    parRMA*       prma;

    Time          tc;
    double        wallTime;
    double        cpuTime;

  };

} // end namespace rma

#endif // RMA_h
