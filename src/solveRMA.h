/**********************************************************
* File name:   solveRMA.h
* Author:      Ai Kagawa
* Description: a header file for RMA driver class
***********************************************************/

#ifndef RMA_h
#define RMA_h

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <pebbl_config.h>
#include <pebbl/utilib/CommonIO.h>
#include <pebbl/bb/branching.h>

#include "baseRMA.h"
#include "dataRMA.h"
#include "serRMA.h"
#include "greedyRMA.h"
#include "Time.h"
#include "utilRMA.h"

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

  class SolveRMA : virtual public base::BaseRMA {

  public:

    SolveRMA(): isParallel(false), rma(NULL), prma(NULL) {}

    virtual ~SolveRMA() {
#ifdef ACRO_HAVE_MPI
      if (isParallel) { CommonIO::end(); uMPI::done(); }
#endif // ACRO_HAVE_MPI
    }

    void setupSolveRMA(int& argc, char**& argv);

    virtual void setData(int& argc, char**& argv) {
      data = new data::DataRMA(argc, argv, (ArgRMA *) this);
    }

    void setupPebblRMA(int& argc, char**& argv);

    void updateWt();

    void resetExactRMA();

    void solveRMA();
    void solveGreedyRMA();
    void solveExactRMA();

    void printRMASolutionTime();

  protected:

    bool                  isParallel;

    data::DataRMA*        data;
    greedyRMA::GreedyRMA* grma;

    pebblRMA::RMA*        rma ;
    pebblRMA::parRMA*     prma;

    Time          tc;

  };

} // end namespace rma

#endif // RMA_h
