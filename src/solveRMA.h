/**********************************************************
* Author:      Ai Kagawa
* Description: a header file for RMA solver class
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
#include "utility.h"

#ifdef ACRO_HAVE_MPI
  #include <pebbl/pbb/parBranching.h>
  #include "parRMA.h"
  #define outstream ucout
  #define IO(action) if (uMPI::iDoIO) { CommonIO: end_tagging(); action; }
#else // ACRO_HAVE_MPI
  typedef void parRMA;
  #define outstream cout;
  #define IO(action) action;
#endif // ACRO_HAVE_MPI

#ifdef ACRO_HAVE_MPI
#define ROOTPROC uMPI::rank==0
#else
#define ROOTPROC true
#endif


namespace rma {

  // SolveRMA class calls selected methods (Greedy RMA and/or PEBBL RMA)
  // to solve the RMA problem
  class SolveRMA : virtual public base::BaseRMA {

  public:

    SolveRMA(): isRMAonly(true), isParallel(false),
                 rma(NULL), prma(NULL) { }

    virtual ~SolveRMA() {
#ifdef ACRO_HAVE_MPI
      //if (isParallel)
      { CommonIO::end(); uMPI::done(); } // MPI_Finalize
#endif // ACRO_HAVE_MPI
    }

    void initMPI(int& argc, char**& argv) {
#ifdef ACRO_HAVE_MPI
      uMPI::init(&argc, &argv, MPI_COMM_WORLD);
#endif // ACRO_HAVE_MPI
    }

    void setupSolveRMA(int& argc, char**& argv);  // setup to sovle RMA

    // set Data RMA class object
    virtual void setDataRMA(int& argc, char**& argv) {
      data = new data::DataRMA(argc, argv, (ArgRMA *) this);
    }

    void setupPebblRMA(int& argc, char**& argv);  // setup PebblRMA

    void resetPebblRMA();    // reset PEBBL RMA variables

    void solveRMA();         // solve RMA using the selected methods
    void solveGreedyRMA();   // solve RMA using the greedy algorithm
    void solvePebblRMA();    // solve RMA using PEBBL

    // check objective value
    void checkObjValue(vector<DataXw> dataInt,
                       vector<unsigned int> a, vector<unsigned int> b);

  protected:

    bool                  isRMAonly;   // it is true for solve RMA, not Boosting

    bool                  isParallel;  // whether or not this program will run in parallel

    data::DataRMA*        data;        // RMA data object

    greedyRMA::GreedyRMA* grma;        // Greedy RMA object

    pebblRMA::RMA*        rma ;        // Serial   PEBBL RMA object
    pebblRMA::parRMA*     prma;        // Parallel PEBBL RMA object

    Time                  tc;          // Time object

  };

} // end namespace rma

#endif // RMA_h
