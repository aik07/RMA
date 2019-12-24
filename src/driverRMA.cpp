/**********************************************************
* File name:   driverRMA.cpp
* Author:      Ai Kagawa
* Description: a source file for RMA driver class
***********************************************************/

#include "driverRMA.h"


namespace rma {

  DriverRMA::DriverRMA(int& argc, char**& argv): rma(NULL), prma(NULL), parallel(false) {

#ifdef ACRO_HAVE_MPI
    uMPI::init(&argc, &argv, MPI_COMM_WORLD);
#endif // ACRO_HAVE_M

    //cout << setprecision(6) << fixed;

    setup(argc, argv);     // setup all paramaters

#ifdef ACRO_HAVE_MPI
if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
    setData(argc, argv);   // set data
#ifdef ACRO_HAVE_MPI
}
#endif //  ACRO_HAVE_MPI

    setupRMA(argc, argv);  // (setup) RMA

  }

  void DriverRMA::setupRMA(int& argc, char**& argv) {

#ifdef ACRO_HAVE_MPI
    int nprocessors = uMPI::size;
    /// Do parallel optimization if MPI indicates that we're using more than one processor
    if (parallel_exec_test<parallelBranching>(argc, argv, nprocessors)) {
      /// Manage parallel I/O explicitly with the utilib::CommonIO tools
      CommonIO::begin();
      CommonIO::setIOFlush(1);
      parallel = true;
      prma     = new parRMA(MPI_COMM_WORLD);
      rma      = prma;
    } else {
#endif // ACRO_HAVE_MPI
      rma = new RMA;
#ifdef ACRO_HAVE_MPI
    }
#endif // ACRO_HAVE_MPI

    rma->setParameters(this); // passing arguments

#ifdef ACRO_HAVE_MPI
if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
    rma->setData(data);
#ifdef ACRO_HAVE_MPI
}
#endif //  ACRO_HAVE_MPI

    exception_mngr::set_stack_trace(false);
    rma->setup(argc,argv);
    exception_mngr::set_stack_trace(true);

  }


  void DriverRMA::solveRMA() {
    if (exactRMA()) {

      resetExactRMA();

      if (initGuess()) {
	       solveGreedyRMA();
         rma->setInitialGuess(grma->isPosIncumb, grma->maxObjValue,
                              grma->L, grma->U);
      }
      solveExactRMA();
    } else {
      solveGreedyRMA();
    }
  }


  void DriverRMA::solveGreedyRMA() {
    grma = new GreedyRMA(this, data);
    grma->runGreedyRangeSearch();
  }

  void DriverRMA::resetExactRMA() {

#ifdef ACRO_HAVE_MPI
    if (parallel) {
      prma->reset();
      if (printBBdetails()) prma->printConfiguration();
      CommonIO::begin_tagging();
    } else {
#endif //  ACRO_HAVE_MPI
      rma->reset();
#ifdef ACRO_HAVE_MPI
    }
#endif //  ACRO_HAVE_MPI

    rma->mmapCachedCutPts.clear();
    rma->workingSol.value = -inf;
    rma->numDistObs       = data->numTrainObs;	    // only use training data
    rma->setSortObsNum(data->vecTrainData);

  }

  // solve RMA
  void DriverRMA::solveExactRMA() {

    rma->resetTimers();
    InitializeTiming();
    if (printBBdetails()) rma->solve();  // print out B&B details
    else                  rma->search();

#ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
      rma->printSolutionTime();
#ifdef ACRO_HAVE_MPI
    }
#endif //  ACRO_HAVE_MPI

    // CommonIO::end();
    // uMPI::done();

  } // end function solveExactRMA()

} // end namespace rma
