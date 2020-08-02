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

    setup(argc, argv);     // setup all paramaters

    setData(argc, argv);   // set data

    setupRMA(argc, argv);

  }


  void DriverRMA::setData(int& argc, char**& argv) {
    data = new data::DataRMA(argc, argv, (ArgRMA *) this);
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
      prma     = new pebblRMA::parRMA(MPI_COMM_WORLD);
      rma      = prma;
    } else {
#endif // ACRO_HAVE_MPI
      rma = new pebblRMA::RMA;
#ifdef ACRO_HAVE_MPI
    }
#endif // ACRO_HAVE_MPI

    rma->setParameters(this); // passing arguments
    rma->setData(data);

#ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
      rma->setSortedObsIdx(data->vecTrainData);
#ifdef ACRO_HAVE_MPI
    }
#endif //  ACRO_HAVE_MPI

    //exception_mngr::set_stack_trace(false);
    rma->setup(argc,argv);
    //exception_mngr::set_stack_trace(true);

  }


  void DriverRMA::solveRMA() {
    if (exactRMA()) {

      resetExactRMA();

      if (initGuess()) {
	/*
#ifdef ACRO_HAVE_MPI
	if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
	*/
	  solveGreedyRMA();
	  rma->setInitialGuess(grma->isPosIncumb, grma->maxObjValue,
			       grma->L, grma->U);
	  /*
#ifdef ACRO_HAVE_MPI
	}
#endif //  ACRO_HAVE_MPI
	  */

      }

      solveExactRMA();

    } else {

#ifdef ACRO_HAVE_MPI
      if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
	solveGreedyRMA();
#ifdef ACRO_HAVE_MPI
      }
#endif //  ACRO_HAVE_MPI

    }
  }


  void DriverRMA::solveGreedyRMA() {
    grma = new greedyRMA::GreedyRMA(this, data);
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

  }


  // solve RMA
  void DriverRMA::solveExactRMA() {

    rma->resetTimers();
    InitializeTiming();

    tc.startTime();

    if (printBBdetails()) rma->solve();  // print out B&B details
    else                  rma->search();

#ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
      tc.getCPUTime();
      tc.getWallTime();
      printSolutionTime();
#ifdef ACRO_HAVE_MPI
    }
#endif //  ACRO_HAVE_MPI

    int global_sum;
    MPI_Reduce(&rma->subCount[2], &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (uMPI::size>1 and uMPI::rank==0) ucout << "\tNum of Nodes: " << global_sum << "\n";

    CommonIO::end();
    uMPI::done();

  } // end function solveExactRMA()


  void DriverRMA::printSolutionTime() {
    ucout << "ERMA Solution: "  << rma->workingSol.value
          << "\tCPU time: "     << tc.getCPUTime();

    if (uMPI::size==1)
      ucout << "\tNum of Nodes: " << rma->subCount[2] << "\n";
  }

} // end namespace rma
