/**********************************************************
* File name:   driverRMA.cpp
* Author:      Ai Kagawa
* Description: a source file for RMA driver class
***********************************************************/

#include "driverRMA.h"


namespace rma {

  DriverRMA::DriverRMA(int& argc, char**& argv):
             isParallel(false), rma(NULL), prma(NULL) {

#ifdef ACRO_HAVE_MPI
    uMPI::init(&argc, &argv, MPI_COMM_WORLD);
#endif // ACRO_HAVE_MPI

    setup(argc, argv);     // setup all paramaters

    setData(argc, argv);   // set DataRMA class

    setupRMA(argc, argv);  // set RMA class

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
      isParallel = true;
      prma       = new pebblRMA::parRMA(MPI_COMM_WORLD);
      rma        = prma;
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

    // exception_mngr::set_stack_trace(false);
    rma->setup(argc,argv);
    // exception_mngr::set_stack_trace(true);

  }


  void DriverRMA::solveRMA() {
    if (exactRMA()) {

      resetExactRMA();

      if (initGuess()) {
	/*
#ifdef ACRO_HAVE_MPI
	if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
	/*/
    	  solveGreedyRMA();
    	  rma->setInitialGuess(grma->isPosIncumb, grma->maxObjValue,
    			       grma->L, grma->U);
	  /*
#ifdef ACRO_HAVE_MPI
	}
#endif //  ACRO_HAVE_MPI
	  /*/

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
    if (isParallel) {
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
    rma->workingSol.value = getInf();

  }


  // solve RMA
  void DriverRMA::solveExactRMA() {

    rma->resetTimers();
    InitializeTiming();

    tc.startTime();

    if (printBBdetails()) rma->solve();  // print out B&B details
    else                  rma->search();

    tc.getCPUTime();
    tc.getWallTime();
    printSolutionTime();

    CommonIO::end();
    uMPI::done();

  } // end function solveExactRMA()


  void DriverRMA::printSolutionTime() {

    double global_solution = rma->workingSol.value;
    int    total_nodes     = rma->subCount[2];

    if (uMPI::size>1) {
      // ucout << "reduce " << " " << total_nodes << "\n";
      MPI_Reduce(&rma->workingSol.value, &global_solution,
                 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce(&rma->subCount[2],      &total_nodes,
                 1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    }

#ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
    std::cout << std::fixed << std::setprecision(4)
              << "ERMA Solution: "  << global_solution;
    std::cout << std::fixed << std::setprecision(2)
              << " \tCPU time: "     << tc.getCPUTime()
              << " \tNum of Nodes: " << total_nodes << "\n";
#ifdef ACRO_HAVE_MPI
    }
#endif //  ACRO_HAVE_MPI

  }

} // end namespace rma
