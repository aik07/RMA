/**********************************************************
* File name:   solveRMA.cpp
* Author:      Ai Kagawa
* Description: a source file for RMA driver class
***********************************************************/

#include "solveRMA.h"


namespace rma {

  // setup to sovle RMA
  void SolveRMA::setupSolveRMA(int& argc, char**& argv) {

    setup(argc, argv);           // setup all paramaters

    setData(argc, argv);         // set DataRMA class

    setupPebblRMA(argc, argv);   // set RMA class

  }


  void SolveRMA::setupPebblRMA(int& argc, char**& argv) {

#ifdef ACRO_HAVE_MPI
    uMPI::init(&argc, &argv, MPI_COMM_WORLD);
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
      rma = new pebblRMA::RMA(this);
#ifdef ACRO_HAVE_MPI
    }
#endif // ACRO_HAVE_MPI

    rma->setParameters(this);    // passing arguments
    rma->setData(data);          // set data

#ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
      rma->setSortedObsIdx(data->vecTrainObsIdx);
#ifdef ACRO_HAVE_MPI
    }
#endif //  ACRO_HAVE_MPI

    // exception_mngr::set_stack_trace(false);
    rma->setup(argc,argv);
    // exception_mngr::set_stack_trace(true);

  }


  // solve RMA using the chosen methods
  void SolveRMA::solveRMA() {

    if (isPebblRMA()) { // if RMA is solved using PEBBL

      resetPebblRMA();  // reset PEEBL RMA variables

      if (isInitGuess()) {  // if the PEBBL get initial guess by solving the greedy RMA

// TODO: the greedy RMA can be solved using only one process

  /*
#ifdef ACRO_HAVE_MPI
	if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
	/*/
    	  solveGreedyRMA();  // solve the greedy RMA
        // set the initial guess solution using the greedy RMA solution
        // (positive or negative solution, initial objective value,
        //  lower and upper bounds)
    	  rma->setInitialGuess(grma->isPostObjVal(),   grma->getObjVal(),
    			                   grma->getLowerBounds(), grma->getUpperBounds());
	  /*
#ifdef ACRO_HAVE_MPI
	}
#endif //  ACRO_HAVE_MPI
	  /*/

      }

      solvePebblRMA();  // solve RMA using PEBBL

    } else { // if RMA is solved using PEBBL

#ifdef ACRO_HAVE_MPI
      if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
	      solveGreedyRMA();  // solve the greedy RMA
#ifdef ACRO_HAVE_MPI
      }
#endif //  ACRO_HAVE_MPI

    }
  }


  // solve Greedy RMA
  void SolveRMA::solveGreedyRMA() {
    grma = new greedyRMA::GreedyRMA(this, data);
    grma->runGreedyRangeSearch();
  }


  // reset PEBBL RMA variables
  void SolveRMA::resetPebblRMA() {

#ifdef ACRO_HAVE_MPI
    if (isParallel) { // if in parallel
      prma->reset();                 // reset the prallel RMA class
      if (isPrintBBdetails())        // if B&B details should be shown
        prma->printConfiguration();  // print the PEBBL configuration
      CommonIO::begin_tagging();
    } else { // if in serial
#endif //  ACRO_HAVE_MPI
      rma->reset(); // reset the serial PEBBL RMA settings
#ifdef ACRO_HAVE_MPI
    }
#endif //  ACRO_HAVE_MPI

    rma->mmapCachedCutPts.clear();     // clean up the cached cut points storage
    rma->workingSol.value = -getInf(); // set the working solution value to be negative infinity

  }


  // solve RMA using PEBBL
  void SolveRMA::solvePebblRMA() {

    rma->resetTimers();  // reset PEBBL timers
    InitializeTiming();

    tc.startTime();      // start the timer

    if (isPrintBBdetails()) rma->solve();  // print out B&B details
    else                    rma->search(); // only get the RMA optimal value (no B&B details)

    printPebblRMASolutionTime();

  } // end function solvePebblRMA()


  void SolveRMA::printPebblRMASolutionTime() {

    double global_solution = rma->workingSol.value; // set the current node's RMA solution value
    int    total_nodes     = rma->subCount[2];      // set the current node's subproblems which

    if (uMPI::size>1) { // if the PEBBL RMA is solved in parllel

      // ucout << "reduce " << " " << total_nodes << "\n";

      // set the global solution value to be the maximum among all nodes's solutions
      MPI_Reduce(&rma->workingSol.value, &global_solution,
                 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

      // set the total node to be the sum of the nodes among all nodes
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
