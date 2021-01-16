/**********************************************************
* File name:   solveRMA.cpp
* Author:      Ai Kagawa
* Description: a source file for RMA solver class
***********************************************************/


#include "solveRMA.h"


namespace rma {


  // setup to sovle RMA
  void SolveRMA::setupSolveRMA(int& argc, char**& argv) {

    setup(argc, argv);           // setup all paramaters

    setData(argc, argv);         // set DataRMA class

    setupPebblRMA(argc, argv);   // set PEBBL RMA class

  } // end setupSolveRMA function


  // setup for PEBBL RMA
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
      prma       = new pebblRMA::parRMA(this, MPI_COMM_WORLD);
      rma        = prma;
    } else {
#endif // ACRO_HAVE_MPI
      rma = new pebblRMA::RMA(this);
#ifdef ACRO_HAVE_MPI
    }
#endif // ACRO_HAVE_MPI

    rma->setParameters(this);    // passing arguments
    rma->setData(data);          // set data class in RMA

    // exception_mngr::set_stack_trace(false);
    rma->setup(argc,argv);
    // exception_mngr::set_stack_trace(true);

  } // end setupPebblRMA function


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

    rma->setSortedObsIdx(data->vecNonZeroWtObsIdx);

  }  // end resetPebblRMA function


  // solve RMA using the chosen methods
  void SolveRMA::solveRMA() {

    if (isRMAonly) { // if RMA only, not Boosting
      data->setDataIntWeight();   // set weights for dataIntTrain
      data->setNumPosNegObs();    // set # of positive and negative observations
    }

    data->removeZeroWtObs();  // remot zero weight observations

    if (isPebblRMA()) { // if RMA is solved using PEBBL
      resetPebblRMA();  // reset PEEBL RMA variables

      if (isInitGuess()) {  // if the PEBBL get initial guess by solving the greedy RMA

        // TODO: the greedy RMA can be solved using only one process
        // if (ROOTPROC) { // if root process

           solveGreedyRMA();  // solve the greedy RMA

          // set the initial guess solution using the greedy RMA solution
          // (positive or negative solution, initial objective value,
          //  lower and upper bounds)
          rma->setInitialGuess(grma->isPostObjVal(),   grma->getObjVal(),
                               grma->getLowerBounds(), grma->getUpperBounds());

        // } // end if root process

      } // end if isInitialGuess()

      solvePebblRMA();  // solve RMA using PEBBL

    } else { // if RMA is solved using PEBBL

      if (ROOTPROC) { // if root process
	      solveGreedyRMA();  // solve the greedy RMA
      } // end if root process

    } // end if RMA is solved using PEBBL

  }  // end solveRMA function


  // solve Greedy RMA
  void SolveRMA::solveGreedyRMA() {
    grma = new greedyRMA::GreedyRMA(this, data);
    grma->runGreedyRangeSearch();
  }


  // solve RMA using PEBBL
  void SolveRMA::solvePebblRMA() {

    rma->resetTimers();  // reset PEBBL timers
    InitializeTiming();

    tc.startTime();      // start the timer

    if (isPrintBBdetails()) rma->solve();  // print out B&B details
    else                    rma->search(); // only get the RMA optimal value (no B&B details)

    rma->printSolutionTime(tc.getCPUTime());

  } // end function solvePebblRMA()


} // end namespace rma


// #include <pybind11/pybind11.h>
// namespace py = pybind11;
//
//
// PYBIND11_MODULE(SolveRMA, m) {
//
//     // optional module docstring
//     m.doc() = "pybind11 plugin for RMA";
//
//     // bindings to Pet class
//     py::class_<SolveRMA>(m, "SolveRMA")
//         .def(py::init<>())
//         .def("setupSolveRMA", &SolveRMA::setupSolveRMA)
//         .def("solveRMA",      &SolveRMA::solveRMA);
// }
