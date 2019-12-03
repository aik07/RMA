/**********************************************************
* File name:   driverRMA.cpp
* Author:      Ai Kagawa
* Description: a source file for RMA driver class
***********************************************************/

#include "driverRMA.h"


namespace rma {

DriverRMA::DriverRMA(int& argc, char**& argv): rma(NULL), prma(NULL), parallel(false) {

  cout << setprecision(6) << fixed;

  setup(argc, argv);

  setData(argc, argv);
  setupRMA(argc, argv);

}


void DriverRMA::setData(int& argc, char**& argv) {
  args = this;
  data = new Data(argc, argv, args);
}

void DriverRMA::setupRMA(int& argc, char**& argv) {

  #ifdef ACRO_HAVE_MPI
    uMPI::init(&argc, &argv, MPI_COMM_WORLD);
    //uMPI::init(MPI_COMM_WORLD);
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
  rma->setData(data);

}

// solve RMA
void DriverRMA::solveRMA() {

#ifdef ACRO_HAVE_MPI
  if (parallel) {
    prma->reset();
    prma->printConfiguration();
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
  //setDataWts();

  rma->resetTimers();
  InitializeTiming();
  rma->solve();
  //if (args->printBBdetail) rma->solve();  // print out B&B details
  //else                     rma->search();

#ifdef ACRO_HAVE_MPI
  if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
    rma->printSolutionTime();
#ifdef ACRO_HAVE_MPI
  }
#endif //  ACRO_HAVE_MPI

} // end function solveRMA()

} // end namespace rma
