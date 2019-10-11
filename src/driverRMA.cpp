/**********************************************************
* File name:   driverRMA.cpp
* Author:      Ai Kagawa
* Description: a source file for RMA driver class
***********************************************************/

#include "driverRMA.h"


namespace rma {

  DriverRMA::DriverRMA(ArgRMA* args_, Data* data_):
        args(args_), data(data_), rma(NULL), prma(NULL), parallel(false) {

    cout << setprecision(6) << fixed;

    #ifdef ACRO_HAVE_MPI
      //uMPI::init(&data->argc, &data->argv, MPI_COMM_WORLD);
      int nprocessors = uMPI::size;
      /// Do parallel optimization if MPI indicates that we're using more than one processor
      if (parallel_exec_test<parallelBranching>(data->argc, data->argv, nprocessors)) {
        /// Manage parallel I/O explicitly with the utilib::CommonIO tools
        CommonIO::begin();
        CommonIO::setIOFlush(1);
        parallel = true;
        prma     = new parRMA;
        rma      = prma;
      } else {
    #endif // ACRO_HAVE_MPI
        rma = new RMA;
    #ifdef ACRO_HAVE_MPI
      }
    #endif // ACRO_HAVE_MPI

    rma->setParameters(args);
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

  rma->workingSol.value=-inf;
  rma->mmapCachedCutPts.clear();
  rma->numDistObs = data->numTrainObs;	    // only use training data
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
