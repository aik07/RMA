#include "rma.h"

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
  rma->numDistObs = data->numTrain;	    // only use training data
  rma->setSortObsNum(data->vecTrainData);
  //setDataWts();

  rma->resetTimers();
  InitializeTiming();
  if (args->printBBdetail) rma->solve();  // print out B&B details
  else                     rma->search();

#ifdef ACRO_HAVE_MPI
  if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
    rma->printSolutionTime();
#ifdef ACRO_HAVE_MPI
  }
#endif //  ACRO_HAVE_MPI

} // end function solveRMA()
