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

  setData(argc, argv);   // set data
/*
#ifdef ACRO_HAVE_MPI
  if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI+
  setData(argc, argv);   // set data
#ifdef ACRO_HAVE_MPI
  }
#endif //  ACRO_HAVE_MPI
/*/
  setupRMA(argc, argv);  // setup RMA

    setupRMA(argc, argv);

    // cout << "test" << testWeight();

    if (testWeight()!="") {
      vector<double> vecNonUniformWt;
      vector<int>    vecObsIdx;
      vecNonUniformWt.resize(data->numTrainObs);
      vecObsIdx.resize(data->numTrainObs);
      for (int i=0; i<data->numTrainObs; ++i)
	       vecObsIdx[i] = i;
         ifstream inFile(testWeight());
         if (inFile.is_open()) {
      	string line;
      	string tmp;
      	while( getline(inFile,line) ) {
      	  stringstream ss(line);
      	  for (int i=0; i<data->numTrainObs; ++i) {
      	    getline(ss,tmp,',');
      	    // cout << "tmp " << tmp << "\n";
      	    vecNonUniformWt[i] = stod(tmp);
      	    // cout << "vec " <<vecNonUniformWt[i] << "\n";
      	  }
        }
      }
      rma->setWeight(vecNonUniformWt, vecObsIdx);
    }

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

} // end function solveExactRMA()

} // end namespace rma
