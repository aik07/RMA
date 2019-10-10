/**********************************************************
* File name:   driverRMA.cpp
* Author:      Ai Kagawa
* Description: a source file for RMA driver class
***********************************************************/

#include "driverRMA.h"


namespace pebblRMA {

  DriverRMA::DriverRMA(int argc, char** argv):
        parameters_registered(false), min_num_required_args(0),
        rma(NULL), prma(NULL), parallel(false) {

    cout << setprecision(6) << fixed;

    //args = new ArgRMA();
    setup(argc, argv);
    data = new Data(this);
    data->readData(argc, argv);

    #ifdef ACRO_HAVE_MPI
      //uMPI::init(&argc, &argv, MPI_COMM_WORLD);
      int nprocessors = uMPI::size;
      /// Do parallel optimization if MPI indicates that we're using more than one processor
      if (parallel_exec_test<parallelBranching>(argc, argv, nprocessors)) {
        /// Manage parallel I/O explicitly with the utilib::CommonIO tools
        CommonIO::begin();
        CommonIO::setIOFlush(1);
        parallel = true;
        prma     = new parRMA;
        rma      = prma;
        //prma->setParameter(data, data->debug);
        //prma->setParameters(args);
        //prma->setData(data);
      } else {
    #endif // ACRO_HAVE_MPI
        rma = new RMA;
    #ifdef ACRO_HAVE_MPI
      }
    #endif // ACRO_HAVE_MPI

    rma->setParameters(this);
    rma->setData(data);

    //rma->setParameters(data, data->debug);

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




////////////////////// Base class methods ////////////////////////

// Standard serial read-in code.  Returns true if we can continue, false if
// we have to bail out.
bool DriverRMA::setup(int& argc, char**& argv) {

  if (!processParameters(argc, argv, min_num_required_args))
    return false;

  if (plist.size() == 0) {
      ucout << "Using default values for all solver options" << std::endl;
  } else {
    ucout << "User-specified solver options: " << std::endl;
    plist.write_parameters(ucout);
    ucout << std::endl;
  }

  set_parameters(plist,false);

  if ((argc > 0) && !checkParameters(argv[0]))
    return false;

  if (!setupProblem(argc,argv))
    return false;

  if (plist.unused() > 0) {
    ucout << "\nERROR: unused parameters: " << std::endl;
    plist.write_unused_parameters(ucout);
    ucout << utilib::Flush;
    return false;
  }

  return true;

}


bool DriverRMA::processParameters(int& argc, char**& argv,
                unsigned int min_num_required_args) {

  if (argc > 0) solver_name = argv[0];
  else          solver_name = "unknown";

  if (!parameters_registered) {
    register_parameters();
    parameters_registered=true;
  }

  plist.process_parameters(argc, argv, min_num_required_args);

  // Set the name of the problem to be the last thing on the command
  // line. setName will extract the filename root. The setupProblem
  // method can overwrite this later.
  if ((argc > 1) && (argv[argc-1] != NULL))
    setName(argv[1]);

  return true;
}


bool DriverRMA::checkParameters(char const* progName) {

  if (help_parameter) {
    write_usage_info(progName,cout);
    return false;
  }

  if (debug_solver_params) {
    ucout << "---- LPBoost Parameters ----" << endl;
    write_parameter_values(ucout);
    ucout << endl << utilib::Flush;
  }

  return true;
}


void DriverRMA::write_usage_info(char const* progName,std::ostream& os) const {
  writeCommandUsage(progName,os);
  os << endl;
  plist.write_registered_parameters(os);
  os << endl;
}


void DriverRMA::writeCommandUsage(char const* progName,std::ostream& os) const {
  os << "\nUsage: " << progName << " { --parameter=value ... }";
  if (min_num_required_args == 1)
    os << " <problem data file>";
  else if (min_num_required_args == 1)
    os << " <" << min_num_required_args << " problem data files>";
  os << endl;
}


// This sets the official name of the problem by chewing up the
// filename.  It can be overridden.  This version just finds the last
// "/" or "\" in the name and removes it and everything before it.

void DriverRMA::setName(const char* cname) {
#if defined (TFLOPS)
  problemName = cname;
  int i=problemName.size();
  while (i >= 0) {
    if (cname[i] == '/') break;
    i--;
  }
  if (i >= 0)
     problemName.erase(0,i+1);
  // TODO: remove the .extension part for this case
#else
  problemName = cname;
  size_type i = problemName.rfind("/");
  if (i == string::npos)
    i = problemName.rfind("\\");
  if (i != string::npos)
    problemName.erase(0,i+1);

  size_type n = problemName.length();

  if (n < 4)
    return;

  string endOfName(problemName,n-4,4);
  if ((endOfName == ".dat") || (endOfName == ".DAT"))
    problemName.erase(n-4,n);
  if ((endOfName == ".data") || (endOfName == ".DATA"))
      problemName.erase(n-5,n);
#endif
}

} // end namespace pebblRMA
