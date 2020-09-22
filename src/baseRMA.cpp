#include "baseRMA.h"


namespace base {

  ////////////////////// Base class methods ////////////////////////

  // Standard serial read-in code.  Returns true if we can continue, false if
  // we have to bail out.
  bool BaseRMA::setup(int& argc, char**& argv) {

    if (!processParameters(argc,argv,min_num_required_args))
      return false;

#ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
      if (plist.size() == 0) {
        ucout << "Using default values for all solver options" << std::endl;
      } else {
        ucout << "User-specified solver options: " << std::endl;
	      plist.write_parameters(ucout);
	      ucout << std::endl;
      }
#ifdef ACRO_HAVE_MPI
    }
#endif //  ACRO_HAVE_MPI

    set_parameters(plist,false);

    if ((argc > 0) && !checkParameters(argv[0]))
      return false;

    if (!setupProblem(argc,argv))
      return false;

#ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
      if (plist.unused() > 0) {
      	ucout << "\nERROR: unused parameters: " << std::endl;
      	plist.write_unused_parameters(ucout);
      	ucout << utilib::Flush;
      	return false;
      }
#ifdef ACRO_HAVE_MPI
    }
#endif //  ACRO_HAVE_MPI

    return true;

  }


  bool BaseRMA::processParameters(int& argc, char**& argv,
  				  unsigned int min_num_required_args__) {

    if (argc > 0)
      solver_name = argv[0];
    else
      solver_name = "unknown";
    if (!parameters_registered) {
      register_parameters();
      parameters_registered=true;
    }

    plist.process_parameters(argc,argv,min_num_required_args__);

    // Set the name of the problem to be the last thing on the command
    // line. setName will extract the filename root. The setupProblem
    // method can overwrite this later.
    if ((argc > 1) && (argv[argc-1] != NULL))
      setName(argv[1]);

    return true;
  }


  bool BaseRMA::checkParameters(char const* progName) {

    if (help_parameter) {
      write_usage_info(progName,cout);
      return false;
    }

#ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
      if (debug_solver_params) {
	ucout << "---- Parameters ----" << endl;
	write_parameter_values(ucout);
	ucout << endl << utilib::Flush;
      }
#ifdef ACRO_HAVE_MPI
    }
#endif //  ACRO_HAVE_MPI

    return true;
  }


  void BaseRMA::write_usage_info(char const* progName,std::ostream& os) const {
    writeCommandUsage(progName,os);
    os << endl;
    plist.write_registered_parameters(os);
    os << endl;
  }


  void BaseRMA::writeCommandUsage(char const* progName,std::ostream& os) const {
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

  void BaseRMA::setName(const char* cname) {
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

} // end namespacs base
