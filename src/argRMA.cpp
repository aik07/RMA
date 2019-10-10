/*
 *  File name: rmaParams.h
 *  Author:    Ai Kagawa
 */


#include "argRMA.h"


using utilib::ParameterLowerBound;
using utilib::ParameterBounds;
using utilib::ParameterNonnegative;

namespace arg {

ArgRMA::ArgRMA():

  _binarySearchCutVal(false),
  _perCachedCutPts(0.000001),
  _perLimitAttrib(1.0),

  _randSeed(true),
  _initGuess(true),
  _branchSelection(0),
  _countingSort(false),

  _testWt(false),

  _checkObjVal(false),
  _bruteForceEC(false),
  _bruteForceIncumb(false),

  _writeInstances(false),
  _writeNodeTime(false),
  _writeCutPts(false),

  _delta(-1),
  _shrinkDelta(.95),
  _maxInterval(inf),

  _fixedSizeBin(-1),

  _rampUpSizeFact(1.0)

  {
    create_categorized_parameter("binarySearchCutVal", _binarySearchCutVal,
      "<bool>", "false", "binary search cut values in each feature", "RMA");

    create_categorized_parameter("perCachedCutPts", _perCachedCutPts,
      "<double>", "false", "check only cut-points from the cache"
      "if the cache has at least x% of live cut-points out of total cut points",
      "RMA");

    create_categorized_parameter("perLimitAttrib", _perLimitAttrib, "<double>",
        "1.00", "limit number of attributes to check ", "RMA");

    create_categorized_parameter("randSeed", _randSeed, "<bool>",
        "true", "random seed for tied solutions", "RMA");

    create_categorized_parameter("initGuess", _initGuess, "<bool>",
        "true", "enable the initial guess computation", "RMA");

    create_categorized_parameter("branchSelection", _branchSelection, "<int>",
      "0", "Among tied cutpoints, 0: randomize cutpoint to select, "
      "1: always select the first one, 2: always slect the last one", "RMA");

    create_categorized_parameter("countingSort", _countingSort, "<bool>",
      "false", "Use counting sort instead of bucket sort", "RMA");

    create_categorized_parameter("checkObjVal", _checkObjVal, "<bool>",
      "false",	"check the optimal solution in the end ", "RMA");

    create_categorized_parameter("bruteForceEC", _bruteForceEC, "<bool>",
      "false",	"brute force algorithm to create equivalence classes ", "RMA");

    create_categorized_parameter("bruteForceIncumb", _bruteForceIncumb, "<bool>",
      "false",	"brute force algorithm to to compute incumbent in each attribute ",
      "RMA");

    create_categorized_parameter("testWt", _testWt, "<bool>", "false",
      "testing with specified test weights data, testWt.data", "RMA");

    create_categorized_parameter("writeCutPts", _writeCutPts, "<bool>",
      "false", "Write cut points chosen in the solution file ", "RMA");

    create_categorized_parameter("writeInstances", _writeInstances, "<bool>",
      "false", "Write an input file for each weighted problem solved", "RMA");

    create_categorized_parameter("writeNodeTime", _writeNodeTime, "<bool>",
      "false", "Write an input file for the number of B&B node and "
      "CPU time for each iteration", "RMA");

    create_categorized_parameter("delta", _delta, "<double>",
      "0", "delta for recursive discretization", "RMA");

    create_categorized_parameter("shrinkDelta", _shrinkDelta, "<double>",
      ".95", "shrink delta for recursive discretization", "RMA");

    create_categorized_parameter("maxInterval", _maxInterval, "<double>",
      "inf", "set the maximum interval length for recursive integerization", "RMA");

    create_categorized_parameter("rampUpSizeFact", _rampUpSizeFact, "<double>",
      "1.00", "if (#storedCutPts) <= rampUpSizeFact * (#processors),"
      "get out the ramp-up", "RMA");

  }
  ///////////////////// Arguments class methods ////////////////////////

  // Standard serial read-in code.  Returns true if we can continue, false if
  // we have to bail out.
  bool Arguments::setup(int argc, char** argv) {

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


  bool Arguments::processParameters(int& argc, char**& argv,
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


  bool Arguments::checkParameters(char const* progName) {

    if (help_parameter) {
      write_usage_info(progName,cout);
      return false;
    }

    if (debug_solver_params) {
      ucout << "---- Parameters ----" << endl;
      write_parameter_values(ucout);
      ucout << endl << utilib::Flush;
    }

    return true;
  }


  void Arguments::write_usage_info(char const* progName,std::ostream& os) const {
    writeCommandUsage(progName,os);
    os << endl;
    plist.write_registered_parameters(os);
    os << endl;
  }


  void Arguments::writeCommandUsage(char const* progName,std::ostream& os) const {
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

  void Arguments::setName(const char* cname) {
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

} // namespace arguments
