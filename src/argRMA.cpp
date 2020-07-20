/**********************************************************
 *  File name: argRMA.h
 *  Author:    Ai Kagawa
**********************************************************/

#include "argRMA.h"

using utilib::ParameterLowerBound;
using utilib::ParameterBounds;
using utilib::ParameterNonnegative;

namespace arg {

  ArgRMA::ArgRMA():

    _exactRMA(true),

    _binarySearchCutVal(false),
    _perCachedCutPts(0.000001),
    _perLimitAttrib(1.0),

    _randSeed(true),
    _initGuess(true),  
    _branchSelection(0),
    _countingSort(false),

    _testWt(""),

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

    _printBBdetails(false),

    _rampUpSizeFact(1.0)

  {

    create_categorized_parameter("exactRMA", _exactRMA,
      "<bool>", "true", "solve RMA exactly", "RMA");

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

    create_categorized_parameter("testWt", _testWt, "<string>", "",
      "testing with specified test weights data", "RMA");

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


    create_categorized_parameter("printBBdetails", _printBBdetails, "<bool>",
      "false", "print the complete output of the PEBBL branch-and bound", "RMA");

    create_categorized_parameter("rampUpSizeFact", _rampUpSizeFact, "<double>",
      "1.00", "if (#storedCutPts) <= rampUpSizeFact * (#processors),"
      "get out the ramp-up", "RMA");

  }

} // namespace BaseRMA
