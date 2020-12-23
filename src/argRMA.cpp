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

    _isPebblRMA(true),     // using PEBBL RMA or not

    _isBinarySearchCutVal(false),  // using the binary search cut point value or not
    _perCachedCutPts(0.000001),    // % of the cached cut point threshold
    _perLimitAttrib(1.0),          // % of attribute to select

    _isRandSeed(true),       // using the random seed or not
    _isInitGuess(true),      // using the initial guess or not
    _branchSelection(0),     // the branching selection
    _isCountingSort(false),  // is using the counting sort

    _isCheckObjVal(false),       // wheather or not to run the check objective run proecude
    _isBruteForceEC(false),      // whether or not to use a brute force equivalence class construction method
    _isBruteForceIncumb(false),  // weather or not to use a brute force incumbent computation method

    _nonUniformWt(""),           // specify the non-uniform weight

    _isSaveInstances(false),
    _isSaveNodeTime(false),
    _isSaveCutPts(false),

    // _isWriteInstances(false),
    // _isWriteNodeTime(false),
    // _isWriteCutPts(false),

    _delta(-1),
    _shrinkDelta(.95),
    _maxInterval(getInf()),

    _fixedSizeBin(-1),

    _isPrintBBdetails(false),

    _rampUpSizeFact(1.0)

  {

    create_categorized_parameter("isPebblRMA", _isPebblRMA,
      "<bool>", "true", "solve RMA using PEBBL", "RMA");

    // PEBBL options

    create_categorized_parameter("isBinarySearchCutVal", _isBinarySearchCutVal,
      "<bool>", "false", "binary search cut values in each feature", "RMA");

    create_categorized_parameter("perCachedCutPts", _perCachedCutPts,
      "<double>", "false", "check only cut-points from the cache"
      "if the cache has at least x% of live cut-points out of total cut points",
      "RMA");

    create_categorized_parameter("perLimitAttrib", _perLimitAttrib, "<double>",
        "1.00", "limit number of attributes to check ", "RMA");

    create_categorized_parameter("isRandSeed", _isRandSeed, "<bool>",
        "true", "random seed for tied solutions", "RMA");

    create_categorized_parameter("isInitGuess", _isInitGuess, "<bool>",
        "true", "enable the initial guess computation", "RMA");

    create_categorized_parameter("branchSelection", _branchSelection, "<int>",
      "0", "Among tied cutpoints, 0: randomize cutpoint to select, "
      "1: always select the first one, 2: always slect the last one", "RMA");

    create_categorized_parameter("isCountingSort", _isCountingSort, "<bool>",
      "false", "Use counting sort instead of bucket sort", "RMA");

    create_categorized_parameter("isCheckObjVal", _isCheckObjVal, "<bool>",
      "false",	"check the optimal solution in the end ", "RMA");

    create_categorized_parameter("isBruteForceEC", _isBruteForceEC, "<bool>",
      "false",	"brute force algorithm to create equivalence classes ", "RMA");

    create_categorized_parameter("isBruteForceIncumb", _isBruteForceIncumb, "<bool>",
      "false",	"brute force algorithm to to compute incumbent in each attribute ",
      "RMA");

    create_categorized_parameter("nonUniformWt", _nonUniformWt, "<string>", "",
      "read non-uniform weights from a file", "RMA");

    create_categorized_parameter("isSaveCutPts", _isSaveCutPts, "<bool>",
      "false", "Write cut points chosen in the solution file ", "RMA");

    create_categorized_parameter("isSaveInstances", _isSaveInstances, "<bool>",
      "false", "Write an input file for each weighted problem solved", "RMA");

    create_categorized_parameter("isSaveNodeTime", _isSaveNodeTime, "<bool>",
      "false", "Write an input file for the number of B&B node and "
      "CPU time for each iteration", "RMA");

    create_categorized_parameter("delta", _delta, "<double>",
      "0", "delta for recursive discretization", "RMA");

    create_categorized_parameter("shrinkDelta", _shrinkDelta, "<double>",
      ".95", "shrink delta for recursive discretization", "RMA");

    create_categorized_parameter("maxInterval", _maxInterval, "<double>",
      "inf", "set the maximum interval length for recursive integerization", "RMA");


    create_categorized_parameter("isPrintBBdetails", _isPrintBBdetails, "<bool>",
      "true", "print the complete output of the PEBBL branch-and bound", "RMA");

    create_categorized_parameter("rampUpSizeFact", _rampUpSizeFact, "<double>",
      "1.00", "if (#storedCutPts) <= rampUpSizeFact * (#processors),"
      "get out the ramp-up", "RMA");

  }

} // namespace BaseRMA
