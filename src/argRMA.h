/**********************************************************
 * File name:   argRMA.h
 * Author:      Ai Kagawa
 * Description: a header file for RMA argument class
**********************************************************/

#ifndef ARG_RMA_h
#define ARG_RMA_h

#include <limits>
#include <string>

#include <pebbl_config.h>
#include <pebbl/utilib/ParameterSet.h>
#include <pebbl/utilib/ParameterList.h>
#include <pebbl/utilib/CommonIO.h>
#include <pebbl/utilib/memdebug.h>
#include <pebbl/utilib/seconds.h>
#include <pebbl/bb/branching.h>
#include <pebbl/bb/pebblParams.h>
#include <pebbl/pbb/parPebblParams.h>

#include "utilRMA.h"

using namespace std;
using namespace utilib;
using namespace pebbl;

namespace arg {


  /////////////////// Parameters for RMA class  ///////////////////
  class ArgRMA : virtual public ParameterSet, virtual public CommonIO {

  public:

    ArgRMA();
    virtual ~ArgRMA(){};

    ////////////////////// parameters //////////////////////////////

    bool   isPebblRMA()           const {return _isPebblRMA;}

    bool   isBinarySearchCutVal() const {return _isBinarySearchCutVal;}
    double perCachedCutPts()      const {return _perCachedCutPts;}
    double perLimitAttrib()       const {return _perLimitAttrib;}

    bool   isRandSeed()           const {return _isRandSeed;}
    bool   isInitGuess()          const {return _isInitGuess;}
    int    branchSelection()      const {return _branchSelection;}
    bool   isCountingSort()       const {return _isCountingSort;}

    string nonUniformWt()         const {return _nonUniformWt;}

    bool   isCheckObjVal()        const {return _isCheckObjVal;}
    bool   isBruteForceEC()       const {return _isBruteForceEC;}
    bool   isBruteForceIncumb()   const {return _isBruteForceIncumb;}

    bool   isSaveInstances()   const {return _isSaveInstances;}
    bool   isSaveNodeTime()    const {return _isSaveNodeTime;}
    bool   isSaveCutPts()      const {return _isSaveCutPts;}

    double rampUpSizeFact()     const {return _rampUpSizeFact;}

    double delta()              const {return _delta;}
    double shrinkDelta()        const {return _shrinkDelta;}
    double maxInterval()        const {return _maxInterval;}

    bool   isPrintBBdetails()     const {return _isPrintBBdetails;}

    int    fixedSizeBin()       const {return _fixedSizeBin;}

    // void   setPrintBBdetails(bool isPrint) {_printBBdetails = isPrint;}

  private:

    bool   _isPebblRMA;            // solve RMA using PEBBL
    //bool   _isGreedyRMA;           // solve the greedy RMA

    // for non-strong branching ...
    bool   _isBinarySearchCutVal;  // an option for binary-sarching cutpoint
    double _perCachedCutPts;     // check only stored cuts points which is x % of total cut points
    double _perLimitAttrib;      // percentages of features to check

    bool   _isRandSeed;            // random seed for tied solution or bound
    bool   _isInitGuess;	         // compute an initial incumbent
    int    _branchSelection;     // random, first, or last one for tied branch
    bool   _isCountingSort;        // use counting sourt (default is bucket sort)

    // for validation
    bool   _isCheckObjVal;	 // check the solution is right or not
    bool   _isBruteForceEC;        // brute force way to create equivalence classes
    bool   _isBruteForceIncumb;    // brute force way to check incumbent in each atrribute

    string  _nonUniformWt;       // use test weight

    // for saving information
    bool   _isSaveInstances;
    bool   _isSaveNodeTime;       // make an output file containing BoundedSP and run time
    bool   _isSaveCutPts;

    // for recursive integerization
    double _delta;               // if the continuous value is less than delta, aggregate to the same integer
    double _shrinkDelta;         // shrink delta
    double _maxInterval;         // the maximum Interval length

    // for fixed size bin integerization
    int    _fixedSizeBin;

    // for printing more details
    bool   _isPrintBBdetails;

    double _rampUpSizeFact;      // TODO: what is this?

  };

} // namespace arguments

#endif // ARG_RMA_h
