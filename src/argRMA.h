/**********************************************************
 * Author:      Ai Kagawa
 * Description: a source file for ArgRMA class which contains
 *              all parameters for RMA
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
// #include <pebbl/bb/branching.h>
#include <pebbl/bb/pebblParams.h>
#include <pebbl/pbb/parPebblParams.h>

#include "utility.h"

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
    bool   isGreedyRMA()          const {return _isGreedyRMA;}

    bool   isBinarySearchCutVal() const {return _isBinarySearchCutVal;}
    double fracCachedCutPts()      const {return _fracCachedCutPts;}
    double fracLimitAttrib()       const {return _fracLimitAttrib;}

    bool   isRandSeed()           const {return _isRandSeed;}
    bool   isInitGuess()          const {return _isInitGuess;}
    int    branchSelection()      const {return _branchSelection;}
    bool   isCountingSort()       const {return _isCountingSort;}

    string nonUniformWt()         const {return _nonUniformWt;}

    bool   isMinMaxGreedy()       const {return _isMinMaxGreedy;}

    bool   isCheckObjVal()        const {return _isCheckObjVal;}
    bool   isBruteForceEC()       const {return _isBruteForceEC;}
    bool   isBruteForceIncumb()   const {return _isBruteForceIncumb;}

    bool   isSaveInstances()      const {return _isSaveInstances;}
    bool   isSaveNodeTime()       const {return _isSaveNodeTime;}
    bool   isSaveCutPts()         const {return _isSaveCutPts;}

    double rampUpSizeFact()       const {return _rampUpSizeFact;}

    double delta()                const {return _delta;}
    double shrinkDelta()          const {return _shrinkDelta;}
    double maxInterval()          const {return _maxInterval;}

    bool   isPrintBBdetails()     const {return _isPrintBBdetails;}

    unsigned int getNumFixedSizeBins() const {return _numFixedSizeBins;}

  protected:

    bool   _isPebblRMA;            // solve RMA using PEBBL
    bool   _isGreedyRMA;           // solve the greedy RMA

    // for non-strong branching ...
    bool   _isBinarySearchCutVal;  // an option for binary-sarching cutpoint

    // set the threthold of cached cut points to be the fraction of
    // the total cutpoints
    double _fracCachedCutPts;

    // the faction of attributes to check
    double _fracLimitAttrib;

    bool   _isRandSeed;            // random seed for tied solution or bound
    bool   _isInitGuess;	         // compute an initial incumbent
    int    _branchSelection;       // random, first, or last one for tied branch
    bool   _isCountingSort;        // use counting sourt (default is bucket sort)

    // for validation
    bool   _isCheckObjVal;	       // check the solution is right or not
    bool   _isBruteForceEC;        // brute force way to create equivalence classes
    bool   _isBruteForceIncumb;    // brute force way to check incumbent in each atrribute

    string  _nonUniformWt;         // use test weight

    // for recursive integerization
    double _delta;               // if the continuous value is less than delta, aggregate to the same integer
    double _shrinkDelta;         // shrink delta
    double _maxInterval;         // the maximum Interval length

    // for fixed size bin integerization
    unsigned int _numFixedSizeBins;

    // whether or not to sovel min-then-max in the greedy RMA
    bool    _isMinMaxGreedy;

    // for printing more details
    bool   _isPrintBBdetails;

    // for saving information
    bool   _isSaveInstances;
    bool   _isSaveNodeTime;    // make an output file containing BoundedSP and run time
    bool   _isSaveCutPts;

    double _rampUpSizeFact;    // TODO: what is this?

  }; // end ArgRMA class

} // namespace arguments

#endif // ARG_RMA_h
