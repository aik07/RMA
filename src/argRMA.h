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

  // static double inf = numeric_limits<double>::infinity();
  // static int intInf = numeric_limits<int>::max();

  /////////////////// Parameters for RMA class  ///////////////////
  class ArgRMA : virtual public ParameterSet, virtual public CommonIO {

  public:

    ArgRMA();
    virtual ~ArgRMA(){};

    ////////////////////// parameters //////////////////////////////

    bool   exactRMA()           const {return _exactRMA;}

    bool   binarySearchCutVal() const {return _binarySearchCutVal;}
    double perCachedCutPts()    const {return _perCachedCutPts;}
    double perLimitAttrib()     const {return _perLimitAttrib;}

    bool   randSeed()           const {return _randSeed;}
    bool   initGuess()          const {return _initGuess;}
    int    branchSelection()    const {return _branchSelection;}
    bool   countingSort()       const {return _countingSort;}

    string nonUniformWt()       const {return _nonUniformWt;}

    bool   checkObjVal()        const {return _checkObjVal;}
    bool   bruteForceEC()       const {return _bruteForceEC;}
    bool   bruteForceIncumb()   const {return _bruteForceIncumb;}

    bool   writingInstances()   const {return _writeInstances;}
    bool   writingNodeTime()    const {return _writeNodeTime;}
    bool   writingCutPts()      const {return _writeCutPts;}

    double rampUpSizeFact()     const {return _rampUpSizeFact;}

    double delta()              const {return _delta;}
    double shrinkDelta()        const {return _shrinkDelta;}
    double maxInterval()        const {return _maxInterval;}

    bool   printBBdetails()     const {return _printBBdetails;}

    int    fixedSizeBin()       const {return _fixedSizeBin;}

    // void   setPrintBBdetails(bool isPrint) {_printBBdetails = isPrint;}

  private:

    bool   _exactRMA;            // solve exactRMA
    bool   _greedyRMA;            // solve exactRMA

    // for non-strong branching ...
    bool   _binarySearchCutVal;  // an option for binary-sarching cutpoint
    double _perCachedCutPts;     // check only stored cuts points which is x % of total cut points
    double _perLimitAttrib;      // percentages of features to check

    bool   _randSeed;            // random seed for tied solution or bound
    bool   _initGuess;	         // compute an initial incumbent
    int    _branchSelection;     // random, first, or last one for tied branch
    bool   _countingSort;        // use counting sourt (default is bucket sort)

    // for validation
    bool   _checkObjVal;	 // check the solution is right or not
    bool   _bruteForceEC;        // brute force way to create equivalence classes
    bool   _bruteForceIncumb;    // brute force way to check incumbent in each atrribute

    string  _nonUniformWt;       // use test weight

    // for saving information
    bool   _writeInstances;
    bool   _writeNodeTime;       // make an output file containing BoundedSP and run time
    bool   _writeCutPts;

    // for recursive integerization
    double _delta;               // if the continuous value is less than delta, aggregate to the same integer
    double _shrinkDelta;         // shrink delta
    double _maxInterval;         // the maximum Interval length

    // for fixed size bin integerization
    int    _fixedSizeBin;

    // for printing more details
    bool   _printBBdetails;

    double _rampUpSizeFact;      // TODO: what is this?

  };

} // namespace arguments

#endif // ARG_RMA_h
