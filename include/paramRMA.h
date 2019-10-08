/*
 *  File name: paramRMA.h
 *  Author:    Ai Kagawa
 */

#ifndef PARAM_RMA_h
#define PARAM_RMA_h

#include <pebbl/utilib/ParameterSet.h>
#include <pebbl_config.h>
#include <limits>

using namespace std;


namespace pebblRMA {

 //static double inf = numeric_limits<double>::infinity();
 //static int intInf = numeric_limits<int>::max();

 /////////////////// Parameters for RMA class  ///////////////////
class rmaParams :
  virtual public utilib::ParameterSet,
  virtual public utilib::CommonIO {

public:

  rmaParams();
  ~rmaParams(){};

  bool binarySearchCutVal() const {return _binarySearchCutVal;}
  double perCachedCutPts()  const {return _perCachedCutPts;}
  double perLimitAttrib()   const {return _perLimitAttrib;}

  bool getInitialGuess()    const {return _getInitialGuess;}

  bool checkObjVal()        const {return _checkObjVal;}
  bool bruteForceEC()       const {return _bruteForceEC;}
  bool bruteForceIncumb()   const {return _bruteForceIncumb;}

  bool writingInstances()   const {return _writeInstances;}
  bool writingNodeTime()    const {return _writeNodeTime;}
  bool writingCutPts()      const {return _writeCutPts;}

  bool testWeight()         const {return _testWt;}

  double rampUpSizeFact()   const {return _rampUpSizeFact;}

  bool countingSort()       const {return _countingSort;}
  int branchSelection()     const {return _branchSelection;}

  double delta()            const {return _delta;}
	double shrinkDelta()      const {return _shrinkDelta;}
	double limitInterval()    const {return _limitInterval;}

  int fixedSizeBin()        const {return _fixedSizeBin;}

protected:

  // for non-strong branching ...
  bool _binarySearchCutVal;	// an option for binary-sarching cutpoint
  double _perCachedCutPts;	// check only stored cuts points which is x % of total cut points
  double _perLimitAttrib;		// percentages of features to check

  bool _getInitialGuess;		// compute an initial incumbent

  // for validation
  bool _checkObjVal;				// check the solution is right or not
  bool _bruteForceEC; 			// brute force way to create equivalence classes
  bool _bruteForceIncumb;		// brute force way to check incumbent in each atrribute

  // for saving information
  bool _writeInstances;
  bool _writeNodeTime;			// make an output file containing BoundedSP and run time
  bool _writeCutPts;

  bool _testWt;             // use test weight

  double _rampUpSizeFact;   // TODO: what is this?

  bool _countingSort;       // use counting sourt (default is bucket sort)
  int _branchSelection;

  // for recursive integerization
  double _delta;            // if the continuous value is less than delta, aggregate to the same integer
  double _shrinkDelta;      // shrink delta
  double _limitInterval;    // the maximum Interval length

  // for fixed size bin integerization
	int _fixedSizeBin;

 };

} // namespace PARAM_RMA_h

 #endif
