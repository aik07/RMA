/**********************************************************
 * File name:   argRMA.h
 * Author:      Ai Kagawa
 * Description: a header file for RMA argument class
**********************************************************/

#ifndef ARG_RMA_h
#define ARG_RMA_h

#include <limits>

#include <pebbl_config.h>
#include <pebbl/utilib/ParameterSet.h>
#include <pebbl/utilib/ParameterList.h>
#include <pebbl/utilib/CommonIO.h>
#include <pebbl/utilib/memdebug.h>
#include <pebbl/utilib/seconds.h>
#include <pebbl/bb/pebblParams.h>
#include <pebbl/pbb/parPebblParams.h>

using namespace std;


namespace arguments {

 static double inf = numeric_limits<double>::infinity();
 //static int intInf = numeric_limits<int>::max();

 /////////////////// Parameters for RMA class  ///////////////////
class ArgRMA :
  virtual public utilib::ParameterSet,
  virtual public utilib::CommonIO
  {

public:

  ArgRMA();
  virtual ~ArgRMA(){};

  ////////////////////// parameters //////////////////////////////

  bool   binarySearchCutVal() const {return _binarySearchCutVal;}
  double perCachedCutPts()    const {return _perCachedCutPts;}
  double perLimitAttrib()     const {return _perLimitAttrib;}

  bool   randSeed()           const {return _randSeed;}
  bool   initGuess()          const {return _initGuess;}
  int    branchSelection()    const {return _branchSelection;}
  bool   countingSort()       const {return _countingSort;}

  bool   testWeight()         const {return _testWt;}

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

  int    fixedSizeBin()       const {return _fixedSizeBin;}

protected:

  // for non-strong branching ...
  bool   _binarySearchCutVal;	// an option for binary-sarching cutpoint
  double _perCachedCutPts;	  // check only stored cuts points which is x % of total cut points
  double _perLimitAttrib;	  	// percentages of features to check

  bool   _randSeed;     // random seed for tied solution or bound
  bool   _initGuess;	       	// compute an initial incumbent
  int    _branchSelection;    // random, first, or last one for tied branch
  bool   _countingSort;       // use counting sourt (default is bucket sort)

  // for validation
  bool   _checkObjVal;				// check the solution is right or not
  bool   _bruteForceEC; 			// brute force way to create equivalence classes
  bool   _bruteForceIncumb;		// brute force way to check incumbent in each atrribute

  bool   _testWt;             // use test weight

  // for saving information
  bool   _writeInstances;
  bool   _writeNodeTime;			// make an output file containing BoundedSP and run time
  bool   _writeCutPts;

  // for recursive integerization
  double _delta;              // if the continuous value is less than delta, aggregate to the same integer
  double _shrinkDelta;        // shrink delta
  double _maxInterval;        // the maximum Interval length

  // for fixed size bin integerization
	int    _fixedSizeBin;

  double _rampUpSizeFact;     // TODO: what is this?

 };


class Arguments : virtual public ArgRMA,
                  virtual public pebbl::pebblParams,
                  virtual public pebbl::parallelPebblParams
 {

  Arguments(int& argc, char**& argv) :
    parameters_registered(false),
    min_num_required_args(0) {
      setup(argc, argv);
    }

public:

  virtual bool   setup(int& argc, char**& argv);

  // Parameter-related methods
  virtual void   write_usage_info(char const* progName, std::ostream& os) const;
  virtual void   writeCommandUsage(char const* progName, std::ostream& os) const;
  virtual bool   processParameters(int& argc, char**& argv,
                          unsigned int min_num_required_args);

  /// Register the parameters into a ParameterList object
  virtual void   register_parameters() { plist.register_parameters(*this); }

  /// Check parameters for setup problems and perform debugging I/O
  virtual bool   checkParameters(char const* progName = "");

  virtual bool   setupProblem(int argc, char** argv) { true; }
  virtual void   setName(const char* cname);


 //////////////////////////////////////////////////////////////////
 ParameterList plist;
 bool          parameters_registered;
 string        problemName;
 string        solver_name;
 unsigned int  min_num_required_args;

};

} // namespace arguments

 #endif // ARG_RMA_h
