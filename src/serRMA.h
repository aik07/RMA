/*
* file name:  serRMA.h
* author:     Ai Kagawa
*
* Solves Rectanglar Maximum Agreement.
*/

#ifndef pebbl_rma_h
#define pebbl_rma_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <stdlib.h>
#include <limits>

#include <pebbl_config.h>
#include <pebbl/bb/branching.h>
#include <pebbl/misc/chunkAlloc.h>
#include <pebbl/utilib/ParameterSet.h>
#include <pebbl/utilib/seconds.h>
#include <pebbl/utilib/exception_mngr.h>
#include <pebbl/utilib/comments.h>

#ifdef ACRO_HAVE_MPI
#include <pebbl/pbb/parBranching.h>
#include <pebbl/utilib/PackBuf.h>
#endif

using namespace utilib;
using namespace std;
using namespace pebbl;

static double inf = numeric_limits<double>::infinity();
static int intInf = numeric_limits<int>::max();

struct IntMinMax { double minOrigVal, maxOrigVal; };

struct Feature {vector<IntMinMax> vecIntMinMax;};

namespace pebblRMA {

  // CutPt class stores chosen chached cut-point (attribute, cut-value)
  struct CutPt {
    int j, v;
    bool operator<(const CutPt& rhs) const {
      return j < rhs.j;
    };
  };
  
  
  // RMAObs class contains info about each observation
  class Data {

  public:
    // Constructors
  Data() : w(0.0) {};
  Data(  const vector<double>& X_, const double & y_ ) :
    X(X_), y(y_) { };

    // Assignment operator
    Data& operator==(const Data& other) {
      X = other.X;
      y = other.y;
      return *this;
    };

    int read(istream& is) { is >> X >> y; return 0; };
    int write(ostream& os) const { os << X << " " << y;	return 0; };

    // Basic properties of the observation
    vector<double> X;	// independent variables
    double y;					// dependent variable
    double w;					// weight of each observation
  };


  // auxiliary classes for choosing branching variables
  class branchItem {

  public:
    double roundedBound, exactBound;
    int whichChild; //, arrayPosition;

  branchItem() : roundedBound(1.0), exactBound(1.0), whichChild(-1) { }; // , arrayPosition(-1)

  branchItem(branchItem& toCopy) :
    roundedBound(toCopy.roundedBound),
      exactBound(toCopy.exactBound),
      whichChild(toCopy.whichChild)
	{ };

    void set(double bound, double roundQuantum);

  };


  //********************************************************************************
  //  The branching choice class...
  class branchChoice {
  public:
    branchItem branch[3];
    int branchVar, cutVal;
    int numTiedSols;
    branchChoice() ;
    branchChoice(double a, double b, double c, int cut, int j);
    void setBounds(double a, double b, double c, int cut, int j);
    void sortBounds();
    bool operator<(const branchChoice& other) const;
    bool operator==(const branchChoice& other) const;

#ifdef ACRO_HAVE_MPI
    static void setupMPI();
    static void freeMPI();
    static MPI_Datatype mpiType;
    static MPI_Op       mpiCombiner;
    static MPI_Op				mpiBranchSelection;
  protected:
    static void setupMPIDatum(void* address, MPI_Datatype thisType,
			      MPI_Datatype* type, MPI_Aint base, MPI_Aint* disp,
			      int* blocklen, int j);
#endif

  protected:
    void possibleSwap(size_type i1, size_type i2);

  };


#ifdef ACRO_HAVE_MPI
  void branchChoiceCombiner(void* invec, void* inoutvec, int* len,
			    MPI_Datatype* datatype) ;
  void branchChoiceRand(branchChoice *inBranch, branchChoice *outBranch,
			int* len, MPI_Datatype *datatype);
#endif




  //********************************************************************************
  // Equivalcnece Class
  class EquivClass {
  public:
  EquivClass():wt(0),obsIdx(-1) {}
  EquivClass(  const int& obs, const double & wt_ ) :
    obsIdx(obs), wt(wt_ ) { }
    ~EquivClass(){}

    void addEC(const EquivClass& ec) { wt += ec.wt;	}

    void addObsWt(const int& obs, const double& weight) {
      if (obsIdx==-1) obsIdx = obs;
      wt+=weight;
    }

    int getObs() const { return obsIdx; } // returns the obervation index
    double getWt() const { return wt;	}   // get the weight of this equiv class

    int write(ostream& os) const { os << obsIdx  << " : " << wt;	return 0; }

  private:
    int obsIdx;
    double wt;
  };


  // Shortcut operators for reading writing RMA items to/from streams
  // Forward declarations...
  class RMA;
  class RMASub;

  //********************************************************************************
  //  The solution class...
  class rmaSolution : public solution
  {
  public:

    rmaSolution(RMA* global_);
    rmaSolution(rmaSolution* toCopy);
    virtual	~rmaSolution() {}

    solution* blankClone() {return new rmaSolution(this);}

    void foundSolution(syncType sync = notSynchronous);
    void fileCutPts(RMA* global_);
    void copy(rmaSolution* toCopy);

    virtual void printContents(ostream& s);
    void const printSolution();
    void checkObjValue();

#ifdef ACRO_HAVE_MPI
    void packContents(PackBuffer& outBuf);
    void unpackContents(UnPackBuffer& inBuf);
    int maxContentsBufSize();
#endif

    vector<int> a, b;
  protected:
    RMA* global;
    virtual double sequenceData();
    size_type sequenceLength() { return a.size()+b.size(); }
  };


  //******************************************************************************
  //  RMA branching class
  class RMA : virtual public branching {

  public:

    RMA();					// constructor
    virtual ~RMA(); // {workingSol.decrementRefs(); }		// Destructor
    bool setupProblem(int& argc,char**& argv);
    branchSub* blankSub();
    solution* initialGuess() ;
    bool haveIncumbentHeuristic() { return true; }

    bool setupIntData(int& argc,char**& argv);
    bool setupOrigData(int& argc,char**& argv);

    void setStdDevX();
    void integerizeData() ;
    void integerizeFixedLengthData();

    void bucketSort(const int& attrib);
    void removeDuplicateObs();

    void setCachedCutPts(const int& j, const int& v);

    // greedy method
    double getMaxRange(const int& j) ;
    double getMinRange(const int& j) ;
    void setObjVec(const int &j) ;
    void dropObsNotCovered(const int &j, const int& lower, const int& upper);

    // get positive or negative observations
    void getPosCovg(set<int> & output, rmaSolution*);
    void getNegCovg(set<int> & output, rmaSolution*);

    // write data to a file, including weights, to a file that
    // can be read by setupProble (added by JE)
    void writeWeightedData(ostream& os);
    void writeInstanceToFile(const int& iterNum);

    // write data to a file, including B&B notes and CPU time,
    // to a file that can be read by setupProble
    void writeStatData(ostream& os);
    void writeStatDataToFile(const int&  iterNum);

    virtual bool verifyLog() {return _verifyLog;}
    ostream& verifyLogFile() { return *_vlFile; };

    void startTime();
    double endTime();

    /************************ get parameters ************************/
    double perCachedCutPts()    const {return _perCachedCutPts;}
    bool   binarySearchCutVal() const {return _binarySearchCutVal;}
    double perLimitAttrib()     const {return _perLimitAttrib;}
    bool   getInitialGuess()    const {return _getInitialGuess;}
    bool   checkObjVal()        const {return _checkObjVal;}
    bool   bruteForceEC()       const {return _bruteForceEC;}
    bool   bruteForceIncumb()   const {return _bruteForceIncumb;}
    bool   writingInstances()   const {return _writeInstances;}
    bool   writingNodeTime()    const {return _writeNodeTime;}
    bool   writingCutPts()      const {return _writeCutPts;}
    bool   testWeight()         const {return _testWt;}
    int    maxBoundedSP()       const {return _maxBoundedSP;}
    double rampUpSizeFact()     const {return _rampUpSizeFact;}
    bool   countingSort()       const {return _countingSort;}
    int    branchSelection()    const {return _branchSelection;}
    double delta()              const {return _delta;}
    double shrinkDelta()        const {return _shrinkDelta;}
    double limitInterval()      const {return _limitInterval;}
    int    fixedSizeBin()       const {return _fixedSizeBin;}

  protected:
    /**************** RMA optional parameters ****************/

    // for non-strong branching ...
    double _perCachedCutPts;		// check only stored cuts points which is x % of total cut points
    bool _binarySearchCutVal;	// binarySearchCutVal
    double _perLimitAttrib;			// percentages of features to check

    bool _getInitialGuess;		// compute an initial incumbent

    bool _checkObjVal;				// check the solution is right in the end
    bool _bruteForceEC; 			// brute force way to create equivalence classes
    bool _bruteForceIncumb;		// brute force way to check incumbent in each atrribute

    bool _writeInstances;
    bool _writeNodeTime;			// make an output file containing BoundedSP and run time
    bool _writeCutPts;

    bool _testWt;

    double _rampUpSizeFact;
    int _maxBoundedSP;				// set a maximum number of bounded subproblems to check

    bool _countingSort;
    int _branchSelection;

    // for recursive integerization
    double _delta;
    double _shrinkDelta;
    double _limitInterval;

    // for fixed size bin integerization
    int _fixedSizeBin;

  public:

    // store observations in original order and observations in sorted order
    vector<Data> intData, origData;

    // contains l_j-1 = (# of distinct value observed in the feature)
    vector<int> distFeat;
    vector<int> sortedObsIdx; 	// store sorted observations

    size_type numObs;				// # of observations
    size_type numAttrib;		// # of attribute
    size_type numDistObs;		// # of distinct observations

    rmaSolution workingSol;
    rmaSolution* guess;

    // for cut-point caching
    int numCC_SP;				// # of subproblems using cutpoint caching
    int numTotalCutPts;    // # of total cutpoints
    multimap<int, int> mmapCachedCutPts; // map to store already chosen cut points in another branches

    bool _verifyLog;
    ostream* _vlFile;

  private:

    bool foundBox;
    vector<int> vecCoveredObs;
    vector<int> L, U, Lmax, Umax, Lmin, Umin;
    vector<double> W;
    int maxL, tmpL, tmpU; double tmpObj;

    vector<double> minX, maxX, sdX;
    vector<Feature> vecFeature;		// contains features original and integeried values

    clock_t timeStart, timeEnd, clockTicksTaken;
    double timeInSeconds;

  }; // end class RMA ************************************************************************


  inline void rmaSolution::foundSolution(syncType sync) {
    global->foundSolution(new rmaSolution(this),sync);
    //fileCutPts(global);
  };


  //******************************************************************************
  //  RMA branchSub class
  class RMASub : virtual public branchSub {

  public:

  RMASub() : globalPtr(NULL) {};	// A constructor for a subproblem
    virtual ~RMASub() {};  // A virtual destructor for a subproblem

    /// Return a pointer to the base class of the global branching object
    branching* bGlobal() const { return global(); };

    rmaSolution* workingSol() { return &(globalPtr->workingSol); };

    /// Return a pointer to the global branching object
    inline RMA* global() const { return globalPtr; };

    void setGlobalInfo(RMA* glbl) {globalPtr = glbl;}

    void RMASubFromRMA(RMA* master);
    void RMASubAsChildOf(RMASub* parent, int whichChild) ;

    // Initialize this subproblem to be the root of the branching tree
    void setRootComputation();

    void boundComputation(double* controlParam);

    virtual int splitComputation() ;

    /// Create a child subproblem of the current subproblem
    virtual branchSub* makeChild(int whichChild) ;

    // if it returns true, the computed bound is exact and don't need to separate
    // terminal node of Branch and Bound tree
    bool candidateSolution();

    solution* extractSolution() { return new rmaSolution(workingSol()); }

    //void incumbentHeuristic();
    void foundSolution(syncType sync = notSynchronous) {
      workingSol()->foundSolution(sync);
    }

    //************************* helper functions (start) *******************************************

    // different ways to branch
    void strongBranching();
    void cachedBranching();
    void binaryBranching();
    void hybridBranching();

    void branchingProcess(const int& j, const int& v) ;

    void setNumLiveCutPts();
    int getNumLiveCachedCutPts();
    void setLiveCachedCutPts();
    void sortCachedCutPtByAttrib() ;

    //void countingSortObs(const int& j) ;
    void countingSortEC(const int& j) ;
    void bucketSortObs(const int& j);
    void bucketSortEC(const int& j);

    void checkIncumbent(const int& j);
    double getMaxRange1(const int& j) ;
    double getMinRange1(const int& j) ;
    double getObjValue(const int& j, const int& v);

    double getBoundMerge() const;
    double getBoundDrop() const;
    void setInitialEquivClass();
    void mergeEquivClass(const int& j, const int& al_, const int& au_,
			 const int& bl_, const int& bu_) ;
    void dropEquivClass(const int& j, const int& al_, const int& bu_);
    bool isInSameClass(const int& obs1, const int& obs2,
		       const int& j, const int& au_, const int& bl_);
    void setEquivClassBF(const int& j, const int& au_, const int& bl_);

    void printSP(const int& j, const int& al, const int& au,
		 const int& bl, const int& bu) const;
    void printCurrentBounds() ;
    void printBounds(vector<double> Bounds, vector<int> Order,
		     const int& j) const;

    //**************************  helper functions (end) ******************************

  protected:
    RMA* globalPtr;  // A pointer to the global branching object
    inline int numObs() { return global()->numObs; };
    inline int numDistObs() { return global()->numDistObs; };
    inline int numAttrib() { return global()->numAttrib; };
    inline vector<int> distFeat() { return global()->distFeat; };

  public:

    vector<int> al, au, bl, bu; // lower and upper bound for a and b vectors size of N features
    int curObs, aj, bj;
    int NumTiedSols;

    vector<int> coveredObs;		// observations which are covered in this subproblem (al<= feat <= bu)
    vector<int> coveredObs1;	// (only used for bruteForceEC) covered observations for each child
    vector<int> sortedECidx;	// equivalcnece class which are covered in this subproblem
    vector<int> sortedECidx1; // equivalence class which are covered in child
    vector< EquivClass > vecEquivClass;	// initial equivalence class
    vector< EquivClass > vecEquivClass1;	// merged equivalence class

    vector<double> vecBounds;	// bounds for 2 or 3 childrens' bounds
    branchChoice _branchChoice;

    vector<CutPt> cachedCutPts,  sortedCachedCutPts;

    //deque<bool> vecParentInd;

  };//******************** class RMASub (end) ********************************

  // Now we have enough information to define...
  inline branchSub* RMA::blankSub() {
    RMASub *temp = new RMASub();
    temp->RMASubFromRMA(this);
    return temp;
  };

} //********************* namespace pebbl ********************************

ostream& operator<<(ostream& os, pebblRMA::branchChoice& bc);
ostream& operator<<(ostream& os, const deque<bool>& v);
ostream& operator<<(ostream& os, pebblRMA::Data& obj);
istream& operator>>(istream& is, pebblRMA::Data& obj);
ostream& operator<<(ostream& os, const vector<int>& v);
ostream& operator<<(ostream& os, pebblRMA::EquivClass& obj);

#endif
