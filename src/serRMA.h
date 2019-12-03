/**********************************************************
*  File name:   serRMA.cpp
*  Author:      Ai Kagawa
*  Description: a header file for the serial RMA solver using PEBBL
**********************************************************/

#ifndef pebbl_rma_h
#define pebbl_rma_h

#include <iostream>
#include <deque>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <mpi.h>

#include <pebbl_config.h>
#include <pebbl/bb/branching.h>
#include <pebbl/misc/chunkAlloc.h>
#include <pebbl/utilib/memdebug.h>
#include <pebbl/utilib/ParameterSet.h>
#include <pebbl/utilib/seconds.h>
//#include <pebbl/misc/fundamentals.h>
//#include <SimpleHashTable.h>

#ifdef ACRO_HAVE_MPI
  #include <pebbl/pbb/parBranching.h>
  #include <pebbl/utilib/PackBuf.h>
#endif

#include "Time.h"
#include "argRMA.h"
#include "dataRMA.h"

using namespace std;
using namespace pebbl;
using namespace utilib;
using namespace arg;
using namespace data;
//using namespace greedyRMA;


namespace pebblRMA {


// struct CutPt stores chosen chached cut-point (attribute, cut-value)
struct CutPt {
  int j, v;
  bool operator<(const CutPt& rhs) const { return j < rhs.j; }
};


// auxiliary classes for choosing branching variables
 class branchItem {

public:
   double roundedBound, exactBound;
   int    whichChild;

   branchItem() : roundedBound(1.0), exactBound(1.0), whichChild(-1) { }; // , arrayPosition(-1)

   branchItem(branchItem& toCopy) : roundedBound(toCopy.roundedBound),
     exactBound(toCopy.exactBound), whichChild(toCopy.whichChild) { };

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
  static MPI_Op	      mpiBranchSelection;
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


// to plot cut point in order
class CutPtOrder {
public:
  CutPtOrder(){};
  CutPtOrder(int _order, int _j, int _v) : order(_order), j(_j), v(_v) {};
  ~CutPtOrder(){};
  void setCutPt(CutPtOrder cp) { order = cp.order; j = cp.j; v = cp.v; }
  int order, j, v;
};


//********************************************************************************
// Equivalcnece Class
class EquivClass {
public:
  EquivClass():wt(0),obsIdx(-1) {}
  EquivClass(const int& obs, const double & wt_): obsIdx(obs), wt(wt_) { }
  ~EquivClass(){}

  void addEC(const EquivClass& ec) { wt += ec.wt; }

  void addObsWt(const int& obs, const double& weight) {
    if (obsIdx==-1) obsIdx = obs;
    wt+=weight;
  }

  int getObs() const { return obsIdx; } // returns the obervation index
  double getWt() const { return wt; }   // get the weight of this equiv class

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
class rmaSolution : virtual public solution {

public:

  rmaSolution(){};
  rmaSolution(RMA* global_);
  rmaSolution(rmaSolution* toCopy);
  virtual ~rmaSolution() {}

  solution* blankClone() {return new rmaSolution(this);}

  void foundSolution(syncType sync = notSynchronous);
  void fileCutPts(RMA* global_);
  void copy(rmaSolution* toCopy);

  virtual void printContents(ostream& s);
  void const printSolution();
  void checkObjValue();
  void checkObjValue1(vector<int> &A, vector<int> &B,
		      vector<int> &coveredObs, vector<int> &sortedECidx);

#ifdef ACRO_HAVE_MPI
  void packContents(PackBuffer& outBuf);
  void unpackContents(UnPackBuffer& inBuf);
  int maxContentsBufSize();
#endif

  vector<int> a, b;
  bool isPosIncumb;

protected:
  RMA*      global;
  double    sequenceData();
  size_type sequenceLength() { return a.size()+b.size()+sizeof(isPosIncumb); }

};


//******************************************************************************
//  RMA branching class
class RMA : virtual public branching, public ArgRMA {

friend class LPB;

public:

  RMA();					// constructor
  virtual ~RMA(); // {workingSol.decrementRefs(); }		// Destructor

  //void setParameters(ArgPMA* args, Data* data, const int& deb_int);
  bool setParameters(ArgRMA* args_) { args=args_; }
  bool setData(Data* data_);

  bool       setupProblem(int& argc,char**& argv) { return true; }
  branchSub* blankSub();
  solution*  initialGuess();
  bool       haveIncumbentHeuristic() { return true; }

  //void setsortedObsIdx();
  void setSortObsNum(vector<int> & train) { sortedObsIdx = train; }
  void setCachedCutPts(const int& j, const int& v);

  //double getWeight(double pred, set<int> CovgIdx);
  void setWeight( vector<double> wt, vector<int> train);
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

  bool     verifyLog()     { return _verifyLog; }
  ostream& verifyLogFile() { return *_vlFile; }

  void   startTime();
  double endTime();

  virtual void printSolutionTime() const {
    ucout << "ERMA Solution: " << incumbentValue
	  << "\tCPU time: "    << searchTime << "\n";
  }

  // contains l_j-1 = (# of distinct value observed in the feature)
  vector<int> distFeat;
  vector<int> sortedObsIdx; 	// store sorted observations

  vector<vector<CutPtOrder> > CutPtOrders;				// to plot cut points

  size_type numObs;				// # of observations
  size_type numAttrib;		// # of attribute
  size_type numDistObs;		// # of distinct observations

  rmaSolution workingSol;
  rmaSolution guess;

  // for cut-point caching
  int                numCC_SP;	       // # of subproblems using cutpoint caching
  int                numTotalCutPts;   // # of total cutpoints
  multimap<int, int> mmapCachedCutPts; // map to store already chosen cut points in another branches

  bool     _verifyLog;
  ostream* _vlFile;

  clock_t timeStart, timeEnd, clockTicksTaken;
  double timeInSeconds;

  Data*      data;
  ArgRMA*    args;
  //GreedyRMA* grma;

}; // end class RMA ************************************************************************


inline void rmaSolution::foundSolution(syncType sync) {
  global->foundSolution(new rmaSolution(this), sync);
  //fileCutPts(global);
};


//******************************************************************************
//  RMA branchSub class
class RMASub : virtual public branchSub {

public:

  RMASub() : globalPtr(NULL) {};	// A constructor for a subproblem
  virtual ~RMASub() {};  // A virtual destructor for a subproblem

  /// Return a pointer to the base class of the global branching object
  branching* bGlobal() const { return global(); }

  rmaSolution* workingSol() { return &(globalPtr->workingSol); }

  /// Return a pointer to the global branching object
  inline RMA* global() const { return globalPtr; }

  void setGlobalInfo(RMA* glbl) {globalPtr = glbl;}

  void RMASubFromRMA(RMA* master);
  void RMASubAsChildOf(RMASub* parent, int whichChild);

  // Initialize this subproblem to be the root of the branching tree
  void setRootComputation();

  void boundComputation(double* controlParam);

  virtual int splitComputation();

  /// Create a child subproblem of the current subproblem
  virtual branchSub* makeChild(int whichChild);

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

  void branchingProcess(const int& j, const int& v);

  void setNumLiveCutPts();
  int  getNumLiveCachedCutPts();
  void cutpointCaching();
  void sortCachedCutPtByAttrib();

  //void countingSortObs(const int& j) ;
  void countingSortEC(const int& j);
  void bucketSortObs(const int& j);
  void bucketSortEC(const int& j);

  // functions to copute incumbent value
  void   compIncumbent(const int& j);
  void   chooseMinOrMaxRange();
  double runMinKadane(const int& j);
  double runMaxKadane(const int& j);
  void   setOptMin(const int& j);
  void   setOptMax(const int& j);
  double getObjValue(const int& j, const int& v);

  // fuctions for tree rotations to compute bounds
  double getBoundMerge() const;
  double getBoundDrop() const;
  void   setInitialEquivClass();
  void   mergeEquivClass(const int& j, const int& al_, const int& au_,
			 const int& bl_, const int& bu_) ;
  void   dropEquivClass(const int& j, const int& al_, const int& bu_);
  bool   isInSameClass(const int& obs1, const int& obs2,
	                     const int& j, const int& au_, const int& bl_);
  void   setEquivClassBF(const int& j, const int& au_, const int& bl_);

  // functions for printing
  void   printSP(const int& j, const int& al, const int& au,
                 const int& bl, const int& bu) const;
  void   printCurrentBounds() ;
  void   printBounds(vector<double> Bounds, vector<int> Order,
		     const int& j) const;

  void setCutPts(); // TODO: What is this for?????

  //**************************  helper functions (end) ******************************

protected:
  RMA* globalPtr;  // A pointer to the global branching object
  //inline double getObjectiveVal() const {return abs(posCovgWt-negCovgWt); };
  inline int         numObs()     { return global()->numObs; };
  inline int         numDistObs() { return global()->numDistObs; };
  inline int         numAttrib()  { return global()->numAttrib; };
  inline vector<int> distFeat()   { return global()->distFeat; };

public:

  vector<int> al, au, bl, bu; // lower and upper bound for a and b vectors size of N features

  int curObs, aj, bj;
  int NumTiedSols;

  vector<int> coveredObs;     // observations which are covered in this subproblem (al<= feat <= bu)
  vector<int> sortedECidx;    // equivalcnece class which are covered in this subproblem
  vector<int> sortedECidx1;   // equivalence class which are covered in child

  vector< EquivClass > vecEquivClass;	// initial equivalence class
  vector< EquivClass > vecEquivClass1;	// merged equivalence class

  vector<double> vecBounds;	// bounds for 2 or 3 childrens' bounds
  branchChoice   _branchChoice;

  vector<CutPt> cachedCutPts, sortedCachedCutPts;

  // variables for incumbent computations
  int    NumPosTiedSols, NumNegTiedSols;
  double tmpMin, tmpMax, minVal, maxVal;
  double optMinLower, optMinUpper, optMaxLower, optMaxUpper;
  int    optMinAttrib, optMaxAttrib;
  double rand_num;

  int         numRestAttrib;
  deque<bool> deqRestAttrib;

  // TODO: what are these ????
  vector<int> listExcluded, excCutFeat, excCutVal; 	// store excluded cut-points

};//******************** class RMASub (end) ********************************

// Now we have enough information to define...
inline branchSub* RMA::blankSub() {
  RMASub *temp = new RMASub();
  temp->RMASubFromRMA(this);
  return temp;
};

} //********************* namespace pebbl ********************************

ostream& operator<<(ostream& os, pebblRMA::branchChoice& bc);
ostream& operator<<(ostream& os, pebblRMA::EquivClass&   obj);

#endif
