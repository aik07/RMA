/**********************************************************
 *  File name:   serRMA.cpp
 *  Author:      Ai Kagawa
 *  Description: a header file for the serial RMA solver using PEBBL
 **********************************************************/


#ifndef pebbl_rma_h
#define pebbl_rma_h

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <iostream>
#include <mpi.h>
#include <vector>

#include <pebbl/bb/branching.h>
#include <pebbl/misc/chunkAlloc.h>
#include <pebbl/utilib/ParameterSet.h>
#include <pebbl/utilib/memdebug.h>
#include <pebbl/utilib/seconds.h>
#include <pebbl_config.h>

#ifdef ACRO_HAVE_MPI
#include <pebbl/pbb/parBranching.h>
#include <pebbl/utilib/PackBuf.h>
#endif

#include "Time.h"
#include "argRMA.h"
#include "dataRMA.h"
#include "utilRMA.h"

using namespace std;
using namespace pebbl;
using namespace utilib;


namespace pebblRMA {

// struct CutPt stores chosen chached cut-point (attribute, cut-value)
struct CutPt {
  unsigned int j, v;
  bool operator<(const CutPt &rhs) const { return j < rhs.j; }
};

/******************************************************************************/
// auxiliary classes for choosing branching variables
class branchItem {

public:
  double roundedBound, exactBound;
  int whichChild;

  branchItem()
      : roundedBound(1.0), exactBound(1.0),
        whichChild(-1){}; // , arrayPosition(-1)

  branchItem(branchItem &toCopy)
      : roundedBound(toCopy.roundedBound), exactBound(toCopy.exactBound),
        whichChild(toCopy.whichChild){};

  void set(double bound, double roundQuantum);
};

/******************************************************************************/
//  The branching choice class...
class branchChoice {

public:

  branchChoice();
  branchChoice(double a, double b, double c, int cut, int j);

  void setBounds(double a, double b, double c, int cut, int j);
  void sortBounds();
  bool operator<(const branchChoice &other) const;
  bool operator==(const branchChoice &other) const;

#ifdef ACRO_HAVE_MPI
  static void setupMPI();
  static void freeMPI();
  static MPI_Datatype mpiType;
  static MPI_Op       mpiCombiner;
  static MPI_Op       mpiBranchSelection;

protected:
  static void setupMPIDatum(void *address, MPI_Datatype thisType,
                            MPI_Datatype *type, MPI_Aint base, MPI_Aint *disp,
                            int *blocklen, unsigned int j);
#endif // ACRO_HAVE_MPI

protected:
  void possibleSwap(size_type i1, size_type i2);

public:
  branchItem branch[3];
  unsigned int branchVar, cutVal;
  unsigned int numTiedSols;
};

#ifdef ACRO_HAVE_MPI
void branchChoiceCombiner(void *invec, void *inoutvec, int *len,
                          MPI_Datatype *datatype);
void branchChoiceRand(branchChoice *inBranch, branchChoice *outBranch,
                      int *len, MPI_Datatype *datatype);
#endif

/******************************************************************************/
// CutPt class to plot cut point in order
class CutPtOrder {
public:
  CutPtOrder(){};
  CutPtOrder(unsigned int _order, unsigned int _j, unsigned int _v) : order(_order), j(_j), v(_v){};
  ~CutPtOrder(){};
  void setCutPt(CutPtOrder cp) {
    order = cp.order;
    j = cp.j;
    v = cp.v;
  }
  unsigned int order, j, v;
};

/******************************************************************************/
// Equivalcnece Class
class EquivClass {
public:
  EquivClass() : obsIdx(-1), wt(0.0) {}
  EquivClass(const int &obs, const double &wt_) : obsIdx(obs), wt(wt_) {}
  ~EquivClass() {}

  // merging aother equivalence class to this one,
  // just add the other equivalence class's weight to this one
  void addEC(const EquivClass &ec) { wt += ec.wt; }

  // add this observation index (obs), and weight to this equivalence class
  void addObsWt(const int &obs, const double &weight) {
    if (obsIdx == -1)  // if this equivalence class has no representative observation
      obsIdx = obs;    // save this observation as the representative observation
    wt += weight;  // adding up all weight in the same class
  }

  int    getObs() const { return obsIdx; } // returns the obervation index
  double getWt() const { return wt; }   // get the weight of this equiv class

  int write(ostream &os) const {
    os << obsIdx << " : " << wt;
    return 0;
  }

private:
  int    obsIdx;  // representive observation idex
  double wt;      // the total weight of thie equivalence class
}; // end EquivClass


// Shortcut operators for reading writing RMA items to/from streams
// Forward declarations...
class RMA;
class RMASub;


/******************************************************************************/
//  RMA solution class
// (Note: read PEBBL user guide to understand how to use this class)
class rmaSolution : virtual public solution {

public:
  rmaSolution(){};
  rmaSolution(RMA *global_);
  rmaSolution(rmaSolution *toCopy);
  virtual ~rmaSolution() {}

  void reset(const unsigned int numAttrib, vector<unsigned int> &vecNumDistVals);
  void setData(const bool isPosIncumb, const double objVal,
               vector<unsigned int> &a, vector<unsigned int> &b);

  solution *blankClone() { return new rmaSolution(this); }

  void foundRMASolution(syncType sync = notSynchronous);
  void fileCutPts(RMA *global_);
  void copy(rmaSolution *toCopy);

  virtual void printContents(ostream &s);
  void const printSolution();

  void checkObjValue();

  // check the objective value computed by PEBBL RMA is right
  // using the lower and upper bounds (A and B),
  // the covered observation index list (coveredObs),
  // and the sorted equivalence class index list (sortedECidx)
  void checkObjValue1(vector<unsigned int> &A, vector<unsigned int> &B,
                      vector<unsigned int> &coveredObs,
                      vector<unsigned int> &sortedECidx);

#ifdef ACRO_HAVE_MPI
  void packContents(PackBuffer &outBuf) const;
  void unpackContents(UnPackBuffer &inBuf);
  int  maxContentsBufSize();
#endif

  vector<unsigned int> a, b;   // a and b are the lower and upper bound vectors
  bool isPosIncumb;            // whether or not it is positive incumbent

protected:
  RMA *global;
  double sequenceData();
  size_type sequenceLength() {
    return a.size() + b.size() + sizeof(isPosIncumb);
  }
};


/******************************************************************************/
//  RMA branching class
class RMA : virtual public branching, public ArgRMA {

  // friend class LPB;

public:
  RMA();          // constructor
  RMA(pebblParams *param) ;
  virtual ~RMA(); // { workingSol.decrementRefs(); }		// Destructor

  // copy the RMA argument classe, ArgRMA clsss
  void setParameters(arg::ArgRMA *args_) { args = args_; }

  // copy the RMA data class, DataRMA class
  void setData(data::DataRMA *data_)     { data = data_; }

  // PEBBL functions
  // (Note: read PEBBL user guide to understand how to use these functions)
  bool setupProblem(int &argc, char **&argv) { return true; } // read data file
  branchSub *blankSub();
  solution *initialGuess();
  bool haveIncumbentHeuristic() { return true; }

  // set the initial guess
  void setInitialGuess(bool isPosIncumb, double maxObjValue,
                       vector<unsigned int> L, vector<unsigned int> U);

  // set sorted observation idex
  void setSortedObsIdx(vector<unsigned int> &train) { sortedObsIdx = train; }

  // set cached cut points
  void setCachedCutPts(const unsigned int &j, const unsigned int &v);

  // double getWeight(double pred, set<int> CovgIdx);

  void setWeight(vector<double> wt, vector<unsigned int> train);
  void getPosCovg(set<int> &output, rmaSolution *);
  void getNegCovg(set<int> &output, rmaSolution *);

  // write data to a file, including weights, to a file that
  // can be read by setupProble (added by JE)
  void writeWeightedData(ostream &os);
  void writeInstanceToFile(const int &iterNum);

  // write data to a file, including B&B notes and CPU time,
  // to a file that can be read by setupProble
  void writeStatData(ostream &os);
  void writeStatDataToFile(const int &iterNum);

  bool verifyLog() { return _verifyLog; }
  ostream &verifyLogFile() { return *_vlFile; }

  /*
    void printSolutionTime() const {
    ucout << "ERMA Solution: " << incumbentValue
    << "\tCPU time: "    << searchTime << "\n";
    }
  */
  vector<unsigned int> sortedObsIdx;      // store sorted observations

  vector<vector<CutPtOrder>> CutPtOrders; // to plot cut points

  rmaSolution workingSol;  // RMA solution class object for the working solution
  rmaSolution *guess;      // RMA solution class object for the initial guess

  // for cut-point caching
  unsigned int numCC_SP; // # of subproblems using cutpoint caching

  // map to store already chosen cut points in another branches
  multimap<unsigned int, unsigned int> mmapCachedCutPts;

  bool    _verifyLog;
  ostream *_vlFile;

  data::DataRMA *data;
  arg::ArgRMA   *args;

}; // end RMA class
// ************************************************************************

inline void rmaSolution::foundRMASolution(syncType sync) {
  global->foundSolution(new rmaSolution(this), sync);
  // fileCutPts(global);
};

//******************************************************************************
//  RMA branchSub class
class RMASub : virtual public branchSub {

public:
  RMASub() : globalPtr(NULL){}; // A constructor for a subproblem
  virtual ~RMASub(){};          // A virtual destructor for a subproblem


  // PEBBL functions
  // (Note: read PEBBL user guide to understand how to use these functions)

  /// Return a pointer to the base class of the global branching object
  branching *bGlobal() const { return global(); }

  rmaSolution *workingSol() { return &(globalPtr->workingSol); }

  /// Return a pointer to the global branching object
  inline RMA *global() const { return globalPtr; }

  void setGlobalInfo(RMA *glbl) { globalPtr = glbl; }

  void RMASubFromRMA(RMA *master);
  void RMASubAsChildOf(RMASub *parent, int whichChild);

  // Initialize this subproblem to be the root of the branching tree
  void setRootComputation();

  void boundComputation(double *controlParam);

  /// Create a child subproblem of the current subproblem
  virtual branchSub *makeChild(int whichChild);

  // if it returns true, the computed bound is exact and don't need to separate
  // terminal node of Branch and Bound tree
  bool candidateSolution();

  solution *extractSolution() { return new rmaSolution(workingSol()); }

  virtual int splitComputation();

  // void incumbentHeuristic();
  void foundRMASolution(syncType sync = notSynchronous) {
    workingSol()->foundRMASolution(sync);
  }

  //**************** helper functions (start) ******************************

  // different ways to branch
  void strongBranching();  // strong branching
  void cachedBranching();  // branching by the cutpoint caching method
  void binaryBranching();  // branching by the binary search method
  void hybridBranching();  // hybrid branching method

  void branchingProcess(const unsigned int &j, const unsigned int &v);

  void setNumLiveCutPts();         // set # of cached cut points
  int  getNumLiveCachedCutPts();   // get # of cached cutpoints
  void cutpointCaching();          // cutpoint caching algorithm
  void sortCachedCutPtByAttrib();  // sort cached cutpoints by attribute

  void countingSortObs(const int& j);         // counting sort observations
  void countingSortEC(const unsigned int &j); // counting sort equivalnce classes
  void bucketSortObs(const unsigned int &j);  // bucket sort observations
  void bucketSortEC(const unsigned int &j);   // bucekt sort equivalence classes

  // functions to copute incumbent value
  void   compIncumbent(const unsigned int &j);  // compute the incumbent
  void   chooseMinOrMaxRange();                 // choose the min or max sol
  double runMinKadane(const unsigned int &j);   // run min Kadane's algorithm
  double runMaxKadane(const unsigned int &j);   // run max Kadane's algorithm
  void   setOptMin(const unsigned int &j);
  void   setOptMax(const unsigned int &j);
  double getObjValue(const unsigned int &j, const unsigned int &v);

  // fuctions for tree rotations to compute bounds
  double getBoundMerge() const;   // get the bound when merging equivalence classes
  double getBoundDrop()  const;   // get the bound when dropping equivalence classes
  void   setInitialEquivClass();  // set the initial equivalence class

  // merge the equivalence class for for attribute j,
  // al, au, bl, bu bounds
  void   mergeEquivClass(const unsigned int &j,
                         const unsigned int &al_, const unsigned int &au_,
                         const unsigned int &bl_, const unsigned int &bu_);

  // drop the equivalence class which are not covered
  // by al and bu bounds of attribute j
  void   dropEquivClass(const unsigned int &j,
                        const unsigned int &al_, const unsigned int &bu_);

  // check whether or not the observation 1 and observation 2 are in the same classe
  // based on au and bl bounds of attribte j
  bool   isInSameClass(const unsigned int &obs1, const unsigned int &obs2,
                       const unsigned int &j,
                       const unsigned int &au_, const unsigned int &bl_);

  // construct the equivalence class for for attribute j,
  // au and bu bounds using a brute force algorithm (for compareing)
  void   setEquivClassBF(const unsigned int &j,
                         const unsigned int &au_, const unsigned int &bl_);

  // functions for printing
  // print the current subproblems by chainging current attribute's bounds
  void   printSP(const unsigned int &j,
                 const unsigned int &al, const unsigned int &au,
                 const unsigned int &bl, const unsigned int &bu) const;

  // print the current bounds
  void   printCurrentBounds();

  // TODO: what is this for?
  void   printBounds(vector<double> Bounds, vector<int> Order,
                     const int &j) const;

  void setCutPts(); // TODO: What is this for?????

  //**************************  helper functions (end) ******************************

protected:
  RMA *globalPtr; // A pointer to the global branching object
  // inline double getObjectiveVal() const {return abs(posCovgWt-negCovgWt); };
  // inline unsigned int numObs()           { return global()->data->numTrainObs; };
  // inline unsigned int numDistObs()       { return global()->data->numTrainObs; };
  // TODO: why there are numObs numDistObs?

  inline unsigned int numAttrib()        { return global()->data->numAttrib; };

  // inline unsigned int idxTrain(const int &i) {
  //         return global()->data->vecTrainObsIdx[i];
  // };
  //
  // inline unsigned int idxObsEC(const int &i) {
  //         return vecEquivClass[sortedECidx[i]].getObs();
  // };

  inline vector<unsigned int> vecNumDistVals() {
    return global()->data->vecNumDistVals;
  };

public:

  vector<unsigned int> al, au, bl, bu; // lower and upper bound for a and b vectors size of N features

  unsigned int curObs;      // current observation indx
  unsigned int aj, bj;      // lower and upper bound (a, b) for attribute j
  unsigned int numTiedSols; // # of tied solutions

  vector<unsigned int> coveredObs; // observations which are covered in this subproblem
                          // (al<= feat <= bu)
  vector<unsigned int> sortedECidx; // equivalcnece class which are covered in this subproblem
  vector<unsigned int> sortedECidx1; // equivalence class which are covered in child

  vector<EquivClass> vecEquivClass;  // initial equivalence class
  vector<EquivClass> vecEquivClass1; // merged equivalence class

  vector<double> vecBounds;          // bounds for 2 or 3 childrens' bounds
  branchChoice   _branchChoice;

  vector<CutPt> cachedCutPts, sortedCachedCutPts;

  // variables for incumbent computations
  unsigned int numPosTiedSols, numNegTiedSols;
  double       tmpMin, tmpMax;
  double       minVal, maxVal;
  unsigned int optMinLower, optMinUpper;
  unsigned int optMaxLower, optMaxUpper;
  unsigned int optMinAttrib, optMaxAttrib;
  double       rand_num;

  // TODO: What are those?
  unsigned int numRestAttrib;
  deque<bool>  deqRestAttrib;

  // TODO: what are these ????
  vector<unsigned int> listExcluded, excCutFeat, excCutVal; // store excluded cut-points

}; // end RMA class


// Now we have enough information to define...
inline branchSub *RMA::blankSub() {
  RMASub *temp = new RMASub();
  temp->RMASubFromRMA(this);
  return temp;
} // end brancchSub function

} // namespace pebblRMA


// operator overloading functions to output branching choice info
ostream &operator<<(ostream &os, pebblRMA::branchChoice &bc);

// operator overloading functions to output euivalence class info
ostream &operator<<(ostream &os, pebblRMA::EquivClass &obj);

#endif
