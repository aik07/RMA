/*
*  Author:      Ai Kagawa
*  Description: a serial greedy rectangular maximum agreement problem solver
*/


#ifndef GREEDY_RMA_h
#define GREEDY_RMA_h

#include <ostream>
#include <vector>

#include "Time.h"
#include "argRMA.h"
#include "dataRMA.h"
#include "utility.h"

using namespace std;
using namespace arg;
using namespace data;


namespace greedyRMA {

  enum OptimizationDirection {MIN, MAX};


  // greedy RMA class compute the greedy RMA solution using Kadane's algorithm
  class GreedyRMA {

  public:

    // constructor
    GreedyRMA(ArgRMA *args_, DataRMA *data_) : args(args_), data(data_) { }

    // reset Greedy RMA solutions
    void reset();

    // run the greedy range search alogrithm
    void   runGreedyRangeSearch();

    // search oprimal range for minimum or maximum objective value
    void   searchGreedyRange(const bool &isMax);

    // reset variables for the optimal range search
    void   resetBestRange();

    // set the current objective value, lower and upper bounds
    // by running the min  or max ver. of Kadane's algorithm for attribute j
    void   runKadaneAlgo(const bool &isMax, const unsigned int &j);

    // set current optimal objective value, lower bound, upper bound for attribute j
    // for the max ver.
    void   setBestRange(const unsigned int &j);

    // whether or not to update the optimal solution
    // when there are tied solution for min or max versions
    // (break the tie using a fair probability based on # of tied solutions)
    bool  isUpdateBestSol(const int &numTiedSol);

    // drop observations which are not covered for attribute j's lower and upper bound
    void   dropObsNotCovered(const unsigned int &j,
                             const unsigned int &lower,
                             const unsigned int &upper);

    // set total weights for each value for attribut j
    void   setVecValueWeight(const unsigned int &j);

    void   updateBestBounds() {
      copy(vecLowerWorking.begin(), vecLowerWorking.end(), vecLower.begin());
      copy(vecUpperWorking.begin(), vecUpperWorking.end(), vecUpper.begin());
    }

    // print Greedy RMA solution
    void   printSolution();

    /****************** get functions ************************/
    bool                 isPostObjVal()   const { return isPosObjVal; }
    double               getObjVal()      const { return bestObjVal; }
    vector<unsigned int> getLowerBounds() const { return vecLower; }
    vector<unsigned int> getUpperBounds() const { return vecUpper; }

    // void   setInit1DRules();
    // void   set1DOptRange(const int& j);

private:

    // whether or not the objective value is coming from the max or positive value
    bool         isPosObjVal;

    bool         isImprovedSol;

    // vectors of lower and upper bounds of the box
    vector<unsigned int> vecLower, vecUpper;

    // # of tied solutions for the current best objective value
    unsigned int numTiedSols;

    // current lower and upper bounds for a given attribute
    unsigned int curLower, curUpper;

    // current objective value
    double       curObjVal;

    // optimal objective value, min ver. and max ver. of the objective values
    double       bestObjVal;

    // current optimal attribute
    unsigned int bestAttrib;

    // previously selected or restricted attribute
    unsigned int prevBestAttrib;

    // optimal lower or upper bound for the current optimal attribute
    unsigned int bestLower, bestUpper;

    // vector of lower and upper bounds of the box for the min or max version
    vector<unsigned int> vecLowerWorking, vecUpperWorking;

    // a vector of total weights for each value of a given attribute
    vector<double>       vecValueWeight;

    // a vector of the covered observation indices
    vector<unsigned int> vecCvdObsIdx;

    ArgRMA*   args;  // argument RMA class
    DataRMA*  data;  // data RMA class
    Time      ts;    // Time class

  }; // end GreedyRMA class


} // namespace greedyRMA

#endif  // GREEDY_RMA_h
