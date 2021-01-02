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


  // greedy RMA class compute the greedy RMA solution using Kadane's algorithm
  class GreedyRMA {

  public:

    // constructor
    GreedyRMA(ArgRMA *args_, DataRMA *data_) : args(args_), data(data_) { }

    // run the greedy range search alogrithm
    void   runGreedyRangeSearch();

    // search oprimal range for minimum objective value
    void   searchMinOptRange();

    // search oprimal range for maximum objective value
    void   searchMaxOptRange();

    // choose the best solution by comparing the min and max versions
    void   chooseMinOrMaxRange();

    //TODO: combine min max range search together

    // set the current objective value, lower and upper bounds
    // by running the min ver. of Kadane's algorithm for attribute j
    void   runMinKadane(const unsigned int &j);

    // set the current objective value, lower and upper bounds
    // by running the max ver. of Kadane's algorithm for attribute j
    void   runMaxKadane(const unsigned int &j);

    // set current optimal objective value, lower bound, upper bound for attribute j
    // for the min ver.
    void   setOptMin(const unsigned int &j);

    // set current optimal objective value, lower bound, upper bound for attribute j
    // for the max ver.
    void   setOptMax(const unsigned int &j);

    // set total weights for each value for attribut j
    void   setVecValueWeight(const unsigned int &j);

    // drop observations which are not covered for attribute j's lower and upper bound
    void   dropObsNotCovered(const unsigned int &j,
                             const unsigned int &lower,
                             const unsigned int &upper);

    // print Greedy RMA solution
    void   printSolution();

    /****************** utility functions ********************/

    // reset variables for the minimum optimal range search
    void resetMinOptRange();

    // reset variables for the maximum optimal range search
    void resetMaxOptRange();

    // whether or not to choose the max solution
    // when min and max objective value are the same
    // (break the tie using a fair probability)
    bool  isChooseMaxWhenTiedMinMax();

    // whether or not to update the optimal solution
    // when there are tied solution for min or max versions
    // (break the tie using a fair probability based on # of tied solutions)
    bool  isUpdateOptSol(const int &numTiedSol);

    /****************** get functions ************************/
    bool   isPostObjVal() const { return isPosIncumb; }
    double getObjVal()    const { return optObjVal; }
    vector<unsigned int> getLowerBounds() const { return vecLower; }
    vector<unsigned int> getUpperBounds() const { return vecUpper; }

    // void   setInit1DRules();
    // void   set1DOptRange(const int& j);

private:

    // whether or not the objective value is coming from the max or positive value
    bool         isPosIncumb;

    // whether or not the algorithm found a new restructing range
    // which can improve the objective value
    bool         isFondNewBox;

    // # of negative and positive solution which are tied
    unsigned int numNegTiedSols, numPosTiedSols;

    // current lower and upper bounds for a given attribute
    unsigned int curLower, curUpper;

    // current objective value
    double       curMinObjVal, curMaxObjVal;

    // optimal objective value, min ver. and max ver. of the objective values
    double       optObjVal, minObjVal, maxObjVal;

    // current optimal attribute
    unsigned int optAttrib;

    // previously selected or restricted attribute
    unsigned int prevAttrib;

    // optimal lower or upper bound for the current optimal attribute
    unsigned int optLower, optUpper;

    // a vector of the covered observation indices
    vector<unsigned int> vecCvdObsIdx;

    // vectors of lower and upper bounds of the box
    vector<unsigned int> vecLower, vecUpper;

    // vector of lower and upper bounds of the box for the max version
    vector<unsigned int> vecLowerMax, vecUpperMax;

    // vector of lower and upper bounds of the box for the min version
    vector<unsigned int> vecLowerMin, vecUpperMin;

    // a vector of total weights for each value of a given attribute
    vector<double>       vecValueWeight;

    ArgRMA*   args;  // argument RMA class
    DataRMA*  data;  // data RMA class
    Time      ts;    // Time class

  }; // end GreedyRMA class


} // namespace greedyRMA

#endif  // GREEDY_RMA_h
