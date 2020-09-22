/*
*  File name:   greedyRMA.h
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

  class GreedyRMA {

  public:

    GreedyRMA(ArgRMA *args_, DataRMA *data_) : args(args_), data(data_) { }

    void   runGreedyRangeSearch();

    void   chooseMinOrMaxRange();

    void   getMinOptRange();
    void   getMaxOptRange();

    void   setOptMin(const unsigned int &j);
    void   setOptMax(const unsigned int &j);

    //TODO: combine min max range search together
    double runMinKadane(const unsigned int& j) ;
    double runMaxKadane(const unsigned int& j) ;

    void   setObjVec(const unsigned int &j);

    void   dropObsNotCovered(const unsigned  int &j, const unsigned int& lower,
                             const unsigned int& upper);

    double getObjCovered(const unsigned int& j, const unsigned int& v);

    void   printSolution();

    // void   setInit1DRules();
    // void   set1DOptRange(const int& j);

    // private:

    bool   isPosIncumb;
    bool   foundBox;
    unsigned int    NumNegTiedSols;
    unsigned int    NumPosTiedSols;

    unsigned int    tmpL, tmpU;
    double tmpObj, tmpMin, tmpMax;
    double maxObjValue, minVal, maxVal;

    unsigned int    optAttrib;
    unsigned int    prevAttrib;
    unsigned int    optLower, optUpper;

    unsigned int    obs;
    bool   fondNewBox;

    vector<unsigned int>    vecCoveredObs;
    vector<unsigned int>    L, U;
    vector<unsigned int>    Lmax, Umax;
    vector<unsigned int>    Lmin, Umin;
    vector<double> vecWeight;

    ArgRMA*   args;
    DataRMA*  data;
    Time      ts;

  };


} // namespace greedyRMA

#endif
