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

    void   setOptMin(const int &j);
    void   setOptMax(const int &j);

    //TODO: combine min max range search together
    double runMinKadane(const int& j) ;
    double runMaxKadane(const int& j) ;

    void   setObjVec(const int &j);

    void   dropObsNotCovered(const int &j, const int& lower, const int& upper);

    double getObjCovered(const int& j, const int& v);

    void   printSolution();

    // void   setInit1DRules();
    // void   set1DOptRange(const int& j);

    // private:

    bool   isPosIncumb;
    bool   foundBox;
    unsigned int    NumNegTiedSols;
    unsigned int    NumPosTiedSols;

    int    tmpL, tmpU;
    double tmpObj, tmpMin, tmpMax;
    double maxObjValue, minVal, maxVal;

    unsigned int    optAttrib;
    unsigned int    prevAttrib;
    unsigned int    optLower, optUpper;

    int    obs;
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
