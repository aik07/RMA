/*
*  File name:   greedyRMA.h
*  Author:      Ai Kagawa
*  Description: a serial greedy rectangular maximum agreement problem solver
*/

#ifndef greedyRMA_h
#define greedyRMA_h

#include <ostream>
#include <vector>

#include "Time.h"
#include "base.h"

using namespace base;
using namespace std;


namespace greedyRMA {

class GreedyRMA {

public:

  GreedyRMA(Data *d) : data(d) {}

  void runGreedyRangeSearch();

  void chooseMinOrMaxRange();

  void getMinOptRange();
  void getMaxOptRange();

  void setOptMin(const int &j);
  void setOptMax(const int &j);

  //TODO: combine min max range search together
  double runMinKadane(const int& j) ;
  double runMaxKadane(const int& j) ;

  void setObjVec(const int &j);

  void dropObsNotCovered(const int &j, const int& lower, const int& upper);

  double getObjCovered(const int& j, const int& v);

  void printSolution();

  void setInit1DRules();
  void set1DOptRange(const int& j);

// private:

  vector<int> vecCoveredObs;
  vector<int> L, U;
  vector<int> Lmax, Umax;
  vector<int> Lmin, Umin;
  vector<double> W;

  bool isPosIncumb;

  int NumNegTiedSols;
  int NumPosTiedSols;

  bool foundBox;

  int tmpL;
  int tmpU;

  double tmpObj;
  double tmpMax;
  double tmpMin;

  double minVal;
  double maxVal;
  double maxObjValue;

  int optAttrib;
  int optLower;
  int optUpper;
  int prevAttrib;

  int obs;
  bool fondNewBox;

  Data* data;

  Time ts;

};


} // namespace greedyRMA

#endif
