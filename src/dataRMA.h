/**********************************************************
 *  File name:   dataRMA.h
 *  Author:      Ai Kagawa
 *  Description: a header file for RMA data class
/**********************************************************/

#ifndef DATA_h
#define DATA_h

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Time.h"
#include "argRMA.h"

using namespace std;
using namespace arg;

namespace data {

  struct IntMinMax { double minOrigVal, maxOrigVal; };
  struct Feature   { vector<IntMinMax> vecIntMinMax; };

//////////////////// a clsss for integerized dataset ////////////////////
class DataXw {

public:

  DataXw() : w(0.0) {}
  DataXw( const vector<int>& X_ ) : X(X_) {}

  int read(istream& is)        { is >> X >> w; return 0; }
  int write(ostream& os) const { os << X << w; return 0; }

//private:
  vector<int> X;  // integerized explanatory variables
  double      w;	// weight of each observation

  friend class DataRMA;

};


//////////////////// a clsss for original dataset ////////////////////
class DataXy {

public:

  DataXy() : y(0.0) {}
  DataXy( const vector<double>& X_, const int & y_ ) : X(X_), y(y_) { }

  int read(istream& is)        { is >> X >> y;        return 0; }
  int write(ostream& os) const { os << X << " " << y; return 0; }

//private:
  vector<double> X;  // explanatory variables
  double         y;	 // dependent variable

friend class DataRMA;

};


/////////////////////////// DataRMA class ///////////////////////////
class DataRMA {

public:

  DataRMA() {}
  DataRMA(ArgRMA *args_): args(args_) {}
  DataRMA(int& argc, char**& argv, ArgRMA *args_);

  bool readData(int& argc, char**& argv);
  bool readRandObs(int argc, char** argv);
  void readNonUniformWt();

  void setDataDimensions();
  void setPosNegObs();
  void setIntTrainData();
  void setWeight();
  void setNumMaxDistVal();

  void setXStat(vector<DataXy> &origData);
  void setYStat(vector<DataXy> &origData);

  void setStandDataX (vector<DataXy>& origData, vector<DataXy> &standData);
  void setStandDataY (vector<DataXy>& origData, vector<DataXy> &standData);
  void integerizeData(vector<DataXy>& origData, vector<DataXw> &intData);
  void integerizeFixedLengthData(vector<DataXy> &origData, vector<DataXw> &stdandData);

  template <class T> void writeObs(T vecData);

//protected:

  double         avgY, sdY;
  vector<double> avgX, sdX;
  vector<double> minX, maxX;

  int numOrigObs;                // # of observations in original data
  int numTrainObs;               // # of distinct observation after discretization

  int numAttrib;                 // # of attributes
  int numPosTrainObs;
  int numNegTrainObs;
  int numTotalCutPts;            // # of cutpoints for RMA
  int numMaxDistVal;             // maximum distinct value among attributes

  vector<int>     distFeat;	     // distinct features after discretization
  vector<int>     vecRandObs;    // contains randomize all observations
  vector<int>     vecTrainData;  // contains only training dataset observations
  vector<int>     vecTestData;   // contains only training dataset observations

  vector<DataXy>  origTrainData;      // original datasets X and y
  vector<DataXw>  intTrainData;       // discretized data X abd w (weight)
  vector<DataXy>  standTrainData;

  vector<Feature> vecFeature;    // contains features original and integeried values

  Time     tc;
  double   wallTime;
  double   cpuTime;

  ArgRMA   *args;

};

} // end namespace data


ostream& operator<<(ostream& os, const deque<bool>& v);
ostream& operator<<(ostream& os, const vector<int>& v);
ostream& operator<<(ostream& os, const vector<double>& v);
ostream& operator<<(ostream& os, const vector<vector<int> >& v);
ostream& operator<<(ostream& os, const vector<vector<double> >& v);

ostream& operator<<(ostream& os, data::DataXw& obj);
istream& operator>>(istream& is, data::DataXw& obj);
ostream& operator<<(ostream& os, data::DataXy& obj);
istream& operator>>(istream& is, data::DataXy& obj);

#endif
