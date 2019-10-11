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

#include <pebbl_config.h>
#include <pebbl/utilib/ParameterList.h>
#include <pebbl/utilib/memdebug.h>
#include <pebbl/utilib/seconds.h>
#include <pebbl/utilib/CommonIO.h>
#include <pebbl/bb/pebblParams.h>
#include <pebbl/pbb/parPebblParams.h>

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
  DataXw( const vector<int>& X_ ) : X(X_) { }

  int read(istream& is)        { is >> X >> w; return 0; }
  int write(ostream& os) const { os << X << w; return 0; }

//private:
  vector<int> X;  // integerized explanatory variables
  double      w;	// weight of each observation

  friend class Data;

};


//////////////////// a clsss for original dataset ////////////////////
class DataXy {

public:

  DataXy() {}
  DataXy( const vector<double>& X_, const int & y_ ) : X(X_), y(y_) { }

  int read(istream& is)        { is >> X >> y; return 0; }
  int write(ostream& os) const { os << X << " " << y; return 0; }

//private:
  vector<double> X;  // explanatory variables
  double         y;	 // dependent variable

friend class Data;

};


/////////////////////////// Data class ///////////////////////////
class Data {

public:

  Data() {}
  Data(int argc_, char** argv_, ArgRMA *args_):
      argc(argc_), argv(argv_), args(args_) {
    if (args->debug>=10) cout << "Data::readData\n";
    readData();
    setDataDimensions();
    numTrainObs = numOrigObs;
    vecTrainData.resize(numTrainObs);
    standData.resize(numTrainObs);
    for (int i=0; i<numTrainObs; ++i) vecTrainData[i]=i;
    if (args->debug>=10) cout << "numTrainObs: " << numTrainObs << "\n";
    setStandData();
    if (args->delta() != -1) integerizeData();
    else {
      intData.resize(numOrigObs);
      for (int i=0; i<numOrigObs; ++i) { // for each observation
        intData[i].X.resize(numAttrib);
        for (int j=0; j<numAttrib; j++) // for each attribute
          intData[i].X[j] = origData[i].X[j];
        intData[i].w = origData[i].y;
      } // end while
    }
    //setPosNegObs();
  }

  //bool readData(int argc, char** argv);
  bool readData();
  bool readRandObs(int argc, char** argv);

  void setDataDimensions();

  void integerizeData();
  void setStandData();

  void setPosNegObs();

  void integerizeFixedLengthData();

  void writeIntObs();
  void writeOrigObs();

  void setXStat();

//protected:

  double         avgY, sdY;
  vector<double> avgX, sdX;
  vector<double> minX, maxX;

  int numOrigObs;                // # of observations in original data
  int numTrainObs;               // # of distinct observation after discretization
  int numTestObs;                // # of testing observations
  int numAttrib;                 // # of attributes
  int numPosTrainObs;
  int numNegTrainObs;
  int numTotalCutPts;            // # of cutpoints for RMA
  int maxL;	                     // maximum distinct value among attributes

  vector<int>     distFeat;	     // distinct features after discretization
  vector<int>     vecRandObs;    // contains randomize all observations
  vector<int>     vecTrainData;  // contains only training dataset observations
  vector<int>     vecTestData;   // contains only training dataset observations

  vector<DataXy>  origData;      // original datasets X and y
  vector<DataXw>  intData;       // discretized data X abd w (weight)
  vector<DataXy>  standData;

  vector<Feature> vecFeature;    // contains features original and integeried values

  Time     tc;
	double   wallTime;
  double   cpuTime;

  ArgRMA   *args;
  int      argc;
  char**   argv;

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
