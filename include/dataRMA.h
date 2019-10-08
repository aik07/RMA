/*
 *  File name:   data.h
 *  Author:      Ai Kagawa
 *  Description: a header file for LPBase and Data classes
 */

#ifndef Data_h
#define Data_h

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
#include "allParams.h"


namespace data {

//////////////////// a clsss for integerized dataset ////////////////////
class DataXw {

public:

  DataXw() : w(0.0) {}
  DataXw( const vector<int>& X_ ) : X(X_) { }

  int read(istream& is) { is >> X >> w; return 0; }
  int write(ostream& os) const { os << X << w;	return 0; }

private:
  vector<int> X;  // integerized explanatory variables
  double w;	  // weight of each observation

  friend class Data;

};


//////////////////// a clsss for original dataset ////////////////////
class DataXy {

public:

  DataXy() {}
  DataXy( const vector<double>& X_, const int & y_ ) : X(X_), y(y_) { }

  int read(istream& is) { is >> X >> y; return 0; }
  int write(ostream& os) const { os << X << " " << y; return 0; }

private:
  vector<double> X;  // explanatory variables
  double y;	     // dependent variable

  friend class Data;

};


/////////////////////////// Data class ///////////////////////////
class Data {

public:

  Data() {}

  bool readData(int argc, char** argv);
  bool readRandObs(int argc, char** argv);

  void setDataDimensions();

  void integerizeData() ;
  void setStandData();

  void setPosNegObs();

  void integerizeFixedLengthData();

  void writeIntObs();
  void writeOrigObs();

  void setXStat();

protected:

  double avgY, sdY;
  vector<double> avgX, sdX;
  vector<double> minX, maxX;

  int numOrigObs;             // # of observations in original data
  int numTrainObs;            // # of distinct observation after discretization
  int numTestObs;             // # of testing observations
  int numAttrib;              // # of attributes
  int numPosTrainObs;
  int numNegTrainObs;

  int numTotalCutPts;         // # of cutpoints for RMA
  int maxL;	                  // maximum distinct value among attributes

  vector<DataXy> origData;    // original datasets X and y
  vector<DataXw> intData;     // discretized data X abd w (weight)
  vector<DataXy> standData;

  vector<int> distFeat;	      // distinct features after discretization
  vector<int> vecRandObs;     // contains randomize all observations
  vector<int> vecTrainData;   // contains only training dataset observations
  vector<int> vecTestData;    // contains only training dataset observations
  vector<Feature> vecFeature; // contains features original and integeried values

};

} // data namespace


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
