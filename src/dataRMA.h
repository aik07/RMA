//**********************************************************
//  File name:   dataRMA.h
//  Author:      Ai Kagawa
//  Description: a header file for RMA data class
//*********************************************************

#ifndef DATA_h
#define DATA_h

#include <limits>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Time.h"
#include "argRMA.h"
#include "utilRMA.h"

using namespace std;
using namespace arg;
using namespace utilib;


namespace data {

  struct IntMinMax { double minOrigVal, maxOrigVal; };
  struct Feature   { vector<IntMinMax> vecIntMinMax; };


//////////////////// a clsss for integerized dataset (X and w) ////////////////////
class DataXw {

public:

  DataXw() : w(0.0) {}
  DataXw( const vector<unsigned int>& X_ ) : X(X_) {}

  int read(istream& is)        { is >> X >> w; return 0; }
  int write(ostream& os) const { os << X << w; return 0; }

//private:
  vector<unsigned int> X; // integerized explanatory variables
  double               w;	// weight of each observation

  friend class DataRMA;

};


//////////////////// a clsss for original dataset (X and y)////////////////////
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
  DataRMA(int& argc, char**& argv, ArgRMA *args_);

  // DataRMA& operator=( const DataRMA& other );

  // read data from the data file, and set dataOrigTrain
  bool readData(int& argc, char**& argv);

  // bool readRandObs(int argc, char** argv);

  // read nonuniform weights from a specified file
  void readNonUniformWt();

  void setDataDimensions();

  void setNumPosNegObs();
  void setDataIntTrainX();
  void setDataIntTrainWeight();

  void removeZeroWtObs();

  void setVecNumDistFeats();
  void setMaxNumDistFeats();
  void setNumTotalCutPts();

  void setVecAvgX(vector<DataXy> &origData);  // set vecAvgX
  void setVecSdX(vector<DataXy> &origData);   // set vecSdX

  void setAvgY(vector<DataXy> &origData);     // set avgY
  void setSdY(vector<DataXy>  &origData);     // set sdY

  // set the standerdized data of X and Y
  void setStandDataX (vector<DataXy>& origData, vector<DataXy> &standData);
  void setStandDataY (vector<DataXy>& origData, vector<DataXy> &standData);

  // integerize data by using episolon
  void integerizeData(vector<DataXy>& origData, vector<DataXw> &intData);

  // integerize data by using fixed interval length
  void integerizeFixedLengthData(vector<DataXy> &origData,
                                 vector<DataXw> &stdandData);

  template <class T> void saveXObs(T vecData);

//protected:

  // # of observations in original data
  unsigned int numOrigObs;

  // # of distinct observation after discretization
  unsigned int numTrainObs;
  // unsigned int numTestObs;       // # of testing observations

  unsigned int numAttrib;        // # of attributes
  unsigned int numPosTrainObs;   // # of positive training observations
  unsigned int numNegTrainObs;   // # of negative training observations
  unsigned int numTotalCutPts;   // # of cutpoints for RMA
  unsigned int maxNumDistFeats;  // maximum distinct value among all attributes

  // average and standard deviation of Y
  double         avgY, sdY;

  // average and standard deviation vectors of X for each attribute
  vector<double> vecAvgX, vecSdX;

  // minimum and maximum deviation vectors of X for each attribute
  vector<double> vecMinX, vecMaxX;

  // a vector contains # of distinct featrues for each attribute after discretization
  vector<unsigned int>  vecNumDistFeats;

  // a vector contains only training observation indices
  vector<unsigned int>  vecTrainObsIdx;

  // vector<unsigned int>  vecRandObs;    // contains randomize all observations
  // vector<unsigned int>  vecTestObsIdx;   // contains only training dataset observations

  vector<DataXy>  dataOrigTrain;      // original datasets X and y
  vector<DataXy>  dataStandTrain;     // starndardized datasets of X and y
  vector<DataXw>  dataIntTrain;       // discretized data X abd w (weight)

  // vector<DataXy>  dataOrigTest;      // original datasets X and y
  // vector<DataXw>  dataIntTest;       // discretized data X abd w (weight)
  // vector<DataXy>  dataStandTest;

  vector<Feature> vecFeature;    // contains features original and integeried values






  Time     tc;     // Time class object
  ArgRMA   *args;  // ArgRMA class object

};

} // end namespace data

ostream& operator<<(ostream& os, data::DataXw& obj);
istream& operator>>(istream& is, data::DataXw& obj);
ostream& operator<<(ostream& os, data::DataXy& obj);
istream& operator>>(istream& is, data::DataXy& obj);

#endif
