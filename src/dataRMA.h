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

  // each bin has lower and upper bound of the original values
  struct Bin           { double lowerBound, upperBound; };

  // each attronites have manmy bins where each bin for each integer value
  struct BinsPerAttrib { vector<Bin> vecBins; };

//////////////////// a clsss for integerized dataset (X and w) ////////////////////
class DataXw {

public:

  DataXw() : w(0.0) {}
  DataXw( const vector<unsigned int>& X_ ) : X(X_) {}

  int read(istream& is)        { is >> X >> w; return 0; }
  int write(ostream& os) const { os << X << " : "<< w; return 0; }

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
  int write(ostream& os) const { os << X << " : " << y; return 0; }

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

  // count # of positive and negative observations
  void setNumPosNegObs();

  // set dataIntTrainX
  void setDataIntX();

  // set setDataIntTrainWeight
  void setDataIntWeight();

  // remove onservations with zero weights from the training data
  void removeZeroWtObs();

  void setVecNumDistVals();  // set vecNumDistVals
  void setMaxNumDistVals(); // set maxNumDistVals
  void setNumTotalCutPts();

  void setVecAvgX();  // set vecAvgX
  void setVecSdX();   // set vecSdX

  void setAvgY();   // set avgY
  void setSdY();    // set sdY

  // set the standerdized data of X and Y
  void setDataStandX ();
  void setDataStandY ();

  // integergize data by using the episilon aggregation
  void integerizeEpsData();

  // integerize data by using fixed interval length
  void integerizeFixedData();

  template <class T> void saveXObs(T vecData);

  inline unsigned int idxTrain(const int &i) { return vecTrainObsIdx[i]; };
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
  unsigned int maxNumDistVals;  // maximum distinct value among all attributes

  // average and standard deviation of Y
  double           avgY, sdY;

  // average and standard deviation vectors of X for each attribute
  vector<double>   vecAvgX, vecSdX;

  // a vector contains # of distinct featrues for each attribute after discretization
  vector<unsigned int>  vecNumDistVals;

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

  Time     tc;     // Time class object
  ArgRMA   *args;  // ArgRMA class object

  void setSetDistVals(const int &j);

  // set episilon for attribute j
  void setEpsilon(const int &j);

  void printIntegerizationInfo();

  // assign integer values without the recursive technique
  void assignIntNotRecursively(const unsigned int &j);

  // assign integer values recursively
  void assignIntRecursively(const unsigned int &j);

  void setVecBinsCopy(const unsigned int &j);

  void setMapOrigInt(const unsigned int &j);

  void setDataIntEps(const unsigned int &j);

  void printRecursiveIntInfo(int j, int k, int countExtraBins,
                             int countL, int countR, int q);

  void printLowerUpperInfo(int lower, int upper);

  void printVecAttribIntInfo(const unsigned int &j,
                             const unsigned int &countExtraBins);

  void printAfterEpsIntegerization();

  // fixed length integerization

  void setInitVecMinMaxDevX();

  // set dataIntTrain
  void setDataIntFixed();

  void setVecAttribIntInfoIntFixed();

  // contains info about bins of lower and upper bounds info for all attributes
  vector<BinsPerAttrib> vecAttribIntInfo;

  vector<Bin>      vecBinsCopy; // a vector contins min and max for each integerized value

  // minimum and maximum deviation vectors of X for each attribute
  vector<double>   vecMinDevX, vecMaxDevX;

  // a set continas all distinct values for each attribute
  set<double>      setDistVals;

  // a container maps from an original value to an
  map<double, int> mapOrigInt;

  double           eps; // episilon, aggregation level

  double           interval;  // confidence interval range

};

} // end namespace data

ostream& operator<<(ostream& os, data::DataXw& obj);
istream& operator>>(istream& is, data::DataXw& obj);
ostream& operator<<(ostream& os, data::DataXy& obj);
istream& operator>>(istream& is, data::DataXy& obj);

#endif
