//**********************************************************
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
#include "utility.h"

#ifdef ACRO_HAVE_MPI
#define ROOTRANK (uMPI::rank==0)
#else
#define ROOTRANK true
#endif

using namespace std;
using namespace arg;
using namespace utilib;


namespace data {

  enum TestTrainData {TEST, TRAIN, VALID};

  // each bin has lower and upper bound of the original values
  struct Bin           { double lowerBound, upperBound; };

  // each attronites have manmy bins where each bin for each integer value
  struct BinsPerAttrib { vector<Bin> vecBins; };

  /////////////// a clsss for integerized dataset (X and w) ////////////////
  class DataXw {

  public:

    DataXw() : w(0.0) {}
    DataXw( const vector<unsigned int>& X_ ) : X(X_) {}

    // read X and w
    int read(istream& is)        { is >> X >> w; return 0; }

    // output X and w
    int write(ostream& os) const { os << X << " : "<< w; return 0; }

  //private:
    vector<unsigned int> X; // integerized explanatory variables
    double               w;	// weight of each observation

    friend class DataRMA;

  }; // end DataXw class


  ///////////////// a clsss for original dataset (X and y)/////////////////
  class DataXy {

  public:

    DataXy() : y(0.0) {}
    DataXy( const vector<double>& X_, const int & y_ ) : X(X_), y(y_) { }

    // read X and y
    int read(istream& is)        { is >> X >> y;        return 0; }

    // output X and w
    int write(ostream& os) const { os << X << " : " << y; return 0; }

  //private:
    vector<double> X;  // explanatory variables
    double         y;	 // dependent variable

  friend class DataRMA;

}; // end DataXy class


  /////////////////////////// DataRMA class ///////////////////////////
  class DataRMA {

  public:

    DataRMA() {}
    DataRMA(int& argc, char**& argv, ArgRMA *args_);

    // read data from the data file, and set dataOrigTrain
    void readData(int& argc, char**& argv,
                  const bool &isTest, vector<DataXy> &dataOrig);

    // set dataIntTrainX
    void setDataIntX();

    // remove onservations with zero weights from the training data
    void removeZeroWtObs();

    // set # of distinct values for each attribute
    void setVecNumDistVals();

    // set the maximum # of distinct values among all attributes
    void setMaxNumDistVals();

    // set setDataIntTrainWeight
    void setRMAonlyWeight();

    // read nonuniform weights from a specified file
    void readNonUniformWt();

    // count # of positive and negative observations
    void setNumPosNegObs();

    // set # of total cutpoints
    void setNumTotalCutPts();

    void setVecAvgX();  // set vecAvgX
    void setVecSdX();   // set vecSdX

    void setAvgY();   // set avgY
    void setSdY();    // set sdY

    // set the standerdized data of X and Y
    void setStandardizedX ();
    void setStandardizedY ();

    // integergize data by using the episilon aggregation
    void integerizeEpsData();

    // integerize data by using fixed interval length
    void integerizeFixedData();

    // save X values for all observations
    template <class T> void saveXObs(T vecData);

  //protected:

    unsigned int numAttrib;       // # of attributes
    unsigned int numTrainObs;     // # of training observations
    unsigned int numTestObs;      // # of testing observations
    unsigned int numNonZeroWtObs; // # of non-zero weight observations

    unsigned int maxNumDistVals;  // maximum distinct value among all attributes

    // for RMA only, not needed for Boosting
    unsigned int numPosTrainObs;  // # of positive training observations
    unsigned int numNegTrainObs;  // # of negative training observations
    unsigned int numTotalCutPts;  // # of total cutpoints for PEBBL RMA

    // average and standard deviation of Y
    double               avgY, sdY;

    // average and standard deviation vectors of X for each attribute
    vector<double>       vecAvgX, vecSdX;

    // a vector contains # of distinct featrues for each attribute after discretization
    vector<unsigned int> vecNumDistVals;

    // a vector contains only training observation indices
    vector<unsigned int> vecNonZeroWtObsIdx;

    // a vector contains only testing dataset observations
    // vector<unsigned int>  vecTestObsIdx;

    vector<DataXy>  dataOrigTrain;     // original datasets X and y
    vector<DataXy>  dataStandTrain;    // starndardized datasets of X and y
    vector<DataXw>  dataIntTrain;      // discretized data X abd w (weight)

    vector<DataXy>  dataOrigTest;      // original datasets X and y

    Time     tc;     // Time class object
    ArgRMA   *args;  // ArgRMA class object

    // set setDistVals, distinct values for attribute j
    void setSetDistVals(const int &j);

    // set episilon for attribute j
    void setEpsilon(const int &j);

    // print integerization information
    void printIntegerizationInfo();

    // assign integer values without the recursive technique
    void assignIntNotRecursively(const unsigned int &j);

    // assign integer values recursively
    void assignIntRecursively(const unsigned int &j);

    // set vecBinsCopy
    void setVecBinsCopy(const unsigned int &j);

    // set matOrigLower and matOrigUpper
    void setMapOrigInt(const unsigned int &j);

    // set data for episilon integerization
    void setDataIntEps(const unsigned int &j);

    // print recursive integerization info
    void printRecursiveIntInfo(int j, int k, int countExtraBins,
                               int countL, int countR, int q);

    // print lower and upper bounds
    void printLowerUpperInfo(int lower, int upper);

    // print vecAttribIntInfo
    void printVecAttribIntInfo(const unsigned int &j,
                               const unsigned int &countExtraBins);

    // print info after episilon integerization
    void printAfterEpsIntegerization();

    /***************** fixed length integerization ******************/

    // set the initial values for vecMinValX and vecMaxValX
    void setInitVecMinMaxValX();

    // set dataIntTrain
    void setDataIntFixed();

    // set vecAttribIntInfo while using fixed length integerization
    void setVecAttribIntInfoIntFixed();

    // contains info about bins of lower and upper bounds info for all attributes
    vector<BinsPerAttrib> vecAttribIntInfo;

    /********************** for integerization *****************/

    double           interval;    // confidence interval range

    /********************** for episilon integerization *****************/

    double           eps;         // episilon, aggregation level

    // a set continas all distinct values for each attribute
    set<double>      setDistVals;

    // a container maps from an original value to an
    map<double, int> mapOrigInt;

    // a vector contins min and max for each integerized value
    vector<Bin>      vecBinsCopy;

    /********************** for integerization by fixed bins *****************/
    // minimum and maximum deviation vectors of X for each attribute
    vector<double>   vecMinValX, vecMaxValX;

  }; // end DataRMA class

} // end namespace data

// operator overloading to write and read DataXw class object info
ostream& operator<<(ostream& os, data::DataXw& obj);
istream& operator>>(istream& is, data::DataXw& obj);

// operator overloading to write and read DataXy class object info
ostream& operator<<(ostream& os, data::DataXy& obj);
istream& operator>>(istream& is, data::DataXy& obj);

#endif // end namespace data
