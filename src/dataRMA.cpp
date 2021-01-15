/**********************************************************
 * Author:      Ai Kagawa
 * Description: a source file for Data class
 ***********************************************************/


#include "dataRMA.h"


namespace data {


  DataRMA::DataRMA(int &argc, char **&argv, ArgRMA *args_) : args(args_) {

    readData(argc, argv); // read the data and set dataOrigTrain

    // for now, numOrigObs = numTrain
    numTrainObs = numOrigObs;
    vecNonZeroWtObsIdx.resize(numTrainObs);
    for (unsigned int i=0; i < numTrainObs; ++i)  vecNonZeroWtObsIdx[i] = i;

    // Note: It is more efficient to remove observations with zero weights first,
    //     then integerized data. However, assuming that RMA is used for Boosting,
    //     I am integerizing all datasets for now

    setDataIntX();        // set dataIntTrain X (integerization)

    setDataIntWeight();   // set weights for dataIntTrain

    if (args->nonUniformWt() != "")  // if the nonUniform weight file is given
      removeZeroWtObs();             // remove observations with zero weights

    setVecNumDistVals();  // set a vector of # of distinct values for each attribute
    setMaxNumDistVals();  // set the maximum # of the distinct values among attributes

    if (args->isPebblRMA()) setNumTotalCutPts();  // set # of total cut points for B&B

    setNumPosNegObs();    // set # of positive and negative observations

  } // end constructor DataRMA


  // read data file and set dataOrigTrain
  bool DataRMA::readData(int &argc, char **&argv) {

    unsigned int i, j;
    double tmp;
    string line;

    if (args->debug >= 10)
      cout << "Data::readData\n";

    // read data from the data file
    if (argc <= 1) {
      cerr << "No filename specified\n";
      return false;
    }

    if (args->debug>=2) tc.startTime();  // start the timer

    numOrigObs = 0;
    numAttrib  = 0;

    ifstream s(argv[1]); // open the data file

    // check whether or not the file is opened correctly
    if (!s) {
      cerr << "Could not open file \"" << argv[1] << "\"\n";
      return false;
    }

    // read how many columns and rows
    while (getline(s, line)) { // for each lline
      if (numOrigObs == 0) {
        istringstream streamCol(line);
        while (streamCol >> tmp)
          ++numAttrib;
      }
      ++numOrigObs;
    }
    --numAttrib; // last line is response value

  #ifdef ACRO_HAVE_MPI
    if (uMPI::rank == 0) {
  #endif //  ACRO_HAVE_MPI
      cout << "(mxn): " << numOrigObs << "\t" << numAttrib << "\n";
  #ifdef ACRO_HAVE_MPI
    }
  #endif //  ACRO_HAVE_MPI

    s.clear();
    s.seekg(0, ios::beg);

    dataOrigTrain.resize(numOrigObs);
    for (i = 0; i < numOrigObs; ++i) { // for each observation
      dataOrigTrain[i].X.resize(numAttrib);
      for (j = 0; j < numAttrib; j++) // for each attribute
        s >> dataOrigTrain[i].X[j];
      s >> dataOrigTrain[i].y;
    } // end for each observation

    // if the original data has 0 as -1 class, change from 0 to -1
    for (i = 0; i < numOrigObs; ++i) // for each observation
      if (dataOrigTrain[i].y == 0)   // if y=0
        dataOrigTrain[i].y = -1.0;

    s.close(); // close the data file

    if (args->debug>=2)
      cout << "setupProblem CPU time: " << tc.getCPUTime() ; // << tc.getWallTime();

    return true;

  } // end readData function


  // set X for the dataIntTrain
  void DataRMA::setDataIntX() {

    // TODO: add fixed bins integerization

    // set dataIntTrainX dimensions
    dataIntTrain.resize(numTrainObs);
    for (unsigned int i = 0; i < numTrainObs; ++i)
      dataIntTrain[idxTrain(i)].X.resize(numAttrib);

    // set integerized X values
    if (args->delta() != -1) { // integerize data

      vecAttribIntInfo.resize(numAttrib);
      integerizeEpsData();

    } else { // assuming that the original data is already integerized

      for (unsigned int i = 0; i < numTrainObs; ++i)  // for each observation
        for (unsigned int j = 0; j < numAttrib; ++j)  // for each attribute
          // assign the original X as integerized X
          dataIntTrain[idxTrain(i)].X[j] = dataOrigTrain[idxTrain(i)].X[j];
    }

    if (args->debug>=10) saveXObs(dataIntTrain);

  } // end setDataIntX function


  // set weights of dataIntTrain
  void DataRMA::setDataIntWeight() {

    if (args->nonUniformWt() != "") { // if the nonUniform weight file is given
      readNonUniformWt();  // set the weights from the file

    } else {
      // give equal weight for each observation (1/numTrainObs)
      for (unsigned int i = 0; i < numTrainObs; ++i) // for each observation
        dataIntTrain[idxTrain(i)].w
              = dataOrigTrain[idxTrain(i)].y * 1.0 / (double)numTrainObs;
    } // end if

  } // end setDataIntWeight function


  // remove observations with zero weight
  void DataRMA::removeZeroWtObs() {

    unsigned int numNonZeroObs = -1;

    vecNonZeroWtObsIdx.resize(numOrigObs);

    for (unsigned int i = 0; i < numOrigObs; ++i) // for each training observation
      if (dataIntTrain[i].w != 0)  // if objservation i has non-zero wiehgt
        vecNonZeroWtObsIdx[++numNonZeroObs] = i;

    numTrainObs = numNonZeroObs+1;
    vecNonZeroWtObsIdx.resize(numTrainObs);

    if (args->debug >= 10)
      cout << "numTrainObs: " << numTrainObs << "\n";

  } // end removeZeroWtObs function


  // read non uniform weight and set the weight for each observation
  void DataRMA::readNonUniformWt() {

    // TODO: for now, each process is reading non-uniform weight file
    /*
  #ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
  #endif //  ACRO_HAVE_MPI
    */

     // Input stream class to operate on non-uniform weight file
    ifstream inFile(args->nonUniformWt());

    if (inFile.is_open()) { // if the file is open

      string line, tmp;

      // read a string or a line from an input stream
      // inFile is an object of istream class
      // line is a string object that input is stored in this object after being read from the stream
      while (getline(inFile, line)) {

        stringstream ss(line);  // Used for breaking words

        for (unsigned int i = 0; i < numOrigObs; ++i) { // for each training observation
          getline(ss, tmp, ',');
          // cout << "tmp " << tmp << "\n";
          dataIntTrain[i].w = stod(tmp);  // stod: convert string to double
          // cout << "vec " <<vecNonUniformWt[i] << "\n";
        } // end for each train observation

      } // end while

    } else { // if the program could not open the file
      cerr << "error: cannot read nonuniform wt";
    } // end if the file is open

    // rma->setDataIntWeight(vecNonUniformWt, vecObsIdx);

    /*
  #ifdef ACRO_HAVE_MPI
  }
  #endif //  ACRO_HAVE_MPI
    */

    if (args->debug >= 10) {
      ucout << "rank: " << uMPI::rank << " wt: ";
      for (unsigned int i = 0; i < numOrigObs; ++i) // for each observation
        ucout << dataIntTrain[i].w << ", ";          // print out weight
      ucout << "\n";
    } // end if

  } // end readNonUniformWt function


  // set vecNumDistVals (# of distinc values for each attributes)
  void DataRMA:: setVecNumDistVals() {

    vecNumDistVals.resize(numAttrib);  // resize vecNumDistVals

    for (unsigned int j = 0; j < numAttrib; ++j) { // for each attribute

      vecNumDistVals[j] = 1;  // the minimum distinct # is 1

      // if there are more distinct values, increment the # of distinct values
      for (unsigned int i = 0; i < numTrainObs; ++i)  // for each observation
        if (vecNumDistVals[j] < dataIntTrain[idxTrain(i)].X[j] + 1)
          vecNumDistVals[j] = dataIntTrain[idxTrain(i)].X[j] + 1;

    } // end for each attribute

  } // end setVecNumDistVals function


  // set the maximum distinct value of all attributes
  void DataRMA::setMaxNumDistVals() {

    maxNumDistVals     = 0;

    for (unsigned int j = 0; j < numAttrib; ++j)  // for each attribute
      if (maxNumDistVals < vecNumDistVals[j])
        maxNumDistVals = vecNumDistVals[j];

  } // end setMaxNumDistVals function


  // set the maximum distinct value of all attributes
  void DataRMA::setNumTotalCutPts() {

    numTotalCutPts = 0;

    for (unsigned int j = 0; j < numAttrib; ++j)  // for each attribute
      numTotalCutPts += vecNumDistVals[j]-1; // sum up the tatoal cut point

  } // end setNumTotalCutPts function


  // set positive and negative observations
  void DataRMA::setNumPosNegObs() {

    numPosTrainObs = 0;
    numNegTrainObs = 0;

    for (unsigned int i = 0; i < numTrainObs; ++i) { // for each training observation
      if (dataOrigTrain[idxTrain(i)].y == 1) // if the response value is 1
        ++numPosTrainObs;  // count observation i as postive
      else
        ++numNegTrainObs;  // count observation i as negative
    } // end for each observation

  #ifdef ACRO_HAVE_MPI
    if (uMPI::rank == 0) {
  #endif //  ACRO_HAVE_MPI
      ucout << "# of Positive : Negative Observations: "
            << numPosTrainObs << " : " << numNegTrainObs << "\n";
  #ifdef ACRO_HAVE_MPI
    }
  #endif //  ACRO_HAVE_MPI

  } // end setNumPosNegObs function


  // set vecAvgX, a vector of average X for each attribute
  void DataRMA::setVecAvgX() {

    vecAvgX.resize(numAttrib);

    for (unsigned int j = 0; j < numAttrib; ++j)
      vecAvgX[j] = 0;

    for (unsigned int j = 0; j < numAttrib; ++j) {   // for each attribute j
      for (unsigned int i = 0; i < numTrainObs; ++i) // for each observation i
        vecAvgX[j] += dataOrigTrain[idxTrain(i)].X[j]; // set the sum of X value for attribute j
      vecAvgX[j] /= numTrainObs;        // set the average of X value for attribute j
    } // for each attribute j

    if (args->debug >= 2) cout << "vecAvgX: " << vecAvgX;

  } // end setVecAvgX function


  // set vecSdX, a vector of standard deviation of X for each attribute
  void DataRMA::setVecSdX() {

    vecSdX.resize(numAttrib);

    for (unsigned int  j = 0; j < numAttrib; ++j)
      vecSdX[j]  = 0;

    // set the vecAvgX and vecSdX for all attributes
    for (unsigned int  j = 0; j < numAttrib; ++j) {   // for each attribute j

      for (unsigned int i = 0; i < numTrainObs; ++i)  // for each observation
        vecSdX[j] += pow(dataOrigTrain[idxTrain(i)].X[j] - vecAvgX[j], 2);

      vecSdX[j] /= numTrainObs;     // divied by # of observations
      vecSdX[j] = sqrt(vecSdX[j]);  // square root

    } // end for each attribute j

    if (args->debug >= 2) cout << "vecSdX: " << vecSdX;

  } // end setVecSdX function


  // set avgY, the average of Y values
  void DataRMA::setAvgY() {

    avgY = 0;

    for (unsigned int i = 0; i < numTrainObs; ++i)
      avgY += dataOrigTrain[idxTrain(i)].y; // get avg of y

    avgY /= numTrainObs; // get average response value

  } // end setAvgY function


  // set sdY, the standard devication of Y values
  void DataRMA::setSdY() {

    sdY = 0;

    for (unsigned int i = 0; i < numTrainObs; ++i)
      sdY += pow(dataOrigTrain[idxTrain(i)].y - avgY, 2); // get std dev of y

    sdY /= numTrainObs;
    sdY = sqrt(sdY);

  } // end setSdY function


  // set Standadize Data
  void DataRMA::setDataStandX() {

    dataStandTrain.resize(numTrainObs);

    if (vecAvgX.size()==0) setVecAvgX();
    if (vecSdX.size()==0)  setVecSdX();

    // standardize X in each attribute
    for (unsigned int j = 0; j < numAttrib; ++j) { // for each attribute

      for (unsigned int i = 0; i < numTrainObs; ++i) { // for each observation

        // resize each observation's X
        dataStandTrain[i].X.resize(numAttrib);

        if (vecSdX[j]!=0)
          dataStandTrain[idxTrain(i)].X[j]
            = (dataOrigTrain[idxTrain(i)].X[j] - vecAvgX[j]) / vecSdX[j];
        else  // if the standard deviation is 0, do not divide
          dataStandTrain[idxTrain(i)].X[j]
             = dataOrigTrain[idxTrain(i)].X[j] - vecAvgX[j];
      } // end for each observation

    } // end for each attribute

    if (args->debug >= 10) {
      // print the dataStandTrain
      cout << "dataStandTrain: \n";
      for (unsigned int i = 0; i < numTrainObs; ++i)
        cout << "obs: " << idxTrain(i)
             << ": "    << dataStandTrain[idxTrain(i)] << "\n";

    }  // end debug

  }  // end setDataStandX function


  // set the standardized data for Y-value
  void DataRMA::setDataStandY() {

    dataStandTrain.resize(numTrainObs); // TODO: resized twice

    setAvgY(); // set avgY
    setSdY();  // set sdY

    // standardize y
    for (unsigned int i = 0; i < numTrainObs; ++i) // for each obseration i
      dataStandTrain[idxTrain(i)].y
        = (dataOrigTrain[idxTrain(i)].y - avgY) / sdY;

  } // end setDataStandY function

  /////////////////////// integerization ///////////////////////////////

  // integergize data by using the episilon aggregation
  void DataRMA::integerizeEpsData() {

    if (args->debug >= 1) tc.startTime(); // start the timer

    vecNumDistVals.resize(numAttrib);

    for (unsigned int j = 0; j < numAttrib; ++j) { // for each attribute

      setSetDistVals(j);           // set setDostVals

      setEpsilon(j);               // set epsilon

      assignIntNotRecursively(j);  // non-recursive integerization

      assignIntRecursively(j);     // recursive integerization

      setDataIntEps(j);            // set dataInt for Epsilon Integerization

    } // end for each attribute

    if (args->debug >= 1)  printAfterEpsIntegerization();

  } // end integerizeEpsData


  // set setDistVals for attribute j
  void DataRMA::setSetDistVals(const int &j) {

    setDistVals.clear(); // Clear content

    for (unsigned int i = 0; i < numTrainObs; ++i) { // for each training observation
      if (args->debug >= 10)
        cout << "dataOrigTrain: " << dataOrigTrain[idxTrain(i)].X[j] << "\n";
      setDistVals.insert(dataOrigTrain[idxTrain(i)].X[j]); // TODO: This is not efficient
    } // end for each observation

    if (args->debug >= 10) {
      cout << "size setDistVals: " << setDistVals.size() << "\n";
      cout << "setDistVals: "      << setDistVals;
    }

  } // end setSetDistVals function


  // set episilon for attribute j
  void DataRMA::setEpsilon(const int &j) {

    if (vecAvgX.size()==0) setVecAvgX();
    if (vecSdX.size()==0)  setVecSdX();

    // TODO: avoid hard-cording of 4.0
    // get 95% confidence interval range
    // minimm of 4 standard deviations or the entire rnage

    interval = min(4.0 * vecSdX[j],
                   *setDistVals.rbegin() - *setDistVals.begin());

    // minimum of delta or max interval limit, and multiple the interval
    eps = min(args->delta(), args->maxInterval()) * interval;

    if (args->debug >= 2) printIntegerizationInfo();

  } // end setEpsilon function


  /************ assign integer without recursive integerization ************/
  void DataRMA::assignIntNotRecursively(const unsigned int &j) {

    unsigned int k = 0; // counter for each distinct value
    set<double>::iterator it, itp; // iterator for the set

    mapOrigInt.clear();

    vecAttribIntInfo[j].vecBins.resize(setDistVals.size());

    // The first bin's lower bound is eqal to the smallest distinct value
    itp = setDistVals.begin();
    vecAttribIntInfo[j].vecBins[0].lowerBound = *itp;
    vecAttribIntInfo[j].vecBins[0].upperBound = *itp;

    // walk thorugh the set of distincet values
    // some value can be aggregated by the level of the episilon
    for (it = setDistVals.begin(); it != setDistVals.end(); ++it) {

      if (args->debug >= 2) printLowerUpperInfo(*itp, *it);

      // if the distance between the current lower and upper is
      // greater than episilon, assign the next number for the current upper
      if ( (*it - *itp) > eps ) {

        // the previous bin's upper bound is the last value
       vecAttribIntInfo[j].vecBins[k].upperBound  = *itp; // *(--it);

        // the current bin's lower bound is the current value
       vecAttribIntInfo[j].vecBins[++k].lowerBound = *it; // *(++it);

      } // end if creating a new bin

      // set the previous distinct value as the current distinct value
      itp             = it;

      // set the current original values to map the "k" integer value
      mapOrigInt[*it] = k;

    } // for each setDistVals

    // set the final bin's upper bound
    vecAttribIntInfo[j].vecBins[k].upperBound = *it; // TODO: check this
    // vecAttribIntInfo[j].vecBins[k].upperBound = *(--it);

    if (args->debug >= 2)  cout << "mapOrigInt contains: " << mapOrigInt;

    vecNumDistVals[j] = k+1; // get distinct # of feature

    vecAttribIntInfo[j].vecBins.resize(vecNumDistVals[j]);

  } // end assignIntNotRecursively function


  // assign the integer recursively
  void DataRMA::assignIntRecursively(const unsigned int &j) {

    bool         isSplit;
    unsigned int countL, countR, countExtraBins, countIn;
    double       curLower, curUpper, curLowerIn, curUpperIn, curUpperIn_updated;

    set<double>::iterator it;

    // if there is interval limit
    // and the size of distince value is not same as the number of integers
    // assigned
    if ( args->maxInterval() != getInf()
        || vecNumDistVals[j] != setDistVals.size() ) {

      setVecBinsCopy(j);

      countExtraBins = 0; // TOTO: what is countExtraBins? (It was p)

      for (unsigned int k = 0; k < vecNumDistVals[j]; ++k) { // for each distinct value

        isSplit = true;
        countR  = 0;        // TODO: what is r?

        // set tempL and temU
        curLower = vecBinsCopy[k].lowerBound;
        curUpper = vecBinsCopy[k].upperBound;

        // if each interval violates the limit, splitting, and episilon > threshold
        while ( (curUpper - curLower) > args->maxInterval() * interval
                && isSplit
                && eps > .0001) { // TODO: specify .0001

          isSplit = false;
          eps     *= args->shrinkDelta(); // shrink eps

          if (args->debug >= 2)
            cout << "shrinked episilon: " << eps << '\n';

          for (countIn = 0; countIn <= countR; ++countIn) { // TODO: waht is q

            countL = 0;  // TODO: what is it?

            curLowerIn = vecAttribIntInfo[j].vecBins[k + countExtraBins + countIn]
                       .lowerBound;
            curUpperIn = vecAttribIntInfo[j].vecBins[k + countExtraBins + countIn]
                       .upperBound;

            if (args->debug >= 2) printLowerUpperInfo(curLowerIn, curUpperIn);

            // if the interval violates the limit
            if ((curUpperIn - curLowerIn) > args->maxInterval() * interval) {

              isSplit = true;

              for (it = setDistVals.find(curLowerIn); ; ++it) { // TODO: what is this?

                curUpperIn_updated = *it; // update the curren upper bound

                if (args->debug >= 2) {
                  cout << "updated upper: ";
                  printLowerUpperInfo(curUpperIn_updated, curUpperIn);
                }

                // if the updated interval is greater than episilon
                // create an extra bin
                if ((curUpperIn_updated - curLowerIn) > eps) {

                  ++countR;
                  ++countL; // count for an extra bin

                  // set the previous bin's upper bound
                 vecAttribIntInfo[j].vecBins[k + countExtraBins + countL + countIn - 1]
                               .upperBound = curLowerIn;

                  // set the current bin's lower bound
                 vecAttribIntInfo[j].vecBins[k + countExtraBins + countL + countIn]
                               .lowerBound = curUpperIn_updated;

                  if (args->debug >= 2)
                    printRecursiveIntInfo(j, k, countExtraBins, countL, countR, countIn);

                // error if the interval is greater than the limit
                }  else if ( (curUpper - curLower) < 0 ) {

                  cerr << "!!!!!!!!!!Something Wrong!!!!!!!!!!!\n";

                  if (args->debug >= 2) printVecAttribIntInfo(j, countExtraBins);

                } // end if

                // current lower bound = the current upper bound which is updated
                curLowerIn = curUpperIn_updated;
                vecAttribIntInfo[j].vecBins[k + countExtraBins + countL + countIn]
                         .upperBound = curUpper;

                // if the current upper updated and the old upper values are the same
                if (curUpperIn_updated == curUpper)
                  break;

              } // end for each inner sub interval

            }   // end if each interval is greater than the threthold

          }   // end for (countExtraBins=0; p<=r; ++countExtraBins)

        } // end while ( (curUpper-curLower) > getLimitInterval()*interval && isSplit)

        countExtraBins += countR; // count the extra bins

        // end if the interval does not violate the limit and
        // we need extra bins
        if ( (curUpper - curLower) <= args->maxInterval() * interval
             && countExtraBins > 0) {

         // set the vecAttribIntInfo
         vecAttribIntInfo[j].vecBins[k + countExtraBins].lowerBound =
              vecBinsCopy[k].lowerBound;

         vecAttribIntInfo[j].vecBins[k + countExtraBins].upperBound =
              vecBinsCopy[k].upperBound;

        } // end if

      } // for each distinct value

      if (args->debug >= 2) printVecAttribIntInfo(j, countExtraBins);

      setMapOrigInt(j);

      vecNumDistVals[j] += countExtraBins; // update distinct # of feature

      vecAttribIntInfo[j].vecBins.resize(vecNumDistVals[j]);

    } // end if recursive discretization

  } // end assignIntRecursively function


  // set vecBinsCopy
  void DataRMA::setVecBinsCopy(const unsigned int &j) {

    unsigned int k;

    // set vecBinsCopy; // TODO: why copy?
    vecBinsCopy.resize(vecNumDistVals[j]);

    for (k = 0; k < vecNumDistVals[j]; ++k) { // for each distinct value

      vecBinsCopy[k].lowerBound = vecAttribIntInfo[j].vecBins[k].lowerBound;

      vecBinsCopy[k].upperBound = vecAttribIntInfo[j].vecBins[k].upperBound;

    } // end for each distinct value

    if (args->debug >= 2) {
      cout << "\nBin Lower Bound: ";
      for (k = 0; k < vecNumDistVals[j]; ++k)
        cout << vecBinsCopy[k].lowerBound << ' ';
      cout << "\nBin Upper Bound:";
      for (k = 0; k < vecNumDistVals[j]; ++k)
        cout << vecBinsCopy[k].upperBound << ' ';
      cout << '\n';
    } // end debug

  } // end setVecBinsCopy function


  // set mapOrigInt, map a original value to it integer
  void DataRMA::setMapOrigInt(const unsigned int &j) {

    unsigned int idxBin = 0;
    set<double>::iterator it;

    // for each distinct values in setDistVals
    for (it = setDistVals.begin(); it != setDistVals.end(); ++it) {
      // if current value is greater than the original upper bound
      if (*it >vecAttribIntInfo[j].vecBins[idxBin].upperBound)
        ++idxBin; // go to the next bucket

     mapOrigInt[*it] = idxBin;

   } // end for each distinct value

   if (args->debug >= 2)
     cout << "at recurisve iteration, mapOrigInt: " << mapOrigInt;

  } // end setMapOrigInt function


  // set dataIntTrain for recursively discretization
  void DataRMA::setDataIntEps(const unsigned int &j) {

    for (unsigned int i = 0; i < numTrainObs; ++i) { // for each observation
      dataIntTrain[idxTrain(i)].X.resize(numAttrib);
      dataIntTrain[idxTrain(i)].X[j] = mapOrigInt[dataOrigTrain[idxTrain(i)].X[j]];
    } // end for each observation

  } // end setDataIntEps function


  // print integerization info
  void DataRMA::printIntegerizationInfo() {
    cout << "delta: " << args->delta()
         << ", max: " << *setDistVals.rbegin()
         << ", min: " << *setDistVals.begin()
         << ", eps: " << eps << "\n"
         << ", maxInterval: " << args->maxInterval() * interval << "\n";
  } // end printIntegerizationInfo function


  // print the current lower and upper bounds and its difference
  void DataRMA::printLowerUpperInfo(int lower, int upper) {
    cout << "current lower: "   << lower
         << ", current upper: " << upper
         << ", diff: "          << (upper - lower) << "\n";
  }  // end printLowerUpperInfo function

  // print recursive integerization info
  void DataRMA::printRecursiveIntInfo(int j, int k, int countExtraBins,
                                      int countL, int countR, int countIn) {

    cout << " idx: " << k + countExtraBins + countL + countIn - 1
         << " Prev Bin's Upper: "
         << vecAttribIntInfo[j].vecBins
              [k + countExtraBins + countL + countIn - 1].upperBound
         << " Curr Bin's Lower: "
         << vecAttribIntInfo[j].vecBins
              [k + countExtraBins + countL + countIn].lowerBound
         << "\n";

    cout << " pre-recursive bin k: " << k  // for each distinct value
         << ", countExtraBins: "     << countExtraBins
         << ", countR: "             << countR
         << ", countL: "             << countL
         << ", countIn: "            << countIn
         << "\n";

  } // end printRecursiveIntInfo function


  void DataRMA::printVecAttribIntInfo(const unsigned int &j,
                                      const unsigned int &countExtraBins) {

    cout << "\nvecAttribIntInfo[j].vecBins's lowerBound";
    for (unsigned int k = 0; k < vecNumDistVals[j] + countExtraBins; ++k)
      cout << vecAttribIntInfo[j].vecBins[k].lowerBound << ' ';

    cout << "\nvecAttribIntInfo[j].vecBins's lupperBound ";
    for (unsigned int k = 0; k < vecNumDistVals[j] + countExtraBins; ++k)
      cout << vecAttribIntInfo[j].vecBins[k].upperBound << ' ';

    cout << '\n';

  }  // end printVecAttribIntInfo function


  // print info after episilon integerization
  void DataRMA::printAfterEpsIntegerization() {

    cout << "vecNumDistVals: " << vecNumDistVals << "\n";

    for (unsigned int i = 0; i < numTrainObs; ++i)
      cout << "obs: " << idxTrain(i) << ": " << dataIntTrain[idxTrain(i)] << "\n";

    cout << "integerizeProblem: \t";
    tc.getCPUTime();

    if (args->debug >= 2)
      tc.getWallTime();

    /*
      #ifdef ACRO_HAVE_MPI
        if (uMPI::rank==0) {
      #endif //  ACRO_HAVE_MPI
        if (writePred()) {
          saveXObs();
          saveXObs();
        }
      #ifdef ACRO_HAVE_MPI
        }
      #endif //  ACRO_HAVE_MPI
    */

  } // end printAfterEpsIntegerization function


  /******************* integerize X using a fixed length ****************/

  //  integerize values into the fixed bins
  void DataRMA::integerizeFixedData() {

    setInitVecMinMaxValX();         // initialize vecMinValX and vecMaxValX

    setDataIntFixed();              // setIntData

    setVecNumDistVals();            // set vecNumDistVals

    setVecAttribIntInfoIntFixed();  // setVecAttribIntInfo

  } // end integerizeFixedData function


  // set vecMinValX and vecMaxValX, min and max values of each attribute
  void DataRMA::setInitVecMinMaxValX() {

    vecMaxValX.resize(numAttrib);
    vecMinValX.resize(numAttrib);

    // reset the minimum and maximum values of X for each attribute
    for (unsigned int j = 0; j < numAttrib; ++j) {
      vecMinValX[j] = getInf();
      vecMaxValX[j] = -getInf();
    }

    for (unsigned int j = 0; j < numAttrib; ++j) {  // for each attribute

      for (unsigned int i = 0; i < numTrainObs; ++i) {  // for each observation

        // if current observation feature is less than the minimum for the attribute j
        if (dataOrigTrain[idxTrain(i)].X[j] < vecMinValX[j])
          vecMinValX[j] = dataOrigTrain[idxTrain(i)].X[j]; // set vecMinValX[j]

        // if the current observation feature is more than the maximum for the attribute j
        if (dataStandTrain[idxTrain(i)].X[j] > vecMaxValX[j])
          vecMaxValX[j] = dataStandTrain[idxTrain(i)].X[j]; // set vecMaxValX[j]

      } // end for each observation

    } // for each attribute

  } // end setInitVecMinMaxValX function


  // set integerized data for the integerization by the fixed length
  void DataRMA::setDataIntFixed() {

    unsigned int idxBin;
    vector<bool> vecIsEmptyBins;
    vector<unsigned int> vecAccEmptyIdx;

    for (unsigned int i = 0; i < numTrainObs; ++i)  // for each observation
      dataIntTrain[idxTrain(i)].X.resize(numAttrib);

    for (unsigned int j = 0; j < numAttrib; ++j) { // for each attribute

      vecIsEmptyBins.resize(args->getNumFixedSizeBins());
      fill(vecIsEmptyBins.begin(), vecIsEmptyBins.end(), true);

      interval = ( vecMaxValX[j] - vecMinValX[j] )
                  / (double) args->getNumFixedSizeBins();

      for (unsigned int i = 0; i < numTrainObs; ++i) { // for each observation

        idxBin = floor( (dataOrigTrain[idxTrain(i)].X[j] - vecMinValX[j])
                 / interval );

        vecIsEmptyBins[idxBin] = false;

        dataIntTrain[idxTrain(i)].X[j] = idxBin;

      } // end for each observation

    } // end for attribute

    for (unsigned int j = 0; j < numAttrib; ++j) { // for each attribute

      idxBin = 0;

      for (unsigned int k = 0; k < args->getNumFixedSizeBins(); ++k) { // for each bin

        if (vecIsEmptyBins[k]) { // if this index is not empty
          ++idxBin;
          vecAccEmptyIdx[k] = idxBin;  // accumurate how many empty bins so far
        } // end if this index is not empty

      } // end for each bin

      // adjust integer values considering the empty bins
      for (unsigned int i = 0; i < numTrainObs; ++i) // for each observation
        dataIntTrain[idxTrain(i)].X[j]
          -= vecAccEmptyIdx[dataIntTrain[idxTrain(i)].X[j]];

    } // end for attribute

  } // end setDataIntFixed function


  // set vecAttribIntInfo uning the fixed interval
  void DataRMA::setVecAttribIntInfoIntFixed() {

    for (unsigned int j = 0; j < numAttrib; ++j) { // for each attribute

      interval = ( vecMaxValX[j] - vecMinValX[j] )
                 / (double) args->getNumFixedSizeBins();

      vecAttribIntInfo[j].vecBins.resize(vecNumDistVals[j]);

      for (unsigned int k = 0; k < maxNumDistVals; ++k) { // for each distinct value

        // lower bound in the original value in the fixed interval
        vecAttribIntInfo[j].vecBins[0].lowerBound =
            (double) k       * interval + vecMinValX[j];

        // upper bound in the original value in the fixed interval
        vecAttribIntInfo[j].vecBins[0].upperBound =
            (double) (k + 1) * interval + vecMinValX[j];

      } // for each distinct feature

    } // for each attribute

  } // end setVecAttribIntInfoIntFixedIntFixed function


  // save integerized or original X values of all the training observations
  template <class T> void DataRMA::saveXObs(T vecData) {

    // set the file name
    stringstream s;
    (typeid(T) == typeid(unsigned int)) ? s << "intX.out" : s << "origX.out";
    ofstream os(s.str().c_str());

    // output X values
    for (unsigned i = 0; i < numTrainObs; ++i)  // for each observation
      os << vecData[idxTrain(i)].X << "\n";

    os.close();  // close the file

  } // end saveXObs


} // namespace data


/*********************** operators ************************/

// Operators to read and write RMA Objs to streams
ostream &operator<<(ostream &os, data::DataXy &obj) {
  obj.write(os);
  return os;
}

istream &operator>>(istream &is, data::DataXy &obj) {
  obj.read(is);
  return is;
}

// Operators to read and write RMA Objs to streams
ostream &operator<<(ostream &os, data::DataXw &obj) {
  obj.write(os);
  return os;
}

istream &operator>>(istream &is, data::DataXw &obj) {
  obj.read(is);
  return is;
}
