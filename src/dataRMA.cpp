/**********************************************************
 * File name:   dataRMA.cpp
 * Author:      Ai Kagawa
 * Description: a source file for Data class
 ***********************************************************/
#include "dataRMA.h"


namespace data {


DataRMA::DataRMA(int &argc, char **&argv, ArgRMA *args_) : args(args_) {

  readData(argc, argv); // read the data and set dataOrigTrain
  setDataDimensions();  // set data dimensions

  // setStandData();
  setDataIntTrainX();        // set dataIntTrain
  setDataIntTrainWeight();  // set weights for dataIntTrain

  if (args->nonUniformWt() != "")  // if the nonUniform weight file is given
    removeZeroWtObs();   // remove observations with zero weights

  setVecNumDistFeats(); // set a vector of # of distinct values for each attribute
  setMaxNumDistFeats(); // set the maximum number of the distinct values among attributes
  setNumTotalCutPts();  // set # of total cut points for B&B

  setNumPosNegObs();      // set each observation to be positive or negative

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

  tc.startTime();  // start the timer

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

  // for (int i=0; i<numOrigObs; ++i)
  //   ucout << "obs: " << i << ": " << origData[i] << "\n" ;

  // print out original obs info
  // for (int i=0; i<numOrigObs; ++i)
  //   DEBUGPR(1, ucout << "obs: " << i << ": " << origData[i] << "\n" );

  // if (readShuffledObs())
  //   readRandObs(argc, argv);

  // DEBUGPRX(2, this, "setupProblem: \n");
  // DEBUGPRX(2, this, tc.endCPUTime());
  // DEBUGPRX(2, this, tc.endWallTime());

  return true;
}


// TODO: not using this any more
// read shuffled observation from the data file
// bool DataRMA::readRandObs(int argc, char **argv) {
//
//   ucout << "Use Shuffled Obs\n";
//   ifstream s(argv[2]); // open the data file
//
//   // check whether or not the file is opened correctly
//   if (!s) {
//     cerr << "Could not open file \"" << argv[2] << "\"\n";
//     return false;
//   }
//
//   vecRandObs.resize(numOrigObs);
//
//   // read data
//   for (unsigned int i = 0; i < numOrigObs; ++i)
//     s >> vecRandObs[i];
//
//   s.close(); // close the data file
//
//   if (args->debug >= 2)
//     cout << "vecRandObs: " << vecRandObs;
//   return true;
// }


// remove observations with zero weight
void DataRMA::removeZeroWtObs() {

  int numNonZeroObs = 0;

  for (unsigned int i = 0; i < numTrainObs; ++i) { // for each training observation
    if (dataIntTrain[i].w != 0) {  // if objservation i has non-zero wiehgt
      vecTrainObsIdx[numNonZeroObs] = i;
      ++numNonZeroObs;
    }
  }

  numTrainObs = numNonZeroObs;
  vecTrainObsIdx.resize(numTrainObs);

}


// read non uniform weight and set the weight for each observation
void DataRMA::readNonUniformWt() {

  vector<double> vecNonUniformWt;
  vecNonUniformWt.resize(numTrainObs);

  // TODO: for now, each process is reading non-uniform weight file
  /*
#ifdef ACRO_HAVE_MPI
  if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
  */
  ifstream inFile(args->nonUniformWt()); // Input stream class to operate on non-uniform weight file
  if (inFile.is_open()) { // if the file is open
    string line, tmp;
    // read a string or a line from an input stream
    // inFile is an object of istream class
    // line is a string object that input is stored in this object after being read from the stream
    while (getline(inFile, line)) {
      stringstream ss(line);  // Used for breaking words
      for (unsigned int i = 0; i < numTrainObs; ++i) { // for each training observation
        getline(ss, tmp, ',');
        // cout << "tmp " << tmp << "\n";
        vecNonUniformWt[i] = stod(tmp);  // stod: convert string to double
        // cout << "vec " <<vecNonUniformWt[i] << "\n";
      } // end for each train observation
    } // end while
  } else { // if the program could not open the file
    cout << "error: cannot read nonuniform wt";
  }
  // rma->setDataIntTrainWeight(vecNonUniformWt, vecObsIdx);

  for (unsigned int i = 0; i < numTrainObs; ++i) { // for each training observation
    dataIntTrain[i].w = vecNonUniformWt[i]; // assigned the non-uniform weight
  }
  /*
#ifdef ACRO_HAVE_MPI
}
#endif //  ACRO_HAVE_MPI
  */

  if (args->debug >= 10) {
    ucout << "rank: " << uMPI::rank << " wt: ";
    for (unsigned int i = 0; i < numTrainObs; ++i) // for each observation
      ucout << dataIntTrain[i].w << ", ";          // print out weight
    ucout << "\n";
  } // end if

}


// set data dimensions
void DataRMA::setDataDimensions() {

  numTrainObs = numOrigObs; // TODO: need to be fixed

  if (args->debug >= 10)
    cout << "numTrainObs: " << numTrainObs << "\n";

  dataIntTrain.resize(numTrainObs);
  dataStandTrain.resize(numTrainObs);

  for (unsigned int i = 0; i < numTrainObs; ++i) {
    dataIntTrain[i].X.resize(numAttrib);
    dataStandTrain[i].X.resize(numAttrib);
  }


  vecFeature.resize(numAttrib);

  vecTrainObsIdx.resize(numTrainObs);
  for (unsigned int i = 0; i < numTrainObs; ++i)
    vecTrainObsIdx[i] = i;
}

// set X for the dataIntTrain
void DataRMA::setDataIntTrainX() {

  // set X values
  if (args->delta() != -1) // integerize data
    integerizeData(dataOrigTrain, dataIntTrain);
  else  // assuming that the original data is already integerized
    for (unsigned int i = 0; i < numTrainObs; ++i) { // for each observation
      for (unsigned int j = 0; j < numAttrib; ++j) { // for each attribute
        // assign the original X as integerized X
        dataIntTrain[i].X[j] = dataOrigTrain[i].X[j];
      } // end for each attribute
    } // end for observation

} // end setDataIntTrainX function


// set weights of dataIntTrain
void DataRMA::setDataIntTrainWeight() {
  if (args->nonUniformWt() != "") { // if the nonUniform weight file is given
    readNonUniformWt();  // set the weights from the file
  } else {
    // give equal weight for each observation (1/numTrainObs)
    for (unsigned int i = 0; i < numTrainObs; ++i) // for each observation
      dataIntTrain[i].w = dataOrigTrain[i].y * 1.0 / (double)numTrainObs;
  }
}


// set vecNumDistFeats (# of distinc values for each attributes)
void DataRMA:: setVecNumDistFeats() {
  vecNumDistFeats.resize(numAttrib);
  for (unsigned int j = 0; j < numAttrib; ++j) { // for each attribute
    vecNumDistFeats[j] = 0;
    for (unsigned int i = 0; i < numTrainObs; ++i)  // for each observation
      if (vecNumDistFeats[j] < dataIntTrain[i].X[j] + 1)
        vecNumDistFeats[j] = dataIntTrain[i].X[j] + 1;
  } // end for each attribute
}


// set the maximum distinct value of all attributes
void DataRMA::setMaxNumDistFeats() {
  maxNumDistFeats     = 0;
  for (unsigned int j = 0; j < numAttrib; ++j)  // for each attribute
    if (maxNumDistFeats < vecNumDistFeats[j])
      maxNumDistFeats = vecNumDistFeats[j];
}


// set the maximum distinct value of all attributes
void DataRMA::setNumTotalCutPts() {
  numTotalCutPts = 0;
  for (unsigned int j = 0; j < numAttrib; ++j)  // for each attribute
    numTotalCutPts += vecNumDistFeats[j]-1; // sum up the tatoal cut point
}


// set positive and negative observations
void DataRMA::setNumPosNegObs() {

  numPosTrainObs = 0;
  numNegTrainObs = 0;

  for (unsigned int i = 0; i < numTrainObs; ++i) { // for each training observation
    if (dataOrigTrain[vecTrainObsIdx[i]].y == 1) // if the response value is 1
      ++numPosTrainObs;  // count observation i as postive
    else
      ++numNegTrainObs;  // count observation i as negative
  }

#ifdef ACRO_HAVE_MPI
  if (uMPI::rank == 0) {
#endif //  ACRO_HAVE_MPI
    ucout << "# of Positive : Negative Observations: "
          << numPosTrainObs << " : " << numNegTrainObs << "\n";
#ifdef ACRO_HAVE_MPI
  }
#endif //  ACRO_HAVE_MPI
}


// set vecAvgX, a vector of average X for each attribute
void DataRMA::setVecAvgX(vector<DataXy> &origData) { // TODO: why passing this data

  unsigned int i, j;
  vecAvgX.resize(numAttrib);

  for (j = 0; j < numAttrib; ++j)
    vecAvgX[j] = 0;

  for (j = 0; j < numAttrib; ++j) {   // for each attribute j
    for (i = 0; i < numTrainObs; ++i) // for each training observation i
      vecAvgX[j] += origData[i].X[j]; // set the sum of X value for attribute j
    vecAvgX[j] /= numTrainObs;        // set the average of X value for attribute j
  } // for each attribute j

} // end setVecAvgX function


// set vecSdX, a vector of standard deviation of X for each attribute
void DataRMA::setVecSdX(vector<DataXy> &origData) {

  unsigned int i, j;
  vecSdX.resize(numAttrib);

  for (j = 0; j < numAttrib; ++j)
    vecSdX[j]  = 0;

  // set the vecAvgX and vecSdX for all attributes
  for (j = 0; j < numAttrib; ++j) {    // for each attribute j
    for (i = 0; i < numTrainObs; ++i)  // for each training observation
      vecSdX[j] += pow(origData[i].X[j] - vecAvgX[j], 2);
    vecSdX[j] /= numTrainObs;
    vecSdX[j] = sqrt(vecSdX[j]);
  } // end for each attribute j

} // end setVecSdX function


// set avgY, the average of Y values
void DataRMA::setAvgY(vector<DataXy> &origData) {

  unsigned int obs;
  avgY = 0;

  for (unsigned int i = 0; i < numTrainObs; ++i) {
    obs = vecTrainObsIdx[i];
    avgY += origData[obs].y; // get avg of y
  }

  avgY /= numTrainObs; // get average response value

}


// set sdY, the standard devication of Y values
void DataRMA::setSdY(vector<DataXy> &origData) {

  unsigned int obs;
  sdY = 0;

  for (unsigned int i = 0; i < numTrainObs; ++i){
    obs = vecTrainObsIdx[i];
    sdY += pow(origData[obs].y - avgY, 2); // get std dev of y
  }

  sdY /= numTrainObs;
  sdY = sqrt(sdY);

}


// set Standadize Data
void DataRMA::setStandDataX(vector<DataXy> &origData,
                            vector<DataXy> &standData) {

  unsigned int i, j, obs;

  setVecAvgX(origData);
  setVecSdX(origData);

  // standardize X in each attribute
  for (j = 0; j < numAttrib; ++j) // for each attribute
    for (i = 0; i < numTrainObs; ++i) { // for each training observation
      obs = vecTrainObsIdx[i];
      standData[obs].X[j] = (origData[obs].X[j] - vecAvgX[j]) / vecSdX[j];
    }

  // if (args->debug >= 100)
  //   // print the standData
  //   for (i = 0; i < numTrainObs; ++i) {
  //     obs = vecTrainObsIdx[i];
  //     cout << "obs: " << obs << ": " << standData[obs] << "\n";
  //   }

}


// set the standardized data for Y-value
void DataRMA::setStandDataY(vector<DataXy> &origData,
                            vector<DataXy> &standData) {

  unsigned int obs;

  setAvgY(origData); // set avgY
  setSdY(origData);  // set sdY

  // standardize y
  for (unsigned int i = 0; i < numTrainObs; ++i) { // for each obseration i
    obs = vecTrainObsIdx[i];
    standData[obs].y = (origData[obs].y - avgY) / sdY;
  }

}


void DataRMA::integerizeData(vector<DataXy> &origData,
                             vector<DataXw> &intData) {

  bool isSplit;
  unsigned int i, j, k, l, r, p, q, o, obs;
  double tmpL, tmpU, tmpL1, tmpU1, tmp1U;

  double interval;  // confidence interval range
  double eps, eps0; // episilon, aggregation level

  vector<double> vecTemp(numOrigObs);

  set<double> setDistVal;        // a set continas all distinct values for each attribute
  set<double>::iterator it, itp; // iterator for the set

  map<double, int> mapDblInt; // a container maps from an original value to an
                              // integeried value
  map<double, int>::iterator itm; // iterator for the map

  vector<IntMinMax> copyIntMinMax; // a vector contins min and max for each integerized value

  tc.startTime();

  vecMaxX.resize(numAttrib);
  vecMinX.resize(numAttrib);

  // reset the minimum and maximum values of X for each attribute
  for (j = 0; j < numAttrib; ++j) {
    vecMinX[j] = getInf();
    vecMaxX[j] = -getInf();
  }

  // if (isLPBoost()) setXStat();

  for (j = 0; j < numAttrib; ++j) { // for each attribute

    if (args->debug >= 2)
      cout << "feat: " << j << "\n";
    setDistVal.clear();
    cout << "test: ";
    for (i = 0; i < numTrainObs; ++i) { // for each training observation
      obs = vecTrainObsIdx[i];  // get the observation index
      if (args->debug >= 10)
        cout << "test: " << origData[obs].X[j] << "\n";
      setDistVal.insert(origData[obs].X[j]);
    }

    if (args->debug >= 10)
      cout << "size setDistVal: " << setDistVal.size() << "\n";

    if (args->debug >= 2) {
      cout << "setDistVal: ";
      for (it = setDistVal.begin(); it != setDistVal.end(); ++it)
        cout << *it << " ";
      cout << '\n';
    }

    // get 95% confidence interval range
    interval = min(4.0 * vecSdX[j], *setDistVal.rbegin() - *setDistVal.begin());

    // episiolon, aggregation level, for integerization
    eps = min(args->delta(), args->maxInterval()) * interval;

    eps0 = eps;
    if (args->debug >= 2)
      cout << "delta: " << args->delta() << "\n"
           << "max: "   << *setDistVal.rbegin()
           << ", min: " << *setDistVal.begin() << "\n"
           << "eps: "   << eps << "\n"
           << "maxInterval: " << args->maxInterval() * interval << "\n";

    /************ assign integer without recursive integerization ************/
    k = 0;
    mapDblInt.clear();
    vecFeature[j].vecIntMinMax.resize(setDistVal.size());

    // the min value is equal to the maximum value for each integer assigned
    itp = setDistVal.begin();
    vecFeature[j].vecIntMinMax[0].minOrigVal = *itp;
    vecFeature[j].vecIntMinMax[0].maxOrigVal = *itp;

    // walk thorugh the set of distincet values
    // some value can be aggregated by the level of the episilon
    for (it = setDistVal.begin(); it != setDistVal.end(); ++it) {
      if (args->debug >= 2)
        cout << "tmpL: " << *itp << " tmpU: " << *it
             << " diff: " << (*it - *itp) << "\n";
      if ((*it - *itp) > eps) { // aggregating some value
        vecFeature[j].vecIntMinMax[++k - 1].maxOrigVal = *(--it);
        vecFeature[j].vecIntMinMax[k].minOrigVal = *(++it);
      }
      itp = it;
      mapDblInt[*it] = k;
    }
    vecFeature[j].vecIntMinMax[k].maxOrigVal = *(--it);

    if (args->debug >= 2) {
      cout << "mapDblInt contains:";
      for (itm = mapDblInt.begin(); itm != mapDblInt.end(); ++itm)
        cout << " [" << itm->first << ':' << itm->second << ']';
      cout << '\n';
    }

    vecFeature[j].vecIntMinMax.resize(k + 1);
    vecNumDistFeats[j] = k; // get distinct # of feature

    /************************ recursive integerization ************************/

    // if there is interval limit
    // and the size of distince value is not same as the number of integers
    // assigned
    if (args->maxInterval() != getInf() || k != setDistVal.size() - 1) {

      copyIntMinMax.resize(k + 1);
      for (i = 0; i <= k; ++i) {
        copyIntMinMax[i].minOrigVal = vecFeature[j].vecIntMinMax[i].minOrigVal;
        copyIntMinMax[i].maxOrigVal = vecFeature[j].vecIntMinMax[i].maxOrigVal;
      }

      if (args->debug >= 2) {
        cout << "\nvecIntMin ";
        for (i = 0; i <= k; ++i)
          cout << copyIntMinMax[i].minOrigVal << ' ';
        cout << "\nvecIntMax ";
        for (i = 0; i <= k; ++i)
          cout << copyIntMinMax[i].maxOrigVal << ' ';
        cout << '\n';
      }

      p = 0;
      for (i = 0; i <= k; ++i) {

        isSplit = true;
        eps = eps0;
        r = 0;

        tmpL = copyIntMinMax[i].minOrigVal;
        tmpU = copyIntMinMax[i].maxOrigVal;

        while ((tmpU - tmpL) > args->maxInterval() * interval && isSplit &&
               eps > .0001) {

          isSplit = false;
          eps *= args->shrinkDelta();

          if (args->debug >= 2)
            cout << "new eps: " << eps << '\n';

          for (q = 0; q <= r; ++q) {
            l = 0;
            tmpL1 = vecFeature[j].vecIntMinMax[i + p + q].minOrigVal;
            tmpU1 = vecFeature[j].vecIntMinMax[i + p + q].maxOrigVal;

            if (args->debug >= 2)
              cout << " q: " << q << " tmpL2: " << tmpL1 << " tmpU2: " << tmpU1
                   << " diff: " << tmpU1 - tmpL1 << "\n";

            if ((tmpU1 - tmpL1) < 0) {

              if (args->debug >= 2) {

                cout << "!!!!!!!!!!Something Wrong!!!!!!!!!!!\n";

                cout << "\nvecIntMin2 ";
                for (o = 0; o <= k + p; ++o)
                  cout << vecFeature[j].vecIntMinMax[o].minOrigVal << ' ';
                cout << "\nvecIntMax2 ";
                for (o = 0; o <= k + p; ++o)
                  cout << vecFeature[j].vecIntMinMax[o].maxOrigVal << ' ';
                cout << '\n';
              }

            } else if ((tmpU1 - tmpL1) > args->maxInterval() * interval) {

              isSplit = true;

              for (it = setDistVal.find(tmpL1);; ++it) {
                tmp1U = *it;
                if (args->debug >= 2)
                  cout << "tmpL1: " << tmpL1 << " tmpU1: " << tmp1U
                       << " diff: " << tmp1U - tmpL1 << "\n";

                if ((tmp1U - tmpL1) > eps) {
                  ++l;
                  ++r;
                  vecFeature[j].vecIntMinMax[i + p + l + q - 1].maxOrigVal =
                      tmpL1;
                  vecFeature[j].vecIntMinMax[i + p + l + q].minOrigVal = tmp1U;

                  if (args->debug >= 2) {
                    cout << " idx: " << i + p + l + q - 1 << " tmpL4: "
                         << vecFeature[j]
                                .vecIntMinMax[i + p + l + q - 1]
                                .maxOrigVal
                         << " tmpU4: "
                         << vecFeature[j].vecIntMinMax[i + p + l + q].minOrigVal
                         << "\n";
                    cout << " i: " << i << "p: " << p << " r: " << r
                         << " l: " << l << " q: " << q << "\n";
                  }
                }

                tmpL1 = tmp1U;
                vecFeature[j].vecIntMinMax[i + p + l + q].maxOrigVal = tmpU;

                if (tmp1U == tmpU1)
                  break;

              } // end for each inner sub interval
            }   // end if each interval is less than the threthold
          }     // end for (p=0; p<=r; ++p)
        } // end while ( (tmpU-tmpL) > getLimitInterval()*interval && isSplit)

        p += r;

        if ((tmpU - tmpL) <= args->maxInterval() * interval && p > 0) {
          vecFeature[j].vecIntMinMax[i + p].minOrigVal =
              copyIntMinMax[i].minOrigVal;
          vecFeature[j].vecIntMinMax[i + p].maxOrigVal =
              copyIntMinMax[i].maxOrigVal;
        }

      } // end for (i=0; i<=k; ++i), each original interval

      if (args->debug >= 2) {
        cout << "\nvecIntMin1 ";
        for (i = 0; i <= k + p; ++i)
          cout << vecFeature[j].vecIntMinMax[i].minOrigVal << ' ';
        cout << "\nvecIntMax1 ";
        for (i = 0; i <= k + p; ++i)
          cout << vecFeature[j].vecIntMinMax[i].maxOrigVal << ' ';
        cout << '\n';
      }

      o = 0;
      for (it = setDistVal.begin(); it != setDistVal.end(); ++it) {
        if (*it > vecFeature[j].vecIntMinMax[o].maxOrigVal)
          ++o;
        mapDblInt[*it] = o;
      }

      if (args->debug >= 2) {
        cout << "mapDblInt1 contains:";
        for (itm = mapDblInt.begin(); itm != mapDblInt.end(); ++itm)
          cout << " [" << itm->first << ':' << itm->second << ']';
        cout << '\n';
      }

      vecFeature[j].vecIntMinMax.resize(k + p + 1);
      vecNumDistFeats[j] = k + p; // get distinct # of feature

    } // end if recursive discretization applies

    // set intData sets
    for (i = 0; i < numTrainObs; ++i) {
      obs = vecTrainObsIdx[i];
      intData[obs].X.resize(numAttrib);
      intData[obs].X[j] = mapDblInt[origData[obs].X[j]];
    }

  } // end for each attr// set vecAvgX and vecSdX, average and standard devication vectors for X

  /*
  for (i=0; i<numTrainObs ; ++i) {
    obs = vecTrainObsIdx[i];
    DEBUGPRX(20, this, "IntObs: " << obs << ": "
      << intData[obs] << '\n');
  }*/

  if (args->debug >= 1)
    cout << "vecNumDistFeats: " << vecNumDistFeats << "\n";

  /*
    #ifdef ACRO_HAVE_MPI
      if (uMPI::rank==0) {
    #endif //  ACRO_HAVE_MPI
      if (writePred()) {
        saveXObs(dataOrigTrain);
        saveXObs(dataIntTrain);
      }
    #ifdef ACRO_HAVE_MPI
      }
    #endif //  ACRO_HAVE_MPI
  */
  maxNumDistFeats = 0;
  for (j = 0; j < numAttrib; ++j) {
    numTotalCutPts += vecNumDistFeats[j];
    if (maxNumDistFeats - 1 < vecNumDistFeats[j])
      maxNumDistFeats = vecNumDistFeats[j] + 1;
  }

  ////////////////////////////////////////////////////////////////////////////
  if (args->debug >= 1)
    for (i = 0; i < numTrainObs; ++i) {
      obs = vecTrainObsIdx[i];
      cout << "obs: " << obs << ": " << intData[obs] << "\n";
    }

  if (args->debug >= 1) {
    cout << "integerizeProblem: \t";
    tc.getCPUTime();
  }
  if (args->debug >= 2)
    tc.getWallTime();

} // end integerizeData


//  integerize into fixed bin
void DataRMA::integerizeFixedLengthData(vector<DataXy> &origData,
                                        vector<DataXw> &intData) {

  unsigned int i, j, obs, glMaxL = -1;
  unsigned int sizeBin = args->fixedSizeBin();
  maxNumDistFeats = 0;

  // fix X matrix
  for (i = 0; i < numTrainObs; ++i) {
    obs = vecTrainObsIdx[i];
    for (j = 0; j < numAttrib; ++j) {
      if (dataStandTrain[obs].X[j] < vecMinX[j])
        vecMinX[j] = dataStandTrain[obs].X[j]; // get vecMinX[j]
      if (dataStandTrain[obs].X[j] > vecMaxX[j])
        vecMaxX[j] = dataStandTrain[obs].X[j]; // get vecMaxX[j]
    }
  }

  vecNumDistFeats.resize(numAttrib);
  for (j = 0; j < numAttrib; ++j) {
    maxNumDistFeats = -1;
    for (i = 0; i < numTrainObs; ++i) {
      obs = vecTrainObsIdx[i];
      intData[obs].X.resize(numAttrib);
      intData[obs].X[j] = floor((origData[obs].X[j] - vecMinX[j]) /
                                ((vecMaxX[j] - vecMinX[j]) / (double)sizeBin));
      if (maxNumDistFeats < intData[obs].X[j])
        maxNumDistFeats = intData[obs].X[j];
    }

    vecNumDistFeats[j] = maxNumDistFeats;
    if (glMaxL < maxNumDistFeats)
      glMaxL = maxNumDistFeats;

    vecFeature[j].vecIntMinMax.resize(maxNumDistFeats);
    for (i = 0; i < maxNumDistFeats; ++i) {
      vecFeature[j].vecIntMinMax[0].minOrigVal =
          (double)i * ((vecMaxX[j] - vecMinX[j]) / (double)sizeBin) + vecMinX[j];
      vecFeature[j].vecIntMinMax[0].maxOrigVal =
          (double)(i + 1) * ((vecMaxX[j] - vecMinX[j]) / (double)sizeBin) + vecMinX[j];
    }
  }
}


// TODO: need this?
// save X values of all the training observations
template <class T> void DataRMA::saveXObs(T vecData) {

  unsigned int i, j, obs;
  stringstream s;
  (typeid(T) == typeid(int)) ? s << "int" << '.' : s << "orig" << '.';
  ofstream os(s.str().c_str());

  for (i = 0; i < numTrainObs; ++i) { // for each training observation
    for (j = 0; j < numAttrib; ++j) { // for each attribute
      obs = vecTrainObsIdx[i];
      os << vecData[obs].X[j] << " ";
    }
    os << "\n";
  }
  os.close();
}

} // namespace data


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
