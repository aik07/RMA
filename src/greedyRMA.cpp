/**********************************************************
 * File name:   greedyRMA.cpp
 * Author:      Ai Kagawa
 * Description: a source file for greedy RMA solver
 ***********************************************************/

#include "greedyRMA.h"

namespace greedyRMA {

void GreedyRMA::runGreedyRangeSearch() {

  ts.startTime();

  getMinOptRange(); // get oprimal range for minimum objective value

  getMaxOptRange(); // get oprimal range for maximum objective value

  chooseMinOrMaxRange();

  printSolution();

} // end function greedyRMA

void GreedyRMA::printSolution() {

#ifdef ACRO_HAVE_MPI
  if (uMPI::rank == 0) {
#endif //  ACRO_HAVE_MPI

    std::cout << "GRMA Solution: ";
    isPosIncumb ? ucout << "+" : ucout << "-";
    std::cout << std::fixed << std::setprecision(4)
              << maxObjValue << "\t";
    std::cout << std::fixed << std::setprecision(2)
              << "CPU Time: " << ts.getCPUTime() << "\n";
    if (args->debug >= 2)
      std::cout << ts.getWallTime();
    // if (args->printBoost()) {
    if (args->debug >= 2)
      std::cout << "L: " << L << "U: " << U;
      //}

#ifdef ACRO_HAVE_MPI
  }
#endif //  ACRO_HAVE_MPI
}

/************************** Final Optimal Range **************************/
void GreedyRMA::chooseMinOrMaxRange() {

  L.resize(data->numAttrib);
  U.resize(data->numAttrib);

  (args->isRandSeed())
      ? srand((NumNegTiedSols + NumPosTiedSols) * time(NULL) * 100)
      : srand(1);
  double rand_num = (rand() % 10001) / 10000.0;
  // DEBUGPRX(0, global(), "rand: " << rand_num  << "\n");
  // DEBUGPRX(0, global(), "rand1: " << 1.0 /  NumTiedSols << "\n");

  if (maxVal > -minVal ||
      (maxVal == minVal &&
       rand_num <=
           NumPosTiedSols / (double)(NumNegTiedSols + NumPosTiedSols))) {

    maxObjValue = maxVal;
    isPosIncumb = true;

    copy(Lmax.begin(), Lmax.end(), L.begin());
    copy(Umax.begin(), Umax.end(), U.begin());

    if (args->debug >= 10)
      ucout << "chose max\n";

  } else {

    maxObjValue = -minVal;
    isPosIncumb = false;

    copy(Lmin.begin(), Lmin.end(), L.begin());
    copy(Umin.begin(), Umin.end(), U.begin());

    if (args->debug >= 10)
      ucout << "chose min\n";
  }

  if (args->debug >= 10)
    ucout << "Lmin: " << Lmin << "Umin: " << Umin;
  if (args->debug >= 10)
    ucout << "Lmax: " << Lmax << "Umax: " << Umax;
}

/********************* Optimal Range for Min Objective *********************/
void GreedyRMA::getMinOptRange() {

  tmpMin = getInf();
  minVal = getInf();

  optAttrib = -1;
  prevAttrib = -1;

  Lmin.clear();
  Lmin.resize(data->numAttrib);

  Umin.resize(data->distFeat.size());
  copy(data->distFeat.begin(), data->distFeat.end(), Umin.begin());

  vecCoveredObs.resize(data->numTrainObs);
  copy(data->vecTrainData.begin(), data->vecTrainData.end(),
       vecCoveredObs.begin());

  do {

    fondNewBox = false;

    for (unsigned int j = 0; j < data->numAttrib; ++j) { // for each feature

      if (j != prevAttrib) { // if this attribute is not restricted

        setObjVec(j);
        tmpMin = runMinKadane(j);

        if (tmpMin == minVal) {

          NumNegTiedSols++;
          (args->isRandSeed()) ? srand(NumNegTiedSols * time(NULL) * 100)
                             : srand(1);
          double rand_num = (rand() % 10001) / 10000.0;
          // DEBUGPRX(0, global(), "rand: " << rand_num  << "\n");
          // DEBUGPRX(0, global(), "rand1: " << 1.0 /  NumTiedSols << "\n");

          if (rand_num <= 1.0 / NumNegTiedSols)
            setOptMin(j);

        } else if (tmpMin < minVal) {
          NumNegTiedSols = 1;
          setOptMin(j);
        }

      } // end if this attribute is not restricted

    } // end each feature

    if (fondNewBox) {

      dropObsNotCovered(optAttrib, optLower, optUpper);

      prevAttrib = optAttrib;
      Lmin[optAttrib] = optLower;
      Umin[optAttrib] = optUpper;

      if (args->debug >= 10)
        ucout << "vecCoveredObs: " << vecCoveredObs
              << "final optAttrib: (a,b): " << optAttrib << ": (" << optLower
              << ", " << optUpper << "), min: " << minVal << "\n";
    }

  } while (fondNewBox);
}

/********************* Optimal Range for Max Objective *********************/
void GreedyRMA::getMaxOptRange() {

  tmpMax = -getInf();
  maxVal = -getInf();

  optAttrib = -1;
  prevAttrib = -1;

  Lmax.clear();
  Lmax.resize(data->numAttrib);

  Umax.resize(data->distFeat.size());
  copy(data->distFeat.begin(), data->distFeat.end(), Umax.begin());

  vecCoveredObs.resize(data->vecTrainData.size());
  copy(data->vecTrainData.begin(), data->vecTrainData.end(),
       vecCoveredObs.begin());

  do {

    fondNewBox = false;

    for (unsigned int j = 0; j < data->numAttrib; ++j) { // for each feature

      if (j != prevAttrib) { // if this attribute is not restricted

        setObjVec(j);
        tmpMax = runMaxKadane(j);

        if (tmpMax == maxVal) {

          NumPosTiedSols++;
          (args->isRandSeed()) ? srand(NumPosTiedSols * time(NULL) * 100)
                             : srand(1);
          double rand_num = (rand() % 10001) / 10000.0;
          // DEBUGPRX(0, global(), "rand: " << rand_num  << "\n");
          // DEBUGPRX(0, global(), "rand1: " << 1.0 /  NumTiedSols << "\n");

          if (rand_num <= 1.0 / NumPosTiedSols)
            setOptMax(j);

        } else if (tmpMax > maxVal) {
          NumPosTiedSols = 1;
          setOptMax(j);
        }

      } // end if this attribute is not restricted
    }   // end for each attribute

    if (fondNewBox) { // if better solution found

      dropObsNotCovered(optAttrib, optLower, optUpper);

      prevAttrib = optAttrib;
      Lmax[optAttrib] = optLower;
      Umax[optAttrib] = optUpper;

      if (args->debug >= 10)
        ucout << "vecCoveredObs: " << vecCoveredObs
              << "final optAttrib: (a,b): " << optAttrib << ": (" << optLower
              << ", " << optUpper << "), max: " << maxVal << "\n";
    }

  } while (fondNewBox);
}

/********************* Optimal Range for 1D rules  *********************/
/*
void GreedyRMA::setInit1DRules() {

  Lmin.clear();
  Lmin.resize(data->numAttrib);

  Umin.resize(data->distFeat.size());
  copy(data->distFeat.begin(), data->distFeat.end(), Umin.begin());

  Lmax.clear();
  Lmax.resize(data->numAttrib);

  Umax.resize(data->distFeat.size());
  copy(data->distFeat.begin(), data->distFeat.end(), Umax.begin());

  vecCoveredObs.resize(data->numTrainObs);
  copy(data->vecTrainData.begin(), data->vecTrainData.end(),
vecCoveredObs.begin());
}


void GreedyRMA::set1DOptRange(const int& j) {
  setObjVec(j);
  tmpMin = runMinKadane(j);
  maxObjValue = tmpMin;
  optLower = tmpL;
  optUpper = tmpU;
  isPosIncumb = false;
  tmpMax = runMaxKadane(j);
  if (-tmpMin <= tmpMax) {
    maxObjValue = tmpMax;
    optLower = tmpL;
    optUpper = tmpU;
    isPosIncumb = true;
  }
}
*/

void GreedyRMA::setOptMin(const unsigned int &j) {

  minVal = tmpMin;
  fondNewBox = true;

  optAttrib = j;
  optLower = tmpL;
  optUpper = tmpU;

  if (args->debug >= 10)
    ucout << "optAttrib: (a,b): " << optAttrib << ": (" << optLower << ", "
          << optUpper << "), min: " << minVal << "\n";
}

void GreedyRMA::setOptMax(const unsigned int &j) {

  maxVal = tmpMax;
  fondNewBox = true;

  optAttrib = j;
  optLower = tmpL;
  optUpper = tmpU;

  if (args->debug >= 10)
    ucout << "optAttrib: (a,b): " << optAttrib << ": (" << optLower << ", "
          << optUpper << "), max: " << maxVal << "\n";
}

// get Miniumum range for the feature
double GreedyRMA::runMinKadane(const unsigned int &j) {

  int s = Lmin[j];
  tmpL = Lmin[j];
  tmpU = Umin[j];

  double minEndHere = 0;
  double minSoFar = getInf();

  for (unsigned int i = Lmin[j]; i <= Umin[j]; ++i) {
    minEndHere += vecWeight[i]; // getObjCovered(j, i);
    if (minEndHere < minSoFar) {
      minSoFar = minEndHere;
      tmpL = s;
      tmpU = i;
    }
    if (minEndHere > 0) {
      minEndHere = 0;
      s = i + 1;
    }
  }

  if (args->debug >= 10)
    ucout << "Minimum contiguous sum is " << minSoFar << " "
          << " feat (L,U): " << j << " (" << tmpL << ", " << tmpU << ")\n";
  return minSoFar;
}

// get Maximum range for the feature
double GreedyRMA::runMaxKadane(const unsigned int &j) {

  int s = Lmax[j];
  tmpL = Lmax[j];
  tmpU = Umax[j];

  double maxEndHere = 0;
  double maxSoFar = -getInf(); // min so far

  for (unsigned int i = Lmax[j]; i <= Umax[j]; ++i) {
    maxEndHere += vecWeight[i]; // getObjCovered(j, i);
    if (maxEndHere > maxSoFar) {
      maxSoFar = maxEndHere;
      tmpL = s;
      tmpU = i;
    }
    if (maxEndHere < 0) {
      maxEndHere = 0;
      s = i + 1;
    }
  }

  if (args->debug >= 10)
    ucout << "Maximum contiguous sum is " << maxSoFar << " feat (L,U): " << j
          << " (" << tmpL << ", " << tmpU << ")\n";
  return maxSoFar;
}

void GreedyRMA::dropObsNotCovered(const unsigned int &j, const unsigned int &lower,
                                  const unsigned int &upper) {

  unsigned int obs;
  int l = -1;
  if (args->debug >= 10)
    ucout << "Before drop: " << vecCoveredObs;

  for (unsigned int i = 0; i < vecCoveredObs.size(); ++i) {

    obs = vecCoveredObs[i];

    if (lower <= data->intTrainData[obs].X[j] &&
        data->intTrainData[obs].X[j] <= upper)
      //&& data->intTrainData[obs].w!=0)
      vecCoveredObs[++l] = obs; // store covered observations
  }

  vecCoveredObs.resize(l + 1); // shrink the size of vecCoveredObs
  if (args->debug >= 10)
    ucout << "After drop: " << vecCoveredObs;
}

void GreedyRMA::setObjVec(const unsigned int &j) {

  unsigned int i, v, obs;

  vecWeight.resize(data->numMaxDistVal);
  for (i = 0; i < data->numMaxDistVal; ++i)
    vecWeight[i] = 0;

  for (i = 0; i < vecCoveredObs.size(); ++i) {
    obs = vecCoveredObs[i];
    v = data->intTrainData[obs].X[j];
    vecWeight[v] += data->intTrainData[obs].w;
  }//#include <acro_config.h>


  if (args->debug >= 10)
    ucout << "vecWeight: " << vecWeight;
}

} // namespace greedyRMA
