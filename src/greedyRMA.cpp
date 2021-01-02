/**********************************************************
 * Author:      Ai Kagawa
 * Description: a source file for greedy RMA solver
 ***********************************************************/


#include "greedyRMA.h"


namespace greedyRMA {


  void GreedyRMA::runGreedyRangeSearch() {

    ts.startTime();         // start the timer

    vecLowerMin.resize(data->numAttrib);
    vecUpperMin.resize(data->numAttrib);
    vecLowerMax.resize(data->numAttrib);
    vecUpperMax.resize(data->numAttrib);
    vecLower.resize(data->numAttrib);
    vecUpper.resize(data->numAttrib);
    vecValueWeight.resize(data->maxNumDistVals);

    searchMinOptRange();    // search oprimal range for minimum objective value

    searchMaxOptRange();    // search oprimal range for maximum objective value

    chooseMinOrMaxRange();  // choose the optimal solution by compareing the min and max versions

    printSolution();        // print out the solution

  } // end function greedyRMA


  /************************** Final Optimal Range **************************/
  void GreedyRMA::chooseMinOrMaxRange() {

    // if the absulute value of maximum objective value is larger than its mimum objective value
    //   or selecting maximum objective value when they are tied
    if (maxObjVal > -minObjVal || isChooseMaxWhenTiedMinMax() ) {

      optObjVal   = maxObjVal;  // set the maximum objective value from the max version of the range serach
      isPosIncumb = true;    // incumbent is the positive value

      // set the optimal lower and upper bounds of the box
      copy(vecLowerMax.begin(), vecLowerMax.end(), vecLower.begin());
      copy(vecUpperMax.begin(), vecUpperMax.end(), vecUpper.begin());

      if (args->debug >= 10)
        ucout << "chose max\n";

    } else {  // the optimal solutiona is from the min ver.

      optObjVal   = -minObjVal;  // set the minimun objective value from the min version of the range serach
      isPosIncumb = false;       // incumbent is not positive value

      // set the optimal lower and upper bounds of the box
      copy(vecLowerMin.begin(), vecLowerMin.end(), vecLower.begin());
      copy(vecUpperMin.begin(), vecUpperMin.end(), vecUpper.begin());

      if (args->debug >= 10)
        ucout << "chose min\n";

    }

    if (args->debug >= 10) {
      ucout << "Lower Min: " << vecLowerMin << "Upper Min: " << vecUpperMin;
      ucout << "Lower Max: " << vecLowerMax << "Upper Max: " << vecUpperMax;
    }

  }  // end function chooseMinOrMaxRange


  /********************* Optimal Range for Min Objective *********************/
  void GreedyRMA::searchMinOptRange() {

    resetMinOptRange();

    do {  // repeat while finding a new solution

      isFondNewBox = false;

      for (unsigned int j = 0; j < data->numAttrib; ++j) { // for each feature

        if (j != prevAttrib) { // if this attribute is not previously chosen

          // set the vector of total weights for each value of attribute j
          setVecValueWeight(j);

          // run objective value for attribut j
          runMinKadane(j);

          // if the current min objective value is equal to
          // the min objective value found so far
          if (curMinObjVal == minObjVal) {

            numNegTiedSols++; // count the negative tied solutions

            if (isUpdateOptSol(numNegTiedSols))
              setOptMin(j); // update the current optimal solutiona for min ver

          // if the current min objective value is less than
          // the min objective value found so far
          } else if (curMinObjVal < minObjVal) {
            numNegTiedSols = 1;
            setOptMin(j); // update the current optimal solutiona for min ver
          }

        } // end if this attribute is not restricted

      } // end each feature

      if (isFondNewBox) { // if a new solution is discovered

        // drop observations which are not covered
        // by the optimal attribute, and its optimal lower and upper bounds
        dropObsNotCovered(optAttrib, optLower, optUpper);

        // set the previous chosen attribute
        prevAttrib             = optAttrib;

        // update the lower and upper bounds for the optimal attribute
        vecLowerMin[optAttrib] = optLower;
        vecUpperMin[optAttrib] = optUpper;

        if (args->debug >= 10)
          ucout << "; final optAttrib: " << optAttrib
                << ";  (a,b): " << ": (" << optLower << ", " << optUpper
                << "), minObjVal: " << minObjVal << "\n";

      } // end if

    } while (isFondNewBox);

  } // end searchMinOptRange function


  /********************* Optimal Range for Max Objective *********************/
  void GreedyRMA::searchMaxOptRange() {

    resetMaxOptRange();

    do { // repeat while finding a new solution

      isFondNewBox = false;

      for (unsigned int j = 0; j < data->numAttrib; ++j) { // for each attribute j

        if (j != prevAttrib) { // if this attribute is not restricted

          setVecValueWeight(j);
          runMaxKadane(j);

          // if the current max objective value is equal to
          // the max objective value found so far
          if (curMaxObjVal == maxObjVal) {

            numPosTiedSols++;  // count the tied solutions

            if (isUpdateOptSol(numPosTiedSols))
              setOptMax(j); // update the current optimal solutiona for max ver

          // if the current max objective value is greater than
          // the max objective value found so far
          } else if (curMaxObjVal > maxObjVal) {
            numPosTiedSols = 1;
            setOptMax(j); // update the current optimal solutiona for max ver
          }

        } // end if this attribute is not restricted
      }   // end for each attribute

      if (isFondNewBox) { // if a better solution was discovered

        // drop observations which are not covered
        // by the optimal attribute, and its optimal lower and upper bounds
        dropObsNotCovered(optAttrib, optLower, optUpper);

        // set the previous chosen attribute
        prevAttrib             = optAttrib;

        // update the lower and upper bounds for the optimal attribute
        vecLowerMax[optAttrib] = optLower;
        vecUpperMax[optAttrib] = optUpper;

        if (args->debug >= 10)
          ucout << "; final optAttrib: " << optAttrib
                << ";  (a,b): " << ": (" << optLower << ", " << optUpper
                << "); maxObjVal: " << maxObjVal << "\n";
      }

    } while (isFondNewBox);

  } // enf function searchMaxOptRange


  // update the optimal attribute (j), objective value, and lower and upper bounds
  // for the minimization version
  void GreedyRMA::setOptMin(const unsigned int &j) {

    optAttrib    = j;
    minObjVal    = curMinObjVal;
    optLower     = curLower;
    optUpper     = curUpper;

    isFondNewBox = true;  // a new improved box found

    if (args->debug >= 10)
      ucout << "optAttrib: " << optAttrib
            << "; (a,b): (" << optLower << ", " << optUpper << ")"
            << "; minObjVal: " << minObjVal << "\n";
  }


  // update the optimal attribute (j), objective value, and lower and upper bounds
  // for the maximization version
  void GreedyRMA::setOptMax(const unsigned int &j) {

    optAttrib    = j;
    maxObjVal    = curMaxObjVal;
    optLower     = curLower;
    optUpper     = curUpper;

    isFondNewBox = true;  // a new improved solution found

    if (args->debug >= 10)
      ucout << "optAttrib: " << optAttrib
            << "; (a,b): (" << optLower << ", " << optUpper << ")"
            << "; maxObjVal: " << maxObjVal << "\n";

  }


  // set the current objective value, lower and upper bounds
  // by running the min ver. of Kadane's algorithm for attribute j
  void GreedyRMA::runMinKadane(const unsigned int &j) {

    int s    = vecLowerMin[j];
    curLower = vecLowerMin[j];
    curUpper = vecUpperMin[j];

    double minEndHere = 0;
    curMinObjVal = getInf();

    for (unsigned int i = vecLowerMin[j]; i <= vecUpperMin[j]; ++i) {
      minEndHere += vecValueWeight[i];
      if (minEndHere < curMinObjVal) {
        curMinObjVal = minEndHere;
        curLower = s;
        curUpper = i;
      }
      if (minEndHere > 0) {
        minEndHere = 0;
        s = i + 1;
      }
    }

    if (args->debug >= 10)
      ucout << "Minimum contiguous sum is " << curMinObjVal << " "
            << "; feat: " << j
            << "; (L,U): (" << curLower << ", " << curUpper << ")\n";

  } // enf runMinKadane function


  // set the current objective value, lower and upper bounds
  // by running the max ver. of Kadane's algorithm for attribute j
  void GreedyRMA::runMaxKadane(const unsigned int &j) {

    int s    = vecLowerMax[j];
    curLower = vecLowerMax[j];
    curUpper = vecUpperMax[j];

    double maxEndHere = 0;
    curMaxObjVal      = -getInf(); // min so far

    for (unsigned int i = vecLowerMax[j]; i <= vecUpperMax[j]; ++i) {
      maxEndHere += vecValueWeight[i];
      if (maxEndHere > curMaxObjVal) {
        curMaxObjVal = maxEndHere;
        curLower = s;
        curUpper = i;
      }
      if (maxEndHere < 0) {
        maxEndHere = 0;
        s = i + 1;
      }
    }

    if (args->debug >= 100)
      ucout << "Maximum contiguous sum is " << curMaxObjVal
            << "; feat: " << j
            << "; (L,U): (" << curLower << ", " << curUpper << ")\n";

  } // end runMaxKadane function


  // drop observations which are not covered for attribute j's lower and upper bound
  void GreedyRMA::dropObsNotCovered(const unsigned int &j, const unsigned int &lower,
                                    const unsigned int &upper) {

    unsigned int obs;
    int l = -1;

    if (args->debug >= 100)
      ucout << "Before drop vecCvdObsIdx: " << vecCvdObsIdx;

    for (unsigned int i = 0; i < vecCvdObsIdx.size(); ++i) {

      obs = vecCvdObsIdx[i];

      // if this observation is covered by the j'th attribute's lower and upper bounds
      if (lower <= data->dataIntTrain[obs].X[j] &&
          data->dataIntTrain[obs].X[j] <= upper)
        //&& data->dataIntTrain[obs].w!=0)
        vecCvdObsIdx[++l] = obs; // store covered observations
    }

    vecCvdObsIdx.resize(l + 1); // shrink the size of vecCvdObsIdx

    if (args->debug >= 100)
      ucout << "After drop: ";

    if (args->debug >= 50)
      ucout << "vecCvdObsIdx: " << vecCvdObsIdx << "\n";

  } // end dropObsNotCovered function


  // set the vector of total weights for each value of attribute j
  void GreedyRMA::setVecValueWeight(const unsigned int &j) {

    unsigned int i, v, obs;

    // initialize the vecValueWeight
    // vecValueWeight.resize(data->vecNumDistVals[j]);
    // fill(vecValueWeight.begin(), vecValueWeight.end(), 0);
    for (i = 0; i < data->vecNumDistVals[j]; ++i)
      vecValueWeight[i] = 0;

    for (i = 0; i < vecCvdObsIdx.size(); ++i) {  // for each covered observation index
      obs = vecCvdObsIdx[i];              // observation index
      v   = data->dataIntTrain[obs].X[j]; // X value of this observation for attribute j
      // add weights of observations for each X-value of attribute j
      vecValueWeight[v] += data->dataIntTrain[obs].w;
    }

    if (args->debug >= 10)
      ucout << "vecValueWeight: " << vecValueWeight;

  } // end setVecValueWeight function


  // reset variables for the minimum optimal range search
  void GreedyRMA::resetMinOptRange() {

    curMinObjVal = getInf();
    minObjVal    = getInf();

    optAttrib  = -1;
    prevAttrib = -1;

    // set lower vector for the minimim version
    // vecLowerMin.clear();
    // vecLowerMin.resize(data->numAttrib);
    fill(vecLowerMin.begin(), vecLowerMin.end(), 0);

    // set upper vector for the minimim version
    // vecUpperMin.resize(data->vecNumDistVals.size());
    copy(data->vecNumDistVals.begin(), data->vecNumDistVals.end(), vecUpperMin.begin());

    // set the covered observation indices
    vecCvdObsIdx.resize(data->numTrainObs);
    copy(data->vecTrainObsIdx.begin(), data->vecTrainObsIdx.end(),
       vecCvdObsIdx.begin());

  } // end resetMinOptRange function


  // reset variables for the maximum optimal range search
  void GreedyRMA::resetMaxOptRange() {

    curMaxObjVal = -getInf();
    maxObjVal    = -getInf();

    optAttrib  = -1;
    prevAttrib = -1;

    // set lower vector for the maximum version
    // vecLowerMax.clear();
    // vecLowerMax.resize(data->numAttrib);
    fill(vecLowerMax.begin(), vecLowerMax.end(), 0);

    // set upper vector for the maximum version
    // vecUpperMax.resize(data->vecNumDistVals.size());
    copy(data->vecNumDistVals.begin(), data->vecNumDistVals.end(), vecUpperMax.begin());

    vecCvdObsIdx.resize(data->vecTrainObsIdx.size());
    copy(data->vecTrainObsIdx.begin(), data->vecTrainObsIdx.end(),
         vecCvdObsIdx.begin());

  } // end resetMaxOptRange function


  // whether or not to choose the max solution
  // when min and max objective value are the same
  // (break the tie using a fair probability)
  bool  GreedyRMA::isChooseMaxWhenTiedMinMax() {

    // if the max and min objective values are the smae
    if (maxObjVal == minObjVal) {

      // set random seed if specified
      if (args->isRandSeed())
        srand((numNegTiedSols + numPosTiedSols) * time(NULL) * 100);
      else
        srand(1);

      // generate a random numeber
      double rand_num = (rand() % 10001) / 10000.0;

      // DEBUGPRX(0, global(), "rand: " << rand_num  << "\n");
      // DEBUGPRX(0, global(), "rand1: " << 1.0 /  NumTiedSols << "\n");

      if (rand_num <= numPosTiedSols / (double)(numNegTiedSols + numPosTiedSols))
        return true;
    }

    return false;

  } // end isChooseMaxWhenTiedMinMax function


  // whether or not to update the optimal solution
  // when there are tied solution for min or max versions
  // (break the tie using a fair probability based on # of tied solutions)
  bool  GreedyRMA::isUpdateOptSol(const int &numTiedSols) {

    //set random seed if specified
    (args->isRandSeed()) ? srand(numTiedSols * time(NULL) * 100) : srand(1);

    // generate a random numeber
    double rand_num = (rand() % 10001) / 10000.0;

    // DEBUGPRX(0, global(), "rand: " << rand_num  << "\n");
    // DEBUGPRX(0, global(), "rand1: " << 1.0 /  NumTiedSols << "\n");

    if (rand_num <= 1.0 / numTiedSols)
      return true;
    else
      return false;

  }


  // print the greedy RMA solution
  void GreedyRMA::printSolution() {

  #ifdef ACRO_HAVE_MPI
    if (uMPI::rank == 0) {
  #endif //  ACRO_HAVE_MPI

      std::cout << "GRMA Solution: ";
      isPosIncumb ? ucout << "+" : ucout << "-";
      std::cout << std::fixed << std::setprecision(4)
                << optObjVal << "\t";
      std::cout << std::fixed << std::setprecision(2)
                << "CPU Time: " << ts.getCPUTime() << "\n";
      if (args->debug >= 2)
        std::cout << ts.getWallTime();
      // if (args->printBoost()) {
      if (args->debug >= 2)
        std::cout << "Lower: " << vecLower << "\nUpper: " << vecUpper << "\n";
        //}

  #ifdef ACRO_HAVE_MPI
    }
  #endif //  ACRO_HAVE_MPI
  }


  /********************* Optimal Range for 1D rules  *********************/
  /*
  void GreedyRMA::setInit1DRules() {

    Lmin.clear();
    Lmin.resize(data->numAttrib);

    Umin.resize(data->vecNumDistVals.size());
    copy(data->vecNumDistVals.begin(), data->vecNumDistVals.end(), Umin.begin());

    vecLowerMax.clear();
    vecLowerMax.resize(data->numAttrib);

    vecUpperMax.resize(data->vecNumDistVals.size());
    copy(data->vecNumDistVals.begin(), data->vecNumDistVals.end(), vecUpperMax.begin());

    vecCvdObsIdx.resize(data->numTrainObs);
    copy(data->vecTrainObsIdx.begin(), data->vecTrainObsIdx.end(),
  vecCvdObsIdx.begin());
  }


  void GreedyRMA::set1DOptRange(const int& j) {
    setVecValueWeight(j);
    unsigned int curMinObjVal = runMinKadane(j);
    optObjVal = curMinObjVal;
    optLower = curLower;
    optUpper = curUpper;
    isPosIncumb = false;
    curMaxObjVal = runMaxKadane(j);
    if (-curMinObjVal <= curMaxObjVal) {
      optObjVal = curMaxObjVal;
      optLower = curLower;
      optUpper = curUpper;
      isPosIncumb = true;
    }
  }
  */

} // namespace greedyRMA
