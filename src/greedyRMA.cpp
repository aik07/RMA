/**********************************************************
 * Author:      Ai Kagawa
 * Description: a source file for greedy RMA solver
 ***********************************************************/


#include "greedyRMA.h"


namespace greedyRMA {


  void GreedyRMA::reset() {

    bestObjVal = -getInf();

    vecLowerWorking.resize(data->numAttrib);
    vecUpperWorking.resize(data->numAttrib);
    vecLower.resize(data->numAttrib);
    vecUpper.resize(data->numAttrib);

    vecValueWeight.resize(data->maxNumDistVals);

    for (unsigned int k=0; k<data->maxNumDistVals; ++k)
      vecValueWeight[k] = 0;

  } // end reset function


  void GreedyRMA::runGreedyRangeSearch() {

    ts.startTime();         // start the timer

    reset();                // reset variables in the Greedy RMA class

    if (args->isMinMax()) {

      isPosObjVal = false;

      searchGreedyRange(MIN);    // search oprimal range for minimum objective value

      updateBestBounds();

      searchGreedyRange(MAX);    // search oprimal range for maximum objective value

      if (isImprovedSol) {
        isPosObjVal = true;
        updateBestBounds();
      }

    } else {

      isPosObjVal = true;

      searchGreedyRange(MAX);    // search oprimal range for maximum objective value

      updateBestBounds();

      searchGreedyRange(MIN);    // search oprimal range for minimum objective value

      if (isImprovedSol) {
        isPosObjVal = false;
        updateBestBounds();
      }

    } // end if greedy range search

    printSolution();        // print out the solution

  } // end function greedyRMA


  /********************* Greedy Range for Objective *********************/
  void GreedyRMA::searchGreedyRange(const bool &isMax) {

    int  curIter = 0;

    // whether or not the algorithm found an additional restricting range
    // which can improve the objective value in each greedy iteration
    bool isFondNewBox;

    isImprovedSol=false;

    resetBestRange();

    if (args->debug >= 1)
      cout << "************* " << ( isMax ? "MAX" : "MIN") << " Range Search"
           << " *************\n";

    do {  // repeat while finding a new solution

      if (args->debug >= 1)
        cout << "************* current iteration: " << curIter
             << " *************\n";

      isFondNewBox = false;

      for (unsigned int j = 0; j < data->numAttrib; ++j) { // for each feature

        if (j != prevBestAttrib) { // if this attribute is not previously chosen

          // set the vector of total weights for each value of attribute j
          setVecValueWeight(j);

          // run objective value for attribut j
          runKadaneAlgo(isMax, j);

          // if the current min objective value is equal to
          // the min objective value found so far
          if (curObjVal == bestObjVal) {

            ++numTiedSols;  // count the tied solutions

            if (isUpdateBestSol(numTiedSols)) {
              isFondNewBox = true;
              setBestRange(j); // update the current best solution
            }

          // if the current objective value is greater than
          // the best objective value found so far
          } else if (curObjVal > bestObjVal) {

            numTiedSols  = 1 ; // only one tied solution so far
            isFondNewBox = true;
            setBestRange(j); // update the current optimal solutiona for min ver

          } // end if to possibly insert a new restriction

        } // end if this attribute is not restricted in the previous iteration

      } // end each feature

      if (isFondNewBox) { // if a new solution is discovered

        isImprovedSol = isFondNewBox;

        // update the lower and upper bounds for the optimal attribute
        vecLowerWorking[bestAttrib] = bestLower;
        vecUpperWorking[bestAttrib] = bestUpper;

        // drop observations which are not covered
        // by the optimal attribute, and its optimal lower and upper bounds
        dropObsNotCovered(bestAttrib, bestLower, bestUpper);

        // set the previous chosen attribute
        prevBestAttrib = bestAttrib;

        if (args->debug >= 1) {
          cout << "New restrcition: optAttrib: "    << bestAttrib
                << ";  (a,b): ("      << bestLower  << ", " << bestUpper
                << "); bestObjVal: "  << bestObjVal << "\n";
          cout << "vecLowerWorking: " << vecLowerWorking  << "\n";
          cout << "vecUpperWorking: " << vecUpperWorking  << "\n";
        }  // end if debug

      } // end if a new solution is discovered

      ++curIter;

    } while (isFondNewBox);

  } // end searchMinOptRange function


  // reset variables for the best range search
  // reset the lower and upper bounds
  // reset the covered observation list, vecCvdObsIdx
  void GreedyRMA::resetBestRange() {

    // curObjVal  = -getInf();

    bestAttrib     = -1;
    prevBestAttrib = -1;

    // set lower vector for the maximum version
    // vecLowerMax.clear();
    // vecLowerMax.resize(data->numAttrib);
    fill(vecLowerWorking.begin(), vecLowerWorking.end(), 0);

    // set upper vector for the maximum version
    for (unsigned int j=0; j<data->numAttrib; ++j)
      vecUpperWorking[j] = data->vecNumDistVals[j] - 1;

    vecCvdObsIdx.resize(data->numNonZeroWtObs);
    copy(data->vecNonZeroWtObsIdx.begin(), data->vecNonZeroWtObsIdx.end(),
         vecCvdObsIdx.begin());

  } // end resetMaxOptRange function


  // update the optimal attribute (j), objective value, and lower and upper bounds
  // for the maximization version
  void GreedyRMA::setBestRange(const unsigned int &j) {

    bestAttrib    = j;
    bestObjVal    = curObjVal;
    bestLower     = curLower;
    bestUpper     = curUpper;

    if (args->debug >= 2)
      ucout << "Best Attrib: "  << bestAttrib
            << "; (a,b): ("     << bestLower << ", " << bestUpper << ")"
            << "; bestObjVal: " << bestObjVal << "\n";

  } // end setBestRange function


  // set the current objective value, lower and upper bounds
  // by running the max ver. of Kadane's algorithm for attribute j
  void GreedyRMA::runKadaneAlgo(const bool &isMax, const unsigned int &j) {

    int s    = vecLowerWorking[j];
    curLower = vecLowerWorking[j];
    curUpper = vecUpperWorking[j];

    double maxEndHere = 0;
    curObjVal         = -getInf(); // min so far

    for (unsigned int i = vecLowerWorking[j]; i <= vecUpperWorking[j]; ++i) {
      maxEndHere += ( isMax ? vecValueWeight[i] : -vecValueWeight[i] ) ;
      if (maxEndHere > curObjVal) {
        curObjVal = maxEndHere;
        curLower = s;
        curUpper = i;
      }
      if (maxEndHere < 0) {
        maxEndHere = 0;
        s = i + 1;
      }
    }

    if (args->debug >= 100)
      ucout << "Maximum contiguous sum is "    << curObjVal
            << "; feat: " << j
            << "; (L,U): (" << curLower << ", " << curUpper << ")\n";

  } // end runKadaneAlgo function


  // whether or not to update the best solution
  // when there are tied solution
  // (break the tie using a fair probability based on # of tied solutions)
  bool  GreedyRMA::isUpdateBestSol(const int &numTiedSols) {

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


  // drop the observations which are not covered
  // by the attribute j's lower and upper bound
  void GreedyRMA::dropObsNotCovered(const unsigned int &j,
                                    const unsigned int &lower,
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

    // for each covered observation index
    for (i = 0; i < vecCvdObsIdx.size(); ++i) {

      obs = vecCvdObsIdx[i];              // observation index

      // X value of this observation for attribute j
      v   = data->dataIntTrain[obs].X[j];

      // add weights of observations for each X-value of attribute j
      vecValueWeight[v] += data->dataIntTrain[obs].w;

    }

    if (args->debug >= 10)
      ucout << "vecValueWeight: " << vecValueWeight;

  } // end setVecValueWeight function


  // print the greedy RMA solution
  void GreedyRMA::printSolution() {

  #ifdef ACRO_HAVE_MPI
    if (uMPI::rank == 0) {
  #endif //  ACRO_HAVE_MPI

      std::cout << "GRMA Solution: ";
      isPosObjVal ? cout << "+" : ucout << "-";
      std::cout << std::fixed << std::setprecision(4)
                << bestObjVal << "\t";
      std::cout << std::fixed << std::setprecision(2)
                << "CPU Time: " << ts.getCPUTime();
      if (args->debug >= 2)
        std::cout << "Wall Time: " << ts.getWallTime() << "\n";
      std::cout << "\n";
      // if (args->printBoost()) {
      if (args->debug >= 2)
        std::cout << "Greedy vecLower: "   << vecLower
                  << "\nGreedy vecUpper: " << vecUpper << "\n";
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

    vecCvdObsIdx.resize(data->numNonZeroWtObs);
    copy(data->vecNonZeroWtObsIdx.begin(), data->vecNonZeroWtObsIdx.end(),
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
