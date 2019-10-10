/**********************************************************
* File name:   greedyRMA.cpp
* Author:      Ai Kagawa
* Description: a source file for greedy RMA solver
***********************************************************/

#include "greedyRMA.h"

namespace greedyRMA {

	void GreedyRMA::runGreedyRangeSearch() {

		ts.startTime();

		getMinOptRange();		// get oprimal range for minimum objective value

		getMaxOptRange();		// get oprimal range for maximum objective value

    chooseMinOrMaxRange();

		printSolution();

	} // end function greedyRMA


	void GreedyRMA::printSolution() {

#ifdef ACRO_HAVE_MPI
if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI

    ucout << "GRMA Solution: " << maxObjValue << "\t";
		ts.endCPUTime();
		DEBUGPRX(2, data, ts.endWallTime());
    //if (args->printBoost()) {
		DEBUGPRX(2, data, "L: " << L << "U: " << U);
		//}

#ifdef ACRO_HAVE_MPI
}
#endif //  ACRO_HAVE_MPI

	}


	/************************** Final Optimal Range **************************/
	void GreedyRMA::chooseMinOrMaxRange() {

		L.resize(data->numAttrib);
		U.resize(data->numAttrib);

		(args->randSeed()) ? srand((NumNegTiedSols+NumPosTiedSols)*time(NULL)*100) : srand(1);
		double rand_num = (rand() % 10001 ) / 10000.0 ;
		//DEBUGPRX(0, global(), "rand: " << rand_num  << "\n");
		//DEBUGPRX(0, global(), "rand1: " << 1.0 /  NumTiedSols << "\n");

		if ( maxVal>-minVal ||
			 ( maxVal==minVal && rand_num <= NumPosTiedSols/(double)(NumNegTiedSols+NumPosTiedSols) ) ) {

			maxObjValue=maxVal;
			isPosIncumb = true;

			copy(Lmax.begin(), Lmax.end(), L.begin());
			copy(Umax.begin(), Umax.end(), U.begin());

			DEBUGPRX(10, data, "chose max\n");

		} else {

			maxObjValue=-minVal;
			isPosIncumb = false;

			copy(Lmin.begin(), Lmin.end(), L.begin());
			copy(Umin.begin(), Umin.end(), U.begin());

			DEBUGPRX(10, data, "chose min\n");

		}

		DEBUGPRX(10, data, "Lmin: " << Lmin << "Umin: " << Umin);
		DEBUGPRX(10, data, "Lmax: " << Lmax << "Umax: " << Umax);

	}


	/********************* Optimal Range for Min Objective *********************/
	void GreedyRMA::getMinOptRange() {

		tmpMin=inf;
		minVal=inf;

		optAttrib=-1;
		prevAttrib=-1;

		Lmin.clear();
		Lmin.resize(data->numAttrib);

    Umin.resize(data->distFeat.size());
    copy(data->distFeat.begin(), data->distFeat.end(), Umin.begin());

    vecCoveredObs.resize(data->numTrainObs);
    copy(data->vecTrainData.begin(), data->vecTrainData.end(), vecCoveredObs.begin());

    do {

      fondNewBox = false;

      for (int j=0; j<data->numAttrib; ++j) { // for each feature

        if ( j != prevAttrib ) { // if this attribute is not restricted

					setObjVec(j);
          tmpMin = runMinKadane(j);

					if (tmpMin==minVal) {

						NumNegTiedSols++;
						(args->randSeed()) ? srand(NumNegTiedSols*time(NULL)*100) : srand(1);
						double rand_num = (rand() % 10001 ) / 10000.0 ;
						//DEBUGPRX(0, global(), "rand: " << rand_num  << "\n");
						//DEBUGPRX(0, global(), "rand1: " << 1.0 /  NumTiedSols << "\n");

						if ( rand_num  <= 1.0 /  NumNegTiedSols )
							setOptMin(j);

          } else if (tmpMin<minVal) {
						NumNegTiedSols=1;
						setOptMin(j);
					}

        } // end if this attribute is not restricted

      } // end each feature

      if ( fondNewBox ) {

				dropObsNotCovered(optAttrib, optLower, optUpper);

				prevAttrib = optAttrib;
				Lmin[optAttrib] = optLower;
        Umin[optAttrib] = optUpper;

				DEBUGPRX(20, data, "vecCoveredObs: " << vecCoveredObs);
        DEBUGPRX(10, data, "final optAttrib: (a,b): " << optAttrib << ": ("
					<< optLower << ", " << optUpper << "), min: " << minVal << "\n" );

			}

    } while( fondNewBox );

	}


	/********************* Optimal Range for Max Objective *********************/
	void GreedyRMA::getMaxOptRange() {

		tmpMax=-inf;
		maxVal=-inf;

		optAttrib=-1;
		prevAttrib=-1;

    Lmax.clear();
		Lmax.resize(data->numAttrib);

    Umax.resize(data->distFeat.size());
    copy(data->distFeat.begin(), data->distFeat.end(), Umax.begin());

    vecCoveredObs.resize(data->vecTrainData.size());
    copy(data->vecTrainData.begin(), data->vecTrainData.end(), vecCoveredObs.begin());

    do {

      fondNewBox = false;

      for (int j=0; j<data->numAttrib; ++j) { // for each feature

        if ( j != prevAttrib ) { // if this attribute is not restricted

					setObjVec(j);
          tmpMax = runMaxKadane(j);

					if (tmpMax==maxVal) {

						NumPosTiedSols++;
						(args->randSeed()) ? srand(NumPosTiedSols*time(NULL)*100) : srand(1);
						double rand_num = (rand() % 10001 ) / 10000.0 ;
						//DEBUGPRX(0, global(), "rand: " << rand_num  << "\n");
						//DEBUGPRX(0, global(), "rand1: " << 1.0 /  NumTiedSols << "\n");

						if ( rand_num  <= 1.0 /  NumPosTiedSols )
							setOptMax(j);

					} else if (tmpMax>maxVal) {
						NumPosTiedSols=1;
						setOptMax(j);
					}

        } // end if this attribute is not restricted
      } // end for each attribute

      if ( fondNewBox ) {	// if better solution found

        dropObsNotCovered(optAttrib, optLower, optUpper);

				prevAttrib = optAttrib;
        Lmax[optAttrib] = optLower;
        Umax[optAttrib] = optUpper;

        DEBUGPRX(10, data, "vecCoveredObs: " << vecCoveredObs);
        DEBUGPRX(10, data, "final optAttrib: (a,b): " << optAttrib << ": ("
          << optLower << ", " << optUpper << "), max: " << maxVal << "\n" );

      }

    } while( fondNewBox );

	}


	/********************* Optimal Range for 1D rules  *********************/

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
		copy(data->vecTrainData.begin(), data->vecTrainData.end(), vecCoveredObs.begin());
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


	void GreedyRMA::setOptMin(const int &j) {

		minVal = tmpMin;
		fondNewBox = true;

		optAttrib = j;
		optLower = tmpL;
		optUpper = tmpU;

		DEBUGPRX(10, data, "optAttrib: (a,b): " << optAttrib << ": ("
			<< optLower << ", " << optUpper << "), min: " << minVal << "\n") ;

	}


	void GreedyRMA::setOptMax(const int &j) {

		maxVal = tmpMax;
		fondNewBox = true;

		optAttrib = j;
		optLower = tmpL;
		optUpper = tmpU;

		DEBUGPRX(10, data, "optAttrib: (a,b): " << optAttrib << ": ("
			<< optLower << ", " << optUpper << "), max: " << maxVal << "\n" );

	}


	// get Miniumum range for the feature
	double GreedyRMA::runMinKadane(const int &j) {

		int s=Lmin[j];
		tmpL=Lmin[j];
		tmpU=Umin[j];

		double minEndHere=0;
		double minSoFar=inf;

		for (int i=Lmin[j]; i <= Umin[j]; ++i) {
			minEndHere += W[i] ; // getObjCovered(j, i);
			if (minEndHere < minSoFar) {
				minSoFar=minEndHere;
				tmpL=s;
				tmpU=i;
			}
			if ( minEndHere>0 ) {
				minEndHere=0;
				s=i+1;
			}
		}

		DEBUGPRX(10, data, "Minimum contiguous sum is " << minSoFar << " ");
		DEBUGPRX(10, data, " feat (L,U): " << j << " ("
														<< tmpL << ", " << tmpU << ")\n" );
		return minSoFar;

	}


	// get Maximum range for the feature
	double GreedyRMA::runMaxKadane(const int& j) {

		int s=Lmax[j];
		tmpL=Lmax[j];
		tmpU=Umax[j];

		double maxEndHere=0;
		double maxSoFar=-inf;    // min so far

		for (int i=Lmax[j]; i <= Umax[j]; ++i) {
			maxEndHere += W[i] ; // getObjCovered(j, i);
			if ( maxEndHere > maxSoFar ) {
				maxSoFar=maxEndHere;
				tmpL=s;
				tmpU=i;
			}
			if ( maxEndHere<0 ) {
				maxEndHere=0;
				s=i+1;
			}
		}

		DEBUGPRX(10, data, "Maximum contiguous sum is " << maxSoFar );
		DEBUGPRX(10, data, " feat (L,U): " << j << " ("
														<< tmpL << ", " << tmpU << ")\n" );
		return maxSoFar;

	}


  void GreedyRMA::dropObsNotCovered(const int &j, const int &lower, const int &upper) {

	  int obs, l=-1;
    DEBUGPRX(10, data, "Before drop: " << vecCoveredObs;);

		for (int i=0 ; i<vecCoveredObs.size(); ++i) {

		  obs = vecCoveredObs[i];

		  if ( lower <= data->intData[obs].X[j]
				  && data->intData[obs].X[j] <= upper)
				  //&& data->intData[obs].w!=0)
        vecCoveredObs[++l] = obs;   // store covered observations
    }

		vecCoveredObs.resize(l+1);  // shrink the size of vecCoveredObs
    DEBUGPRX(10, data, "After drop: " << vecCoveredObs;);

	}


	void GreedyRMA::setObjVec(const int &j) {

		int i, v, obs;

		W.resize(data->maxL);
		for (i=0; i<data->maxL; ++i)	W[i] = 0;

		for (i=0; i<vecCoveredObs.size(); ++i) {
			obs = vecCoveredObs[i];
			v = data->intData[obs].X[j];
			W[v] += data->intData[obs].w;
		}

		DEBUGPRX(10, data, "W: " << W );

	}

} // namespace greedyRMA


/*
	// bucket sort each attribute
	void GreedyRMA::bucketSortObs(const int& j) {

		int tmpNumObs;
		int bucketIndex;
		int obs, l=-1;

		int bucketSize = data->distFeat[j]+1;
		vector<vector<int> > buckets(bucketSize);

		// for each observation
		for (int i=0; i< data->vecTrainData.size(); ++i) {

      obs = vecCoveredObs[i];
			bucketIndex = data->intData[obs].X[j];

			if (bucketIndex<0)
				cerr<<"Error @ RMA: bucketIndex<0!";
			else if (bucketIndex>=bucketSize)
				cerr<<"Error @ RMA: bucketIndex>=bucketSize!";
			else
				buckets[bucketIndex].push_back(obs);
		}

		// walk buckets to get sorted observation list on this attribute
		for (int v=0; v<bucketSize; ++v)
			for (int k=0; k<buckets[v].size(); ++k)
				data->vecTrainData[++l] = buckets[v][k];

    DEBUGPRX(10, data, "after bucketSort: " << data->vecTrainData;);
	}
*/

	/*
  double GreedyRMAgetObjCovered(const int &j, const int &v) {
    double sumWeight=0.0; int obs, i;

    for (i=curObs; i<vecCoveredObs.size(); ++i) {  // for each observation
      obs = vecCoveredObs[i];
      // if each observation has cutValue in the feat
      if ( data->intData[obs].X[j] == v )
				sumWeight += data->intData[obs].w;
      else { curObs = i; break; }
    } // end for each observation

    return sumWeight;
  }
*/
