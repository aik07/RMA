// rma.cpp
//
//  Implements larger methods in example of how to use object-oriented
//  branching framework (for RMA problems).
//
// Ai Kagawa
//

#include <acro_config.h>
#ifdef ACRO_HAVE_MPI

#include <stdlib.h>
#include <iostream>
#include <map>
#include <vector>
#include <stack>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>    // std::in
#include <mpi.h>
#include <utility>
#include <utilib/logEvent.h>
#include <utilib/_math.h>
#include <utilib/stl_auxiliary.h>
#include <utilib/exception_mngr.h>
#include <utilib/comments.h>
#include <utilib/mpiUtil.h>
#include <utilib/std_headers.h>
#include <utilib/PackBuf.h>
#include <utilib/BitArray.h>
#include <pebbl/fundamentals.h>
#include <pebbl/packedSolution.h>
#include <pebbl/parPebblBase.h>
#include <pebbl/ThreadObj.h>
#include <pebbl/SelfAdjustThread.h>
#include <pebbl/coTree.h>
#include <pebbl/outBufferQ.h>
#include <pebbl/parBranching.h>
#include "parRMA.h"
#include "serRMA.h"


using namespace utilib;
using namespace std;
using namespace pebbl;


namespace pebblRMA {

	////////////////////// cutPtRecThd methods (Begining) //////////////////////////////////

	// constructor of cutPtThd class
	CutPtThd::CutPtThd( parRMA* global_,  MessageID msgID) :
			broadcastPBThread(global_,
							 "Cut Point Receiver",
							 "cutPtRec",
							 "PaleTurquoise3",
							 5,                // logLevel
							 35,               // debug level
							 3*sizeof(int),
							 msgID,
							 3),               // Tree radix -- increase later
			j(-1),
			v(-1),
			ptrParRMA(global_)
			{ };


	// Run method.  This is only invoked if a message arrives.
	bool CutPtThd::unloadBuffer() {

		inBuf >> j >> v >> originator;
		DEBUGPR(10, ucout << "cutPtThd message received from "
						<< status.MPI_SOURCE
						<< "(j, v)=("	<< j << ", " << v << ")"
						<< ", originator=" << originator << '\n');

		bool seenAlready = false;
    multimap<int,int>::iterator it,itlow,itup;

    itlow = ptrParRMA->mmapCachedCutPts.lower_bound(j);  // itlow points to
    itup = ptrParRMA->mmapCachedCutPts.upper_bound(j);   // itup points to

    // print range [itlow,itup):
    for (it=itlow; it!=itup; ++it) {
      if ( (*it).first==j && (*it).second==v ) seenAlready=true;
      DEBUGPR(10, ucout << (*it).first << " => " << (*it).second << '\n');
    }

    DEBUGPR(10, ucout << "cut point (" << j << ", " << v << ") " );
    DEBUGPR(10, ucout << (seenAlready ? "is already in cache\n" : "is new\n"));

		if (originator < 0) {
	  	if (seenAlready)
	  		return false;
	  	originator = uMPI::rank;
	  }

    // if not in the hash table, insert the cut point into the hash table.
	  if (!seenAlready)
			ptrParRMA->mmapCachedCutPts.insert( make_pair<int, int>(j, v) ) ;

		/*
		CutPt tempCutPt = {j, v};
		//tempCutPt.j = j;
		//tempCutPt.v = v;

		multimap<int, int> iterator it = ptrParRMA->mmapCachedCutPts.find(j);
		bool seenAlready =  (it!=ptrParRMA->mmapCachedCutPts.end());

	  DEBUGPR(25, ucout << (seenAlready ? "Seen already\n" : "Looks new\n"));
	  if (originator < 0) {
	  	if (seenAlready)
	  		return false;
	  	originator = uMPI::rank;
	  }
	  if (!seenAlready)
	  	ptrParRMA->mmapCachedCutPts[j] = v;
			//ptrParRMA->mapCachedCutPts.add(hash_bj(tempVecCutPts)%ptrParRMA->totalCutPts, tempCutPt);
		*/
	  return true;
	}


	void CutPtThd::setCutPtThd(const int& _j, const int& _v) {
		j = _j;		v = _v;
	}


	// Logic to relay information to other nodes.
	void CutPtThd::relayLoadBuffer(PackBuffer* buf) {
		*buf << j << v << originator;
		DEBUGPR(20,ucout << "cutPtThd writing (feat, cutVal)=(" << j << ", " << v
					     << "), originator=" << originator << "\n");
	}


	// Special method to send initial message to owning processor
	void CutPtThd::preBroadcastMessage(const int& owningProc) {
		DEBUGPR(25, ucout << "CutPtThd root send to " << owningProc << '\n');
		// A negative value for 'originator' indicates special root message
		originator = -1;
		// Grab a buffer from the same pool used for broadcasts
		PackBuffer* outBuf = outQueue.getFree();
		// Fill it in the usual way
		relayLoadBuffer(outBuf);
		// Send it to the owning processor
		outQueue.send(outBuf,owningProc,tag);
		// Make sure message counters are updated correctly
		global->recordMessageSent(this);
	}


	///////////////////////////////////// parRMA methods //////////////////////////////////////

	parRMA::parRMA() : RMA(), cutPtCaster(NULL) {

		// Default is not to spend time on a dumb ramp up
		rampUpPoolLimitFac = 1.0;
		Parameter& p = get_parameter_object("rampUpPoolLimitFac");
		p.default_value = "1.0";

		rampUpFeatureFac = 1.0;
		create_categorized_parameter("rampUpFeatureFac",
					   rampUpFeatureFac,
					   "<double>",
					   "1.0",
					   "Maximum number of subproblems "
					   "in pool to end ramp-up phase,\n\t"
					   "as a fraction of the total number "
					   "of features.",
					   "Maximum Monomial",
					   utilib::ParameterNonnegative<double>());

		branchChoice::setupMPI();
	}

	// Destructor.
	parRMA::~parRMA() {
		if (cutPtCaster == 0)
			delete cutPtCaster;
		branchChoice::freeMPI();
	}

	/// Note: use VB flag?
	void parRMA::reset(const bool& VBflag) {
		RMA::reset();
		registerFirstSolution(new rmaSolution(this));
		if (cutPtCaster) {
			delete cutPtCaster;
			cutPtCaster = NULL;
		}
		parallelBranching::reset(VBflag);
	}


	void parRMA::placeTasks() {
		parallelBranching::placeTasks();
		cutPtCaster = new CutPtThd(this,cutPtBroadcastTag);
		placeTask(cutPtCaster,true,highPriorityGroup);
	}


	parallelBranchSub* parRMA::blankParallelSub() {
		parRMASub *newSP = new parRMASub();
		newSP->setGlobalInfo(this);
		//newSP->setGlobalInfo(this);
		return newSP;
	};


	// Pack a description of the problem.
	void parRMA::pack(PackBuffer& outBuf) {

	  DEBUGPR(20,ucout << "parRMA::pack invoked..." << '\n');
	  outBuf << numObs << numDistObs << numAttrib ;

		for (int i=0; i<numDistObs; ++i)  outBuf << sortedObsIdx[i];

	  for (int i=0; i<numObs; ++i) {
	   	outBuf << intData[i].X << intData[i].w << intData[i].y;
	  }
	  outBuf << distFeat << _branchSelection
		       << _perLimitAttrib << _perCachedCutPts << numTotalCutPts;

	} // end function parRMA::pack


	// unpack
	void parRMA::unpack(UnPackBuffer& inBuf) {

		DEBUGPR(20,ucout << "parRMA::unpack invoked... " << '\n');

		inBuf >> numObs >> numDistObs >> numAttrib ;

		sortedObsIdx.resize(numDistObs);
		for (int i=0; i<numDistObs; ++i)  inBuf >> sortedObsIdx[i];

		intData.resize(numObs);
		for (int i=0; i<numObs; ++i) {
			intData[i].X.resize(numAttrib);
			inBuf >> intData[i].X >> intData[i].w >> intData[i].y;
		}
		inBuf >> distFeat >> _branchSelection
		      >> _perLimitAttrib >> _perCachedCutPts >> numTotalCutPts;

		DEBUGPR(20,ucout << "parRMA::unpack done." << '\n');

		/*
		DEBUGPR(20,ucout <<" distFeat: ";
				for(size_type i=0; i<numAttrib; i++) {
				ucout << distFeat[i] << ", ";
		});

		DEBUGPR(20,for(size_type i=0; i<numDistObs; i++) {
				ucout <<" wt: " << intData[i].w << '\n';
		});
		*/
	} // end function parRMA::unpack


	int parRMA::spPackSize() {
		return 5*(numAttrib+2)*sizeof(int)	//  al << au << bl << bu << subState
				+  2*sizeof(int)+3*(sizeof(double)+ sizeof(int)) // size of branchChoice
				+  (numAttrib+2)*sizeof(bool);		// vecCheckedFeat
	} // end function parRMA::spPackSize


	// using virtual function
	void parRMA::setCachedCutPts(const int& j, const int& v) {

    bool isAlreadyInCache = false;
    multimap<int,int>::iterator it,itlow,itup;

    itlow = mmapCachedCutPts.lower_bound(j);  // itlow points to
    itup = mmapCachedCutPts.upper_bound(j);   // itup points to

    // print range [itlow,itup):
    for (it=itlow; it!=itup; ++it) {
      if ( (*it).first==j && (*it).second==v ) isAlreadyInCache=true;
      DEBUGPR(10, ucout << (*it).first << " => " << (*it).second << '\n');
    }

    DEBUGPR(10, ucout << "cut point (" << j << ", " << v << ") " );
    DEBUGPR(10, ucout << (isAlreadyInCache ? "is already in cache\n" : "is new\n"));

    // if not in the hash table, insert the cut point into the hash table.
	  if (!isAlreadyInCache)
			mmapCachedCutPts.insert( make_pair<int, int>(j, v) ) ;

		int owningProc = mmapCachedCutPts.find(j)->second % uMPI::size;
		DEBUGPR(20, ucout << "owningProc: " << owningProc << '\n');

		cutPtCaster->setCutPtThd(j, v);

		if (uMPI::rank==owningProc) {
			// This processor is the owning processor
			DEBUGPR(20, ucout << "I am the owner\n");
			cutPtCaster->initiateBroadcast();
		} else {
			DEBUGPR(20, ucout << "Not owner\n");
			cutPtCaster->preBroadcastMessage(owningProc);
		}

	/*
		//int key = rand() % 100000;

		//CutPt tempCutPt = {j,v};
		//tempCutPt.j = j;
		//tempCutPt.v = v;

		//if this cut point is already found before
		//if ( mmapCachedCutPts.count(tempCutPt) > 0 )
		//	return;

		// Store stuff anyway
		mmapCachedCutPts[v]=j;

		int owningProc = mmapCachedCutPts.find(j)->second % uMPI::size;
		DEBUGPR(20, ucout << "owningProc: " << owningProc << '\n');

		cutPtCaster->setCutPtThd(j, v);

		if (uMPI::rank==owningProc) {
			// This processor is the owning processor
			DEBUGPR(20, ucout << "I am the owner\n");
			cutPtCaster->initiateBroadcast();
		} else {
			DEBUGPR(20, ucout << "Not owner\n");
			cutPtCaster->preBroadcastMessage(owningProc);
		}
*/
return;
	} // end function parRMA::setCachedCutPts


	//////////////////////// parRMASub methods /////////////////////////////////

	void parRMASub::pack(utilib::PackBuffer &outBuffer) {

		DEBUGPRXP(20,pGlobal(), "parRMASub::pack invoked...\n");

		outBuffer << al << au << bl << bu ;
		outBuffer << _branchChoice.branchVar << _branchChoice.cutVal;

		for (int i=0; i<3; ++i) {
			outBuffer << _branchChoice.branch[i].roundedBound
					  << _branchChoice.branch[i].whichChild;
			//DEBUGPRX(20,　pGlobal(),"_branchChoice.branch[i].roundedBound. :"
				//<< _branchChoice.branch[i].roundedBound << "\n");
		}

		//outBuffer << vecCheckedFeat;

	  DEBUGPRXP(20, pGlobal(), "parRMASub::pack done. " << " bound: " << bound << "\n");
	} // end function parRMASub::pack


	void parRMASub::unpack(utilib::UnPackBuffer &inBuffer) {

  	//DEBUGPRX(20, pGlobal(),"parRMASub::unpack invoked...\n");

    inBuffer >> al >> au >> bl >> bu;
    inBuffer >> _branchChoice.branchVar >> _branchChoice.cutVal;

    for (int i=0; i<3; ++i) {
    	inBuffer >> _branchChoice.branch[i].roundedBound
    					 >> _branchChoice.branch[i].whichChild;
    	//DEBUGPRX(20,pGlobal(),"branchChoice.branch[i].roundedBound. :"
				//<< _branchChoice.branch[i].roundedBound << "\n");
    }

		//inBuffer >> vecCheckedFeat;

    DEBUGPRX(20,pGlobal(),"parRMASub::unpack done. :" << " bound: " << bound << '\n');

	} // end function parRMASub::unpack


	// makeParallelChild
	parallelBranchSub * parRMASub::makeParallelChild(int whichChild) {

	  DEBUGPRX(20, global(), "parRMASub::makeParallelChild invoked for: "
	    						<< ", whichChild: " << whichChild
	    						<< ", ramp-up flag: " << rampingUp() << '\n');

		if (whichChild==-1) {
			cerr << "which child cannot be -1!";
			return NULL;
		}

#ifdef ACRO_VALIDATING
	    if (whichChild < 0 || whichChild > 2) {
	    	ucout << "parRMASub::makeParallelChild: invalid request "
				  << "for child " << whichChild << '\n';
	    	return NULL;
	    }

	    if ((_branchChoice.branchVar < 0) || (_branchChoice.branchVar >= global()->numAttrib)) {
	    	ucout << "parRMASub::makeParallelChild: invalid branching variable\n";
	    	return NULL;
	    }
#endif

	    // If there are no cached children (because this subproblem was
	    // sent from somewhere else), recreate a child, not necessarily in
	    // bound-sorted order.  Otherwise, grab it from the cache.

		if (_branchChoice.branchVar>numAttrib()) {
			DEBUGPR(20, ucout <<  "ERROR in parallel! "
					 << "_branchChoice.branchVar: " << _branchChoice.branchVar << '\n');
			cerr <<  " ERROR in parallel! "
								 << "_branchChoice.branchVar: " << _branchChoice.branchVar << '\n';
			//return NULL;
		}

		if (whichChild < 0)
			cerr << "whichChild=" << whichChild << '\n';

		DEBUGPR(20, ucout << "whichChild=" << whichChild << '\n');

		parRMASub* temp = new parRMASub();
		temp->setGlobalInfo(globalPtr);

		DEBUGPR(20, ucout << "_branchChoice.branch[whichChild].whichChild="
				<< _branchChoice.branch[whichChild].whichChild << '\n');
		temp->RMASubAsChildOf(this, whichChild);

		DEBUGPR(10,ucout << "Parallel MakeChild produces " << temp << '\n');

		DEBUGPRX(10,global(),"Out of parRMASub::makeParallelChild, "
							"whichChild: " << whichChild << " bound: " << bound << '\n');

		return temp;
	} // end function parRMASub::makeParallelChild


	// split subproblems
	void parRMASub::parStrongBranching(const int& firstIdx, const int& lastIdx) {

		int size = uMPI::size;
		size_type rank = uMPI::rank;
		int quotient  = numAttrib() / size;
		int remainder = numAttrib() % size;
		int firstAttrib = rank*quotient + min((int)rank,remainder);
		int lastAttrib= firstAttrib + quotient + (rank < remainder) - 1;

    int numCutPtsInAttrib, countCutPts=0;

    for (int j=0; j<numAttrib(); ++j ) { // for each attribute

			numCutPtsInAttrib = bu[j] - al[j] - max(0, bl[j]-au[j]);

			if ( countCutPts + numCutPtsInAttrib < firstIdx ) {
				countCutPts += numCutPtsInAttrib;
				if (!global()->bruteForceEC()) bucketSortEC(j);
				if ( firstAttrib<=j && j<=lastAttrib )
					checkIncumbent(j);
				continue;
			}

      DEBUGPR(10,ucout << "original: ");
      printSP(j, al[j], au[j], bl[j], bu[j]);

      for (int v=al[j]; v<bu[j]; ++v ) { // for each cut-value in attribute j

				// if a cut-value is in [au, bl-1]
				if ( au[j]<=v && v<bl[j] ) { v=bl[j]-1; continue; }

				// if this cut point is not assinged in this processor
				if (countCutPts<firstIdx) { ++countCutPts; continue; }
				if (countCutPts>lastIdx) break;

        branchingProcess(j, v);
				++countCutPts;

      } // end for each cut-value in attribute j

      if (!global()->bruteForceEC()) bucketSortEC(j);
			if ( firstAttrib<=j && j<=lastAttrib )
				checkIncumbent(j);

			if ( countCutPts > lastIdx ) break;

    } // end for each attribute

	} // end function parRMASub::splitSP


	void parRMASub::setLiveCachedCutPts() {

    // numLiveCachedCutPts = (# of live cut points from the cache)
    int numLiveCachedCutPts=getNumLiveCachedCutPts();

    // if numCachedCutPts is less than the percentage, check all cut points
    if ( numLiveCachedCutPts
          < globalPtr->numTotalCutPts * globalPtr->perCachedCutPts() )
      return;

    // if not, only check the storedCutPts

    // count number of subproblems only discovering cut-points from the chache
		isCachedCutPts = true;
    ++global()->numCC_SP;
    int j, v, l=-1;
    multimap<int, int>::iterator curr = global()->mmapCachedCutPts.begin();
    multimap<int, int>::iterator end = global()->mmapCachedCutPts.end();
    cachedCutPts.resize(numLiveCachedCutPts);
    while (curr!=end) {
    	j = curr->first;
      v = curr->second;
      //if (j>numAttrib() || v<0) error;
      curr++;
      if (al[j]<=v && v<bu[j]) // if v in [al, bu)
        if ( !( au[j]<bl[j] && au[j]<=v && v<bl[j] ) ) { // if not overlapping
          cachedCutPts[++l].j = j;
          cachedCutPts[l].v = v;
        }
		}

		int size = uMPI::size;
		size_type rank = uMPI::rank;

		if ( (size>cachedCutPts.size() && rank>=cachedCutPts.size())
				|| cachedCutPts.size()==0 ){
			DEBUGPRX(20, global(), "ramp-up, this processor won't compute." << '\n');
			return;
		}

		int quotient  = cachedCutPts.size() / size;
		int remainder = cachedCutPts.size() % size;

		int firstIndex = rank*quotient + min((int)rank,remainder);
		int lastIndex = firstIndex + quotient + (rank < remainder) - 1;

    parCachedBranching(firstIndex, lastIndex);

  }


	// branching using cut-point caching methods
	void parRMASub::parCachedBranching(int firstIdx, int lastIdx) {

		int size = uMPI::size;
		size_type rank = uMPI::rank;
		int quotient  = numAttrib() / size;
		int remainder = numAttrib() % size;
		int firstAttrib = rank*quotient + min((int)rank,remainder);
		int lastAttrib= firstAttrib + quotient + (rank < remainder) - 1;

		//int firstJ = cachedCutPts[firstIdx].j;
		//int lastJ = cachedCutPts[lastIdx].j;

    int k=0;

    sortCachedCutPtByAttrib();
		cachedCutPts = sortedCachedCutPts;

    for (int j=0; j<numAttrib(); ++j) {
			while ( k<cachedCutPts.size() ) {
				if ( j == cachedCutPts[k].j) {
					if (firstIdx<=k && k<=lastIdx)
						branchingProcess(cachedCutPts[k].j, cachedCutPts[k].v);
					++k;
				} else break;
			}
			//if (j==numAttrib()-1) break;
      if (!global()->bruteForceEC()) bucketSortEC(j);
			if ( firstAttrib<=j && j<=lastAttrib ) checkIncumbent(j);
    }

	} // end RMASub::cachedBranching


  // Bound computation -- unless we're in ramp-up, just do the same
  // thing as the serial three-way code.  If we're in ramp-up, try to
  // parallelize the strong branching procedure
	void parRMASub::boundComputation(double* controlParam)  {

		DEBUGPRX(20, global(),
			"In parRMASub::boundComputation, ramp-up flag=" << rampingUp() << '\n');

		if (!rampingUp()) {
			RMASub::boundComputation(controlParam);
			return;
		}

		NumTiedSols=1;

		// Special handling of ramp-up
		DEBUGPRX(20, global(), "Ramp-up bound computation\n");

		//************************************************************************
		coveredObs = global()->sortedObsIdx;
		// sort each feature based on [au, bl]
    for (int j=0; j<numAttrib(); j++)  bucketSortObs(j);
		setInitialEquivClass();	// set vecEquivClass

		//**************************************************************
		// Figure which variables go on which processor.  Make them as even as possible
		// -- the first (remainder) processors have one more variable

		isCachedCutPts = false;
		int size = uMPI::size;
		size_type rank = uMPI::rank;

		// if there are enough discoverd cut points (storedCutPts) check only the list
    /*if ( global()->perCachedCutPts() < 1.0 ) {
			setLiveCachedCutPts();
    }*/

		if (!isCachedCutPts) {  // check all cut points

			setNumLiveCutPts();

			int quotient  = numLiveCutPts / size;
			int remainder = numLiveCutPts % size;

			int firstIndex = rank*quotient + min((int)rank,remainder);
			int lastIndex = firstIndex + quotient + (rank < remainder) - 1;

			DEBUGPRX(10,global(), "numLiveCutPts = " << numLiveCutPts
							<< ", quotient  = " << quotient
							<< ", remainder = " << remainder
							<< ", firstIndex = " << firstIndex
							<< ", lastIndex = " << lastIndex << '\n');

      parStrongBranching(firstIndex, lastIndex);

    }

		printCurrentBounds();

		DEBUGPRX(20, global(), "Best local choice is " <<  _branchChoice <<
                           " NumTiedSols: " << NumTiedSols << '\n');

		/******************* rampUpIncumbentSync *******************/

		DEBUGPRX(1, global(), rank << ": BEFORE rampUpIncumbentSync():"
													 << pGlobal()->rampUpMessages << '\n');

		// Better incumbents may have been found along the way
		pGlobal()->rampUpIncumbentSync();

		DEBUGPRX(1, global(), rank << ": AFTER rampUpIncumbentSync():"
													<< pGlobal()->rampUpMessages << '\n');

		/******************* Global Choice *******************/
		// Now determine the globally best branching choice by global reduction.
		// Use the special MPI type and combiner for branch choices.

		branchChoice bestBranch;

		DEBUGPRX(1, global(), rank << ": before reduceCast: "
										<< pGlobal()->rampUpMessages << '\n');

		if (global()->branchSelection()==0) {
			MPI_Scan(&_branchChoice, &bestBranch, 1,
							 branchChoice::mpiType, branchChoice::mpiBranchSelection,
							 MPI_COMM_WORLD);

		 	MPI_Bcast(&bestBranch, 1, branchChoice::mpiType,
			          size-1, MPI_COMM_WORLD);
		} else {
			uMPI::reduceCast(&_branchChoice, &bestBranch, 1,
		      branchChoice::mpiType, branchChoice::mpiCombiner);
		}

		pGlobal()->rampUpMessages += 2*(uMPI::rank > 0);

		DEBUGPRX(1, global(), rank << ": after reduceCast:"
										<< pGlobal()->rampUpMessages << '\n');

		DEBUGPRX(10, global(), "Best global choice is " << bestBranch << '\n');

		/******************* Cache cut-point *******************/
		if (global()->perCachedCutPts()<1.0)
			globalPtr->setCachedCutPts(bestBranch.branchVar, bestBranch.cutVal);

		/************************************************************/

		// If this processor has the best choice,  there is nothing to do.
		// Otherwise, adjust everything so it looks like we made the globally best choice.
		if (bestBranch.branchVar != _branchChoice.branchVar ||
				bestBranch.cutVal != _branchChoice.cutVal) {
			DEBUGPRX(20, global(),  "Adjusting local choice\n");
			_branchChoice = bestBranch;
		}

		bound = _branchChoice.branch[0].roundedBound;	// look ahead bound
		setState(bounded);

		// If objValue >= bound, then we found a solution.
		if ( workingSol()->value >= _branchChoice.branch[0].roundedBound ) {
		        //workingSol()->printSolution();
			foundSolution();
			setState(dead);
		}

		// If roundedBound=-1, then we found a solution.
		if (_branchChoice.branch[0].roundedBound==-1) {
		        //workingSol()->printSolution();
			foundSolution();
			setState(dead);
		}

		// if the stored cut point is
		if (globalPtr->mmapCachedCutPts.size() >= size * global()->rampUpSizeFact())
			pGlobal()->rampUpPoolLimitFac = 0;

		//DEBUGPRX(50, global(), " bound: " << bound << ", sol val=" << getObjectiveVal() << '\n');
		DEBUGPRX(10, global(), "Ending ramp-up bound computation for bound: " << bound << '\n');

	} // end function parRMASub::boundComputation


  void parRMASub::setNumLiveCutPts()  {
    numLiveCutPts=0;
    // compute the total cut points
    for (int j=0; j<numAttrib(); ++j) {
      // calculate total numbers of cut points
      if ( global()->perLimitAttrib()==1 ) {
        numLiveCutPts += bu[j]-al[j] -  max(0, bl[j]-au[j]);
      }
    }
    if (numLiveCutPts==0) {
      DEBUGPR(20, cout << "No cut points to check!\n");
      setState(dead);
    }
  }

} //************************************ namespace pebbl (end) ************************************

#endif // ACRO_HAVE_MPI
