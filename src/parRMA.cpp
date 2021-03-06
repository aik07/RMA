/**********************************************************
 *  Author:      Ai Kagawa
 *  Description: a source file for the parallel RMA solver using PEBBL
 **********************************************************/


#include "parRMA.h"
#include "serRMA.h"
#include <pebbl/utilib/hash_fn.h>

//#ifdef ACRO_HAVE_MPI

namespace pebblRMA {

////////////// cutPtRecThd methods (Begining) ///////////////////////////

  // constructor of cutPtThd class
  CutPtThd::CutPtThd(parRMA *global_, MessageID msgID)
      : broadcastPBThread(global_, "Cut Point Receiver", "cutPtRec",
                          "PaleTurquoise3",
                          5,  // logLevel
                          35, // debug level
                          3 * sizeof(int), msgID,
                          3), // Tree radix -- increase later
        j(-1), v(-1), ptrParRMA(global_){};


  // Method to read the input buffer and decide what action to take
  bool CutPtThd::unloadBuffer() {

    inBuf >> j >> v >> originator;
    if (ptrParRMA->args->debug >= 10)
      ucout << "cutPtThd message received from " << status.MPI_SOURCE
            << " (j, v)=(" << j << ", " << v << ")"
            << "(j, v)=(" << j << ", " << v << ")"
            << ", originator=" << originator << '\n';

    // See whether this cutpoint is already in the cache (and put it in if not)

    bool seenAlready = ptrParRMA->putInCache(j,v);

    // If this the original send from cutpoint discovery to the owning processor
    // and the cutpoint has been seen already, we want to ignore it.  In this cawse
    // return false to indicate no further action.  Otherwise, we will initiate a 
    // broacast from this processor, so change the originator processor to this one.

    if (originator < 0) {
      if (seenAlready)
        return false;
      originator = uMPI::rank;
    }

    // Return true to indicate that further action may be necessary, either starting 
    // a broadcast or possibly relaying it.

    return true;
  }


  void CutPtThd::setCutPtThd(const unsigned int &_j, const unsigned int &_v) {
    j = _j;
    v = _v;
  }


  // Logic to relay information to other nodes.
  void CutPtThd::relayLoadBuffer(PackBuffer *buf) {
    *buf << j << v << originator;
    if (ptrParRMA->args->debug >= 20)
      ucout << "cutPtThd writing (feat, cutVal)=(" << j << ", " << v
            << "), originator=" << originator << "\n";
  }


  // Special method to send initial message to owning processor
  void CutPtThd::preBroadcastMessage(const int &owningProc) {
    if (ptrParRMA->args->debug >= 25)
      ucout << "CutPtThd root send to " << owningProc << '\n';
    // A negative value for 'originator' indicates special root message
    originator = -1;
    // Grab a buffer from the same pool used for broadcasts
    PackBuffer *outBuf = outQueue.getFree();
    // Fill it in the usual way
    relayLoadBuffer(outBuf);
    // Send it to the owning processor
    outQueue.send(outBuf, owningProc, tag);
    // Make sure message counters are updated correctly
    global->recordMessageSent(this);
  }

  //////////////////////////// parRMA methods ////////////////////////////

  parRMA::parRMA(pebblParams *param, MPI_Comm comm_) : RMA(), cutPtCaster(NULL) { // , mpiComm(comm_)

    // Default is not to spend time on a dumb ramp up
    rampUpPoolLimitFac = 1.0;

    setPebblParameters(param);

    branchChoice::setupMPI();
  }


  // Destructor.
  parRMA::~parRMA() {
    if (cutPtCaster == 0)
      delete cutPtCaster;
    branchChoice::freeMPI();
  }


  /// Note: use VB flag?
  void parRMA::reset(bool VBflag) {
    RMA::reset();
    registerFirstSolution(new rmaSolution(this));
    if (cutPtCaster) {
      delete cutPtCaster;
      cutPtCaster = NULL;
    }
    parallelBranching::reset();
    // parallelBranching::reset(VBflag);
  }


  void parRMA::placeTasks() {
    parallelBranching::placeTasks();
    cutPtCaster = new CutPtThd(this, cutPtBroadcastTag);
    placeTask(cutPtCaster, true, highPriorityGroup);
  }


  parallelBranchSub *parRMA::blankParallelSub() {
    parRMASub *newSP = new parRMASub();
    newSP->setGlobalInfo(this);
    return newSP;
  }


  // Pack a description of the problem.
  void parRMA::pack(PackBuffer &outBuf) {

    if (args->debug >= 10)
      ucout << "parRMA::pack invoked..." << '\n';

    // outBuf << data->numTrainObs << data->numAttrib;
    //
    // for (unsigned int i = 0; i < data->numNonZeroWtObs; ++i)
    //   outBuf << sortedObsIdx[i];
    //
    // for (unsigned int i = 0; i < data->numTrainObs; ++i) {
    //   outBuf << data->dataIntTrain[i].X << data->dataIntTrain[i].w;
    //   // outBuf << data->origTrainData[i].y;
    // }
    //
    // outBuf << data->vecNumDistVals << data->numTotalCutPts;

  } // end function parRMA::pack


  // unpack
  void parRMA::unpack(UnPackBuffer &inBuf) {

    if (args->debug >= 10)
      ucout << "parRMA::unpack invoked... " << '\n';

    // inBuf >> data->numTrainObs >> data->numAttrib;
    //
    // sortedObsIdx.resize(data->numNonZeroWtObs);
    // for (unsigned int i = 0; i < data->numNonZeroWtObs; ++i)
    //   inBuf >> sortedObsIdx[i];
    //
    // data->dataIntTrain.resize(data->numTrainObs);
    // for (unsigned int i = 0; i < data->numTrainObs; ++i) {
    //   data->dataIntTrain[i].X.resize(data->numAttrib);
    //   inBuf >> data->dataIntTrain[i].X >> data->dataIntTrain[i].w;
    //   // inBuf >> data->origTrainData[i].y;
    // }
    //
    // inBuf >> data->vecNumDistVals >> data->numTotalCutPts;

    if (args->debug >= 10)
      ucout << "parRMA::unpack done." << '\n';

    if (args->debug >= 20) {

      ucout << "data->vecNumDistVals: " << data->vecNumDistVals << "\n";

      ucout << "wt: ";
      for (unsigned int i = 0; i < data->numTrainObs; ++i)
        ucout << data->dataIntTrain[i].w << ' ';
      ucout << "\n";

    }

  } // end function parRMA::unpack


  int parRMA::spPackSize() {
    int sizePack =  4 * (data->numAttrib) * sizeof(unsigned int) //  al << au << bl << bu
            + 2 * sizeof(unsigned int) // branchVar, cutVal
            + 3 * sizeof(unsigned int) // which Child for 3 subproblems
            + 3 * sizeof(double)       // roundedBound for 3 subproblems
            + (data->numAttrib) * sizeof(double); // deqRestAttrib
    if (args->debug >= 20)
      cout << "spPackSize: " << sizePack << "\n";

    return sizePack;

  } // end function parRMA::spPackSize


  // Virtual function overriding the simpler one in the serial case
  void parRMA::setCachedCutPts(unsigned int j, unsigned int v) 
  {
    // See if the cutpoint is already in the cache (and put it in if not)
    bool isAlreadyInCache = putInCache(j,v);

    // If this looks like a new cutpoint and we are not inside ramp-up, then
    // try to let the other processors know about it

    if (!isAlreadyInCache && !rampingUp())
    {
      if (args->debug >= 20)
        ucout << "Cut point may need to be broadcast\n";

      // Set the cut point in the broadcaster
      cutPtCaster->setCutPtThd(j,v);

      // Hash the cutpoint to its "owning" processor (much fancier than before)
      int hashTmp[2];
      hashTmp[0] = (int) j;
      hashTmp[1] = (int) v;
      int owningProc = hash_bj(hashTmp, 2) % uMPI::size;

      if (args->debug >= 20)
        ucout << "owningProc: " << owningProc 
              << ((owningProc == uMPI::rank) ? " (local)" : "") << '\n';

      // If we own the cutpoint, initiate a broadcast; otherwise send it to its
      // owning processor to see if it is really new (in which case the owning)
      // processor will broadcast it.

      if (owningProc == uMPI::rank)
        cutPtCaster->initiateBroadcast();
      else
        cutPtCaster->preBroadcastMessage(owningProc);
    }

  } // end function parRMA::setCachedCutPts

  //////////////////////// parRMASub methods /////////////////////////////////

  void parRMASub::pack(utilib::PackBuffer &outBuffer) {

    if (globalPtr->args->debug >= 20)
      ucout << "parRMASub::pack invoked...\n";

    outBuffer << al << au << bl << bu;
    outBuffer << _branchChoice.branchVar << _branchChoice.cutVal;

    for (unsigned int i = 0; i < 3; ++i) {
      outBuffer << _branchChoice.branch[i].roundedBound
                << _branchChoice.branch[i].whichChild;
      // DEBUGPRX(20,　pGlobal(),"_branchChoice.branch[i].roundedBound. :"
      //<< _branchChoice.branch[i].roundedBound << "\n");
    }

    // outBuffer << vecCheckedFeat;
    for (unsigned int j = 0; j < numAttrib(); ++j)
      outBuffer << deqRestAttrib[j];

    if (globalPtr->args->debug >= 20)
      ucout << "parRMASub::pack done. "
           << " bound: " << bound << "\n";
  } // end function parRMASub::pack


  void parRMASub::unpack(utilib::UnPackBuffer &inBuffer) {

    // DEBUGPRX(20, pGlobal(),"parRMASub::unpack invoked...\n");

    inBuffer >> al >> au >> bl >> bu;
    inBuffer >> _branchChoice.branchVar >> _branchChoice.cutVal;

    for (unsigned int i = 0; i < 3; ++i) {
      inBuffer >> _branchChoice.branch[i].roundedBound >>
          _branchChoice.branch[i].whichChild;
      // DEBUGPRX(20,pGlobal(),"branchChoice.branch[i].roundedBound. :"
      //<< _branchChoice.branch[i].roundedBound << "\n");
    }

    deqRestAttrib.resize(numAttrib());
    for (unsigned int j = 0; j < numAttrib(); ++j)
      inBuffer >> deqRestAttrib[j];
    // inBuffer >> vecCheckedFeat;

    if (globalPtr->args->debug >= 20)
      ucout << "parRMASub::unpack done. :"
           << " bound: " << bound << '\n';

  } // end function parRMASub::unpack


  // makeParallelChild
  parallelBranchSub *parRMASub::makeParallelChild(int whichChild) {

    if (global()->args->debug >= 20)
      ucout << "parRMASub::makeParallelChild invoked for: "
           << ", whichChild: " << whichChild << ", ramp-up flag: " << rampingUp()
           << '\n';

    if (whichChild == -1) {
      ucout << "which child cannot be -1!";
      return NULL;
    }

  #ifdef ACRO_VALIDATING
    if (whichChild < 0 || whichChild > 2) {
      ucout << "parRMASub::makeParallelChild: invalid request "
            << "for child " << whichChild << '\n';
      return NULL;
    }

    if ((_branchChoice.branchVar < 0) ||
        (_branchChoice.branchVar >= numAttrib())) {
      ucout << "parRMASub::makeParallelChild: invalid branching variable\n";
      return NULL;
    }
  #endif

    // If there are no cached children (because this subproblem was
    // sent from somewhere else), recreate a child, not necessarily in
    // bound-sorted order.  Otherwise, grab it from the cache.

    if (_branchChoice.branchVar > numAttrib()) {
      if (global()->args->debug >= 20)
        ucout << "ERROR in parallel! "
              << "_branchChoice.branchVar: " << _branchChoice.branchVar << '\n';
      cerr << " ERROR in parallel! "
           << "_branchChoice.branchVar: " << _branchChoice.branchVar << '\n';
      // return NULL;
    }

    if (whichChild < 0)
      ucout << "whichChild=" << whichChild << '\n';

    if (global()->args->debug >= 20)
      ucout << "whichChild=" << whichChild << '\n';

    parRMASub *temp = new parRMASub();
    temp->setGlobalInfo(globalPtr);

    if (global()->args->debug >= 20)
      ucout << "_branchChoice.branch[whichChild].whichChild="
            << _branchChoice.branch[whichChild].whichChild << '\n';
    temp->RMASubAsChildOf(this, whichChild);

    if (global()->args->debug >= 10)
      ucout << "Parallel MakeChild produces " << temp << '\n';

    if (global()->args->debug >= 10)
      ucout << "Out of parRMASub::makeParallelChild, "
              "whichChild: "
           << whichChild << " bound: " << bound << '\n';

    return temp;
  } // end function parRMASub::makeParallelChild


  // split subproblems
  void parRMASub::parStrongBranching(const unsigned int &firstIdx,
                                     const unsigned int &lastIdx) {

    bool isCheckIncumb = false;
    unsigned int numCutPtsInAttrib, countCutPts = 0;

    for (unsigned int j = 0; j < numAttrib(); ++j) { // for each attribute

      numCutPtsInAttrib = bu[j] - al[j];
      if (bl[j] > au[j])
        numCutPtsInAttrib -= (bl[j] + au[j]);

      if (countCutPts + numCutPtsInAttrib <= firstIdx) {
        countCutPts += numCutPtsInAttrib;
        (global()->args->isCountingSort()) ? countingSortEC(j) : bucketSortEC(j);
        // if ( firstAttrib<=j && j<=lastAttrib )
        //	compIncumbent(j);
        continue;
      }

      if (global()->args->debug >= 10)
        ucout << "original: ";
      printSP(j, al[j], au[j], bl[j], bu[j]);

      for (unsigned int v = al[j]; v < bu[j];
           ++v) { // for each cut-value in attribute j

        // if a cut-value is in [au, bl-1]
        if (au[j] <= v && v < bl[j]) {
          v = bl[j] - 1;
          continue;
        }

        // if this cut point is not assinged in this processor
        if (countCutPts < firstIdx) {
          ++countCutPts;
          continue;
        }
        if (countCutPts > lastIdx)
          break;

        branchingProcess(j, v);
        ++countCutPts;

        if (v == al[j])
          isCheckIncumb = true;

      } // end for each cut-value in attribute j

      (global()->args->isCountingSort()) ? countingSortEC(j) : bucketSortEC(j);

      if (isCheckIncumb)
        compIncumbent(j);

      if (j == numAttrib() - 1)
        break;

    } // end for each attribute

  } // end function parRMASub::splitSP


  void parRMASub::setLiveCachedCutPts() {

    // numLiveCachedCutPts = (# of live cut points from the cache)
    int numLiveCachedCutPts = getNumLiveCachedCutPts();

    // if numCachedCutPts is less than the percentage, check all cut points
    if (numLiveCachedCutPts <
        globalPtr->data->numTotalCutPts * globalPtr->args->fracCachedCutPts())
      return;

    // if not, only check the storedCutPts

    // count number of subproblems only discovering cut-points from the chache
    isCachedCutPts = true;
    ++global()->numCC_SP;

    unsigned int j, v;
    int l = -1;
    multimap<unsigned int, unsigned int>::iterator curr = global()->mmapCachedCutPts.begin();
    multimap<unsigned int, unsigned int>::iterator end = global()->mmapCachedCutPts.end();

    cachedCutPts.resize(numLiveCachedCutPts);

    while (curr != end) {
      j = curr->first;
      v = curr->second;
      // if (j>numAttrib() || v<0) error;
      curr++;
      if (al[j] <= v && v < bu[j])                         // if v in [al, bu)
        if (!(au[j] < bl[j] && au[j] <= v && v < bl[j])) { // if not overlapping
          cachedCutPts[++l].j = j;
          cachedCutPts[l].v = v;
        }
    }

    unsigned int size = uMPI::size;
    unsigned int rank = uMPI::rank;

    if ((size > cachedCutPts.size() && rank >= cachedCutPts.size()) ||
        cachedCutPts.size() == 0) {
      DEBUGPRX(20, global(), "ramp-up, this processor won't compute." << '\n');
      return;
    }

    unsigned int quotient = cachedCutPts.size() / size;
    unsigned int remainder = cachedCutPts.size() % size;

    unsigned int firstIndex = rank * quotient + min((int) rank, (int) remainder);
    unsigned int lastIndex = firstIndex + quotient + (rank < remainder) - 1;

    parCachedBranching(firstIndex, lastIndex);

  } // end setLiveCachedCutPts function


  // branching using cut-point caching methods
  void parRMASub::parCachedBranching(unsigned int firstIdx,
                                     unsigned int lastIdx) {

    unsigned int size = uMPI::size;
    unsigned int rank = uMPI::rank;
    unsigned int quotient = numAttrib() / size;
    unsigned int remainder = numAttrib() % size;
    unsigned int firstAttrib = rank * quotient + min((int) rank, (int) remainder);
    unsigned int lastAttrib = firstAttrib + quotient + (rank < remainder) - 1;

    // int firstJ = cachedCutPts[firstIdx].j;
    // int lastJ = cachedCutPts[lastIdx].j;

    unsigned int k = 0;

    sortCachedCutPtByAttrib();
    cachedCutPts = sortedCachedCutPts;

    for (unsigned int j = 0; j < numAttrib(); ++j) {
      while (k < cachedCutPts.size()) {
        if (j == cachedCutPts[k].j) {
          if (firstIdx <= k && k <= lastIdx)
            branchingProcess(cachedCutPts[k].j, cachedCutPts[k].v);
          ++k;
        } else
          break;
      } // end while
      // if (j==numAttrib()-1) break;
      if (firstAttrib <= j && j <= lastAttrib)
        compIncumbent(j);
    } // for each attribute

  } // end RMASub::cachedBranching

  // Bound computation -- unless we're in ramp-up, just do the same
  // thing as the serial three-way code.  If we're in ramp-up, try to
  // parallelize the strong branching procedure
  void parRMASub::boundComputation(double *controlParam) {

    if (global()->args->debug >= 20)
      ucout << "In parRMASub::boundComputation, ramp-up flag=" << rampingUp()
           << '\n';

    if (!rampingUp()) {
      RMASub::boundComputation(controlParam);
      return;
    }

    numTiedSols = 1;

    // Special handling of ramp-up
    if (global()->args->debug >= 20)
      ucout << "Ramp-up bound computation\n";

    //************************************************************************
    coveredObs = global()->sortedObsIdx;
    // sort each feature based on [au, bl]
    for (unsigned int j = 0; j < numAttrib(); j++)
      bucketSortObs(j);
    setInitialEquivClass(); // set vecEquivClass
    /*
      #ifdef ACRO_HAVE_MPI
      if (uMPI::rank==0) {
      #endif //  ACRO_HAVE_MPI
      //if (global()->incumbentValue < globalPtr->guess->value) {
      if (workingSol()->value < globalPtr->guess->value) {
      ucout << "coveredObs3: " << coveredObs;
      //global()->incumbentValue = globalPtr->guess->value;
      workingSol()->value = globalPtr->guess->value;
      workingSol()->a = globalPtr->guess->a;
      workingSol()->b = globalPtr->guess->b;
      foundSolution();
      DEBUGPR(5, workingSol()->printSolution());
      if (global()->args->debug>=10) workingSol()->checkObjValue1(workingSol()->a,
      workingSol()->b, coveredObs,sortedECidx ));
      }
      #ifdef ACRO_HAVE_MPI
      }
      #endif //  ACRO_HAVE_MPI
    */
    // Better incumbents may have been found along the way
    // pGlobal()->rampUpIncumbentSync();

    //**************************************************************
    // Figure which variables go on which processor.  Make them as even as
    // possible
    // -- the first (remainder) processors have one more variable

    isCachedCutPts = false;
    unsigned int size = uMPI::size;
    unsigned int rank = uMPI::rank;

    // if there are enough discoverd cut points (storedCutPts) check only the list
    /*if ( global()->fracCachedCutPts() < 1.0 ) {
      setLiveCachedCutPts();
      }*/

    if (!isCachedCutPts) { // check all cut points

      setNumLiveCutPts();

      unsigned int quotient = numLiveCutPts / size;
      unsigned int remainder = numLiveCutPts % size;

      unsigned int firstIndex = rank * quotient + min((int)rank, (int) remainder);
      unsigned int lastIndex = firstIndex + quotient + (rank < remainder) - 1;

      if (global()->args->debug >= 20)
        ucout << "numLiveCutPts = " << numLiveCutPts
             << ", quotient  = " << quotient << ", remainder = " << remainder
             << ", firstIndex = " << firstIndex << ", lastIndex = " << lastIndex
             << '\n';

      parStrongBranching(firstIndex, lastIndex);
    }

    printCurrentBounds();

    if (global()->args->debug >= 1)
      ucout << "Best local choice is " << _branchChoice
           << " numTiedSols: " << numTiedSols << '\n';

    /******************* rampUpIncumbentSync *******************/

    if (global()->args->debug >= 1)
      ucout << rank
           << ": BEFORE rampUpIncumbentSync():" << pGlobal()->rampUpMessages
           << '\n';

    // Better incumbents may have been found along the way
    pGlobal()->rampUpIncumbentSync();

    if (global()->args->debug >= 1)
      ucout << rank
           << ": AFTER rampUpIncumbentSync():" << pGlobal()->rampUpMessages
           << '\n';

    /******************* Global Choice *******************/
    // Now determine the globally best branching choice by global reduction.
    // Use the special MPI type and combiner for branch choices.

    branchChoice bestBranch;

    if (global()->args->debug >= 1)
      ucout << rank << ": before reduceCast: " << pGlobal()->rampUpMessages
           << '\n';

    if (global()->args->branchSelection() == 0) {
      MPI_Scan(&_branchChoice, &bestBranch, 1, branchChoice::mpiType,
               branchChoice::mpiBranchSelection, MPI_COMM_WORLD);

      MPI_Bcast(&bestBranch, 1, branchChoice::mpiType, size - 1, MPI_COMM_WORLD);
    } else {
      uMPI::reduceCast(&_branchChoice, &bestBranch, 1, branchChoice::mpiType,
                       branchChoice::mpiCombiner);
    }

    pGlobal()->rampUpMessages += 2 * (uMPI::rank > 0);

    if (global()->args->debug >= 20)
      ucout << rank << ": after reduceCast:" << pGlobal()->rampUpMessages << '\n';

    if (global()->args->debug >= 20)
      ucout << "Best global choice is " << bestBranch << '\n';

    /******************* Cache cut-point *******************/
    if (global()->args->fracCachedCutPts() < 1.0)
      globalPtr->setCachedCutPts(bestBranch.branchVar, bestBranch.cutVal);

    /************************************************************/

    // If this processor has the best choice,  there is nothing to do.
    // Otherwise, adjust everything so it looks like we made the globally best
    // choice.
    if (bestBranch.branchVar != _branchChoice.branchVar ||
        bestBranch.cutVal != _branchChoice.cutVal) {
      if (global()->args->debug >= 10)
        ucout << "Adjusting local choice\n";
      _branchChoice = bestBranch;
    }

    bound = _branchChoice.branch[0].roundedBound; // look ahead bound
    setState(bounded);

    // If objValue >= bound, then we found a solution.
    if (workingSol()->value >= _branchChoice.branch[0].roundedBound) {
      // workingSol()->printSolution();
      foundRMASolution(notSynchronous);
      setState(dead);
    }

    // If roundedBound=-1, then we found a solution.
    if (_branchChoice.branch[0].roundedBound == -1) {
      // workingSol()->printSolution();
      foundRMASolution(notSynchronous);
      setState(dead);
    }

    // if the stored cut point is
    if (globalPtr->mmapCachedCutPts.size() >=
        size * global()->args->rampUpSizeFact())
      pGlobal()->rampUpPoolLimitFac = 0;

    // DEBUGPRX(50, global(), " bound: " << bound << ", sol val=" <<
    // getObjectiveVal() << '\n');
    if (global()->args->debug >= 10)
      ucout << "Ending ramp-up bound computation for bound: " << bound << '\n';

  } // end function parRMASub::boundComputation


  void parRMASub::setNumLiveCutPts() {

    numLiveCutPts = 0;
    numRestAttrib = 0;

    if (global()->args->debug >= 10)
      ucout << "deqRestAttrib: " << deqRestAttrib << "\n";
    // compute the total cut points

    for (unsigned int j = 0; j < numAttrib(); ++j) { // for each attribute

      if (deqRestAttrib[j])
        numRestAttrib++; // count how many X are restricted

      // calculate total numbers of cut points
      if ((global()->args->fracLimitAttrib() == 1) ||
          (numRestAttrib > global()->args->fracLimitAttrib() * numAttrib() &&
           deqRestAttrib[j])) {

        numLiveCutPts += bu[j] - al[j];

        if (bl[j] > au[j])
          numLiveCutPts -= (bl[j] + au[j]);

      }

    } // end for each attribute

    if (numLiveCutPts == 0) {
      if (global()->args->debug >= 20)
        ucout << "No cut points to check!\n";
      setState(dead);
    }

  } // end setNumLiveCutPts

} // namespace pebblRMA

//#endif // ACRO_HAVE_MPI
