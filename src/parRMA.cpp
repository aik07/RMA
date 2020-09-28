/**********************************************************
 *  File name:   parRMA.cpp
 *  Author:      Ai Kagawa
 *  Description: a source file for the parallel RMA solver using PEBBL
 **********************************************************/

#include "parRMA.h"
#include "serRMA.h"

//#ifdef ACRO_HAVE_MPI

namespace pebblRMA {

////////////////////// cutPtRecThd methods (Begining)
/////////////////////////////////////

// constructor of cutPtThd class
CutPtThd::CutPtThd(parRMA *global_, MessageID msgID)
    : broadcastPBThread(global_, "Cut Point Receiver", "cutPtRec",
                        "PaleTurquoise3",
                        5,  // logLevel
                        35, // debug level
                        3 * sizeof(int), msgID,
                        3), // Tree radix -- increase later
      j(-1), v(-1), ptrParRMA(global_){};

// Run method.  This is only invoked if a message arrives.
bool CutPtThd::unloadBuffer() {

  inBuf >> j >> v >> originator;
  if (ptrParRMA->args->debug >= 10)
    cout << "cutPtThd message received from " << status.MPI_SOURCE
          << "(j, v)=(" << j << ", " << v << ")"
          << ", originator=" << originator << '\n';

  bool seenAlready = false;
  multimap<unsigned int, unsigned int>::iterator it, itlow, itup;

  itlow = ptrParRMA->mmapCachedCutPts.lower_bound(j); // itlow points to
  itup  = ptrParRMA->mmapCachedCutPts.upper_bound(j);  // itup points to

  // print range [itlow,itup):
  for (it = itlow; it != itup; ++it) {
    if ((*it).first == (unsigned int) j && (*it).second == (unsigned int) v)
      seenAlready = true;
    if (ptrParRMA->args->debug >= 10)
      cout << (*it).first << " => " << (*it).second << '\n';
  }

  if (ptrParRMA->args->debug >= 10)
    cout << "cut point (" << j << ", " << v << ") ";
  if (ptrParRMA->args->debug >= 10)
    cout << (seenAlready ? "is already in cache\n" : "is new\n");

  if (originator < 0) {
    if (seenAlready)
      return false;
    originator = uMPI::rank;
  }

  // if not in the hash table, insert the cut point into the hash table.
  if (!seenAlready)
    ptrParRMA->mmapCachedCutPts.insert(make_pair(j, v));

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
    cout << "cutPtThd writing (feat, cutVal)=(" << j << ", " << v
          << "), originator=" << originator << "\n";
}

// Special method to send initial message to owning processor
void CutPtThd::preBroadcastMessage(const int &owningProc) {
  if (ptrParRMA->args->debug >= 25)
    cout << "CutPtThd root send to " << owningProc << '\n';
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

///////////////////////////////////// parRMA methods
/////////////////////////////////////////

parRMA::parRMA(MPI_Comm comm_) : RMA(), cutPtCaster(NULL) { // , mpiComm(comm_)

  // Default is not to spend time on a dumb ramp up
  rampUpPoolLimitFac = 1.0;

  /*
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
  */
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
};

// Pack a description of the problem.
void parRMA::pack(PackBuffer &outBuf) {

  if (args->debug >= 20)
    cout << "parRMA::pack invoked..." << '\n';

  outBuf << data->numOrigObs << data->numAttrib;

  for (unsigned int i = 0; i < data->numOrigObs; ++i)
    outBuf << sortedObsIdx[i];

  for (unsigned int i = 0; i < data->numOrigObs; ++i) {
    outBuf << data->intTrainData[i].X << data->intTrainData[i].w;
    // outBuf << data->origTrainData[i].y;
  }

  outBuf << data->distFeat << data->numTotalCutPts;

} // end function parRMA::pack

// unpack
void parRMA::unpack(UnPackBuffer &inBuf) {

  if (args->debug >= 20)
    cout << "parRMA::unpack invoked... " << '\n';

  inBuf >> data->numOrigObs >> data->numAttrib;

  sortedObsIdx.resize(data->numOrigObs);
  for (unsigned int i = 0; i < data->numOrigObs; ++i)
    inBuf >> sortedObsIdx[i];

  data->intTrainData.resize(data->numOrigObs);
  for (unsigned int i = 0; i < data->numOrigObs; ++i) {
    data->intTrainData[i].X.resize(data->numAttrib);
    inBuf >> data->intTrainData[i].X >> data->intTrainData[i].w;
    // inBuf >> data->origTrainData[i].y;
  }

  inBuf >> data->distFeat >> data->numTotalCutPts;

  if (args->debug >= 20)
    cout << "parRMA::unpack done." << '\n';

  if (args->debug >= 20) {

    cout << " data->distFeat: ";

    for (unsigned int i = 0; i < data->numAttrib; ++i)
      cout << data->distFeat[i] << ", ";

    for (unsigned int i = 0; i < data->numOrigObs; ++i)
      cout << " wt: " << data->intTrainData[i].w << '\n';
  }

} // end function parRMA::unpack

int parRMA::spPackSize() {
  return 5 * (data->numAttrib + 2) *
             sizeof(int) //  al << au << bl << bu << subState
         + 2 * sizeof(int) +
         3 * (sizeof(double) + sizeof(int))      // size of branchChoice
         + (data->numAttrib + 2) * sizeof(bool); // vecCheckedFeat
} // end function parRMA::spPackSize

// using virtual function
void parRMA::setCachedCutPts(const unsigned int &j, const unsigned int &v) {

  bool isAlreadyInCache = false;
  multimap<unsigned int, unsigned int>::iterator it, itlow, itup;

  itlow = mmapCachedCutPts.lower_bound(j); // itlow points to
  itup  = mmapCachedCutPts.upper_bound(j);  // itup points to

  // print range [itlow,itup):
  for (it = itlow; it != itup; ++it) {
    if ((*it).first == j && (*it).second == v)
      isAlreadyInCache = true;
    if (args->debug >= 10)
      cout << (*it).first << " => " << (*it).second << '\n';
  }

  if (args->debug >= 10)
    cout << "cut point (" << j << ", " << v << ") ";
  if (args->debug >= 10)
    cout << (isAlreadyInCache ? "is already in cache\n" : "is new\n");

  // if not in the hash table, insert the cut point into the hash table.
  if (!isAlreadyInCache)
    mmapCachedCutPts.insert(make_pair(j, v));

  int owningProc = mmapCachedCutPts.find(j)->second % uMPI::size;
  if (args->debug >= 20)
    cout << "owningProc: " << owningProc << '\n';

  cutPtCaster->setCutPtThd(j, v);

  if (uMPI::rank == owningProc) {
    // This processor is the owning processor
    if (args->debug >= 20)
      cout << "I am the owner\n";
    cutPtCaster->initiateBroadcast();
  } else {
    if (args->debug >= 20)
      cout << "Not owner\n";
    cutPtCaster->preBroadcastMessage(owningProc);
  }

  return;
} // end function parRMA::setCachedCutPts

//////////////////////// parRMASub methods /////////////////////////////////

void parRMASub::pack(utilib::PackBuffer &outBuffer) {

  if (globalPtr->args->debug >= 20)
    cout << "parRMASub::pack invoked...\n";

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
    cout << "parRMASub::pack done. "
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
    cout << "parRMASub::unpack done. :"
         << " bound: " << bound << '\n';

} // end function parRMASub::unpack

// makeParallelChild
parallelBranchSub *parRMASub::makeParallelChild(int whichChild) {

  if (global()->args->debug >= 20)
    cout << "parRMASub::makeParallelChild invoked for: "
         << ", whichChild: " << whichChild << ", ramp-up flag: " << rampingUp()
         << '\n';

  if (whichChild == -1) {
    cerr << "which child cannot be -1!";
    return NULL;
  }

#ifdef ACRO_VALIDATING
  if (whichChild < 0 || whichChild > 2) {
    cout << "parRMASub::makeParallelChild: invalid request "
          << "for child " << whichChild << '\n';
    return NULL;
  }

  if ((_branchChoice.branchVar < 0) ||
      (_branchChoice.branchVar >= numAttrib())) {
    cout << "parRMASub::makeParallelChild: invalid branching variable\n";
    return NULL;
  }
#endif

  // If there are no cached children (because this subproblem was
  // sent from somewhere else), recreate a child, not necessarily in
  // bound-sorted order.  Otherwise, grab it from the cache.

  if (_branchChoice.branchVar > numAttrib()) {
    if (global()->args->debug >= 20)
      cout << "ERROR in parallel! "
            << "_branchChoice.branchVar: " << _branchChoice.branchVar << '\n';
    cerr << " ERROR in parallel! "
         << "_branchChoice.branchVar: " << _branchChoice.branchVar << '\n';
    // return NULL;
  }

  if (whichChild < 0)
    cerr << "whichChild=" << whichChild << '\n';

  if (global()->args->debug >= 20)
    cout << "whichChild=" << whichChild << '\n';

  parRMASub *temp = new parRMASub();
  temp->setGlobalInfo(globalPtr);

  if (global()->args->debug >= 20)
    cout << "_branchChoice.branch[whichChild].whichChild="
          << _branchChoice.branch[whichChild].whichChild << '\n';
  temp->RMASubAsChildOf(this, whichChild);

  if (global()->args->debug >= 10)
    cout << "Parallel MakeChild produces " << temp << '\n';

  if (global()->args->debug >= 10)
    cout << "Out of parRMASub::makeParallelChild, "
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
      (global()->args->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
      // if ( firstAttrib<=j && j<=lastAttrib )
      //	compIncumbent(j);
      continue;
    }

    if (global()->args->debug >= 10)
      cout << "original: ";
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

    (global()->args->countingSort()) ? countingSortEC(j) : bucketSortEC(j);

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
      globalPtr->data->numTotalCutPts * globalPtr->args->perCachedCutPts())
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
}

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
    }
    // if (j==numAttrib()-1) break;
    if (firstAttrib <= j && j <= lastAttrib)
      compIncumbent(j);
  }

} // end RMASub::cachedBranching

// Bound computation -- unless we're in ramp-up, just do the same
// thing as the serial three-way code.  If we're in ramp-up, try to
// parallelize the strong branching procedure
void parRMASub::boundComputation(double *controlParam) {

  if (global()->args->debug >= 20)
    cout << "In parRMASub::boundComputation, ramp-up flag=" << rampingUp()
         << '\n';

  if (!rampingUp()) {
    RMASub::boundComputation(controlParam);
    return;
  }

  NumTiedSols = 1;

  // Special handling of ramp-up
  if (global()->args->debug >= 20)
    cout << "Ramp-up bound computation\n";

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
    cout << "coveredObs3: " << coveredObs;
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
  /*if ( global()->perCachedCutPts() < 1.0 ) {
    setLiveCachedCutPts();
    }*/

  if (!isCachedCutPts) { // check all cut points

    setNumLiveCutPts();

    unsigned int quotient = numLiveCutPts / size;
    unsigned int remainder = numLiveCutPts % size;

    unsigned int firstIndex = rank * quotient + min((int)rank, (int) remainder);
    unsigned int lastIndex = firstIndex + quotient + (rank < remainder) - 1;

    if (global()->args->debug >= 20)
      cout << "numLiveCutPts = " << numLiveCutPts
           << ", quotient  = " << quotient << ", remainder = " << remainder
           << ", firstIndex = " << firstIndex << ", lastIndex = " << lastIndex
           << '\n';

    parStrongBranching(firstIndex, lastIndex);
  }

  printCurrentBounds();

  if (global()->args->debug >= 1)
    cout << "Best local choice is " << _branchChoice
         << " NumTiedSols: " << NumTiedSols << '\n';

  /******************* rampUpIncumbentSync *******************/

  if (global()->args->debug >= 1)
    cout << rank
         << ": BEFORE rampUpIncumbentSync():" << pGlobal()->rampUpMessages
         << '\n';

  // Better incumbents may have been found along the way
  pGlobal()->rampUpIncumbentSync();

  if (global()->args->debug >= 1)
    cout << rank
         << ": AFTER rampUpIncumbentSync():" << pGlobal()->rampUpMessages
         << '\n';

  /******************* Global Choice *******************/
  // Now determine the globally best branching choice by global reduction.
  // Use the special MPI type and combiner for branch choices.

  branchChoice bestBranch;

  if (global()->args->debug >= 1)
    cout << rank << ": before reduceCast: " << pGlobal()->rampUpMessages
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
    cout << rank << ": after reduceCast:" << pGlobal()->rampUpMessages << '\n';

  if (global()->args->debug >= 20)
    cout << "Best global choice is " << bestBranch << '\n';

  /******************* Cache cut-point *******************/
  if (global()->args->perCachedCutPts() < 1.0)
    globalPtr->setCachedCutPts(bestBranch.branchVar, bestBranch.cutVal);

  /************************************************************/

  // If this processor has the best choice,  there is nothing to do.
  // Otherwise, adjust everything so it looks like we made the globally best
  // choice.
  if (bestBranch.branchVar != _branchChoice.branchVar ||
      bestBranch.cutVal != _branchChoice.cutVal) {
    if (global()->args->debug >= 10)
      cout << "Adjusting local choice\n";
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
    cout << "Ending ramp-up bound computation for bound: " << bound << '\n';

} // end function parRMASub::boundComputation

void parRMASub::setNumLiveCutPts() {
  numLiveCutPts = 0;
  numRestAttrib = 0;
  if (global()->args->debug >= 10)
    cout << "deqRestAttrib: " << deqRestAttrib << "\n";
  // compute the total cut points
  for (unsigned int j = 0; j < numAttrib(); ++j) {
    if (deqRestAttrib[j])
      numRestAttrib++; // count how many X are restricted
    // calculate total numbers of cut points
    if ((global()->args->perLimitAttrib() == 1) ||
        (numRestAttrib > global()->args->perLimitAttrib() * numAttrib() &&
         deqRestAttrib[j])) {
      numLiveCutPts += bu[j] - al[j];
      if (bl[j] > au[j])
        numLiveCutPts -= (bl[j] + au[j]);
    }
  }
  if (numLiveCutPts == 0) {
    if (global()->args->debug >= 20)
      cout << "No cut points to check!\n";
    setState(dead);
  }
}

} // namespace pebblRMA

//#endif // ACRO_HAVE_MPI
