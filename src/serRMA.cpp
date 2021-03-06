// **********************************************************
// Author:      Ai Kagawa
// Description: a source file for serial RMA solver using PEBBL
// **********************************************************


#include "serRMA.h"


namespace pebblRMA {


#ifdef ACRO_HAVE_MPI


  void branchChoiceCombiner(void *invec, void *inoutvec, int *len,
                            MPI_Datatype *datatype) {

  #ifdef ACRO_VALIDATING
    if (*datatype != branchChoice::mpiType) {
      cerr << "Datatype error in branchChoiceCombiner\n";
      exit(1);
    }
  #endif // ACRO_VALIDATING

    branchChoice *inPtr    = (branchChoice *)invec;
    branchChoice *inOutPtr = (branchChoice *)inoutvec;
    unsigned int n         = *len;

    for (unsigned int i = 0; i < n; ++i)
      if (inPtr[i] < inOutPtr[i])
        inOutPtr[i] = inPtr[i];
  } // end branchChoiceCombiner function


  void branchChoiceRand(branchChoice *in, branchChoice *inout, int *len,
                        MPI_Datatype *datatype) {

  #ifdef ACRO_VALIDATING
    if (*datatype != branchChoice::mpiType) {
      cerr << "Datatype error in branchChoiceRand\n";
      exit(1);
    }
  #endif

    branchChoice c;
    unsigned int n = *len;

    for (unsigned int i = 0; i < n; ++i) {
      if (in < inout)
        c = *in;
      else if (in > inout)
        c = *inout;
      else {
        unsigned int n1 = in->numTiedSols;
        unsigned int n2 = inout->numTiedSols;

        srand((n1 + n2) * time(NULL) * 100); // : srand(1);
        double rand_num = (rand() % (n1 + n2 + 1)) / (double)(n1 + n2);
        // double rand_num = ( rand() % 1001 ) / (double) 1000 ;
        if (rand_num < (double)n1 / (n1 + n2)) {
          // ucout << "rand_num: " << rand_num << endl;
          c = *in;
        } else {
          c = *inout;
        }
        c.numTiedSols = n1 + n2;
      }
      *inout = c;
      in++;
      inout++;
    }

  } // end branchChoiceRand function

  //////////////////// branchChoice class methods ////////////////////

  void branchChoice::setupMPIDatum(void *address, MPI_Datatype thisType,
                                   MPI_Datatype *type, MPI_Aint base,
                                   MPI_Aint *disp, int *blocklen,
                                   unsigned int j) {
    MPI_Get_address(address, &disp[j]);
    disp[j]     -= base;
    type[j]     = thisType;
    blocklen[j] = 1;
  } // end setupMPIDatum function


  void branchChoice::setupMPI() {

    unsigned int arraySize = 3 * 3 + 3;
    int          j = 0;
    MPI_Datatype type[arraySize];
    int          blocklen[arraySize];
    MPI_Aint     disp[arraySize];
    MPI_Aint     base;
    branchChoice example;
    MPI_Get_address(&example, &base);

    for (unsigned int i = 0; i < 3; ++i) {
      setupMPIDatum(&(example.branch[i].roundedBound), MPI_DOUBLE, type, base,
                    disp, blocklen, j++);
      setupMPIDatum(&(example.branch[i].exactBound), MPI_DOUBLE, type, base, disp,
                    blocklen, j++);
      setupMPIDatum(&(example.branch[i].whichChild), MPI_INT, type, base, disp,
                    blocklen, j++);
    }

    setupMPIDatum(&(example.branchVar),   MPI_INT, type, base, disp, blocklen, j++);
    setupMPIDatum(&(example.cutVal),      MPI_INT, type, base, disp, blocklen, j++);
    setupMPIDatum(&(example.numTiedSols), MPI_INT, type, base, disp, blocklen,
                  j++);

    MPI_Type_create_struct(j, blocklen, disp, type, &mpiType);
    MPI_Type_commit(&mpiType);

    MPI_Op_create(branchChoiceCombiner, true, &mpiCombiner);
    MPI_Op_create((MPI_User_function *)branchChoiceRand, true, &mpiBranchSelection);

  } // end setupMPI function


  void branchChoice::freeMPI() {
    MPI_Op_free(&mpiCombiner);
    MPI_Op_free(&mpiBranchSelection);
    MPI_Type_free(&mpiType);
  }


  MPI_Datatype branchChoice::mpiType      ; //= MPI_UB;
  MPI_Op branchChoice::mpiCombiner        = MPI_OP_NULL;
  MPI_Op branchChoice::mpiBranchSelection = MPI_OP_NULL;

  #endif


  branchChoice::branchChoice() : branchVar(MAXINT), cutVal(-1) {
    for (unsigned int i = 0; i < 3; ++i)
      branch[i].set(MAXDOUBLE, 1e-5);
  }


  branchChoice::branchChoice(double a, double b, double c, int cut, int j) {

    branch[0].set(a, 1e-5);
    branch[1].set(b, 1e-5);
    branch[2].set(c, 1e-5);

    for (unsigned int i = 0; i < 3; ++i)
      branch[i].whichChild = i;

    cutVal    = cut;
    branchVar = j;

  } // end branchChoice function


  void branchChoice::setBounds(double a, double b, double c, int cut, int j) {

    branch[0].set(a, 1e-5);
    branch[1].set(b, 1e-5);
    branch[2].set(c, 1e-5);

    for (unsigned int i = 0; i < 3; ++i)
      branch[i].whichChild = i;

    cutVal    = cut;
    branchVar = j;

  } // end setBounds function


  // Primitive sort, but only three elements
  void branchChoice::sortBounds() {
    possibleSwap(0, 1);
    possibleSwap(0, 2);
    possibleSwap(1, 2);
  }


  bool branchChoice::operator<(const branchChoice &other) const {
    for (unsigned int i = 0; i < 3; ++i) {
      if (branch[i].roundedBound < other.branch[i].roundedBound)
        return true;
      else if (branch[i].roundedBound > other.branch[i].roundedBound)
        return false;
    }
    return branchVar < other.branchVar;
  }


  bool branchChoice::operator==(const branchChoice &other) const {
    for (unsigned int i = 0; i < 3; ++i) {
      if (branch[i].roundedBound == other.branch[i].roundedBound)
        continue;
      else
        return false;
    }
    return true;
  }


  void branchChoice::possibleSwap(size_type i1, size_type i2) {

    double roundedBound1 = branch[i1].roundedBound;
    double roundedBound2 = branch[i2].roundedBound;

    if (roundedBound1 < roundedBound2) {
      branchItem tempItem(branch[i1]);
      branch[i1] = branch[i2];
      branch[i2] = tempItem;
    }

  } // end possibleSwap function

  //////////////////// branchItem class methods ////////////////////

  void branchItem::set(double bound, double roundQuantum) {
    exactBound = bound;
    if (roundQuantum == 0)
      roundedBound = bound;
    else
      roundedBound = floor(bound / roundQuantum + 0.5) * roundQuantum;
    whichChild = -1;
  }

  //////////////////// RMA class methods ////////////////////

  // RMA constructor
  RMA::RMA() : workingSol(this), numCC_SP(0) { //, numTotalCutPts(0)

    min_num_required_args = 1;
    branchingInit(maximization, relTolerance, absTolerance);

    workingSol.serial = 0;
    workingSol.sense  = maximization;

  } // end RMA class constructor


  // RMA constructor
  RMA::RMA(pebblParams *param) : workingSol(this), numCC_SP(0) { //, numTotalCutPts(0)

    // version_info += ", RMA example 1.1";
    min_num_required_args = 1;
    branchingInit(maximization, relTolerance, absTolerance);

    workingSol.serial = 0;
    workingSol.sense  = maximization;

    setPebblParameters(param);

  }; //  end RMA class constructor


  RMA::~RMA() {

    // if % of cached cutpoints is less than 100
    if (args->fracCachedCutPts() < 1) {

      int recvbuf = numCC_SP;

      DEBUGPRX(1, this, "Local non-strong branching SP is: " << numCC_SP << "\n");

      uMPI::reduceCast(&numCC_SP, &recvbuf, 1, MPI_INT, MPI_SUM);
      // MPI_Reduce(&sendbuf, &recvbuf, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

      // Print the result
      if (uMPI::rank == 0)
        ucout << "Total non-strong branching SP is: " << recvbuf << "\n";

    } // end if % of cached cutpoints is less than 100

    workingSol.decrementRefs();

    //if (isInitGuess()) workingSol.decrementRefs();

    /*
      if (verifyLog()) {
      verifyLogFile() << endl; //<< "result " << fathomValue() << endl;
      delete _vlFile;    // Doesn't delete file; actually closes it
      }//*/
    // workingSol.decrementRefs();

  } // end RMA class destructor


  void RMA::setPebblParameters(pebblParams *param) {

    this->debug = param->debug;

    this->statusPrintCount = param->statusPrintCount;
    this->statusPrintSeconds = param->statusPrintSeconds;

    this->depthFirst = param->depthFirst;
    this->breadthFirst = param->breadthFirst;

    this->initialDive = param->initialDive;
    this->integralityDive = param->integralityDive;

    this->lazyBounding = param->lazyBounding;
    this->eagerBounding = param->eagerBounding;

    this->relTolerance = param->relTolerance;
    this->absTolerance = param->absTolerance;

    this->earlyOutputMinutes = param->earlyOutputMinutes;
    this->startIncumbent = param->startIncumbent;
    this->validateLog = param->validateLog;
    this->heurLog = param->heurLog;
    this->loadLogSeconds = param->loadLogSeconds;
    this->loadLogWriteSeconds = param->loadLogWriteSeconds;

    this->maxSPBounds = param->maxSPBounds;
    this->maxCPUMinutes = param->maxCPUMinutes;
    this->maxWallMinutes = param->maxWallMinutes;

    this->haltOnIncumbent = param->haltOnIncumbent;

    this->printAbortMessage = param->printAbortMessage;
    this->printIntMeasure = param->printIntMeasure;
    this->printDepth = param->printDepth;

    this->debugPrecision = param->debugPrecision;
    this->suppressWarnings = param->suppressWarnings;
    this->loadLogWriteSeconds = param->loadLogWriteSeconds;
    this->loadMeasureDegree = param->loadMeasureDegree;

    this->enumRelTol = param->enumRelTol;
    this->enumAbsTol = param->enumAbsTol;
    this->enumCutoff = param->enumCutoff;
    this->enumCount = param->debugPrecision;
    this->enumHashSize = param->enumHashSize;

    this->debug_solver_params = param->debug_solver_params;
    this->use_abort = param->use_abort;
    this->version_flag = param->version_flag;
    this->printFullSolution = param->printFullSolution;
    this->solFileName = param->solFileName;
    this->printSpTimes = param->printSpTimes;

  } // end setPebblParameters function


  // rmaSolution *RMA::initialGuess() {
  //
  //   //workingSol.reset(data->numAttrib, data->vecNumDistVals);
  //
  //   if (args->isInitGuess())
  //     return guess;
  //   else
  //     return NULL;
  //
  // } // end initialGuess function


  void RMA::setInitialGuess(bool isPosIncumb, double maxObjValue,
                            vector<unsigned int> lowerBound,
                            vector<unsigned int> upperBound) {

    if (args->debug>=1 ) ucout << "setInitialGuess\n";

    guess = new rmaSolution(this);
    guess    ->setSolution(isPosIncumb, maxObjValue, lowerBound, upperBound);

    // workingSol.setSolution(isPosIncumb, maxObjValue, lowerBound, upperBound);
    // workingSol.foundRMASolution(synchronous);

  } // end setInitialGuess function


  // This routine returns a bool that is true if cutpoint (j,v) is already
  // in the cutpoint cache, otherwise false.  If (j,v) is not in the cache, 
  // it inserts it.

  bool RMA::putInCache(unsigned int j, unsigned int v)
  {
    if (args->debug >= 10)
      ucout << "putInCache called for (" << j << ',' << v << ")\n";

    multimap<unsigned int, unsigned int>::iterator it, itlow, itup;

    itlow = mmapCachedCutPts.lower_bound(j);
    itup  = mmapCachedCutPts.upper_bound(j);

    for (it = itlow; it != itup; ++it)
    {
      unsigned int u = (*it).second;
      if (u == v)
      {
        if (args->debug >= 10)
          ucout << '(' << j << ',' << v << ") already in cache\n";
        return true;             // Found in cache already
      }
      if (u > v)
        break;                   // Already higher values, so not in cache
    }

    if (args->debug >= 10)
      ucout << '(' << j << ',' << v << ") not found in cache\n";

    if (0 > j || j > data->numAttrib)
      ucout << "ERROR! j is out of range for setCachedCutPts";
    else if (0 > v || v > data->vecNumDistVals[j]-2)
      ucout << "ERROR! v is out of range for setCachedCutPts";
    else
      mmapCachedCutPts.insert(make_pair(j, v));

    return false;
  }


  // setWeight
  void RMA::setWeight(vector<double> wt, vector<unsigned int> train) {
    for (unsigned int i = 0; i < wt.size(); ++i)
      data->dataIntTrain[train[i]].w = wt[i];
  }


  // print RMA solution, time to solve, and # of bounded subproblems
  void RMA::printSolutionTime(const double &timeCPU) {

    double global_solution = workingSol.value; // set the current node's RMA solution value
    int    total_nodes     = subCount[2];      // set the current node's subproblems which

  #ifdef ACRO_HAVE_MPI

    if (uMPI::size>1) { // if the PEBBL RMA is solved in parllel

      // ucout << "reduce " << " " << total_nodes << "\n";

      // set the global solution value to be the maximum among all nodes's solutions
      MPI_Reduce(&workingSol.value, &global_solution,
                 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

      // set the total node to be the sum of the nodes among all nodes
      MPI_Reduce(&subCount[2],      &total_nodes,
                 1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);

    } // end if the PEBBL RMA is solved in parllel

    if (uMPI::rank==0) {

      // print RMA solution, Time, # of nodes
      ucout << std::fixed << std::setprecision(4)
                << "ERMA Solution: "  << global_solution
                << std::fixed << std::setprecision(2)
                << " \tCPU time: "     << timeCPU     // searchTime
                << " \tNum of Nodes: " << total_nodes << "\n";

    } // end if (uMPI::rank==0)

  #endif //  ACRO_HAVE_MPI

  } // end printSolutionTime function


  // writes data with weights to a file whose name we concoct
  // from the iteration number argument; added by JE
  void RMA::writeInstanceToFile(const int &iterNum) {

    // create a file name
    stringstream s;
    s << 'w' << iterNum << '.' << problemName;
    ofstream instanceOutputFile(s.str().c_str());

    writeWeightedData(instanceOutputFile);  // write weights

    instanceOutputFile.close();  // close the file

  } // end writeInstanceToFile function


  // Routine added by JE to write out data with weights.  Note that
  // the sign of the observation is just the last attribute in the
  // "_dataStore" vector of vectors.
  void RMA::writeWeightedData(ostream &os) {

    // Set high precision and scientific notation for weights, while
    // saving old flags and precision
    int oldPrecision = os.precision(16);
    std::ios_base::fmtflags oldFlags = os.setf(ios::scientific);

    // Write data
    for (unsigned int i = 0; i < data->numTrainObs; ++i) {
      os << data->dataIntTrain[i].w << ';';
      // Restore stream state
      os.precision(oldPrecision);
      os.flags(oldFlags);
    } // end for each observation

  } // end writeWeightedData function


  // Routine added by AK to write out the number of B&B node and CPU time.
  void RMA::writeStatData(ostream &os) {
    os << subCount[2] << ';' << searchTime << '\n';
  }


  // writes the number of B&B node and CPU time; added by AK
  void RMA::writeStatDataToFile(const int &iterNum) {

    stringstream s;
    s << "BBNode_CPUTime" << '_' << problemName;

    ofstream instanceOutputFile(s.str().c_str());

    if (iterNum == 1)
      ofstream instanceOutputFile(s.str().c_str(), ofstream::out);
    else
      ofstream instanceOutputFile(s.str().c_str(), ofstream::app);

    writeStatData(instanceOutputFile);

    instanceOutputFile.close();  // close the file

  } // end writeStatDataToFile function


  // ********************* RMASub methods (start) *******************************


  void rmaSolution::reset(const unsigned int &numAttrib,
                          const vector<unsigned int> &vecNumDistVals) {
    vector<unsigned int> lower(numAttrib, 0);
    vector<unsigned int> upper(numAttrib, 0);
    for (unsigned int j=0; j<numAttrib; ++j) upper[j] = vecNumDistVals[j]-1;
    // set the current solution is positive, objective value = 0,
    // and the initial lower and upper bounds
    this->setSolution(true, 0, lower, upper);
  } // end reset function


  void rmaSolution::setSolution(const bool isPosIncumb, const double objVal,
                            vector<unsigned int> &a, vector<unsigned int> &b) {
    this->isPosIncumb = isPosIncumb;
    this->value       = objVal;
    this->a           = a;
    this->b           = b;
  } // end setSolution function


  void RMASub::setRootComputation() {

    if (global()->debug>=1 ) ucout << "setRootComputation\n";

    al.resize(numAttrib());
    au.resize(numAttrib());

    fill(al.begin(), al.end(), 0);

    for (unsigned int j=0; j<numAttrib(); ++j)
      au[j] = vecNumDistVals()[j]-1;

    bl << al;
    bu << au;

    deqRestAttrib.resize(numAttrib(), false);

    if (global()->guess != NULL) {
      if (global()->debug>=1 )
        ucout << "setInitGuess in setRootComputation\n";
      //workingSol() = global()->guess;
      workingSol()->setSolution(global()->guess->isPosIncumb,
                                global()->guess->value,
                                global()->guess->a,
                                global()->guess->b);

      foundRMASolution(synchronous);
    }

  } // end setRootComputation function


  void RMASub::boundComputation(double *controlParam) {

    // globalPtr->getSolution();

    if (global()->debug >= 10)
      ucout << "\nal: " << al << ", au: " << au
           << ", bl: " << bl << ", bu: " << bu;

    numTiedSols    = 1;
    numPosTiedSols = 0;
    numNegTiedSols = 0;

    // setCoveredObs();	// find covered observation which are in [al. bu]
    // sort each feature based on [au, bl]
    coveredObs.resize(global()->sortedObsIdx.size());
    copy(global()->sortedObsIdx.begin(), global()->sortedObsIdx.end(),
         coveredObs.begin());

    for (unsigned int j = 0; j < numAttrib(); ++j)
      bucketSortObs(j);

    setInitialEquivClass(); // set initial equivalence class, vecEquivClass

    // if there are enough discoverd cut points (storedCutPts) check only the list
    if (global()->args->fracCachedCutPts() < 1.0 &&
        global()->args->isBinarySearchCutVal())
      hybridBranching(); // hybrid branching

    else if (global()->args->isBinarySearchCutVal())
      binaryBranching(); // binary search cut point caching

    else if (global()->args->fracCachedCutPts() < 1.0)
      cutpointCaching(); // cut point caching

    else // check all cut points (strong branching)
      strongBranching();

    if (global()->debug >= 10)
      printCurrentBounds();

    if (global()->args->debug >= 5)
      ucout << "Branch choice: " << _branchChoice << "\n";

    bound = _branchChoice.branch[0].roundedBound; // look ahead bound
    setState(bounded);

    if (_branchChoice.branch[0].roundedBound < 0) {
      if (global()->debug >= 10)
        ucout << "Bound < 0. \n";
      setState(dead);
      return;
    }

    if (_branchChoice.branchVar > numAttrib()) {

      if (global()->debug >= 10) {
        ucout << "al: " << al << "au: " << au << "bl: " << bl << "bu: " << bu;
        ucout << "branchVar > numAttrib. \n";
      }

      setState(dead);

      return;

    }

    deqRestAttrib[_branchChoice.branchVar] = true;

    // If (current objValue) >= (current bound), we found the solution.
    if (workingSol()->value >= _branchChoice.branch[0].exactBound) {

      if (global()->debug >= 2) {
        workingSol()->printSolution();
        ucout << "Bound: " << _branchChoice.branch[0].exactBound << "\n";
      }

      foundRMASolution(synchronous);
      setState(dead);

      return;

    }


    ///////////////////// create listExcluded list (start) ////////////////
  #ifndef ACRO_HAVE_MPI

    sort(listExcluded.begin(), listExcluded.end());

    listExcluded.erase(unique(listExcluded.begin(), listExcluded.end()),
                       listExcluded.end());

    DEBUGPR(150, ucout << "Excluded: " << _branchChoice.branchVar << listExcluded);
    // DEBUGPR(50, ucout << " bound: " << bound << ", sol val=" <<
    // getObjectiveVal() << "\n");
  #endif
    ////////////////////////// check errors (start) ////////////////////////////
    if (_branchChoice.branchVar >= numAttrib()) {
      DEBUGPR(20, ucout << "ERROR: branch feature is invalid! (j="
                       << _branchChoice.branchVar << ")\n");
      cerr << "ERROR: branch feature is invalid! (j=" << _branchChoice.branchVar
           << ")\n";
      exit(EXIT_FAILURE);
    }

    if (_branchChoice.cutVal < 0) {
      if (global()->debug >= 20)
        ucout << "ERROR: cutValue cannot be less than 0! (cutValue="
             << _branchChoice.cutVal << ")\n";
      exit(EXIT_FAILURE);
    }

    if (_branchChoice.cutVal >= bu[_branchChoice.branchVar]) {
      if (global()->debug >= 20)
        ucout << "ERROR: cutValue cannot be >= bu[" << _branchChoice.branchVar
             << "]! (cutValue=" << _branchChoice.cutVal << ")\n";
      exit(EXIT_FAILURE);
    }

    //////////////// check errors (end) ////////////////////////
  #ifndef ACRO_HAVE_MPI
    listExcluded.push_back(_branchChoice.branchVar);
  #endif
    ///////////////// create listExclided list (end) ////////////////////
  } // end boundComputation function


  int RMASub::getNumLiveCachedCutPts() {
    // numLiveCachedCutPts = (# of live cut points from the cache)
    unsigned int j, v, numLiveCachedCutPts = 0;
    multimap<unsigned int, unsigned int>::iterator curr =
        global()->mmapCachedCutPts.begin();
    multimap<unsigned int, unsigned int>::iterator end =
        global()->mmapCachedCutPts.end();

    // count numLiveCachedCutPts and print out cached cut points
    if (global()->args->debug >= 20)
      ucout << "cached cut-points: ";

    while (curr != end) {

      j = curr->first;
      v = curr->second;

      if (global()->args->debug >= 20)
        ucout << j << ", " << v << "\n";
      // if (j>numAttrib() || v<0) break;

      curr++;

      if (al[j] <= v && v < bu[j])                       // if v in [al, bu)
        if (!(au[j] < bl[j] && au[j] <= v && v < bl[j])) // if not overlapping
          ++numLiveCachedCutPts;
    }

    if (global()->args->debug >= 20)
      ucout << "\n";

    return numLiveCachedCutPts;

  } // end getNumLiveCachedCutPts function


  // return how many children to make from current subproblem
  int RMASub::splitComputation() {

    int numChildren = 0;
    for (unsigned int i = 0; i < 3; ++i)
      if (_branchChoice.branch[i].roundedBound >= 0)
        numChildren++;

    setState(separated);

    return numChildren;

  } // end splitComputation function


  // make a child of current subproblem
  branchSub *RMASub::makeChild(int whichChild) {

    if (whichChild == -1) {
      cerr << "ERROR: No children to make!\n";
      return NULL;
    }

    if (this->_branchChoice.branchVar > numAttrib()) {
      cerr << "ERROR: It is not proper attribute!\n";
      return NULL;
    }

    RMASub *temp = new RMASub();
    temp->RMASubAsChildOf(this, whichChild);

    return temp;

  } // end makeChild function


  void RMASub::RMASubAsChildOf(RMASub *parent, int whichChild) {

    al << parent->al;
    au << parent->au;
    bl << parent->bl;
    bu << parent->bu;
    deqRestAttrib = parent->deqRestAttrib;

  #ifndef ACRO_HAVE_MPI
    listExcluded << parent->listExcluded;
  #endif

    excCutFeat << parent->excCutFeat;
    excCutVal  << parent->excCutVal;

    globalPtr = parent->global();
    branchSubAsChildOf(parent);

    // set bound
    bound      = parent->_branchChoice.branch[whichChild].roundedBound;
    whichChild = parent->_branchChoice.branch[whichChild].whichChild;

    if (global()->args->debug >= 10)
      ucout << "Bound: " << bound << "\n";

    unsigned int j = parent->_branchChoice.branchVar;
    unsigned int lowerBound, upperBound;

    if (j < 0) {
      if (global()->args->debug >= 20)
        ucout << "ERROR: feature j cannot be < 0 (j=" << j << ")\n";
      cerr << "ERROR: feature j cannot be < 0 (j=" << j << ")\n";
      return;
    } else if (j > numAttrib()) {
      if (global()->args->debug >= 20)
        ucout << "ERROR: feature j cannot be > numAttrib (j=" << j << ")\n";
      cerr << "ERROR: feature j cannot be > numAttrib (j=" << j << ")\n";
      return;
    }

    // Case 1: this feature is overlaping bl < au
    if (bl[j] < au[j] && bl[j] <= parent->_branchChoice.cutVal &&
        parent->_branchChoice.cutVal < au[j]) {
      // subproblem 1: P(al, v, v, v)
      if (whichChild == 0) {                  // find a value on [al, bl]
        au[j] = parent->_branchChoice.cutVal; // optCutValue;
        bl[j] = parent->_branchChoice.cutVal; // optCutValue;
        bu[j] = parent->_branchChoice.cutVal; // optCutValue;
      }
      // subproblem 3: P(al, v, v+1, bu)
      else if (whichChild == 2) {                 // find a value on [bl, au]
        au[j] = parent->_branchChoice.cutVal;     // optCutValue;
        bl[j] = parent->_branchChoice.cutVal + 1; // optCutValue+1;
      }
      // subproblem 2: P(v+1, v+1, v+1, bu)
      else if (whichChild == 1) {                 // find a value on [au, bu]
        al[j] = parent->_branchChoice.cutVal + 1; // optCutValue+1;
        au[j] = parent->_branchChoice.cutVal + 1; // optCutValue+1;
        bl[j] = parent->_branchChoice.cutVal + 1; // optCutValue+1;
      }

    } else { // Case 2: this feature is not overlaping au <= bl

      if (au[j] <= bl[j]) {
        lowerBound = au[j];
        upperBound = bl[j];
      } else {
        lowerBound = bl[j];
        upperBound = au[j];
      }

      // find a value on [al, au]
      if (al[j] <= parent->_branchChoice.cutVal &&
          parent->_branchChoice.cutVal < lowerBound) {
        if (whichChild == 0)                    // subproblem 1: P(al, v, bl, bu)
          au[j] = parent->_branchChoice.cutVal; // optCutValue;
        else if (whichChild == 1) { // subproblem 2: P(v+1, au, bl, bu)
          al[j] = parent->_branchChoice.cutVal + 1; // optCutValue+1;
        }
        // find new cut value on [bl, bu]
      } else if (upperBound <= parent->_branchChoice.cutVal &&
                 parent->_branchChoice.cutVal < bu[j]) {
        if (whichChild == 0)                    // subproblem 1: P(al, au, bl, v)
          bu[j] = parent->_branchChoice.cutVal; // optCutValue;
        else if (whichChild == 1) // subproblem 2: P(al, au, v+1, bu)
          bl[j] = parent->_branchChoice.cutVal + 1; // optCutValue+1;
      }

    } // end if

    if (al[j] > au[j]) {
      cerr << "Error al[j]>au[j]! since al[" << j << "]=" << al[j] << "and au["
           << j << "]=" << au[j] << endl;
      return;
    }

    if (global()->args->debug >= 10)
      ucout << "al: " << al << "au: " << au << "bl: " << bl << "bu: " << bu;

  } // end RMASubAsChildOf function


  // find a particular subproblem object to the problem description embodied in
  // the object master
  void RMASub::RMASubFromRMA(RMA *master) {
    globalPtr = master; // set a globalPtr
    // workingSol()->value = getObjectiveVal(); // set bound value as current
    // solution
    if (global()->args->debug >= 20)
      ucout << "Created blank problem, out of rmaSub:::RMASubFromRMA\n";
  } // end RMASubFromRMA function


  bool RMASub::candidateSolution() {

    if (global()->args->debug >= 20)
      ucout << "al: " << al << "au: " << au << "bl: " << bl << "bu: " << bu;

    for (unsigned int j = 0; j < numAttrib(); ++j) {
      if (al[j] != au[j])
        return false;
      if (bl[j] != bu[j])
        return false;
    }

    workingSol()->a << al;
    workingSol()->b << bu;

    // ucout << coveredObs << endl;
    // sort(coveredObs.begin(), coveredObs.begin()+coveredObs.size());
    // ucout << coveredObs << endl;

    // workingSol()->isCovered.resize(numDistObs());
    // for (int i=0; i<numDistObs(); ++i)
    //  workingSol()->isCovered[i]=false;
    // for (int i=0; i<vecEquivClass1.size(); ++i)
    // for (int j=0; j<vecEquivClass1[i].size(); ++j)

    /*
      for (int i=0; i<numDistObs(); ++i)
      for (int j=0; j<numAttrib(); ++j)
      if ( workingSol()->a[j] <= global()->dataIntTrain[i].X[j] &&
      global()->dataIntTrain[i].X[j] <= workingSol()->b[j] ) {
      if ( j==numAttrib()-1)
      workingSol()->isCovered[i]= true;
      } else {
      workingSol()->isCovered[i]= false;
      break;
      }
    */
    if (global()->args->debug >= 5)
      workingSol()->printSolution();

    return true;

  } // end candidateSolution function


  // ***************** RMASub helper functions (start) ****************

  void RMASub::branchingProcess(const unsigned int &j, const unsigned int &v) {

    vecBounds.resize(3);
    vecBounds[2] = -1;

    if (al[j] <= v && v < min(au[j], bl[j])) {

      // Middle Child
      printSP(j, al[j], v, bl[j], bu[j]);
      mergeEquivClass(j, al[j], v, bl[j], bu[j]);
      vecBounds[0] = getBoundMerge();

      // Up Child
      printSP(j, v + 1, au[j], bl[j], bu[j]);
      dropEquivClass(j, v + 1, bu[j]); // al, bu
      vecBounds[1] = getBoundDrop();

    } else if (au[j] <= bl[j] && au[j] <= v && v < bl[j]) {
      return;
    } else if (bl[j] < au[j] && min(au[j], bl[j]) <= v && v < max(au[j], bl[j])) {

      // Down Child
      printSP(j, al[j], v, v, v);
      dropEquivClass(j, al[j], v);
      vecBounds[0] = getBoundDrop();

      // Middle Child
      printSP(j, al[j], v, v + 1, bu[j]);
      mergeEquivClass(j, al[j], v, v + 1, bu[j]);
      vecBounds[2] = getBoundMerge();

      // Up Child
      printSP(j, v + 1, v + 1, v + 1, bu[j]);
      dropEquivClass(j, v + 1, bu[j]);
      vecBounds[1] = getBoundDrop();
      // vecBounds[1] = vecBounds[2] - vecBounds[0];

      // double tmp = vecBounds[2] - vecBounds[0];
      // if (tmp != vecBounds[1]) {
      //   ucout << "BOUNDS: " << tmp << " " << vecBounds[1] << "\n";
      //   ucout << "BOUNDS: " << vecBounds[2] << " " << vecBounds[0] << "\n";
      // }

    } else if (max(au[j], bl[j]) <= v && v < bu[j]) {

      // Down Child
      printSP(j, al[j], au[j], bl[j], v);
      dropEquivClass(j, al[j], v);
      vecBounds[0] = getBoundDrop();

      // Middle Child
      printSP(j, al[j], au[j], v + 1, bu[j]);
      mergeEquivClass(j, al[j], au[j], v + 1, bu[j]);
      vecBounds[1] = getBoundMerge();
    }

    branchChoice thisChoice(vecBounds[0], vecBounds[1], vecBounds[2], v, j);

    if (global()->args->debug >= 15)
      ucout << "Evaluating" << thisChoice << "\n";

    // select variable based on minimum of children
    // bounds given in lexicographically decreasing order
    thisChoice.sortBounds();

    for (unsigned int i = 0; i < vecBounds.size(); ++i)
      if (thisChoice.branch[i].exactBound <= workingSol()->value) {
        thisChoice.branch[i].exactBound = -1;
        thisChoice.branch[i].roundedBound = -1;
      }

    if (global()->args->debug >= 5)
      ucout << "Sorted version is " << thisChoice << "\n";

    if (thisChoice < _branchChoice) { // and thisChoice.branch[0].roundedBound!=-1
      _branchChoice = thisChoice;
      if (global()->args->debug >= 50)
        ucout << "Improves best attribute: " << j << "\n";
      if (global()->args->debug >= 3)
        ucout << "Branch choice now: " << _branchChoice << "\n";
      numTiedSols = 1;
      // foundBound=true;
    } else if (thisChoice == _branchChoice) {
      // ucout << "branchBound: " << thisChoice.branch[0].exactBound << " "
      //     << _branchChoice.branch[0].exactBound;
      if (global()->args->branchSelection() == 0) {
        numTiedSols++;
        (globalPtr->args->isRandSeed()) ? srand(numTiedSols * time(NULL) * 100)
                                      : srand(1);
        double rand_num = (rand() % 10001) / 10000.0;
        // DEBUGPRX(0, global(), "rand: " << rand_num  << "\n");
        // DEBUGPRX(0, global(), "rand1: " << 1.0 /  numTiedSols << "\n");
        if (rand_num <= 1.0 / numTiedSols) {
          _branchChoice = thisChoice;
          if (global()->args->debug >= 50)
            ucout << "Improves best attribute: " << j << "\n";
          if (global()->args->debug >= 10)
            ucout << "Branch choice now is: " << _branchChoice << "\n";
        }
      } else if (global()->args->branchSelection() == 2) {
        _branchChoice = thisChoice;
        if (global()->args->debug >= 50)
          ucout << "Improves best attribute: " << j << "\n";
        if (global()->args->debug >= 10)
          ucout << "Branch choice now is: " << _branchChoice << "\n";
      }
    }

  } // end function void RMASub::branchingProcess


  // strong branching
  void RMASub::strongBranching() {

    int numCutPtsInAttrib=0;

    if (global()->args->debug >= 5) {
      ucout << "\nal: " << al << "\nau: " << au
           << "\nbl: " << bl << "\nbu: " << bu;
      ucout << "\nsortedObs: " << coveredObs << "\n";
    }

    compIncumbent(numAttrib() - 1);

    for (unsigned int j = 0; j < numAttrib(); ++j) {

      if (global()->args->debug >= 10)
        ucout << "original: ";
      printSP(j, al[j], au[j], bl[j], bu[j]);

      numCutPtsInAttrib += bu[j] - al[j];
      if (bl[j] > au[j])
        numCutPtsInAttrib -= (bl[j] + au[j]);
      if (numCutPtsInAttrib == 0)
        continue;

      for (unsigned int v = al[j]; v < bu[j]; ++v) {
        if (au[j] <= v && v < bl[j]) {
          v = bl[j] - 1;
          continue;
        }
        // ucout << "Size of sortedObs: " << coveredObs1.size() << "\n";
        // ucout << "sortedObs: " << coveredObs;
        branchingProcess(j, v);
      }
      if (j == numAttrib() - 1)
        break;
      global()->args->isCountingSort() ? countingSortEC(j) : bucketSortEC(j);
      compIncumbent(j);
    } // end for each feature

  } // end RMASub::strongBranching


  // branching using cut-point caching methods
  void RMASub::cachedBranching() {

    if (global()->args->debug >= 5)
      ucout << "cachedBranching\n";

    if (global()->args->debug >= 5) {
      ucout << "\nal: " << al << "\nau: " << au
           << "\nbl: " << bl << "\nbu: " << bu;
      ucout << "\nsortedObs: " << coveredObs;
    }

   unsigned int k = 0;

    sortCachedCutPtByAttrib();
    cachedCutPts = sortedCachedCutPts;
    compIncumbent(numAttrib() - 1);

    for (unsigned int j = 0; j < numAttrib(); ++j) {
      while (k < cachedCutPts.size()) {
        if (j == cachedCutPts[k].j) {
          branchingProcess(cachedCutPts[k].j, cachedCutPts[k].v);
          ++k;
        } else
          break;
      }

      if (j == numAttrib() - 1)
        break;
      (global()->args->isCountingSort()) ? countingSortEC(j) : bucketSortEC(j);
      compIncumbent(j);
    }

  } // end RMASub::cachedBranching


  // split subproblems and choose cut value by binary search
  void RMASub::hybridBranching() {

    bool firstFewCutPts;
    unsigned int l, u, L, U, cutValue, numCutValues;
    vector<bool> vecCheckedCutVal;
    unsigned int k = 0;

    sortCachedCutPtByAttrib();
    cachedCutPts = sortedCachedCutPts;
    compIncumbent(numAttrib() - 1);

    for (unsigned int j = 0; j < numAttrib(); ++j) { // for each attribute

      if (vecNumDistVals()[j] < 30) {
        while (k < cachedCutPts.size()) {
          if (j == cachedCutPts[k].j) {
            branchingProcess(cachedCutPts[k].j, cachedCutPts[k].v);
            ++k;
          } else
            break;
        }
        if (j == numAttrib() - 1)
          break;
        (global()->args->isCountingSort()) ? countingSortEC(j) : bucketSortEC(j);
        compIncumbent(j);

      } else { // binary search

        numCutValues = bu[j] - al[j];
        if (bl[j] > au[j])
          numCutValues -= bl[j] - au[j];

        if (numCutValues == 0) { // if no cutValue in this feature,
          (global()->args->isCountingSort()) ? countingSortEC(j) : bucketSortEC(j);
          continue; // then go to the next attribute.
        }

        cutValue = -1;
        l = 0;
        u = 0;
        L = al[j];
        U = bu[j];
        firstFewCutPts = true;
        vecCheckedCutVal.clear();
        vecCheckedCutVal.resize(vecNumDistVals()[j]);

        while (true) {
          if (numCutValues > 3) {
            cutValue = L + (U - L) / 2; // integer division
            if (au[j] <= cutValue && cutValue < bl[j]) {
              cutValue = au[j] - 1 - l;
              l++;
              if (cutValue < al[j] && bl[j] < bu[j]) {
                cutValue = bl[j] + u;
                u++;
                // ucout << "1 (j, cutValue) " << j << ", " << cutValue << "\n";
              } else if (cutValue >= bu[j]) {
                // ucout << "2 (j, cutValue) " << j << ", " << cutValue << "\n";
                break; // no cut point in this feature
              }
            }
            if (cutValue >= bu[j]) {
              if (global()->args->debug >= 10)
                ucout << "cutValue>=bu[j] .\n";
              break;
            }
            if (cutValue < al[j]) {
              if (global()->args->debug >= 10)
                ucout << "cutValue<al[j] .\n";
              break;
            }
            if (vecCheckedCutVal[cutValue]) {
              if (global()->args->debug >= 10)
                ucout << "break since oldCutValue=cutValue.\n";
              break;
            }
            vecCheckedCutVal[cutValue] = true;

            printSP(j, al[j], au[j], bl[j], bu[j]);
            if (global()->args->debug >= 10)
              ucout << "j: " << j << " L: " << L << " U: " << U
                   << " cutVal: " << cutValue << "\n";

          } else {
            if (firstFewCutPts) {
              cutValue = L;
              firstFewCutPts = false;
            } else {
              cutValue++;
            }
            if (au[j] <= cutValue && cutValue < bl[j])
              cutValue = bl[j];
            if (cutValue >= bu[j])
              break;
          }

          // ucout << "(j, cutValue) " << j << ", " << cutValue << "\n";
          if (global()->args->debug >= 10)
            ucout << "coveredObs: " << coveredObs;
          branchingProcess(j, cutValue);

          // compare objectives instead of bounds
          // if ( vecObjValue[0] > vecObjValue[1] ) L = cutValue;
          // else U = cutValue;

          // Compare bounds
          if (vecBounds[0] < vecBounds[1])
            L = cutValue;
          else
            U = cutValue;

        } // end while for each cut value of feature f

        if (j == numAttrib() - 1)
          break;

        (global()->args->isCountingSort()) ? countingSortEC(j) : bucketSortEC(j);
        compIncumbent(j);

      } // end binary search
    }   // end for each attribute

  } // end RMASub::hybridBranching


  // split subproblems and choose cut value by binary search
  void RMASub::binaryBranching() {

    bool firstFewCutPts;
    unsigned int l, u, L, U, cutValue, numCutValues;
    vector<bool> vecCheckedCutVal;

    compIncumbent(numAttrib() - 1);

    for (unsigned int j = 0; j < numAttrib(); ++j) { // for each attribute

      numCutValues = bu[j] - al[j];
      if (bl[j] > au[j])
        numCutValues -= bl[j] - au[j];

      if (numCutValues == 0) { // if no cutValue in this feature,
        (global()->args->isCountingSort()) ? countingSortEC(j) : bucketSortEC(j);
        continue; // then go to the next attribute.
      }

      cutValue = 0; // TODO: check this later, it was -1;
      l = 0;
      u = 0;
      L = al[j];
      U = bu[j];
      firstFewCutPts = true;
      vecCheckedCutVal.clear();
      vecCheckedCutVal.resize(vecNumDistVals()[j]);

      while (true) {
        if (numCutValues > 3) {
          cutValue = L + (U - L) / 2; // integer division
          if (au[j] <= cutValue && cutValue < bl[j]) {
            cutValue = au[j] - 1 - l;
            l++;
            if (cutValue < al[j] && bl[j] < bu[j]) {
              cutValue = bl[j] + u;
              u++;
              // ucout << "1 (j, cutValue) " << j << ", " << cutValue << "\n";
            } else if (cutValue >= bu[j]) {
              // ucout << "2 (j, cutValue) " << j << ", " << cutValue << "\n";
              break; // no cut point in this feature
            }
          }
          if (cutValue >= bu[j]) {
            if (global()->args->debug >= 10)
              ucout << "cutValue>=bu[j] .\n";
            break;
          }
          if (cutValue < al[j]) {
            if (global()->args->debug >= 10)
              ucout << "cutValue<al[j] .\n";
            break;
          }
          if (vecCheckedCutVal[cutValue]) {
            if (global()->args->debug >= 10)
              ucout << "break since oldCutValue=cutValue.\n";
            break;
          }
          vecCheckedCutVal[cutValue] = true;

          printSP(j, al[j], au[j], bl[j], bu[j]);
          if (global()->args->debug >= 10)
            ucout << "j: " << j << " L: " << L << " U: " << U
                 << " cutVal: " << cutValue << "\n";

        } else {
          if (firstFewCutPts) {
            cutValue = L;
            firstFewCutPts = false;
          } else {
            cutValue++;
          }
          if (au[j] <= cutValue && cutValue < bl[j])
            cutValue = bl[j];
          if (cutValue >= bu[j])
            break;
        }

        // ucout << "(j, cutValue) " << j << ", " << cutValue << "\n";
        if (global()->args->debug >= 10)
          ucout << "coveredObs: " << coveredObs;
        branchingProcess(j, cutValue);

        // compare objectives instead of bounds
        // if ( vecObjValue[0] > vecObjValue[1] ) L = cutValue;
        // else U = cutValue;

        // Compare bounds
        if (vecBounds[0] < vecBounds[1])
          L = cutValue;
        else
          U = cutValue;

      } // end while for each cut value of feature f

      if (j == numAttrib() - 1)
        break;

      bucketSortEC(j);
      compIncumbent(j);
    } // end for each attribute

  } // end RMASub::binarySplitSP


  void RMASub::cutpointCaching() {

    // numLiveCachedCutPts = (# of live cut points from the cache)
    int numLiveCachedCutPts = getNumLiveCachedCutPts();

    // if numCachedCutPts is less than the percentage, check all cut points
    if (numLiveCachedCutPts <
        global()->data->numTotalCutPts * global()->args->fracCachedCutPts())
      strongBranching();

    else { // if not, only check the storedCutPts
      // count number of subproblems only discovering cut-points from the chache
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

      cachedBranching();
    }

    int j = _branchChoice.branchVar;

    // store cached cut-points
    if (_branchChoice.branchVar < 0 ||
        _branchChoice.branchVar > globalPtr->data->numAttrib)
      return; // ucout << "ERROR! j is out of range for setCachedCutPts";

    else if (_branchChoice.cutVal < 0 ||
             _branchChoice.cutVal > globalPtr->data->vecNumDistVals[j]-2)
      return; // ucout << "ERROR! v is out of range for setCachedCutPts";

    else
      globalPtr->setCachedCutPts(_branchChoice.branchVar, _branchChoice.cutVal);

  }


  void RMASub::bucketSortObs(const unsigned int &j) {
    unsigned int v;
    int l = -1;
    unsigned int size;
    size = bu[j] - al[j] + 1;
    if (bl[j] > au[j])
      size -= bl[j] - au[j];
    vector<vector<unsigned int>> buckets;
    buckets.resize(size);

    for (unsigned int i = 0; i < coveredObs.size(); ++i) {
      v = global()->data->dataIntTrain[coveredObs[i]].X[j];
      if (au[j] < bl[j]) { // no overlapping
        if (v < au[j])
          v -= al[j];
        else if (au[j] <= v && v <= bl[j])
          v = au[j] - al[j];
        else if (bl[j] < v)
          v = v - (bl[j] - au[j]) - al[j];
      } else
        v -= al[j]; // overlapping

      if (v < 0) {
        if (global()->args->debug >= 50)
          ucout << "below covered range \n";
        continue;
      } else if (v >= size) {
        if (global()->args->debug >= 50)
          ucout << "above covered range \n";
        continue;
      }

      buckets[v].push_back(coveredObs[i]);
    }

    // walk buckets to get sorted observation list on this attribute
    for (unsigned int v = 0; v < size; v++)
      for (unsigned int k = 0; k < buckets[v].size(); k++)
        coveredObs[++l] = buckets[v][k];

    coveredObs.resize(l + 1);
  }


  void RMASub::bucketSortEC(const unsigned int &j) {
    unsigned int v, obs;
    int l = -1;
    unsigned int size = bu[j] - al[j] + 1;
    if (bl[j] > au[j])
      size -= bl[j] - au[j];

    if (size == 1)
      return; // do not have to sort since all values are inseperable
    vector<vector<int>> buckets(size);

    for (unsigned int i = 0; i < sortedECidx.size(); ++i) {
      obs = vecEquivClass[sortedECidx[i]].getObs();
      v = global()->data->dataIntTrain[obs].X[j];
      if (au[j] < bl[j]) { // no overlapping
        if (v < au[j])
          v -= al[j];
        else if (au[j] <= v && v <= bl[j])
          v = au[j] - al[j];
        else if (bl[j] < v)
          v = v - (bl[j] - au[j]) - al[j];
      } else
        v -= al[j]; // overlapping

      if (v < 0) {
        if (global()->args->debug >= 10)
          ucout << "below covered range \n";
        continue;
      } else if (v >= size) {
        if (global()->args->debug >= 10)
          ucout << "above covered range \n";
        continue;
      }

      buckets[v].push_back(sortedECidx[i]);
    }

    // walk buckets to get sorted observation list on this attribute
    for (unsigned int v = 0; v < size; v++)
      for (unsigned k = 0; k < buckets[v].size(); k++)
        sortedECidx[++l] = buckets[v][k];

    // sortedECidx.resize(l+1);
  }

  /*
    void RMASub::countingSortObs(const int& j) {

    int i, v, obs;
    int numObs = coveredObs.size();
    int bucketSize = bu[j] - al[j] + 1 - max(0, bl[j]-au[j]);
    vector<int> vecCount;
    vecCount.resize(bucketSize);
    sortedECidx1.resize(numObs);

    for ( i=0; i < numObs ; ++i ) {
    obs = coveredObs[i];
    v = data->arrayObs[obs].X[j];
    if ( au[j] < bl[j] ) { 	// no overlapping
    if (v<au[j]) v -= al[j];
    else if ( au[j]<=v && v<=bl[j] ) v = au[j] - al[j];
    else if ( bl[j] < v ) v = v - (bl[j] - au[j]) - al[j];
    } else v -= al[j];	// overlapping
    ++vecCount[v] ;
    }

    for ( i=1; i < bucketSize ; ++i )
    vecCount[i] += vecCount[i-1] ;

    for ( i=numObs-1; i>=0 ; --i ) {
    obs = coveredObs[i];
    v = data->arrayObs[obs].X[j];
    if ( au[j] < bl[j] ) { 	// no overlapping
    if (v<au[j]) v -= al[j];
    else if ( au[j]<=v && v<=bl[j] ) v = au[j] - al[j];
    else if ( bl[j] < v ) v = v - (bl[j] - au[j]) - al[j];
    } else v -= al[j];	// overlapping
    coveredObs1[ vecCount[v]-1 ] = coveredObs[i] ;
    --vecCount[v];
    }

    coveredObs = coveredObs1;

    }
  */


  void RMASub::countingSortEC(const unsigned int &j) {

    unsigned int i, v, obs;
    unsigned int numEC = sortedECidx.size();
    unsigned int bucketSize = bu[j] - al[j] + 1;
    if (bl[j] > al[j])
      bucketSize -= bl[j] - au[j];
    vector<unsigned int> vecCount(bucketSize);
    sortedECidx1.resize(numEC);

    for (i = 0; i < numEC; ++i) {
      obs = vecEquivClass[sortedECidx[i]].getObs();
      v = global()->data->dataIntTrain[obs].X[j];
      if (au[j] < bl[j]) { // no overlapping
        if (v < au[j])
          v -= al[j];
        else if (au[j] <= v && v <= bl[j])
          v = au[j] - al[j];
        else if (bl[j] < v)
          v = v - (bl[j] - au[j]) - al[j];
      } else
        v -= al[j]; // overlapping
      ++vecCount[v];
    }

    for (i = 1; i < bucketSize; ++i)
      vecCount[i] += vecCount[i - 1];

    for (i = numEC - 1; i >= 0; --i) {
      obs = vecEquivClass[sortedECidx[i]].getObs();
      v = global()->data->dataIntTrain[obs].X[j];
      if (au[j] < bl[j]) { // no overlapping
        if (v < au[j])
          v -= al[j];
        else if (au[j] <= v && v <= bl[j])
          v = au[j] - al[j];
        else if (bl[j] < v)
          v = v - (bl[j] - au[j]) - al[j];
      } else
        v -= al[j]; // overlapping
      sortedECidx1[vecCount[v] - 1] = sortedECidx[i];
      --vecCount[v];
    }

    sortedECidx = sortedECidx1;
  }


  void RMASub::compIncumbent(const unsigned int &j) {

    tmpMin = getInf();
    tmpMax = -getInf();
    // minVal = globalPtr->args->isInitGuess() ?   workingSol()->value : inf;
    // maxVal = globalPtr->args->isInitGuess() ?  -workingSol()->value : -inf;
    minVal = getInf();
    maxVal = -getInf();
    optMinAttrib = -1;
    optMaxAttrib = -1;

    curObs = 0;
    tmpMin = runMinKadane(j);
    if (tmpMin == minVal) {
      numNegTiedSols++;
      (globalPtr->args->isRandSeed()) ? srand(numNegTiedSols * time(NULL) * 100)
                                    : srand(1);
      rand_num = (rand() % 10001) / 10000.0;
      if (rand_num <= 1.0 / numNegTiedSols)
        setOptMin(j);
    } else if (tmpMin < minVal) { // if better min incumbent was found
      numNegTiedSols = 1;
      setOptMin(j);
    }

    curObs = 0;
    tmpMax = runMaxKadane(j);
    if (tmpMax == maxVal) {
      numPosTiedSols++;
      (globalPtr->args->isRandSeed()) ? srand(numNegTiedSols * time(NULL) * 100)
                                    : srand(1);
      rand_num = (rand() % 10001) / 10000.0;
      if (rand_num <= 1.0 / numPosTiedSols)
        setOptMax(j);
    } else if (tmpMax > maxVal) {
      numPosTiedSols = 1;
      setOptMax(j);
    }

    chooseMinOrMaxRange();
  }


  void RMASub::chooseMinOrMaxRange() {

    // int numFoundNewSols = 0;

    if (max(maxVal, -minVal) > workingSol()->value + .000001) {

      // numFoundNewSols = 1;

      (globalPtr->args->isRandSeed())
          ? srand((numNegTiedSols + numPosTiedSols) * time(NULL) * 100)
          : srand(1);

      rand_num = (rand() % 10001) / 10000.0;

      workingSol()->a << al;
      workingSol()->b << bu;

      // if max ver is better than min ver
      // or breaking the tied solution, choose max ver
      if (maxVal > -minVal ||
          (maxVal == minVal &&
           rand_num <=
               numPosTiedSols / (double)(numNegTiedSols + numPosTiedSols))) {

        workingSol()->value = maxVal;
        workingSol()->a[optMaxAttrib] = optMaxLower;
        workingSol()->b[optMaxAttrib] = optMaxUpper;
        workingSol()->isPosIncumb = true;

        if (globalPtr->args->debug >= 10)
          ucout << "positive ";

      } else { // else choose the min ver

        workingSol()->value = -minVal;
        workingSol()->a[optMinAttrib] = optMinLower;
        workingSol()->b[optMinAttrib] = optMinUpper;
        workingSol()->isPosIncumb = false;

        if (globalPtr->args->debug >= 10)
          ucout << "negative ";

      } // end if choosing pos or neg ver.

      if (globalPtr->args->debug >= 1)
        ucout << " new incumbent  " << std::fixed << std::setprecision(4)
             << workingSol()->value << '\n';

      /****** temporarily solution until fixing PEBBL solution *****/

      // globalPtr->globalSol.setSolution(workingSol()->a, workingSol()->b,
      //               workingSol()->isPosIncumb, workingSol()->value);

      foundRMASolution(synchronous);

      if (globalPtr->args->debug >= 2)
        workingSol()->printSolution();

      DEBUGPR(10, workingSol()->checkObjValue1(workingSol()->a, workingSol()->b,
                  coveredObs, sortedECidx ));

    }

  } // end RMASub::chooseMinOrMaxRange function


  void RMASub::setOptMin(const unsigned int &j) {

    minVal = tmpMin;

    optMinAttrib = j;
    optMinLower = aj;
    optMinUpper = bj;
    if (au[j] < bl[j] && au[j] <= optMinUpper && optMinUpper <= bl[j])
      optMinUpper = bl[j];

    if (globalPtr->args->debug >= 5)
      ucout << "optAttrib: (a,b): " << optMinAttrib << ": (" << optMinLower << ", "
           << optMinUpper << "), min: " << minVal << "\n";
  }


  void RMASub::setOptMax(const unsigned int &j) {

    maxVal = tmpMax;

    optMaxAttrib = j;
    optMaxLower = aj;
    optMaxUpper = bj;
    if (au[j] < bl[j] && au[j] <= optMaxUpper && optMaxUpper <= bl[j])
      optMaxUpper = bl[j];

    if (globalPtr->args->debug >= 5)
      ucout << "optAttrib: (a,b): " << optMaxAttrib << ": (" << optMaxLower << ", "
           << optMaxUpper << "), max: " << maxVal << "\n";
  }


  // get Maximum range for the feature
  double RMASub::runMaxKadane(const unsigned int &j) {

    double maxEndHere, maxSoFar; //, tmpObj;
    unsigned int s = al[j];
    aj = al[j];
    bj = al[j];
    // bj=bu[j]; // al[j];
    maxEndHere = 0;
    maxSoFar = -getInf();

    for (unsigned int v = al[j]; v <= bu[j]; ++v) {

      if (au[j] < bl[j] && au[j] < v && v <= bl[j]) {
        v = bl[j];
        continue;
      }

      maxEndHere += getObjValue(j, v);

      if (maxEndHere > maxSoFar) {
        maxSoFar = maxEndHere;
        aj = s;
        bj = v;
      }

      if (maxEndHere < 0) {
        maxEndHere = 0;
        s = v + 1;
        if (au[j] < bl[j] && v == au[j])
          s = bl[j];
      }

    } // end for each value in this attribute

    if (globalPtr->args->debug >= 10)
      ucout << "Maximum contiguous sum is " << maxSoFar
           << " attribute (L,U): " << j << " (" << aj << ", " << bj << ")\n";

    return maxSoFar;
  }


  // get Miniumum range for the feature
  double RMASub::runMinKadane(const unsigned int &j) {

    double minEndHere, minSoFar; //, tmpObj;
    unsigned int s = al[j];
    aj = al[j];
    bj = al[j];
    // bj=bu[j]; // al[j];
    // ucout << "aj: " << aj << ", bj: " << bj <<"\n";
    minEndHere = 0;
    minSoFar = getInf();

    for (unsigned int v = al[j]; v <= bu[j]; ++v) {

      if (au[j] < bl[j] && au[j] < v && v <= bl[j]) {
        v = bl[j];
        continue;
      }

      minEndHere += getObjValue(j, v);

      if (minEndHere < minSoFar) {
        minSoFar = minEndHere;
        aj = s;
        bj = v;
      }

      if (minEndHere > 0) {
        minEndHere = 0;
        s = v + 1;
        if (au[j] < bl[j] && v == au[j])
          s = bl[j];
      }
    }

    if (globalPtr->args->debug >= 10)
      ucout << "Minimum contiguous sum is " << minSoFar
           << " attribute (L,U): " << j << " (" << aj << ", " << bj << ")\n";
    return minSoFar;
  }


  double RMASub::getObjValue(const unsigned int &j, const unsigned int &v) {
    unsigned int obs;
    double covgWt = 0.0;
    if (globalPtr->args->debug >= 20)
      ucout << "j: " << j << ", v: " << v;

    for (unsigned int i = curObs; i < sortedECidx.size(); ++i) {

      obs = vecEquivClass[sortedECidx[i]].getObs();

      // if the observation's jth attribute value = cut-value
      if (global()->data->dataIntTrain[obs].X[j] == v) {
        covgWt += vecEquivClass[sortedECidx[i]].getWt();
        // ucout << "vecEquivClass1: " << sortedECidx[i]
        //     << " covgWt: " << covgWt << endl;
      } else if (au[j] < bl[j] && au[j] <= v && v <= bl[j] &&
                 au[j] <= global()->data->dataIntTrain[obs].X[j] &&
                 global()->data->dataIntTrain[obs].X[j] <= bl[j]) {
        covgWt += vecEquivClass[sortedECidx[i]].getWt();
        // ucout << "vecEquivClass2: " << sortedECidx[i]
        //     << " covgWt: " << covgWt << endl;

      } else if (global()->data->dataIntTrain[obs].X[j] < v) {
        if (globalPtr->args->debug >= 0)
          ucout << "X[j] < v! ";
        //*
        if (globalPtr->args->debug >= 20) {
          ucout << "curObs: " << curObs << " attribute: " << j << "; "
               << global()->data->dataIntTrain[obs].X[j] << " < cutVal: " << v
               << "\n";
          for (unsigned int i = 0; i < coveredObs.size(); ++i)
            ucout << global()->data->dataIntTrain[sortedECidx[i]].X[j]
                 << " bound (" << al[j] << ", " << au[j] << ", " << bl[j] << ", "
                 << bu[j] << ")\n )";
        }
      } else {
        curObs = i;
        break;
      }
      if (i == sortedECidx.size() - 1)
        curObs = sortedECidx.size();
    } // end for each covered observation

    if (globalPtr->args->debug >= 20)
      ucout << ", covgWt: " << covgWt << "\n";
    return covgWt;

  } // end function getObjValue


  double RMASub::getBoundMerge() const {

    // unsigned int obs, idxEC; // observation number
    double pBound = 0.0,
           nBound = 0.0; // weight for positive and negative observation

    for (unsigned int i = 0; i < vecEquivClass1.size();
         ++i) {                          // for each equivalence class
      if (vecEquivClass1[i].getWt() > 0) // for positive observation
        pBound += vecEquivClass1[i].getWt();
      else if (vecEquivClass1[i].getWt() < 0) // for negative observation
        nBound -= vecEquivClass1[i].getWt();
    } // end outer for

    return pBound > nBound ? pBound : nBound;
  } // end functino RMASub::getBound


  double RMASub::getBoundDrop() const {

    //unsigned int obs, idxEC; // observation number
    unsigned int idxEC; // observation number
    double pBound = 0.0, nBound = 0.0;
    // weight for positive and negative observation<

    for (unsigned int i = 0; i < sortedECidx1.size();
         ++i) { // for each equivalence class
      idxEC = sortedECidx1[i];
      if (vecEquivClass[idxEC].getWt() > 0) // for positive observation
        pBound += vecEquivClass[idxEC].getWt();
      else if (vecEquivClass[idxEC].getWt() < 0) // for negative observation
        nBound -= vecEquivClass[idxEC].getWt();
    } // end outer for

    return pBound > nBound ? pBound : nBound;
  } // end functino RMASub::getBound


  void RMASub::setInitialEquivClass() {

    if (coveredObs.size() <= 0) {
      if (globalPtr->args->debug >= 0)
        ucout << "coveredObs is empty.\n";
      return;
    } else if (coveredObs.size() == 1) {
      if (globalPtr->args->debug >= 15)
        ucout << "There is only one covered observation"
             << "\n";
      vecEquivClass.resize(1);
      vecEquivClass[0].addObsWt(coveredObs[0],
                                global()->data->dataIntTrain[coveredObs[0]].w);
      return;
    }

    vecEquivClass.resize(coveredObs.size());

    int obs1 = coveredObs[0];
    int obs2 = coveredObs[1];
    int k = 0;
    vecEquivClass[0].addObsWt(obs1, global()->data->dataIntTrain[obs1].w);

    for (unsigned int i = 1; i < coveredObs.size();
         ++i) { // for each sorted, covered observation

      for (unsigned int j = 0; j < numAttrib(); ++j) { // for each attribute

        if (isInSameClass(obs1, obs2, j, au[j], bl[j])) {
          if (j == numAttrib() - 1) { // if it is in the same equivalent class
            vecEquivClass[k].addObsWt(obs2, global()->data->dataIntTrain[obs2].w);
            if (i != coveredObs.size() - 1) // if not the last observation
              obs2 = coveredObs[i + 1];
          }

        } else { // detected obs1 and obs2 are in different equivClass
          vecEquivClass[++k].addObsWt(obs2, global()->data->dataIntTrain[obs2].w);
          if (i != coveredObs.size() - 1) { // if not the last observation
            obs1 = coveredObs[i];
            obs2 = coveredObs[i + 1];
          }
          break; // as soon as we detect obs1 and obs2 are in different equivClass
                 // compare the next observation combinations
        }

      } // end for each attribute j
    }   // end for each obs, i


    // erase extra space
    vecEquivClass.erase(vecEquivClass.begin() + k + 1, vecEquivClass.end());

    sortedECidx.resize(vecEquivClass.size());
    for (unsigned int i = 0; i < vecEquivClass.size(); ++i)
      sortedECidx[i] = i;

    if (globalPtr->args->debug >= 10)
      ucout << "Size of coveredObs: " << sortedECidx.size() << "\n";
    if (globalPtr->args->debug >= 20)
      ucout << "Size of vecEquivClass: " << vecEquivClass.size() << "\n";
    if (globalPtr->args->debug >= 30)
      for (unsigned int i = 0; i < vecEquivClass.size(); ++i)
        ucout << "EC: " << i << ": " << vecEquivClass[i] << "\n";
  }


  void RMASub::mergeEquivClass(const unsigned int &j, const unsigned int &al_,
                               const unsigned int &au_, const unsigned int &bl_,
                               const unsigned int &bu_) {

    if (al_ != al[j] || bu_ != bu[j])
      dropEquivClass(j, al_, bu_);
    else {
      sortedECidx1.resize(sortedECidx.size());
      copy(sortedECidx.begin(), sortedECidx.end(), sortedECidx1.begin());
    }

    vecEquivClass1.clear();

    if (sortedECidx1.size() <= 0) {
      if (globalPtr->args->debug >= 0)
        ucout << "sortedECidx1 is empty. \n";
      return;
    }

    if (sortedECidx1.size() == 1) {
      if (globalPtr->args->debug >= 15)
        ucout << "There is only one equivalence class"
             << "\n";
      vecEquivClass1.resize(1);
      vecEquivClass1[0] = vecEquivClass[sortedECidx1[0]];
      return;
    }

    unsigned int idxEC1 = sortedECidx1[0];
    unsigned int idxEC2 = sortedECidx1[1];

    unsigned int obs1 = vecEquivClass[idxEC1].getObs();
    unsigned int obs2 = vecEquivClass[idxEC2].getObs();
    unsigned int k = 0, J, _au, _bl;

    vecEquivClass1.resize(sortedECidx1.size());
    vecEquivClass1[0] = vecEquivClass[idxEC1];

    // for each sorted, covered observation
    for (unsigned int i = 1; i < sortedECidx1.size(); ++i) {
      for (unsigned int f = 0; f < numAttrib(); ++f) { // for each attribute

        // always search merging equivalence classes from leaves
        J = f + j;
        if (J >= numAttrib())
          J -= numAttrib(); // if J=numAttribute(), go back to J=0

        // if J is the modified attribute in order to check these two observations
        // are in the same equivalence class or not from the leaves of equivalence
        // tree
        if (J == j) {
          _au = au_;
          _bl = bl_;
        } // DEBUGPR(50,ucout << j << _au << _bl << endl);
        else {
          _au = au[J];
          _bl = bl[J];
        }

        // if the two observations are in the same equivalence class
        if (isInSameClass(obs1, obs2, J, _au, _bl)) {
          if (f == numAttrib() - 1) { // if it is in the same equivalent class
            vecEquivClass1[k].addEC(vecEquivClass[idxEC2]);
            if (i != sortedECidx1.size() - 1) // if not the last equivClass
              idxEC2 = sortedECidx1[i + 1];
            obs2 = vecEquivClass[idxEC2].getObs();
            break;
          }
        } else { // detected obs1 and obs2 are in different equivClass
          vecEquivClass1[++k] =
              vecEquivClass[idxEC2]; // push back all obs in the equivClass
          if (i != sortedECidx1.size() - 1) { // if not the last observation
            idxEC1 = sortedECidx1[i];
            idxEC2 = sortedECidx1[i + 1];
            obs1 = vecEquivClass[idxEC1].getObs();
            obs2 = vecEquivClass[idxEC2].getObs();
          }
          break;
        }
      } // end for each attribute j
    }   // end for each obs

    vecEquivClass1.erase(vecEquivClass1.begin() + k + 1, vecEquivClass1.end());
    sortedECidx1.resize(vecEquivClass1.size());
    for (unsigned int i = 0; i < vecEquivClass1.size(); ++i)
      sortedECidx1[i] = i;

    if (globalPtr->args->debug >= 20)
      ucout << "Size of vecEquivClass1: " << vecEquivClass1.size() << "\n";
    if (globalPtr->args->debug >= 25)
      ucout << "vecEquivClass1: \n";
    if (globalPtr->args->debug >= 30)
      for (unsigned int i = 0; i < vecEquivClass1.size(); ++i)
        ucout << "EC: " << i << ": " << vecEquivClass1[i] << "\n";

  } // end function RMASub::mergeEquivClass


  // drop some equivalence class from the initial equivalence class
  void RMASub::dropEquivClass(const unsigned int &j, const unsigned int &al_,
                              const unsigned int &bu_) {

    unsigned int obs, idxEC;
    int k = -1; // observation number
    sortedECidx1.resize(sortedECidx.size());

    for (unsigned int i = 0; i < sortedECidx.size();
         ++i) { // for each equivalence class
      idxEC = sortedECidx[i];
      obs = vecEquivClass[idxEC].getObs();
      // if covered, put the equiv class index to sortedECidx1
      if (global()->data->dataIntTrain[obs].X[j] >= al_ &&
          global()->data->dataIntTrain[obs].X[j] <= bu_)
        sortedECidx1[++k] = idxEC;
    } // end each equivalence class

    sortedECidx1.resize(k + 1); // erase extra space

  } // end function RMASub::dropEquivClass


  bool RMASub::isInSameClass(const unsigned int &obs1, const unsigned int &obs2,
                             const unsigned int &j, const unsigned int &au_,
                             const unsigned int &bl_) {

    // if obs1.feat == obs2.feat
    if (global()->data->dataIntTrain[obs1].X[j] ==
        global()->data->dataIntTrain[obs2].X[j])
      return true;

    //  OR obs1.feat, obs2.feat in [au, bl]
    if (au_ <= bl_ &&
        (au_ <= global()->data->dataIntTrain[obs1].X[j] &&
         global()->data->dataIntTrain[obs1].X[j] <= bl_) &&
        (au_ <= global()->data->dataIntTrain[obs2].X[j] &&
         global()->data->dataIntTrain[obs2].X[j] <= bl_)) {
      return true;
    }
    return false;
  }


  void RMASub::sortCachedCutPtByAttrib() {

    int l = -1;
    vector<vector<unsigned int>> buckets;
    buckets.resize(numAttrib());
    sortedCachedCutPts.resize(cachedCutPts.size());

    for (unsigned int k = 0; k < cachedCutPts.size(); ++k)
      buckets[cachedCutPts[k].j].push_back(k);

    // walk buckets to get sorted observation list on this attribute
    for (unsigned int i = 0; i < numAttrib(); ++i)
      for (unsigned int j = 0; j < buckets[i].size(); ++j)
        sortedCachedCutPts[++l] = cachedCutPts[buckets[i][j]];
  }


  void RMASub::printSP(const unsigned int &j, const unsigned int &al,
                       const unsigned int &au, const unsigned int &bl,
                       const unsigned int &bu) const {
    if (globalPtr->args->debug >= 10)
      ucout << "j: " << j << " (al, au, bl, bu) = (" << al << ", " << au << ", "
           << bl << ", " << bu << ")\n";
  }


  void RMASub::printCurrentBounds() {
    if (globalPtr->args->debug >= 10)
      ucout << "Best local choice is " << _branchChoice << "\n";
    if (globalPtr->args->debug >= 10)
      ucout << " optFeat=" << _branchChoice.branchVar
           << " optCutValue=" << _branchChoice.cutVal
           << " minBound=" << _branchChoice.branch[0].exactBound << endl;
  } // end junction RMASub::printCurrentBounds


  // TODO explain what was this!!!
  void RMASub::setCutPts() {

    static unsigned int count = 0;

    excCutFeat.push_back(_branchChoice.branchVar);
    excCutVal.push_back(_branchChoice.cutVal);

    unsigned int numBranch = global()->CutPtOrders.size(); // number of branches
    if (excCutFeat.size() ==
        1) {                  // if excludedList is 1, always create a new branch
      vector<CutPtOrder> row; // Create an empty row
      global()->CutPtOrders.push_back(row); // Add the row to the vector
      global()->CutPtOrders[numBranch].push_back(
          CutPtOrder(count, _branchChoice.branchVar, _branchChoice.cutVal));
    } else {
      bool thisBranch = false;
      for (unsigned int i = 0; i < numBranch; ++i) {
        unsigned int sizeBranch = global()->CutPtOrders[i].size();
        if ((sizeBranch + 1) == excCutFeat.size()) {
          for (unsigned int k = 0; k < sizeBranch; ++k) {
            if (global()->CutPtOrders[i][k].j == excCutFeat[k] &&
                global()->CutPtOrders[i][k].v == excCutVal[k]) {
              if (k == sizeBranch - 1)
                thisBranch = true;
            } else
              break;
          }
          if (thisBranch) {
            global()->CutPtOrders[i].push_back(
                CutPtOrder(count, _branchChoice.branchVar, _branchChoice.cutVal));
            break;
          }
        }
      }
      if (!thisBranch) {
        unsigned int sizeThisBranch = excCutFeat.size();
        vector<CutPtOrder> row;               // Create an empty row
        global()->CutPtOrders.push_back(row); // Add the row to the vector
        global()->CutPtOrders[numBranch].resize(sizeThisBranch);

        for (unsigned int i = 0; i < numBranch; ++i) {
          unsigned int sizeBranch = global()->CutPtOrders[i].size();
          if (sizeThisBranch <= sizeBranch)
            for (unsigned int k = 0; k < sizeThisBranch - 1; ++k) {
              if (global()->CutPtOrders[i][k].j == excCutFeat[k] &&
                  global()->CutPtOrders[i][k].v == excCutVal[k]) {
                if (k == sizeThisBranch - 2)
                  thisBranch = true;
              } else {
                break;
              }
            }
          if (thisBranch) {
            for (unsigned int k = 0; k < sizeThisBranch - 1; ++k)
              global()->CutPtOrders[numBranch][k].setCutPt(
                  global()->CutPtOrders[i][k]);
            break;
          }
        }
        global()->CutPtOrders[numBranch][sizeThisBranch - 1].setCutPt(
            CutPtOrder(count, _branchChoice.branchVar, _branchChoice.cutVal));
      }
    }

    count++;
  }


  // ******************************************************************************
  // RMASolution methods

  rmaSolution::rmaSolution(RMA *global_) : solution(global_), global(global_) {
    DEBUGPRX(100, global,
             "Creating rmaSolution at " << (void *)this
              << " with global=" << global << "\n";);

    // Only one solution representation in this application,
    // so typeId can just be 0.
    typeId = 0;
  }


  rmaSolution::rmaSolution(rmaSolution *toCopy) {
    // DEBUGPRX(100, global,  "Copy constructing rmaSolution at "
    //   << (void*) this << " from " << toCopy << endl; );
    copy(toCopy);
    serial = ++(global->solSerialCounter);
  }


  void rmaSolution::copy(rmaSolution *toCopy) {
    solution::copy(toCopy);
    global = toCopy->global;
    a << toCopy->a;
    b << toCopy->b;
    isPosIncumb = toCopy->isPosIncumb;
    if (global->args->debug>=5) {
      ucout << "\ncopy a: " << a << "\n";
      ucout << "copy b: " << b << "\n";
    }
  }


  void rmaSolution::printContents(ostream &os) {

    if (global->enumCount > 1)
      return;

    os << "rectangle: a: " << a << "rectangle: b: " << b ;
    ucout << "rectangle: a: " << a << "rectangle: b: " << b ;

    for (unsigned int i = 0; i < global->data->numAttrib; ++i) {
      if (0 < a[i]) // if lower bound changed
        ucout << a[i] << "<=";
      if (0 < a[i] || b[i] < global->data->vecNumDistVals[i]-1)
        ucout << "x" << i;
      if (b[i] < global->data->vecNumDistVals[i]-1)
        ucout << "<=" << b[i];
      if (0 < a[i] || b[i] < global->data->vecNumDistVals[i]-1)
        ucout << ", ";
    }
    ucout << "\n";

    if (global->args->isCheckObjVal())
      checkObjValue();

    if (global->args->isSaveCutPts()) {
      os << "CutPts:\n";
      for (unsigned int i = 0; i < global->CutPtOrders.size(); ++i) {
        unsigned int sizeBranch = global->CutPtOrders[i].size();
        for (unsigned int j = 0; j < sizeBranch - 1; ++j)
          os << global->CutPtOrders[i][j].order << "-"
                    << global->CutPtOrders[i][j].j << "-"
                    << global->CutPtOrders[i][j].v << "; ";
        os << global->CutPtOrders[i][sizeBranch - 1].order << "-"
                  << global->CutPtOrders[i][sizeBranch - 1].j << "-"
                  << global->CutPtOrders[i][sizeBranch - 1].v << "\n";
      } // end for each cut point
    }   // end if writeingCutPts option

  } // end function printContents


  void const rmaSolution::printSolution() {
    ucout << "printSolution: ";
    ucout << ((isPosIncumb) ? "Positive" : "Negative");
    ucout << "\na: " << a << "b: " << b ;
  }


  void rmaSolution::checkObjValue() {

    double wt = 0.0;

    ucout << "Check RMA solution:";
    ucout << "\na: " << a << "\nb: " << b << "\n";

    // for each observation
    for (unsigned int i = 0; i < global->data->numTrainObs; ++i) {

      // for each attribute
      for (unsigned int j = 0; j < global->data->numAttrib; ++j) {

        if (a[j] <= global->data->dataIntTrain[i].X[j] &&
            global->data->dataIntTrain[i].X[j] <= b[j]) {

          // if this observation is covered by this solution
          if (j == global->data->numAttrib - 1)
            wt += global->data->dataIntTrain[i].w;

        } else
          break; // else go to the next observation
      }          // end for each attribute
    }            // end for each observation

    ucout << "Check RMA ObjValue=" << wt << "\n";

  } // end function rmaSolution::checkObjValue


  void rmaSolution::checkObjValue1(vector<unsigned int> &A, vector<unsigned int> &B,
                                   vector<unsigned int> &coveredObs,
                                   vector<unsigned int> &sortedECidx) {

    unsigned int obs;
    double wt = 0.0;

    for (unsigned int i = 0; i < coveredObs.size(); ++i) { // for each observation
      obs = coveredObs[i];
      for (unsigned int j = 0; j < global->data->numAttrib; ++j) { // for each attribute
        if (A[j] <= global->data->dataIntTrain[obs].X[j] &&
            global->data->dataIntTrain[obs].X[j] <= B[j]) {
          // if this observation is covered by this solution
          if (j == global->data->numAttrib - 1)
            wt += global->data->dataIntTrain[obs].w;
        } else
          break; // else go to the next observation
      }          // end for each attribute
    }            // end for each observation

    ucout << "A: " << A;
    ucout << "B: " << B;
    ucout << "RMA ObjValue=" << wt << "\n";

    if (abs(value - abs(wt)) > .00001) {
      ucout << "RMA Obj Not Match! " << wt << " " << value << "\n";
      ucout << "check coveredObs: " << coveredObs;
      ucout << "check sortedECidx: " << sortedECidx;
    }

  } // end function rmaSolution::checkObjValue


  #ifdef ACRO_HAVE_MPI

  void rmaSolution::packContents(PackBuffer &outBuf) const {
    outBuf << a << b << isPosIncumb;
  }

  void rmaSolution::unpackContents(UnPackBuffer &inBuf) {
    inBuf >> a >> b >> isPosIncumb;
  }

  int rmaSolution::maxContentsBufSize() {
    return 2 * (global->data->numAttrib + 1) * sizeof(int) * 1.5 + sizeof(bool);
  }

  #endif // ACRO_HAVE_MPI


  double rmaSolution::sequenceData() {
    if (sequenceCursor < a.size())
      return a[sequenceCursor++];
    else if (sequenceCursor < a.size() + b.size())
      return b[sequenceCursor++ - a.size()];
    else
      return isPosIncumb;
  } // end sequenceData function

} // namespace pebblRMA


ostream &operator<<(ostream &os, pebblRMA::branchChoice &bc) {
  os << '(' << bc.branch[0].exactBound << ',' << bc.branch[1].exactBound << ','
     << bc.branch[2].exactBound << ")-(" << bc.branch[0].whichChild << ','
     << bc.branch[1].whichChild << ',' << bc.branch[2].whichChild << ")-<"
     << bc.branchVar << ">-(" << bc.cutVal << ')';
  return os;
}


ostream &operator<<(ostream &os, pebblRMA::EquivClass &obj) {
  obj.write(os);
  return os;
}
