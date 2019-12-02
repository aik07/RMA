/**********************************************************
* File name:   serRMA.cpp
* Author:      Ai Kagawa
* Description: a source file for serial RMA solver using PEBBL
***********************************************************/

#include "serRMA.h"

namespace pebblRMA {

//////////////////// MPI methods ////////////////////
//#define MPI_UP1 OMPI_PREDEFINED_GLOBAL(MPI_Datatype,	MPI_TYPE_CREATE_RESIZED)

#ifdef ACRO_HAVE_MPI

  void branchChoiceCombiner(void* invec, void* inoutvec,
                            int* len, MPI_Datatype* datatype) {

#ifdef ACRO_VALIDATING
    if (*datatype != branchChoice::mpiType) {
      cerr << "Datatype error in branchChoiceCombiner\n";
      exit(1);
    }
#endif

    branchChoice* inPtr    = (branchChoice*) invec;
    branchChoice* inOutPtr = (branchChoice*) inoutvec;
    int n = *len;
    for (int i=0; i<n; i++)
      if (inPtr[i] < inOutPtr[i])
	      inOutPtr[i] = inPtr[i];

  }


  void branchChoiceRand(branchChoice *in, branchChoice *inout,
                        int *len, MPI_Datatype *datatype) {

#ifdef ACRO_VALIDATING
    if (*datatype != branchChoice::mpiType) {
      cerr << "Datatype error in branchChoiceRand\n";
      exit(1);
    }
#endif

    branchChoice c;

    for (int i=0; i< *len; ++i) {
      if ( in < inout )      c = *in ;
      else if ( in > inout ) c = *inout ;
      else {
        int n1 = in->numTiedSols;
        int n2 = inout->numTiedSols;

        srand((n1+n2)*time(NULL)*100);// : srand(1);
        double rand_num = ( rand() % (n1+n2+1) ) / (double)(n1+n2) ;
        //double rand_num = ( rand() % 1001 ) / (double) 1000 ;
        if ( rand_num < (double) n1 / (n1+n2) ) {
          //cout << "rand_num: " << rand_num << endl;
          c = *in ;
        } else {
          c = *inout ;
        }
        c.numTiedSols = n1+n2;
      }
      *inout = c;
      in++;
      inout++;
    }

  }


//////////////////// branchChoice class methods ////////////////////

  void branchChoice::setupMPIDatum(void* address, MPI_Datatype thisType,
				   MPI_Datatype* type, MPI_Aint base, MPI_Aint* disp,
				   int* blocklen, int j) {
    MPI_Get_address( address, &disp[j] );
    disp[j] -= base;
    type[j] = thisType;
    blocklen[j] = 1;
  }


  void branchChoice::setupMPI() {

    int arraySize = 3*3 + 3;
    int j = 0;
    MPI_Datatype type[arraySize];
    int          blocklen[arraySize];
    MPI_Aint     disp[arraySize];
    MPI_Aint     base;
    branchChoice example;
    MPI_Get_address(&example, &base);

    for (int i=0; i<3; i++) {
      setupMPIDatum(&(example.branch[i].roundedBound), MPI_DOUBLE,
		    type, base, disp, blocklen, j++);
      setupMPIDatum(&(example.branch[i].exactBound),   MPI_DOUBLE,
		    type, base, disp, blocklen, j++);
      setupMPIDatum(&(example.branch[i].whichChild),   MPI_INT,
		    type, base, disp, blocklen, j++);
    }

    setupMPIDatum(&(example.branchVar),   MPI_INT,
		  type, base, disp, blocklen, j++);
    setupMPIDatum(&(example.cutVal),      MPI_INT,
		  type, base, disp, blocklen, j++);
    setupMPIDatum(&(example.numTiedSols), MPI_INT,
		  type, base, disp, blocklen, j++);

    MPI_Type_create_struct(j, blocklen, disp, type, &mpiType);
    MPI_Type_commit(&mpiType);

    MPI_Op_create( branchChoiceCombiner, true, &mpiCombiner );
    MPI_Op_create( (MPI_User_function *) branchChoiceRand,
                   true, &mpiBranchSelection );

  }


  void branchChoice::freeMPI() {
    MPI_Op_free(&mpiCombiner);
    MPI_Op_free(&mpiBranchSelection);
    MPI_Type_free(&mpiType);
  };


  MPI_Datatype branchChoice::mpiType            = MPI_UB;
  MPI_Op       branchChoice::mpiCombiner        = MPI_OP_NULL;
  MPI_Op       branchChoice::mpiBranchSelection = MPI_OP_NULL;

#endif


  branchChoice::branchChoice() : branchVar(MAXINT), cutVal(-1) {
    for (int i=0; i<3; i++)
      branch[i].set(MAXDOUBLE, 1e-5);
  }

  branchChoice::branchChoice(double a, double b,
                             double c, int cut, int j) {
    branch[0].set(a, 1e-5);
    branch[1].set(b, 1e-5);
    branch[2].set(c, 1e-5);

    for (int i=0; i<3; ++i)
      branch[i].whichChild=i;

    cutVal = cut;
    branchVar = j;
  }


  void branchChoice::setBounds(double a, double b, double c, int cut, int j) {

    branch[0].set(a, 1e-5);
    branch[1].set(b, 1e-5);
    branch[2].set(c, 1e-5);

  	for (int i=0; i<3; ++i)
      branch[i].whichChild=i;

  	cutVal = cut;	branchVar = j;
  }


  // Primitive sort, but only three elements
  void branchChoice::sortBounds() {
    possibleSwap(0,1);
    possibleSwap(0,2);
    possibleSwap(1,2);
  }


  bool branchChoice::operator<(const branchChoice& other) const {
    for(int i=0; i<3; i++) {
      if (branch[i].roundedBound < other.branch[i].roundedBound)
	      return true;
      else if (branch[i].roundedBound > other.branch[i].roundedBound)
	      return false;
    }
    return branchVar < other.branchVar;
  }


  bool branchChoice::operator==(const branchChoice& other) const {
    for(int i=0; i<3; i++) {
      if (branch[i].roundedBound == other.branch[i].roundedBound)
        continue;
      else
        return false;
    }
    return true;
  }


  void branchChoice::possibleSwap(size_type i1, size_type i2) {
    register double roundedBound1 = branch[i1].roundedBound;
    register double roundedBound2 = branch[i2].roundedBound;
      if (roundedBound1 < roundedBound2) {
        branchItem tempItem(branch[i1]);
        branch[i1] = branch[i2];
        branch[i2] = tempItem;
      }
  }


//////////////////// branchItem class methods ////////////////////

  void branchItem::set(double bound, double roundQuantum) {
    exactBound = bound;
    if (roundQuantum == 0)
      roundedBound = bound;
    else
      roundedBound = floor(bound/roundQuantum + 0.5)*roundQuantum;
    whichChild    = -1;
  }


  //////////////////// RMA class methods ////////////////////

  // RMA constructor
  RMA::RMA() : workingSol(this), numTotalCutPts(0), numCC_SP(0) {

    //version_info += ", RMA example 1.1";
    min_num_required_args = 1;
    branchingInit(maximization, relTolerance, absTolerance);

    workingSol.serial = 0;
    workingSol.sense  = maximization;

    //workingSol.isPosIncumb=false;

  };   //  Constructor for RMA class


  RMA::~RMA() {
    if ( args->perCachedCutPts() <1 ) {
      int rank, sendbuf, recvbuf;
      //MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      rank = uMPI::rank;
      sendbuf = numCC_SP ;
      DEBUGPRX(0, this,
	       "Local non-stron branching SP is: " << numCC_SP << "\n");

      uMPI::reduceCast(&sendbuf,&recvbuf,1, MPI_INT, MPI_SUM);

      // create new new communicator and then perform collective communications
      //MPI_Reduce(&sendbuf, &recvbuf, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

      // Print the result
      if (rank == 0) {
        cout << "Total non-strong branching SP is: " << recvbuf << "\n";
      }
    }
    /*
      if (verifyLog()) {
      verifyLogFile() << endl; //<< "result " << fathomValue() << endl;
      delete _vlFile;    // Doesn't delete file; actually closes it
      }//*/
    // workingSol.decrementRefs();
  }


  solution* RMA::initialGuess() {

    if (!args->initGuess()) return NULL;
/*
    rmaSolution* guess = new rmaSolution(this);
    grma = new GreedyRMA(data);
    grma->runGreedyRangeSearch();
    guess->isPosIncumb = grma->isPosIncumb;
    guess->value       = grma->maxObjValue;
    guess->a           = grma->L;
    guess->b           = grma->U;
    return guess;
*/
    return NULL;
  }

/*
  void RMA::setParameters(ArgPMA* args, Data* param, const int& deb_int) {

    debug = deb_int;

    _perCachedCutPts = param->_perCachedCutPts;
    _binarySearchCutVal = param->_binarySearchCutVal;
    _perLimitAttrib = param->_perLimitAttrib;

    _checkObjVal = param->_checkObjVal;

    _writeInstances = param->_writeInstances;
    _writeNodeTime = param->_writeNodeTime;
    _writeCutPts = param->_writeCutPts;

    _rampUpSizeFact = param->_rampUpSizeFact;
    _countingSort = param->_countingSort;
    _branchSelection = param->_branchSelection;

  }
*/


  bool RMA::setData(Data* data_) {
    data       = data_;
    numDistObs = data->numTrainObs;
    numAttrib  = data->numAttrib;
    distFeat   = data->distFeat;

    for (int j=0; j<data->numAttrib; ++j)
      numTotalCutPts += data->distFeat[j];

  }


  // writes data with weights to a file whose name we concoct
  // from the iteration number argument; added by JE
  void RMA::writeInstanceToFile(const int& iterNum) 	{
    stringstream s;
    s << 'w' << iterNum << '.' << problemName;
    ofstream instanceOutputFile(s.str().c_str());
    writeWeightedData(instanceOutputFile);
    instanceOutputFile.close();
  }

  // Routine added by JE to write out data with weights.  Note that
  // the sign of the observation is just the last attribute in the
  // "_dataStore" vector of vectors.
  void RMA::writeWeightedData(ostream& os) {

    // Set high precision and scientific notation for weights, while
    // saving old flags and precision
    int oldPrecision = os.precision(16);
    std::ios_base::fmtflags oldFlags = os.setf(ios::scientific);

    // Write data
    for (size_type i=0; i<numDistObs; i++) {
      os << data->intData[i].w << ';';

      // Restore stream state
      os.precision(oldPrecision);
      os.flags(oldFlags);
    }

  }


  void RMA::setCachedCutPts(const int& j, const int& v) {

    bool isAlreadyInCache = false;
    multimap<int,int>::iterator it,itlow,itup;

    itlow = mmapCachedCutPts.lower_bound(j);  // itlow points to
    itup = mmapCachedCutPts.upper_bound(j);   // itup points to

    // print range [itlow,itup):
    for (it=itlow; it!=itup; ++it) {
      if ( (*it).first==j && (*it).second==v ) isAlreadyInCache=true;
      if (args->debug>=10) ucout << (*it).first << " => " << (*it).second << '\n';
    }

    if (args->debug>=10) ucout << "cut point (" << j << ", " << v << ") " ;
    if (args->debug>=10)
      ucout << (isAlreadyInCache ? "is already in cache\n" : "is new\n");

    // if not in the hash table, insert the cut point into the hash table.
    if (!isAlreadyInCache)
      if (0>j || j>data->numAttrib)
        ucout << "ERROR! j is out of range for setCachedCutPts";
      else if (0>v || v>data->distFeat[j])
        ucout << "ERROR! v is out of range for setCachedCutPts";
      else
	mmapCachedCutPts.insert( make_pair(j, v) ) ;

    /*
      CutPt tempCutPt;//  =  {j,v};
      tempCutPt.j = j;
      tempCutPt.v = v;

      int key = rand() % 100000; // TODO: assign unique key

      //map<CutPt,int>::iterator it = mapCachedCutPts.begin();
      //if ( mapCachedCutPts.find(tempCutPt) == mapCachedCutPts.end() ) {

      multimap<int, int>::const_iterator it = mapCachedCutPts.find(tempCutPt);
      bool seenAlready =  (it!=mapCachedCutPts.end());

      DEBUGPR(25, ucout << (seenAlready ? "Seen already\n" : "Looks new\n"));

      // if not in the hash table, insert the cut point into the hash table.
      if (!seenAlready) {
      mapCachedCutPts.insert( make_pair<CutPt, int>(tempCutPt, key) );
      //mapCachedCutPts[tempCutPt] = rand() % 100000;
      DEBUGPR(10, cout << "key:" << key <<", store CutPt: ("
      << j << ", " << v << ")" << "\n");
      //key++;
      } else {
      DEBUGPR(10, cout << "already in stored cut point ("
      << j << ", " << v << ")\n" );
      }
    */

  }

  void RMA::setWeight(vector<double> wt, vector<int> train) {
    for (unsigned int i=0; i<wt.size(); ++i)
      data->intData[train[i]].w = wt[i];
  }

  // Routine added by AK to write out the number of B&B node and CPU time.
  void RMA::writeStatData(ostream& os) {
    os << subCount[2] << ';' << searchTime << '\n';
  }

  // writes the number of B&B node and CPU time; added by AK
  void RMA::writeStatDataToFile(const int& iterNum)  {
    stringstream s;
    s << "BBNode_CPUTime" << '_' << problemName;
    if (iterNum==1) {
      ofstream instanceOutputFile(s.str().c_str(), ofstream::out);
      writeStatData(instanceOutputFile);
      instanceOutputFile.close();
    } else {
      ofstream instanceOutputFile(s.str().c_str(), ofstream::app);
      writeStatData(instanceOutputFile);
      instanceOutputFile.close();
    }
  }


  void RMA::startTime() {

#ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
      timeStart = clock();
#ifdef ACRO_HAVE_MPI
    }
#endif //  ACRO_HAVE_MPI

  }


  double RMA::endTime() {

#ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
#endif //  ACRO_HAVE_MPI
      timeEnd = clock();
      clockTicksTaken = timeEnd - timeStart;
      timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
      cout <<  "Time: " << timeInSeconds <<"\n";
      return timeInSeconds;
#ifdef ACRO_HAVE_MPI
    }
#endif //  ACRO_HAVE_MPI

  }

  // ********************* RMASub methods (start) *******************************************

  void RMASub::setRootComputation() {
    al.resize(numAttrib());
    bl << al;
    au << distFeat();
    bu <<au;
    deqRestAttrib.resize(numAttrib(), false);
    //workingSol() = globalPtr->guess;
  };


  void RMASub::boundComputation(double* controlParam) {
    //globalPtr->getSolution();

    NumTiedSols=1;
    NumPosTiedSols=0;
    NumNegTiedSols=0;

    //setCoveredObs();	// find covered observation which are in [al. bu]
    // sort each feature based on [au, bl]
    coveredObs.resize(global()->sortedObsIdx.size());
    copy(global()->sortedObsIdx.begin(), global()->sortedObsIdx.end(), coveredObs.begin());

    for (int j=0; j<numAttrib(); j++)  bucketSortObs(j);
    setInitialEquivClass();	// set initial equivalence class, vecEquivClass

    // if there are enough discoverd cut points (storedCutPts) check only the list
    if ( global()->args->perCachedCutPts() < 1.0 &&  global()->args->binarySearchCutVal() )
      hybridBranching();
    else if ( global()->args->binarySearchCutVal() )
      binaryBranching();
    else if ( global()->args->perCachedCutPts() < 1.0 )
      cutpointCaching();
    else   //check all cut points
      strongBranching();

    if (global()->debug>=10) printCurrentBounds();

    bound = _branchChoice.branch[0].roundedBound;	// look ahead bound
    setState(bounded);

    if (_branchChoice.branch[0].roundedBound < 0) {
      if (global()->debug>=10) "Bound < 0. \n";
      setState(dead);
      return;
    }

    if ( _branchChoice.branchVar > numAttrib() ) {
      if (global()->debug>=10)
        ucout << "al: " << al << "au: " << au
              << "bl: " << bl << "bu: " << bu ;
      if (global()->debug>=10) cout << "branchVar > numAttrib. \n";
      setState(dead);
      return;
    }

    deqRestAttrib[_branchChoice.branchVar] = true;

    // If (current objValue) >= (current bound), we found the solution.
    if ( workingSol()->value >= _branchChoice.branch[0].exactBound ) {
      if (global()->debug>=10) workingSol()->printSolution();
      if (global()->debug>=10)
        cout << "Bound: " << _branchChoice.branch[0].exactBound << "\n";
      foundSolution();
      setState(dead);
      return;
    }

    //////////////////////////////////// create listExcluded list (start) ////////////////////////
#ifndef ACRO_HAVE_MPI
    sort(listExcluded.begin(), listExcluded.end());
    listExcluded.erase(unique(listExcluded.begin(), listExcluded.end()), listExcluded.end());
    DEBUGPR(150, ucout << "Excluded: " << _branchChoice.branchVar << listExcluded );
    //DEBUGPR(50, ucout << " bound: " << bound << ", sol val=" << getObjectiveVal() << "\n");
#endif
    /////////////////////////////////// check errors (start) /////////////////////////////////////
    if (_branchChoice.branchVar >= numAttrib()) {
      DEBUGPR(20, ucout << "ERROR: branch feature is invalid! (j="
	      << _branchChoice.branchVar << ")\n" );
      cerr << "ERROR: branch feature is invalid! (j="
	   << _branchChoice.branchVar << ")\n" ;
      exit (EXIT_FAILURE);
    }
    if (_branchChoice.cutVal < 0) {
      if (global()->debug>=20)
        ucout << "ERROR: cutValue cannot be less than 0! (cutValue="
	      << _branchChoice.cutVal << ")\n" ;
      exit (EXIT_FAILURE);
    }
    if (_branchChoice.cutVal >= bu[_branchChoice.branchVar]) {
      if (global()->debug>=20)
        ucout << "ERROR: cutValue cannot be >= bu["
	      << _branchChoice.branchVar << "]! (cutValue="
	      << _branchChoice.cutVal << ")\n" ;
      exit (EXIT_FAILURE);
    }
    /////////////////////////////////// check errors (end) /////////////////////////////////////
#ifndef ACRO_HAVE_MPI
    listExcluded.push_back(_branchChoice.branchVar);
#endif
    //////////////////////////////////// create listExclided list (end) ////////////////////////
  }


  int RMASub::getNumLiveCachedCutPts() {
    // numLiveCachedCutPts = (# of live cut points from the cache)
    int j, v, numLiveCachedCutPts=0;
    multimap<int, int>::iterator curr = global()->mmapCachedCutPts.begin();
    multimap<int, int>::iterator end = global()->mmapCachedCutPts.end();

    // count numLiveCachedCutPts and print out cached cut points
    if (global()->args->debug>=20) ucout << "catched cut-points: ";
    while (curr!=end) {
      j = curr->first;
      v = curr->second;
      if (global()->args->debug>=20) ucout << j << ", " << v << "\n";
      //if (j>numAttrib() || v<0) break;
      curr++;
      if (al[j]<=v && v<bu[j]) // if v in [al, bu)
        if ( !( au[j]<bl[j] && au[j]<=v && v<bl[j] ) )   // if not overlapping
          ++numLiveCachedCutPts;
    }
    if (global()->args->debug>=20) ucout << "\n";
    return numLiveCachedCutPts;
  }


  // return how many children to make from current subproblem
  int RMASub::splitComputation() {
    int numChildren=0;
    for (int i=0; i<3; i++)
      if (_branchChoice.branch[i].roundedBound>=0) numChildren++;
    setState(separated);
    return numChildren;
  }


  // make a child of current subproblem
  branchSub* RMASub::makeChild(int whichChild) {
    if (whichChild==-1) {
      cerr << "ERROR: No children to make!\n";
      return NULL;
    }
    if (this->_branchChoice.branchVar>numAttrib()) {
      cerr<< "ERROR: It is not proper attribute!\n";
      return NULL;
    }
    RMASub *temp = new RMASub();
    temp->RMASubAsChildOf(this, whichChild);
    return temp;
  }


  void RMASub::RMASubAsChildOf(RMASub* parent, int whichChild) {

    al << parent->al;
    au << parent->au;
    bl << parent->bl;
    bu << parent->bu;
    deqRestAttrib = parent->deqRestAttrib;

#ifndef ACRO_HAVE_MPI
    listExcluded << parent-> listExcluded;
#endif
    excCutFeat << parent-> excCutFeat;
    excCutVal << parent-> excCutVal;

    globalPtr = parent->global();
    branchSubAsChildOf(parent);

    // set bound
    bound = parent->_branchChoice.branch[whichChild].roundedBound;
    whichChild = parent->_branchChoice.branch[whichChild].whichChild;

    if (global()->args->debug>=10) ucout << "Bound: " << bound << "\n" ;

    int j = parent->_branchChoice.branchVar;
    int lowerBound, upperBound;

    if (j<0) {
      if (global()->args->debug>=20) ucout << "ERROR: feature j cannot be < 0 (j=" << j << ")\n";
      cerr << "ERROR: feature j cannot be < 0 (j=" << j << ")\n";
      return;
    }
    else if (j>numAttrib()) {
      if (global()->args->debug>=20) ucout << "ERROR: feature j cannot be > numAttrib (j=" << j << ")\n";
      cerr << "ERROR: feature j cannot be > numAttrib (j=" << j << ")\n";
      return;
    }

    // Case 1: this feature is overlaping bl < au
    if ( bl[j]<au[j]
	 && bl[j]<=parent->_branchChoice.cutVal
	 && parent->_branchChoice.cutVal<au[j]) {
      // subproblem 1: P(al, v, v, v)
      if (whichChild == 0) { // find a value on [al, bl]
	au[j] = parent->_branchChoice.cutVal;//optCutValue;
	bl[j] = parent->_branchChoice.cutVal;//optCutValue;
	bu[j] = parent->_branchChoice.cutVal;//optCutValue;
      }
      // subproblem 3: P(al, v, v+1, bu)
      else if (whichChild == 2) { // find a value on [bl, au]
	au[j] = parent->_branchChoice.cutVal;//optCutValue;
	bl[j] = parent->_branchChoice.cutVal+1;//optCutValue+1;
      }
      // subproblem 2: P(v+1, v+1, v+1, bu)
      else if (whichChild == 1) { // find a value on [au, bu]
	al[j] = parent->_branchChoice.cutVal+1;//optCutValue+1;
	au[j] = parent->_branchChoice.cutVal+1;//optCutValue+1;
	bl[j] = parent->_branchChoice.cutVal+1;//optCutValue+1;
      }

    } else { // Case 2: this feature is not overlaping au <= bl

      if (au[j]<=bl[j]) { lowerBound = au[j]; upperBound = bl[j]; }
      else { lowerBound = bl[j]; upperBound = au[j]; }

      // find a value on [al, au]
      if ( al[j]<=parent->_branchChoice.cutVal && parent->_branchChoice.cutVal<lowerBound ) {
	if (whichChild == 0) // subproblem 1: P(al, v, bl, bu)
	  au[j] = parent->_branchChoice.cutVal;//optCutValue;
	else if (whichChild == 1) { // subproblem 2: P(v+1, au, bl, bu)
	  al[j] = parent->_branchChoice.cutVal+1;//optCutValue+1;
	}
	// find new cut value on [bl, bu]
      } else if ( upperBound<=parent->_branchChoice.cutVal && parent->_branchChoice.cutVal<bu[j] ) {
	if (whichChild == 0) // subproblem 1: P(al, au, bl, v)
	  bu[j] = parent->_branchChoice.cutVal;//optCutValue;
	else if (whichChild == 1) // subproblem 2: P(al, au, v+1, bu)
	  bl[j] = parent->_branchChoice.cutVal+1;//optCutValue+1;
      }

    } // end if

    if (al[j]>au[j]) {
      cerr << "Error al[j]>au[j]! since al[" << j << "]=" << al[j]
	   << "and au[" << j << "]=" << au[j] << endl;
      return;
    }
    if (global()->args->debug>=10)
      ucout << "al: " << al << "au: " << au
	    << "bl: " << bl << "bu: " << bu ;
  }


  // find a particular subproblem object to the problem description embodied in the object master
  void RMASub::RMASubFromRMA(RMA* master) {
    globalPtr = master;		// set a globalPtr
    //workingSol()->value = getObjectiveVal(); // set bound value as current solution
    if (global()->args->debug>=20)
      ucout << "Created blank problem, out of rmaSub:::RMASubFromRMA" << "\n";
  };


  bool RMASub::candidateSolution() {
    if (global()->args->debug>=20)
      ucout << "al: " << al << "au: " << au
	    << "bl: " << bl << "bu: " << bu ;
    for (int j=0; j<numAttrib() ; ++j) {
      if ( al[j] != au[j] ) return false;
      if ( bl[j] != bu[j] ) return false;
    }
    workingSol()->a << al;
    workingSol()->b << bu;

    //cout << coveredObs << endl;
    //sort(coveredObs.begin(), coveredObs.begin()+coveredObs.size());
    //cout << coveredObs << endl;

    //workingSol()->isCovered.resize(numDistObs());
    //for (int i=0; i<numDistObs(); ++i)
    //  workingSol()->isCovered[i]=false;
    //for (int i=0; i<vecEquivClass1.size(); ++i)
    //for (int j=0; j<vecEquivClass1[i].size(); ++j)

    /*
      for (int i=0; i<numDistObs(); ++i)
      for (int j=0; j<numAttrib(); ++j)
      if ( workingSol()->a[j] <= global()->intData[i].X[j] &&
      global()->intData[i].X[j] <= workingSol()->b[j] ) {
      if ( j==numAttrib()-1)
      workingSol()->isCovered[i]= true;
      } else {
      workingSol()->isCovered[i]= false;
      break;
      }
    */
    if (global()->args->debug>=5) workingSol()->printSolution();
    return true;
  }


  // ******************* RMASub helper functions (start) *******************************

  void RMASub::branchingProcess(const int& j, const int& v) {

    vecBounds.resize(3);
    vecBounds[2]=-1;

    if ( al[j]<=v && v < min(au[j],bl[j]) ) {

      // Middle Child
      printSP(j, al[j], v, bl[j], bu[j]);
      mergeEquivClass(j, al[j], v, bl[j], bu[j]);
      vecBounds[0]=getBoundMerge();

      // Up Child
      printSP(j, v+1, au[j], bl[j], bu[j]);
      dropEquivClass(j, v+1, bu[j]); // al, bu
      vecBounds[1]=getBoundDrop();

    } else if ( au[j]<=bl[j] && au[j]<=v && v<bl[j]  ) {
      return;
    } else if ( bl[j]<au[j] &&
      min(au[j],bl[j])<=v && v<max(au[j],bl[j]) ) {

      // Down Child
      printSP(j, al[j], v, v, v);
      dropEquivClass(j, al[j], v);
      vecBounds[0]=getBoundDrop();

      // Middle Child
      printSP(j, al[j], v, v+1,  bu[j]);
      mergeEquivClass(j, al[j], v, v+1,  bu[j]);
      vecBounds[2]=getBoundMerge();

      vecBounds[1] = vecBounds[2] - vecBounds[0];
      /*
      // Up Child
      printSP(j, v+1, v+1, v+1, bu[j]);
      dropEquivClass(j, v+1, bu[j]);
      vecBounds[1]=getBoundDrop();
      */

    } else if ( max(au[j],bl[j])<=v && v<bu[j]) {

      // Down Child
      printSP(j, al[j], au[j], bl[j], v);
      dropEquivClass(j, al[j], v);
      vecBounds[0]=getBoundDrop();

      // Middle Child
      printSP(j, al[j], au[j], v+1, bu[j]);
      mergeEquivClass(j, al[j], au[j], v+1, bu[j]);
      vecBounds[1]=getBoundMerge();

    }

    branchChoice thisChoice(vecBounds[0], vecBounds[1], vecBounds[2], v, j);

    if (global()->args->debug>=15) ucout << "Evaluating" << thisChoice << "\n";

    // select variable based on minimum of children
    // bounds given in lexicographically decreasing order
    thisChoice.sortBounds();

    for (int i=0; i<vecBounds.size(); i++)
      if (thisChoice.branch[i].exactBound <= workingSol()->value ) {
        thisChoice.branch[i].exactBound=-1;
        thisChoice.branch[i].roundedBound=-1;
      }

    if (global()->args->debug>=15) cout << "Sorted version is " << thisChoice << "\n";

    if (thisChoice < _branchChoice) {
      //cout << "branchBound: " << thisChoice.branch[0].exactBound << " "
      //     << _branchChoice.branch[0].exactBound;
      _branchChoice = thisChoice;
      if (global()->args->debug>=50) ucout << "Improves best attribute: " << j << "\n";
      if (global()->args->debug>=10) ucout << "Branch choice now: " << _branchChoice << "\n";
      NumTiedSols=1;
      //foundBound=true;
    } else if (thisChoice == _branchChoice) {
      //cout << "branchBound: " << thisChoice.branch[0].exactBound << " "
      //     << _branchChoice.branch[0].exactBound;
      if (global()->args->branchSelection()==0) {
        NumTiedSols++;
        (globalPtr->args->randSeed()) ? srand(NumTiedSols*time(NULL)*100) : srand(1);
        double rand_num = (rand() % 10001 ) / 10000.0 ;
        //DEBUGPRX(0, global(), "rand: " << rand_num  << "\n");
        //DEBUGPRX(0, global(), "rand1: " << 1.0 /  NumTiedSols << "\n");
        if ( rand_num  <= 1.0 /  NumTiedSols ) {
          _branchChoice = thisChoice;
          if (global()->args->debug>=50) ucout << "Improves best attribute: " << j << "\n";
          if (global()->args->debug>=10) ucout << "Branch choice now is: " << _branchChoice << "\n";
        }
      } else if (global()->args->branchSelection()==2) {
        _branchChoice = thisChoice;
        if (global()->args->debug>=50) ucout << "Improves best attribute: " << j << "\n";
        if (global()->args->debug>=10) ucout <<"Branch choice now is: " << _branchChoice << "\n";
      }
    }

  } // end function void RMASub::branchingProcess


  // strong branching
  void RMASub::strongBranching() {

    int numCutPtsInAttrib;
    if (global()->args->debug>=10)
      ucout << "al: " << al << "au: " << au << "bl: " << bl << "bu: " << bu ;

    if (global()->args->debug>=10) ucout << "sortedObs: " << coveredObs;
    compIncumbent(numAttrib()-1);

    for (int j=0; j<numAttrib(); ++j ) {

      if (global()->args->debug>=10) ucout << "original: ";
      printSP(j, al[j], au[j], bl[j], bu[j]);

      numCutPtsInAttrib = bu[j]-al[j] - max(0, bl[j]-au[j]);
      if (numCutPtsInAttrib==0) continue;

      for (int v=al[j]; v<bu[j]; ++v ) {
        if (au[j]<=v && v<bl[j] ) { v=bl[j]-1; continue; }
        //cout << "Size of sortedObs: " << coveredObs1.size() << "\n";
        //cout << "sortedObs: " << coveredObs;
        branchingProcess(j, v);
      }
      if (j==numAttrib()-1) break;
      global()->args->countingSort() ? countingSortEC(j) : bucketSortEC(j);
      compIncumbent(j);
    } // end for each feature

  } // end RMASub::strongBranching


  // branching using cut-point caching methods
  void RMASub::cachedBranching() {

    int k=0;

    sortCachedCutPtByAttrib();
    cachedCutPts = sortedCachedCutPts;
    compIncumbent(numAttrib()-1);

    for (int j=0; j<numAttrib(); ++j) {
      while ( k<cachedCutPts.size() ) {
        if ( j == cachedCutPts[k].j ) {
          branchingProcess(cachedCutPts[k].j, cachedCutPts[k].v);
          ++k;
        } else break;
      }

      if (j==numAttrib()-1) break;
      (global()->args->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
      compIncumbent(j);
    }

  } // end RMASub::cachedBranching


  // split subproblems and choose cut value by binary search
  void RMASub::hybridBranching() {

    bool firstFewCutPts;
    int l, u, L, U, cutValue, numCutValues;
    vector<bool> vecCheckedCutVal;
    int k=0;

    sortCachedCutPtByAttrib();
    cachedCutPts = sortedCachedCutPts;
    compIncumbent(numAttrib()-1);

    for (int j=0; j<numAttrib(); ++j ) {  // for each attribute

      if (distFeat()[j]<30) {
        while ( k<cachedCutPts.size() ) {
          if ( j == cachedCutPts[k].j ) {
            branchingProcess(cachedCutPts[k].j, cachedCutPts[k].v);
            ++k;
          } else break;
        }
        if (j==numAttrib()-1) break;
        (global()->args->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
        compIncumbent(j);

      } else {  // binary search

        numCutValues = bu[j]-al[j];
        if (bl[j]>au[j]) numCutValues -= bl[j]-au[j];

        if (numCutValues==0)	{// if no cutValue in this feature,
	        (global()->args->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
          continue;			// then go to the next attribute.
        }

        cutValue=-1; l=0; u=0; L = al[j]; U = bu[j]; firstFewCutPts=true;
        vecCheckedCutVal.clear();
        vecCheckedCutVal.resize(distFeat()[j]+1);

        while ( true ) {
          if (numCutValues>3) {
            cutValue = L+(U-L)/2;   //integer division
            if (au[j]<=cutValue && cutValue<bl[j] ) {
              cutValue = au[j]-1-l;
              l++;
              if ( cutValue<al[j] && bl[j]<bu[j] ) {
                cutValue = bl[j]+u;
                u++;
                //cout << "1 (j, cutValue) " << j << ", " << cutValue << "\n";
              } else if  ( cutValue>=bu[j] ) {
                //cout << "2 (j, cutValue) " << j << ", " << cutValue << "\n";
                break; // no cut point in this feature
              }
            }
            if ( cutValue>=bu[j] ) {
              if (global()->args->debug>=10) ucout << "cutValue>=bu[j] .\n";
              break;
            }
            if ( cutValue<al[j] ) {
              if (global()->args->debug>=10) ucout << "cutValue<al[j] .\n";
              break;
            }
            if (vecCheckedCutVal[cutValue]) {
              if (global()->args->debug>=10) ucout << "break since oldCutValue=cutValue.\n";
              break;
            }
            vecCheckedCutVal[cutValue] = true;

            printSP(j, al[j], au[j], bl[j], bu[j]);
            if (global()->args->debug>=10) ucout << "j: " << j << " L: "<< L << " U: " << U
                << " cutVal: " << cutValue << "\n";

          } else {
            if (firstFewCutPts) {
              cutValue=L;
              firstFewCutPts=false;
            } else {
              cutValue++;
            }
            if ( au[j]<=cutValue && cutValue<bl[j] ) cutValue=bl[j];
            if ( cutValue>=bu[j] ) break;
          }

          //cout << "(j, cutValue) " << j << ", " << cutValue << "\n";
          if (global()->args->debug>=10) cout << "coveredObs: " << coveredObs;
          branchingProcess(j, cutValue);

          // compare objectives instead of bounds
          //if ( vecObjValue[0] > vecObjValue[1] ) L = cutValue;
          //else U = cutValue;

          // Compare bounds
          if ( vecBounds[0] < vecBounds[1] ) L = cutValue; else U = cutValue;

        }  // end while for each cut value of feature f

        if (j==numAttrib()-1) break;

        (global()->args->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
        compIncumbent(j);

      } // end binary search
    }	// end for each attribute

  } // end RMASub::hybridBranching


  // split subproblems and choose cut value by binary search
  void RMASub::binaryBranching() {

    bool firstFewCutPts;
    int l, u, L, U, cutValue, numCutValues;
    vector<bool> vecCheckedCutVal;

    compIncumbent(numAttrib()-1);

    for (int j=0; j<numAttrib(); ++j ) {  // for each attribute

      numCutValues = bu[j]-al[j];
      if (bl[j]>au[j]) numCutValues -= bl[j]-au[j];

      if (numCutValues==0)	{// if no cutValue in this feature,
        (global()->args->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
        continue;			// then go to the next attribute.
      }

      cutValue=-1; l=0; u=0; L = al[j]; U = bu[j]; firstFewCutPts=true;
      vecCheckedCutVal.clear();
      vecCheckedCutVal.resize(distFeat()[j]+1);

      while ( true ) {
        if (numCutValues>3) {
          cutValue = L+(U-L)/2;   //integer division
          if (au[j]<=cutValue && cutValue<bl[j] ) {
            cutValue = au[j]-1-l;
            l++;
            if ( cutValue<al[j] && bl[j]<bu[j] ) {
              cutValue = bl[j]+u;
              u++;
              //cout << "1 (j, cutValue) " << j << ", " << cutValue << "\n";
            } else if  ( cutValue>=bu[j] ) {
              //cout << "2 (j, cutValue) " << j << ", " << cutValue << "\n";
              break; // no cut point in this feature
            }
          }
          if ( cutValue>=bu[j] ) {
            if (global()->args->debug>=10) ucout << "cutValue>=bu[j] .\n";
            break;
          }
          if ( cutValue<al[j] ) {
            if (global()->args->debug>=10) ucout << "cutValue<al[j] .\n";
            break;
          }
          if (vecCheckedCutVal[cutValue]) {
            if (global()->args->debug>=10) ucout << "break since oldCutValue=cutValue.\n";
            break;
          }
          vecCheckedCutVal[cutValue] = true;

          printSP(j, al[j], au[j], bl[j], bu[j]);
          if (global()->args->debug>=10) ucout << "j: " << j << " L: "<< L << " U: " << U
              << " cutVal: " << cutValue << "\n";

        } else {
          if (firstFewCutPts) {
            cutValue=L;
            firstFewCutPts=false;
          } else {
            cutValue++;
          }
          if ( au[j]<=cutValue && cutValue<bl[j] ) cutValue=bl[j];
          if ( cutValue>=bu[j] ) break;
        }

        //cout << "(j, cutValue) " << j << ", " << cutValue << "\n";
        if (global()->args->debug>=10) ucout << "coveredObs: " << coveredObs;
  			branchingProcess(j, cutValue);

        // compare objectives instead of bounds
        //if ( vecObjValue[0] > vecObjValue[1] ) L = cutValue;
        //else U = cutValue;

        // Compare bounds
        if ( vecBounds[0] < vecBounds[1] ) L = cutValue; else U = cutValue;

      }  // end while for each cut value of feature f

      if (j==numAttrib()-1) break;

      bucketSortEC(j);
      compIncumbent(j);
    }	// end for each attribute

  } // end RMASub::binarySplitSP


  void RMASub::cutpointCaching() {

    // numLiveCachedCutPts = (# of live cut points from the cache)
    int numLiveCachedCutPts=getNumLiveCachedCutPts();

    // if numCachedCutPts is less than the percentage, check all cut points
    if ( numLiveCachedCutPts
          < global()->numTotalCutPts * global()->args->perCachedCutPts() )
      strongBranching();

    else { // if not, only check the storedCutPts
      // count number of subproblems only discovering cut-points from the chache
      ++global()->numCC_SP;
      int j, v, l=-1;
      multimap<int, int>::iterator curr = global()->mmapCachedCutPts.begin();
      multimap<int, int>::iterator end  = global()->mmapCachedCutPts.end();
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

      cachedBranching();
    }

    int j = _branchChoice.branchVar;

    // store cached cut-points
    if ( _branchChoice.branchVar < 0
      || _branchChoice.branchVar > globalPtr->data->numAttrib )
      return; //ucout << "ERROR! j is out of range for setCachedCutPts";
    else if ( _branchChoice.cutVal < 0
      || _branchChoice.cutVal > globalPtr->data->distFeat[j] )
      return; // ucout << "ERROR! v is out of range for setCachedCutPts";
    else
      globalPtr->setCachedCutPts(_branchChoice.branchVar, _branchChoice.cutVal);

  }


  void RMASub::bucketSortObs(const int& j) {
    int v, l=-1;
    int size = bu[j] - al[j] + 1 - max(0, bl[j]-au[j]);
    vector<vector<int> > buckets;
    buckets.resize(size);

    for (int i=0; i<coveredObs.size(); i++) {
      v = global()->data->intData[coveredObs[i]].X[j];
      if ( au[j] < bl[j] ) { 	// no overlapping
        if (v<au[j]) v -= al[j];
        else if ( au[j]<=v && v<=bl[j] ) v = au[j] - al[j];
        else if ( bl[j] < v ) v = v - (bl[j] - au[j]) - al[j];
      } else v -= al[j];	// overlapping

      if (v<0) {
	if (global()->args->debug>=10) cout << "below covered range \n"; continue;
      } else if (v>=size) {
	if (global()->args->debug>=10) cout << "above covered range \n"; continue;
      }

      buckets[v].push_back(coveredObs[i]);
    }

    // walk buckets to get sorted observation list on this attribute
    for (int v=0; v<size; v++)
      for (int k=0; k<buckets[v].size(); k++)
	coveredObs[++l] = buckets[v][k];

    coveredObs.resize(l+1);
  }


  void RMASub::bucketSortEC(const int& j) {
    int v, obs, l=-1;
    int size = bu[j] - al[j] + 1 - max(0, bl[j]-au[j]);
    if (size==1) return; // do not have to sort since all values are inseperable
    vector<vector<int> > buckets(size);

    for (int i=0; i<sortedECidx.size(); i++) {
      obs = vecEquivClass[sortedECidx[i]].getObs();
      v = global()->data->intData[obs].X[j];
      if ( au[j] < bl[j] ) { 	// no overlapping
        if (v<au[j]) v -= al[j];
        else if ( au[j]<=v && v<=bl[j] ) v = au[j] - al[j];
        else if ( bl[j] < v ) v = v - (bl[j] - au[j]) - al[j];
      } else v -= al[j];	// overlapping

      if (v<0) {
        if (global()->args->debug>=10) cout << "below covered range \n"; continue;
      } else if (v>=size) {
        if (global()->args->debug>=10) cout << "above covered range \n"; continue;
      }

      buckets[v].push_back(sortedECidx[i]);
    }

    // walk buckets to get sorted observation list on this attribute
    for (int v=0; v<size; v++)
      for (int k=0; k<buckets[v].size(); k++)
        sortedECidx[++l] = buckets[v][k];

    //sortedECidx.resize(l+1);
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


  void RMASub::countingSortEC(const int& j) {

    int i, v, obs;
    int numEC = sortedECidx.size();
    int bucketSize = bu[j] - al[j] + 1 - max(0, bl[j]-au[j]);
    vector<int> vecCount(bucketSize);
    sortedECidx1.resize(numEC);

    for ( i=0; i < numEC ; ++i ) {
      obs = vecEquivClass[sortedECidx[i]].getObs();
      v = global()->data->intData[obs].X[j];
      if ( au[j] < bl[j] ) { 	// no overlapping
        if (v<au[j]) v -= al[j];
        else if ( au[j]<=v && v<=bl[j] ) v = au[j] - al[j];
        else if ( bl[j] < v ) v = v - (bl[j] - au[j]) - al[j];
      } else v -= al[j];	// overlapping
      ++vecCount[v] ;
    }

    for ( i=1; i < bucketSize ; ++i )
      vecCount[i] += vecCount[i-1] ;

    for ( i=numEC-1; i>=0 ; --i ) {
      obs = vecEquivClass[sortedECidx[i]].getObs();
      v = global()->data->intData[obs].X[j];
      if ( au[j] < bl[j] ) { 	// no overlapping
        if (v<au[j]) v -= al[j];
        else if ( au[j]<=v && v<=bl[j] ) v = au[j] - al[j];
        else if ( bl[j] < v ) v = v - (bl[j] - au[j]) - al[j];
      } else v -= al[j];	// overlapping
      sortedECidx1[ vecCount[v]-1 ] = sortedECidx[i] ;
      --vecCount[v];
    }

    sortedECidx = sortedECidx1;

  }


  void RMASub::compIncumbent(const int& j) {

    tmpMin = inf;
    tmpMax = -inf;
    //minVal = inf;
    //maxVal = -inf;
    minVal = globalPtr->args->initGuess() ?   workingSol()->value : inf;
    maxVal = globalPtr->args->initGuess() ?  -workingSol()->value : -inf;
    optMinAttrib=-1;
    optMaxAttrib=-1;

    curObs=0;
    tmpMin = runMinKadane(j);
    if (tmpMin==minVal) {
      NumNegTiedSols++;
      (globalPtr->args->randSeed()) ? srand(NumNegTiedSols*time(NULL)*100) : srand(1);
      rand_num = (rand() % 10001 ) / 10000.0 ;
      if ( rand_num <= 1.0/NumNegTiedSols ) setOptMin(j);
    } else if (tmpMin<minVal) { // if better min incumbent was found
      NumNegTiedSols=1;
      setOptMin(j);
    }

    curObs=0;
    tmpMax = runMaxKadane(j);
    if (tmpMax==maxVal) {
      NumPosTiedSols++;
      (globalPtr->args->randSeed()) ? srand(NumNegTiedSols*time(NULL)*100) : srand(1);
      rand_num = (rand() % 10001 ) / 10000.0 ;
      if ( rand_num <= 1.0/NumPosTiedSols ) setOptMax(j);
    } else if (tmpMax>maxVal) {
      NumPosTiedSols=1;
      setOptMax(j);
    }

    chooseMinOrMaxRange();

  }


  void RMASub::chooseMinOrMaxRange() {
    if ( max(maxVal, -minVal) > workingSol()->value +.000001) {
      (globalPtr->args->randSeed()) ?srand((NumNegTiedSols+NumPosTiedSols)*time(NULL)*100) : srand(1);
      rand_num = (rand() % 10001 ) / 10000.0 ;
      workingSol()->a << al;
      workingSol()->b << bu;
      if (maxVal>-minVal ||
	  ( maxVal==minVal &&
	    rand_num <= NumPosTiedSols/(double)(NumNegTiedSols+NumPosTiedSols) ) ) {
        workingSol()->value = maxVal;
        workingSol()->a[optMaxAttrib]=optMaxLower;
        workingSol()->b[optMaxAttrib]=optMaxUpper;
        workingSol()->isPosIncumb=true;
        if (globalPtr->args->debug>=5) cout << "positive ";
      } else  {
        workingSol()->value =-minVal;
        workingSol()->a[optMinAttrib]=optMinLower;
        workingSol()->b[optMinAttrib]=optMinUpper;
        workingSol()->isPosIncumb=false;
        if (globalPtr->args->debug>=5) cout << "negative ";
      }
      foundSolution();
      if (globalPtr->args->debug>=1) cout << "new incumbent  " << workingSol()->value << '\n';
      if (globalPtr->args->debug>=5) workingSol()->printSolution();
      //DEBUGPR(10, workingSol()->checkObjValue1(workingSol()->a, workingSol()->b,
      //        coveredObs,sortedECidx ));
    }
  }


  void RMASub::setOptMin(const int &j){

    minVal = tmpMin;

    optMinAttrib = j;
    optMinLower = aj;
    optMinUpper = bj;
    if (au[j]<bl[j] && au[j]<=optMinUpper && optMinUpper<=bl[j])
      optMinUpper = bl[j];

    if (globalPtr->args->debug>=10)
      cout << "optAttrib: (a,b): " << optMinAttrib
	   << ": (" << optMinLower << ", " << optMinUpper
	   << "), min: " << minVal << "\n" ;

  }


  void RMASub::setOptMax(const int &j){

    maxVal = tmpMax;

    optMaxAttrib = j;
    optMaxLower = aj;
    optMaxUpper = bj;
    if (au[j]<bl[j] && au[j]<=optMaxUpper && optMaxUpper<=bl[j])
      optMaxUpper = bl[j];

    if (globalPtr->args->debug>=10)
      cout << "optAttrib: (a,b): " << optMaxAttrib
	   << ": (" << optMaxLower  << ", " << optMaxUpper
	   << "), max: " << maxVal << "\n" ;

  }


  // get Maximum range for the feature
  double RMASub::runMaxKadane(const int& j) {

    double maxEndHere, maxSoFar, tmpObj;
    int s=al[j];
    aj=al[j];
    bj=al[j];
    maxEndHere = 0;
    maxSoFar   = -inf;

    for (int v=al[j]; v <= bu[j]; ++v) {

      if (au[j]<bl[j] && au[j]<v && v<=bl[j]) {
        v=bl[j];
        continue;
      }

      maxEndHere += getObjValue(j, v);

      if (maxEndHere > maxSoFar) {
        maxSoFar = maxEndHere;
        aj=s;
        bj=v;
      }

      if ( maxEndHere<0 ) {
        maxEndHere = 0;
        s = v+1;
        if (au[j]<bl[j] && v==au[j])
          s = bl[j];
      }

    } // end for each value in this attribute

    if (globalPtr->args->debug>=10)
      cout << "Maximum contiguous sum is " << maxSoFar ;
    if (globalPtr->args->debug>=10)
      cout << " attribute (L,U): " << j << " (" << aj << ", " << bj << ")\n" ;

    return maxSoFar;

  }


  // get Miniumum range for the feature
  double RMASub::runMinKadane(const int& j) {

    double minEndHere, minSoFar, tmpObj;
    int s=al[j];
    aj=al[j];
    bj=al[j];
    minEndHere = 0;
    minSoFar   = inf;

    for (int v=al[j]; v <= bu[j]; ++v) {

      if (au[j]<bl[j] && au[j]<v && v<=bl[j]) {
        v=bl[j];
        continue;
      }

      minEndHere += getObjValue(j, v);

      if (minEndHere < minSoFar) {
        minSoFar = minEndHere;
        aj=s;
        bj=v;
      }

      if ( minEndHere>0 ) {
        minEndHere = 0;
        s = v+1;
        if (au[j]<bl[j] && v==au[j])
          s = bl[j];
      }

/*
      if ( tmpObj < minEndHere+tmpObj ) {
        minEndHere = tmpObj;
        if (minEndHere < minSoFar) {
          minSoFar=minEndHere;
          aj=v;
          bj=v;
          if (au[j]<bl[j] && au[j]<v && v<=bl[j])
            bj=bl[j];
        }
      } else {
        minEndHere += tmpObj;
        if (minEndHere < minSoFar) {
          minSoFar=minEndHere;
          bj=v;
          if (au[j]<bl[j] && au[j]<v && v<=bl[j])
            bj=bl[j];
        }
      }*/
    }

    if (globalPtr->args->debug>=10) cout << "Minimum contiguous sum is " << minSoFar ;
    if (globalPtr->args->debug>=10) cout << " attribute (L,U): " << j << " ("
                            <<  aj << ", " << bj << ")\n" ;
    return minSoFar;
  }


  double RMASub::getObjValue(const int& j, const int& v) {
    int obs;
    double covgWt = 0.0;
    if (globalPtr->args->debug>=20) ucout << "j: " << j << " v:" << v ;

    for (int i=curObs; i<sortedECidx.size(); ++i) {

      obs = vecEquivClass[sortedECidx[i]].getObs();

      // if the observation's jth attribute value = cut-value
      if ( global()->data->intData[obs].X[j] == v ) {
        covgWt += vecEquivClass[sortedECidx[i]].getWt();
        //cout << "vecEquivClass1: " << sortedECidx[i]
        //     << " covgWt: " << covgWt << endl;
      } else if (au[j]<bl[j] && au[j]<=v && v<=bl[j]
              && au[j]<=global()->data->intData[obs].X[j]
              && global()->data->intData[obs].X[j]<=bl[j]) {
        covgWt += vecEquivClass[sortedECidx[i]].getWt();
        //cout << "vecEquivClass2: " << sortedECidx[i]
        //     << " covgWt: " << covgWt << endl;

      } else if (global()->data->intData[obs].X[j]<v) {
        if (globalPtr->args->debug>=0) cout << "X[j] < v! ";
        //*
        if (globalPtr->args->debug>=20) cout << "curObs: " << curObs << " attribute: " << j << "; "
             << global()->data->intData[obs].X[j]
             << " < cutVal: " << v << "\n";
        for (int i=0; i<coveredObs.size(); ++i)
          if (globalPtr->args->debug>=20) cout << global()->data->intData[sortedECidx[i]].X[j]
            << " bound (" << al[j] << ", " << au[j] << ", "
               << bl[j] << ", " << bu[j] << ")\n )";
        ///
      } else {
        curObs = i;
        break;
      }
      if (i==sortedECidx.size()-1) curObs = sortedECidx.size();
    }  // end for each covered observation

    if (globalPtr->args->debug>=20) ucout << "covgWt: " << covgWt << "\n";
    return covgWt ;

  }  // end function getObjValue


	double RMASub::getBoundMerge() const {

    int obs, idxEC;	// observation number
    double pBound=0.0, nBound=0.0; // weight for positive and negative observation

    for (int i=0; i< vecEquivClass1.size(); i++) { // for each equivalence class
      if (vecEquivClass1[i].getWt() > 0 ) // for positive observation
        pBound +=  vecEquivClass1[i].getWt();
      else if (vecEquivClass1[i].getWt() < 0 )// for negative observation
        nBound -=  vecEquivClass1[i].getWt();
    } // end outer for

    return pBound>nBound ? pBound : nBound;
  } // end functino RMASub::getBound


  double RMASub::getBoundDrop() const {

    int obs, idxEC;	// observation number
    double pBound= 0.0, nBound= 0.0;; // weight for positive and negative observation

    for (int i=0; i< sortedECidx1.size(); i++) { // for each equivalence class
      idxEC = sortedECidx1[i];
      if (vecEquivClass[idxEC].getWt() > 0 ) // for positive observation
        pBound +=  vecEquivClass[idxEC].getWt();
      else if (vecEquivClass[idxEC].getWt() < 0 )// for negative observation
        nBound -=  vecEquivClass[idxEC].getWt();
    } // end outer for

    return pBound>nBound ? pBound : nBound;
  } // end functino RMASub::getBound


  void RMASub::setInitialEquivClass() {

    if (coveredObs.size()<=0) {
      if (globalPtr->args->debug>=0) cout << "coveredObs is empty.\n";
      return;
    } else if (coveredObs.size()==1) {
      if (globalPtr->args->debug>=15) cout << "There is only one covered observation" << "\n";
      vecEquivClass.resize(1);
      vecEquivClass[0].addObsWt(coveredObs[0],
				global()->data->intData[coveredObs[0]].w);
      return;
    }

    vecEquivClass.resize(coveredObs.size());

    int obs1 = coveredObs[0];
    int obs2 = coveredObs[1];
    int k=0;
    vecEquivClass[0].addObsWt(obs1, global()->data->intData[obs1].w);

    for (int i=1; i<coveredObs.size(); ++i) { // for each sorted, covered observation

      for (int j=0; j<numAttrib(); ++j) { // for each attribute

        if ( isInSameClass(obs1, obs2, j, au[j], bl[j]) ) {
          if (j==numAttrib()-1) { // if it is in the same equivalent class
            vecEquivClass[k].addObsWt(obs2, global()->data->intData[obs2].w);
            if (i!=coveredObs.size()-1)  // if not the last observation
              obs2=coveredObs[i+1];
          }

        } else {  // detected obs1 and obs2 are in different equivClass
          vecEquivClass[++k].addObsWt(obs2, global()->data->intData[obs2].w);
          if (i!=coveredObs.size()-1) { // if not the last observation
            obs1 = coveredObs[i];
            obs2 = coveredObs[i+1];
          }
          break; // as soon as we detect obs1 and obs2 are in different equivClass
                 // compare the next observation combinations
        }

      } // end for each attribute j
    } // end for each obs, i

    // erase extra space
    vecEquivClass.erase(vecEquivClass.begin()+k+1, vecEquivClass.end());

    sortedECidx.resize(vecEquivClass.size());
    for (int i=0; i<vecEquivClass.size(); ++i) sortedECidx[i] = i;

    if (globalPtr->args->debug>=10) ucout << "Size of coveredObs: " << sortedECidx.size() << "\n";
    if (globalPtr->args->debug>=20) ucout << "Size of vecEquivClass: " << vecEquivClass.size() << "\n";
    if (globalPtr->args->debug>=30)
      for (int i=0; i<vecEquivClass.size(); ++i)
        cout << "EC: " << i << ": " << vecEquivClass[i] << "\n" ;
  }


  void RMASub::mergeEquivClass(const int& j, const int& al_, const int& au_,
			       const int& bl_, const int& bu_) {

    if ( al_!=al[j] || bu_!=bu[j] )
      dropEquivClass(j, al_, bu_);
    else {
      sortedECidx1.resize(sortedECidx.size());
      copy(sortedECidx.begin(), sortedECidx.end(), sortedECidx1.begin());
    }

    vecEquivClass1.clear();

    if (sortedECidx1.size()<=0){
      if (globalPtr->args->debug>=0) cout << "sortedECidx1 is empty. \n";
      return;
    }

    if (sortedECidx1.size()==1){
      if (globalPtr->args->debug>=15) cout << "There is only one equivalence class" << "\n";
      vecEquivClass1.resize(1);
      vecEquivClass1[0] = vecEquivClass[sortedECidx1[0]];
      return;
    }

    int idxEC1 = sortedECidx1[0];
    int idxEC2 = sortedECidx1[1];

    int obs1 = vecEquivClass[idxEC1].getObs();
    int obs2 = vecEquivClass[idxEC2].getObs();
    int k=0, J, _au, _bl;

    vecEquivClass1.resize(sortedECidx1.size());
    vecEquivClass1[0] = vecEquivClass[idxEC1];

    for (int i=1; i<sortedECidx1.size(); ++i) { // for each sorted, covered observation

      for (int f=0; f<numAttrib(); ++f) { // for each attribute

        // always search merging equivalence classes from leaves
        J = f+j;
        if (J>=numAttrib()) J-=numAttrib();	// if J=numAttribute(), go back to J=0

	// if J is the modified attribute in order to check these two observations
	// are in the same equivalence class or not from the leaves of equivalence tree
	if (J==j) { _au = au_; _bl =  bl_; } //DEBUGPR(50,ucout << j << _au << _bl << endl);
	else {_au = au[J]; _bl = bl[J]; }

	if ( isInSameClass(obs1, obs2, J, _au, _bl) ) {
	  if (f==numAttrib()-1) { // if it is in the same equivalent class
            vecEquivClass1[k].addEC(vecEquivClass[idxEC2]);
	    if (i!=sortedECidx1.size()-1)  // if not the last equivClass
	      idxEC2 = sortedECidx1[i+1];
            obs2 = vecEquivClass[idxEC2].getObs();
            break;
	  }
	} else {  // detected obs1 and obs2 are in different equivClass
	  vecEquivClass1[++k]=vecEquivClass[idxEC2];  // push back all obs in the equivClass
	  if (i!=sortedECidx1.size()-1) { // if not the last observation
            idxEC1 = sortedECidx1[i];
            idxEC2 = sortedECidx1[i+1];
	    obs1 = vecEquivClass[idxEC1].getObs();
	    obs2 = vecEquivClass[idxEC2].getObs();
	  }
	  break;
	}
      } // end for each attribute j
    } // end for each obs

    vecEquivClass1.erase(vecEquivClass1.begin()+k+1, vecEquivClass1.end());
    sortedECidx1.resize(vecEquivClass1.size());
    for (int i=0; i<vecEquivClass1.size(); ++i) sortedECidx1[i] = i;

    if (globalPtr->args->debug>=20) ucout << "Size of vecEquivClass1: " << vecEquivClass1.size() << "\n";
    if (globalPtr->args->debug>=25) ucout << "vecEquivClass1: \n";
    if (globalPtr->args->debug>=30) for (int i=0; i<vecEquivClass1.size(); ++i)
				      cout << "EC: " << i << ": " << vecEquivClass1[i] << "\n";

  }  // end function RMASub::mergeEquivClass


  // drop some equivalence class from the initial equivalence class
  void RMASub::dropEquivClass(const int& j, const int&  al_, const int&  bu_) {

    int obs, idxEC, k=-1;	// observation number
    sortedECidx1.resize(sortedECidx.size());

    for (int i=0; i< sortedECidx.size(); i++) { // for each equivalence class
      idxEC = sortedECidx[i];
      obs = vecEquivClass[idxEC].getObs();
      // if covered, put the equiv class index to sortedECidx1
      if ( global()->data->intData[obs].X[j] >= al_
	   && global()->data->intData[obs].X[j] <= bu_ )
        sortedECidx1[++k] = idxEC;
    } // end each equivalence class

    sortedECidx1.resize(k+1); // erase extra space

  } // end function RMASub::dropEquivClass


  bool RMASub::isInSameClass(const int& obs1, const int& obs2,
			     const int& j, const int& au_, const int& bl_) {

    // if obs1.feat == obs2.feat
    if ( global()->data->intData[obs1].X[j]
	 == global()->data->intData[obs2].X[j] )
      return true;

    //  OR obs1.feat, obs2.feat in [au, bl]
    if ( au_<=bl_
	 && ( au_<= global()->data->intData[obs1].X[j]
	      && global()->data->intData[obs1].X[j]<=bl_)
	 && ( au_<= global()->data->intData[obs2].X[j]
	      && global()->data->intData[obs2].X[j]<=bl_)	)  {
      return true;
    }
    return false;
  }

  void RMASub::sortCachedCutPtByAttrib() {

    int l=-1;
    vector<vector<int> > buckets;
    buckets.resize(numAttrib());
    sortedCachedCutPts.resize(cachedCutPts.size());

    for (int k=0; k<cachedCutPts.size(); ++k)
      buckets[cachedCutPts[k].j].push_back(k);

    // walk buckets to get sorted observation list on this attribute
    for (int i=0; i<numAttrib(); ++i)
      for (int j=0; j<buckets[i].size(); ++j)
	sortedCachedCutPts[++l] = cachedCutPts[buckets[i][j]];

  }


  void RMASub::printSP(const int& j, const int& al, const int& au,
		       const int& bl, const int& bu) const {
    if (globalPtr->args->debug>=10)
      cout << "j: " << j << " (al, au, bl, bu) = ("
	   << al << ", " << au << ", " << bl << ", " << bu << ")\n" ;
  }


  void RMASub::printCurrentBounds() {
    if (globalPtr->args->debug>=10) cout << "Best local choice is " <<  _branchChoice << "\n";
    if (globalPtr->args->debug>=10)
      cout << " optFeat="     << _branchChoice.branchVar
	   << " optCutValue=" << _branchChoice.cutVal
	   << " minBound="    << _branchChoice.branch[0].exactBound << endl;
  } // end junction RMASub::printCurrentBounds


  //TODO explain what was this!!!
  void RMASub::setCutPts() {

    static int count=0;

    excCutFeat.push_back(_branchChoice.branchVar);
    excCutVal.push_back(_branchChoice.cutVal);

    int numBarnch = global()->CutPtOrders.size();	// number of branches
    if ( excCutFeat.size()==1 ) {	// if excludedList is 1, always create a new branch
      vector<CutPtOrder> row; // Create an empty row
      global()->CutPtOrders.push_back(row); // Add the row to the vector
      global()->CutPtOrders[numBarnch].push_back
               (CutPtOrder(count, _branchChoice.branchVar, _branchChoice.cutVal));
    } else {
      bool thisBranch = false;
      for (int i=0; i<numBarnch; ++i) {
	int sizeBranch = global()->CutPtOrders[i].size();
	if ((sizeBranch+1)==excCutFeat.size()) {
	  for (int k=0; k<sizeBranch; ++k) {
	    if (global()->CutPtOrders[i][k].j==excCutFeat[k]
                && global()->CutPtOrders[i][k].v==excCutVal[k]) {
	      if (k==sizeBranch-1)
		thisBranch = true;
	    } else break;
	  }
	  if (thisBranch) {
	    global()->CutPtOrders[i].push_back
                     (CutPtOrder(count, _branchChoice.branchVar, _branchChoice.cutVal));
	    break;
	  }
	}
      }
      if (!thisBranch) {
	int sizeThisBranch = excCutFeat.size();
	vector<CutPtOrder> row; // Create an empty row
	global()->CutPtOrders.push_back(row); // Add the row to the vector
	global()->CutPtOrders[numBarnch].resize(sizeThisBranch);

	for (int i=0; i<numBarnch; ++i) {
	  int sizeBranch = global()->CutPtOrders[i].size();
	  if (sizeThisBranch<=sizeBranch)
	    for (int k=0; k<sizeThisBranch-1; ++k) {
	      if (global()->CutPtOrders[i][k].j==excCutFeat[k]
                  && global()->CutPtOrders[i][k].v==excCutVal[k]) {
		if (k==sizeThisBranch-2)
		  thisBranch = true;
	      } else {
		break;
	      }
	    }
	  if (thisBranch) {
	    for (int k=0; k<sizeThisBranch-1; ++k)
	      global()->CutPtOrders[numBarnch][k].setCutPt(global()->CutPtOrders[i][k]);
	    break;
	  }
	}
	global()->CutPtOrders[numBarnch][sizeThisBranch-1].
          setCutPt(CutPtOrder(count, _branchChoice.branchVar, _branchChoice.cutVal));
      }
    }

    count++;
  }


  // ******************************************************************************
  // RMASolution methods

  rmaSolution::rmaSolution(RMA* global_) : solution(global_), global(global_) {
    DEBUGPRX(100, global, "Creating rmaSolution at " << (void*) this
	     << " with global=" << global << endl; );

    // Only one solution representation in this application,
    // so typeId can just be 0.
    typeId = 0;
  }


  rmaSolution::rmaSolution(rmaSolution* toCopy) {
    DEBUGPRX(100, global,  "Copy constructing rmaSolution at "
	     << (void*) this << " from " << toCopy << endl; );
    copy(toCopy);
    serial = ++(global->solSerialCounter);
  }


  void rmaSolution::copy(rmaSolution* toCopy) {
    solution::copy(toCopy);
    global = toCopy->global;
    a << toCopy->a;
    b << toCopy->b;
    isPosIncumb = toCopy->isPosIncumb;
  }


  void rmaSolution::printContents(ostream& outStream) {

    if(global->enumCount>1) return;

    outStream << "rectangle: a: " << a << "rectangle: b: " << b << "\n";
    cout << "rectangle: a: " << a << "rectangle: b: " << b << "\n";

    for (int i=0; i<global->numAttrib; ++i) {
      if (0<a[i])	// if lower bound changed
	cout << a[i] << "<=";
      if ( 0<a[i] || b[i]<global->distFeat[i] )
	cout << "x" << i ;
      if (b[i]<global->distFeat[i])
	cout << "<=" << b[i];
      if ( 0<a[i] || b[i]<global->distFeat[i] )
	cout << ", ";
    }
    cout << "\n";

    if ( global->args->checkObjVal() ) checkObjValue();

    if (global->args->writingCutPts()) {
      outStream << "CutPts:\n";
      for (int i=0; i<global->CutPtOrders.size(); ++i) {
	int sizeBranch = global->CutPtOrders[i].size();
	for (int j=0; j<sizeBranch-1; ++j)
	  outStream << global->CutPtOrders[i][j].order << "-"
		    << global->CutPtOrders[i][j].j
		    << "-" << global->CutPtOrders[i][j].v << "; ";
	outStream << global->CutPtOrders[i][sizeBranch-1].order << "-"
		  << global->CutPtOrders[i][sizeBranch-1].j << "-"
		  << global->CutPtOrders[i][sizeBranch-1].v << "\n";
      }  // end for each cut point
    }  // end if writeingCutPts option

  }  // end function printContents


  void const rmaSolution::printSolution() {
    ucout << (isPosIncumb) ? "Positive" : "Negative"; ucout << "\n";
    ucout << "printSolution: a: " << a << "printSolution: b: " << b << "\n";
  }


  void rmaSolution::checkObjValue() {

    int obs;
    double wt=0.0;

    for (int i=0; i<global->numDistObs; ++i) { // for each observation
      obs = global->sortedObsIdx[i];
      for (int j=0; j< global->data->numAttrib; ++j) { // for each attribute
      	if ( a[j] <= global->data->intData[obs].X[j]
	     && global->data->intData[obs].X[j] <= b[j] ) {
          // if this observation is covered by this solution
          if (j== global->data->numAttrib-1)
            wt+= global->data->intData[obs].w;
        } else break; // else go to the next observation
      }  // end for each attribute
    }  // end for each observation

    cout << "RMA ObjValue=" << wt << "\n";

  }  // end function rmaSolution::checkObjValue


  void rmaSolution::checkObjValue1(vector<int> &A, vector<int> &B,
				   vector<int> &coveredObs, vector<int> &sortedECidx) {
    int obs;  	double wt=0.0;

    for (int i=0; i<coveredObs.size(); ++i) { // for each observation
      obs = coveredObs[i];
      for (int j=0; j< global->data->numAttrib; ++j) { // for each attribute
      	if ( A[j] <= global->data->intData[obs].X[j]
	     && global->data->intData[obs].X[j] <= B[j] ) {
          // if this observation is covered by this solution
          if (j== global->data->numAttrib-1)
            wt+= global->data->intData[obs].w;
        } else break; // else go to the next observation
      }  // end for each attribute
    }  // end for each observation

    cout << "A: " << A ;
    cout << "B: " << B ;
    cout << "RMA ObjValue=" << wt << "\n";
    if ( abs(value - abs(wt)) >.00001) {
      cout << "RMA Obj Not Match! " << wt << " " << value << "\n";
      cout << "check coveredObs: " << coveredObs;
      cout << "check sortedECidx: " << sortedECidx;
    }
  }  // end function rmaSolution::checkObjValue


#ifdef ACRO_HAVE_MPI

  void rmaSolution::packContents(PackBuffer & outBuf) {
    outBuf << a << b << isPosIncumb;
  }

  void rmaSolution::unpackContents(UnPackBuffer &inBuf) {
    inBuf >> a >> b >> isPosIncumb;
  }

  int rmaSolution::maxContentsBufSize() {
    return 2*(global->data->numAttrib+1) *sizeof(int) * 1.5 ;
  }

#endif

  double rmaSolution::sequenceData() {
    if (sequenceCursor < a.size())
      return a[sequenceCursor++];
    else if (sequenceCursor < a.size() + b.size())
      return b[sequenceCursor++ - a.size()];
    else
      return isPosIncumb;
  }


} // **************** namespace pebblRMA (end) *******************

ostream& operator<<(ostream& os, pebblRMA::branchChoice& bc)  {
  os << '(' << bc.branch[0].exactBound << ',' << bc.branch[1].exactBound << ','
     << bc.branch[2].exactBound << ")-(" << bc.branch[0].whichChild << ','
     << bc.branch[1].whichChild << ',' << bc.branch[2].whichChild
     << ")-<" << bc.branchVar << ">-(" << bc.cutVal << ')';
  return os;
}

ostream& operator<<(ostream& os, pebblRMA::EquivClass& obj)  {
  obj.write(os);
  return os;
}

//*/
