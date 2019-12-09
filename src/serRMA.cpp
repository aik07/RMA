//
// serRMA.cpp
//
//  Implements larger methods in example of how to use object-oriented
//  branching framework (for RMA problems).
//
// Ai Kagawa
//

#include "serRMA.h"

using namespace utilib;
using namespace std;
using namespace pebbl;


namespace pebblRMA {
  
  
  void branchItem::set(double bound, double roundQuantum) {
    exactBound = bound;
    if (roundQuantum == 0)
      roundedBound = bound;
    else
      roundedBound = floor(bound/roundQuantum + 0.5)*roundQuantum;
    whichChild    = -1;
  }
  
  
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
        srand ( (n1+n2)*time(NULL)*100 );
        double rand_num = ( rand() % (n1+n2+1) ) / (double)(n1+n2) ;
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


  void branchChoice::setupMPIDatum(void* address, MPI_Datatype thisType,
				   MPI_Datatype* type, MPI_Aint base, MPI_Aint* disp,
				   int* blocklen, int j) {
    MPI_Address( address, &disp[j] );
    disp[j] -= base;
    type[j] = MPI_DOUBLE;
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
    MPI_Address(&example, &base);

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

    MPI_Type_struct(j, blocklen, disp, type, &mpiType);
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

  /***************** branchChoice methods *****************/

  branchChoice::branchChoice() : branchVar(MAXINT), cutVal(-1) {
    for (int i=0; i<3; i++) {
      branch[i].set(MAXDOUBLE, 1e-5);
      //branch[i].arrayPosition = i;
    }
  }

  branchChoice::branchChoice(double a, double b,
                             double c, int cut, int j) {
    branch[0].set(a, 1e-5);
    branch[1].set(b, 1e-5);
    branch[2].set(c, 1e-5);
    for (int i=0; i<3; ++i)		branch[i].whichChild=i;
    cutVal = cut;
    branchVar = j;
  }

  void branchChoice::setBounds(double a, double b, double c, int cut, int j) {
    branch[0].set(a, 1e-5); branch[1].set(b, 1e-5);	branch[2].set(c, 1e-5);
    for (int i=0; i<3; ++i) branch[i].whichChild=i;
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



  //********************************************************************************
  // RMA methods

  // RMA constructor
  RMA::RMA() :  workingSol(this), numObs(0), numDistObs(0), numAttrib(0),
		numTotalCutPts(0), numCC_SP(0),
		_perCachedCutPts(1.0), _binarySearchCutVal(false), _perLimitAttrib(1.0),
		_writeNodeTime(false), _writeCutPts(false),
		_rampUpSizeFact(1.0), _maxBoundedSP(intInf),
		_bruteForceEC(false), _bruteForceIncumb(false), _checkObjVal(false),
		_getInitialGuess(true), _countingSort(false), _branchSelection(0),
		_delta(-1), _shrinkDelta(.95), _limitInterval(inf)
		//, _vlFile(NULL)
  {

    //version_info += ", RMA example 1.1";
    min_num_required_args = 1;
    branchingInit(maximization, relTolerance, absTolerance);

    workingSol.serial = 0;
    workingSol.sense  = maximization;

    create_categorized_parameter("perCachedCutPts", _perCachedCutPts,
				 "<double>", "false", "check only cut-points from the cache"
				 "if the cache has at least x% of live cut-points out of total cut points",
				 "RMA");

    create_categorized_parameter("binarySearchCutVal", _binarySearchCutVal,
				 "<bool>", "false", "binary search cut values in each feature", "RMA");

    create_categorized_parameter("perLimitAttrib", _perLimitAttrib, "<double>",
				 "1.00", "limit number of attributes to check ", "RMA");

    create_categorized_parameter("getInitialGuess", _getInitialGuess, "<bool>",
				 "true", "enable the initial guess computation", "RMA");

    create_categorized_parameter("checkObjVal", _checkObjVal, "<bool>",
				 "false",	"check the optimal solution in the end ", "RMA");

    create_categorized_parameter("bruteForceEC", _bruteForceEC, "<bool>",
				 "false",	"brute force algorithm to create equivalence classes ", "RMA");

    create_categorized_parameter("bruteForceIncumb", _bruteForceIncumb, "<bool>",
				 "false",	"brute force algorithm to to compute incumbent in each attribute ",
				 "RMA");

    create_categorized_parameter("writeCutPts", _writeCutPts, "<bool>",
				 "false", "Write cut points chosen in the solution file ", "RMA");

    create_categorized_parameter("writeInstances", _writeInstances, "<bool>",
				 "false", "Write an input file for each weighted problem solved", "RMA");

    create_categorized_parameter("writeNodeTime", _writeNodeTime, "<bool>",
				 "false", "Write an input file for the number of B&B node and "
				 "CPU time for each iteration", "RMA");

    create_categorized_parameter("testWt", _testWt, "<bool>", "false",
				 "testing with specified test weights data, testWt.data", "RMA");

    create_categorized_parameter("maxBoundedSP", _maxBoundedSP, "<int>",
				 "intInf", "maximum number of bouneded subproblems", "RMA");

    create_categorized_parameter("rampUpSizeFact", _rampUpSizeFact, "<double>",
				 "1.00", "if (#storedCutPts) <= rampUpSizeFact * (#processors),"
				 "get out the ramp-up", "RMA");

    create_categorized_parameter("countingSort", _countingSort, "<bool>",
				 "false", "Use counting sort instead of bucket sort", "RMA");

    create_categorized_parameter("branchSelection", _branchSelection, "<int>",
				 "0", "Among tied cutpoints, 0: randomize cutpoint to select, "
				 "1: always select the first one, 2: always slect the last one", "RMA");

    create_categorized_parameter("delta", _delta, "<double>",
				 "0", "delta for recursive discretization", "RMA");

    create_categorized_parameter("shrinkDelta", _shrinkDelta, "<double>",
				 ".95", "shrink delta for recursive discretization", "RMA");

    create_categorized_parameter("limitInterval", _limitInterval, "<double>",
				 "inf", "limit Interval length of bouneded subproblems", "RMA");

    if (_bruteForceEC) _bruteForceIncumb=true;

  };   //  Constructor for RMA class


  RMA::~RMA() {
    if ( perCachedCutPts() <1 ) {
      int rank, sendbuf, recvbuf;
      //MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      rank = uMPI::rank;
      sendbuf = numCC_SP ;
      DEBUGPR(10, cout << "Local non-stron branching SP is: " << numCC_SP << "\n");

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


  bool RMA::setupProblem(int& argc,char**& argv) {

    startTime();

    if (delta()==-1) { // read already integerized dataset
      if ( !setupIntData(argc, argv) ) return false;
    } else { // read original (yet integerized) dataset
      if ( setupOrigData(argc, argv) ) {
        if ( fixedSizeBin() != -1 ) { // recursive integerization
          integerizeData();
        } else {  // fixed bin integerization
          integerizeFixedLengthData();
        }
      } else return false;
    }

    // remove duplicate data and adjust weights
    for (int j=0; j<numAttrib; ++j) bucketSort(j);
    removeDuplicateObs();

    cout << "numObs: " << numObs << endl;
    cout << "numDistObs: " << numDistObs << endl;

    for (int i=0; i<numObs; ++i)
      DEBUGPR(30, cout << "obs: " << i << ": " << intData[i] << "\n" );

    //startVerifyLogIfNeeded();

    cout << "setupProblem " ;
    endTime();

    return true;
  };


  bool RMA::setupIntData(int& argc,char**& argv) {

    int tmp, i=0;
    double _y;
    vector<double> _X;

    // read data from the data file
    if (argc <= 1) {
      cerr << "No filename specified\n";
      return false;
    }

    ifstream s(argv[1]); // open the data file

    // check whether or not the file is opened correctly
    if (!s) {
      cerr << "Could not open file \"" << argv[1] << "\"\n";
      return false;
    }

    // read number of observation and attribute in the first line of the data
    s >> numObs >> numAttrib ;

    intData.resize(numObs);
    distFeat.resize(numAttrib);
    _X.resize(numAttrib);

    while(!(s.eof())) { // while this raw past the last row

      for (int j=0; j<numAttrib; ++j) {
	s >> _X[j];
	if (s.eof()) break;
      } // end for

      s >> _y ;
      if (s.eof()) break;

      intData[i].X.resize(numAttrib);

      maxL=0;
      for (int j=0; j<numAttrib+1 ; ++j) { // for each attribute and response
	if (j<numAttrib) {
	  intData[i].X[j] = _X[j] ;
	  if ( _X[j] > distFeat[j] )  // get distinct # of feature
	    distFeat[j] = _X[j] ;
          if (maxL < distFeat[j]) maxL = distFeat[j];
	} // end if j<numAttrib
	else { // last column is the response value
	  intData[i].y = _y ;
	}
      } // end for each attribute

      // assign weight to each observation
      if (intData[i].y==1) intData[i].w = 1.0/numObs;
      else                  intData[i].w = -1.0/numObs;

      ++i; // count each row

    } // end while each row

    s.close();  // close the data file

    DEBUGPR(0, cout << "distFeat :" << distFeat << "\n");

    if (getInitialGuess()) W.resize(maxL+1);

    sortedObsIdx.resize(numObs);
    for (int i=0; i<numObs; ++i)	sortedObsIdx[i]=i;
    // compute how many cut points in this data set
    for (int j=0; j<numAttrib ; ++j) numTotalCutPts += distFeat[j];

    return true;
  }


  void RMA::setStdDevX() {

    int i, j;
    vector<double> avgX(numAttrib);
    sdX.resize(numAttrib);
    maxX.resize(numAttrib);
    minX.resize(numAttrib);

    for (j=0; j<numAttrib; ++j) {
      avgX[j]=0;
      sdX[j]=0;
      minX[j] = inf ;
      maxX[j] = -inf;
    }

    //////////////////////////////////////////////////////////////
    for (i=0; i<numObs; ++i)
      for (j=0; j<numAttrib; ++j)
	avgX[j] += origData[i].X[j];

    //////////////////////////////////////////////////////////////
    // get std dev of X in each attribute
    for (j=0; j<numAttrib; ++j) {
      avgX[j] /= numObs;
      for (i=0; i<numObs; ++i)
	sdX[j] += pow(origData[i].X[j]-avgX[j], 2);
      sdX[j] /= numObs;
      sdX[j] = sqrt(sdX[j]);
    }

  }


  bool RMA::setupOrigData(int& argc,char**& argv) {

    unsigned int i, j;
    double tmp;  string line;

    // read data from the data file
    if (argc <= 1) { cerr << "No filename specified\n"; return false;	}
    ifstream s(argv[1]); // open the data file
    // check whether or not the file is opened correctly
    if (!s) {	cerr << "Could not open file \"" << argv[1] << "\"\n"; return false; }

    // read how many columns and rows
    while (getline(s, line)) {
      if (numObs==0) {
        istringstream streamCol(line);
        while ( streamCol >> tmp ) ++numAttrib;
      }
      ++numObs;
    }
    --numAttrib; // last line is response value

    cout << "(mxn): "<< numObs << "\t" << numAttrib << "\n";

    s.clear();
    s.seekg(0, ios::beg);

    origData.resize(numObs);
    for (i=0; i<numObs; ++i) { // for each observation
      origData[i].X.resize(numAttrib);
      for (j=0; j<numAttrib; j++) // for each attribute
        s >> origData[i].X[j];
      s >> origData[i].y ;
    } // end while

    s.close();  // close the data file

    // print out original obs info
    for (int i=0; i<numObs; ++i)
      DEBUGPR(1, cout << "obs: " << i << ": " << origData[i] << "\n" );

  }


  void RMA::integerizeData() {

    double limitInterval = _limitInterval; // limitInterval();
    double delta = _delta; //delta();
    double shrinkDelta = _shrinkDelta; //.9; //shrinkDelta();

    bool isSplit, flag;
    int i, j, k, l, r, p, q, o, obs;
    double eps, eps0, interval, tmpL, tmpU, tmpL1, tmpU1, tmp1U;
    vector<double> vecTemp;
    set<double> setDistVal;
    set<double>::iterator it, itp;
    map<double, int> mapDblInt;
    map<double, int>::iterator itm;
    vector<IntMinMax> copyIntMinMax;
    distFeat.resize(numAttrib);
    vecTemp.resize(numObs);

    setStdDevX();

    for (j=0; j<numAttrib; ++j) {

      DEBUGPR(2, cout << "feat: " << j << "\n");
      setDistVal.clear();
      for (i=0; i<numObs; ++i)
	setDistVal.insert(origData[i].X[j]);

      DEBUGPR(2, cout << "setDistVal: " ;
	       for (it=setDistVal.begin(); it!=setDistVal.end(); ++it) cout << *it << " ";
	       cout << '\n');

      interval = min(4.0*sdX[j], *setDistVal.rbegin() - *setDistVal.begin()) ; // (avgX[j]+ 2*sdX[j]) - (avgX[j]-2*sdX[j])

      eps = min(delta, limitInterval) * interval ;

      eps0 = eps;
      DEBUGPR(2, cout << "delta: " << delta << "\n";
	       cout << "max: " << *setDistVal.rbegin()
	       << " min: " << *setDistVal.begin() << "\n";
	       cout << "eps: " << eps << "\n";
	       cout << "limitInterval: " << limitInterval*interval << endl);

      /************ assign integer without recursive integerization ************/
      k=0;
      mapDblInt.clear();
      vecFeature[j].vecIntMinMax.resize(setDistVal.size());
      itp = setDistVal.begin();
      vecFeature[j].vecIntMinMax[0].minOrigVal = *itp;
      vecFeature[j].vecIntMinMax[0].maxOrigVal = *itp;

      for (it=setDistVal.begin(); it!=setDistVal.end(); ++it) {
	DEBUGPR(2, cout << "tmpL: " << *itp << " tmpU: " << *it
		 << " diff: " << (*it-*itp) << endl);
        if ( (*it-*itp)>eps ) {
	  vecFeature[j].vecIntMinMax[++k-1].maxOrigVal = *(--it);
	  vecFeature[j].vecIntMinMax[k].minOrigVal     = *(++it);
	}
	itp = it;
        mapDblInt[*it] = k;
      }
      vecFeature[j].vecIntMinMax[k].maxOrigVal = *(--it);

      DEBUGPR(2, cout << "mapDblInt contains:";
	       for (itm = mapDblInt.begin(); itm != mapDblInt.end(); ++itm)
		 cout << " [" << itm->first << ':' << itm->second << ']';
	       cout << '\n');
      vecFeature[j].vecIntMinMax.resize(k+1);
      distFeat[j] = k ; // get distinct # of feature

      /************************ recursive integerization ************************/
      if (limitInterval!=inf || k!=setDistVal.size()-1) { // if there is interval limit

	copyIntMinMax.resize(k+1);
	for (i=0; i<=k; ++i) {
	  copyIntMinMax[i].minOrigVal = vecFeature[j].vecIntMinMax[i].minOrigVal;
	  copyIntMinMax[i].maxOrigVal = vecFeature[j].vecIntMinMax[i].maxOrigVal;
	}
	DEBUGPR(2, cout << endl << "vecIntMin ";
		 for (i=0; i<=k; ++i) cout << copyIntMinMax[i].minOrigVal << ' ';
		 cout << "\nvecIntMax ";
		 for (i=0; i<=k; ++i) cout << copyIntMinMax[i].maxOrigVal << ' ';
		 cout << '\n');

	p=0;
	for (i=0; i<=k; ++i) {
	  isSplit=true; eps=eps0; r=0;
	  tmpL = copyIntMinMax[i].minOrigVal;
	  tmpU = copyIntMinMax[i].maxOrigVal;
	  while ( (tmpU-tmpL) > limitInterval*interval && isSplit && eps>.0001) {
	    isSplit=false;	eps *= shrinkDelta ;
	    DEBUGPR(2, cout << "new eps: " << eps << '\n');
	    for (q=0; q<=r; ++q) {
	      l=0;
	      tmpL1 = vecFeature[j].vecIntMinMax[i+p+q].minOrigVal;
	      tmpU1 = vecFeature[j].vecIntMinMax[i+p+q].maxOrigVal;
	      DEBUGPR(2, cout << " q: " << q
		       << " tmpL2: " << tmpL1 << " tmpU2: " << tmpU1
		       << " diff: " << tmpU1 - tmpL1 << endl);
	      if ( ( tmpU1 - tmpL1 ) < 0 ) {
		DEBUGPR(2, cout << "Something Wrong!!!!!!!!!!!!!!!!!!!!!\n");
		DEBUGPR(2, cout << endl << "vecIntMin2 ";
			 for (o=0; o<=k+p; ++o) cout << vecFeature[j].vecIntMinMax[o].minOrigVal << ' ';
			 cout << "\nvecIntMax2 ";
			 for (o=0; o<=k+p; ++o) cout << vecFeature[j].vecIntMinMax[o].maxOrigVal << ' ';
			 cout << '\n');
	      } else if ( ( tmpU1 - tmpL1 ) > limitInterval*interval ) {
		isSplit=true;
		for (it=setDistVal.find(tmpL1); ; ++it) {
		  tmp1U = *it;
		  DEBUGPR(2, cout << "tmpL1: " << tmpL1 << " tmpU1: " << tmp1U
			   << " diff: " << tmp1U-tmpL1 << endl);
		  if ( ( tmp1U-tmpL1 ) > eps ) {
		    ++l; ++r; flag=true;
		    vecFeature[j].vecIntMinMax[i+p+l+q-1].maxOrigVal = tmpL1;
		    vecFeature[j].vecIntMinMax[i+p+l+q].minOrigVal   = tmp1U;
		    DEBUGPR(2, cout << " idx: " << i+p+l+q-1
			     << " tmpL4: " << vecFeature[j].vecIntMinMax[i+p+l+q-1].maxOrigVal
			     << " tmpU4: " << vecFeature[j].vecIntMinMax[i+p+l+q].minOrigVal
			     << endl);
		    DEBUGPR(2, cout << " i: " << i << "p: " << p << " r: " << r
			     << " l: " << l << " q: " << q    << endl);
		  }
		  tmpL1 = tmp1U ;
		  vecFeature[j].vecIntMinMax[i+p+l+q].maxOrigVal = tmpU;
		  if ( tmp1U==tmpU1 ) break;
		} // end for each inner sub interval
	      }  // end if each interval is less than the threthold
	    } // end for (p=0; p<=r; ++p)
	  } // end while ( (tmpU-tmpL) > limitInterval*interval && isSplit)
	  p+=r;
	  if ( (tmpU-tmpL) <= limitInterval*interval && p>0 ) {
	    vecFeature[j].vecIntMinMax[i+p].minOrigVal = copyIntMinMax[i].minOrigVal;
	    vecFeature[j].vecIntMinMax[i+p].maxOrigVal = copyIntMinMax[i].maxOrigVal;
	  }

	} // end for (i=0; i<=k; ++i), each original interval

	DEBUGPR(2, cout << endl << "vecIntMin1 ";
		 for (i=0; i<=k+p; ++i) cout << vecFeature[j].vecIntMinMax[i].minOrigVal << ' ';
		 cout << "\nvecIntMax1 ";
		 for (i=0; i<=k+p; ++i) cout << vecFeature[j].vecIntMinMax[i].maxOrigVal << ' ';
		 cout << '\n');

	o=0;
	for (it = setDistVal.begin(); it != setDistVal.end(); ++it) {
	  if ( *it > vecFeature[j].vecIntMinMax[o].maxOrigVal ) ++o;
	  mapDblInt[*it] = o;
	}

	DEBUGPR(2, cout << "mapDblInt1 contains:";
		 for (itm = mapDblInt.begin(); itm != mapDblInt.end(); ++itm)
		   cout << " [" << itm->first << ':' << itm->second << ']';
		 cout << '\n');

	vecFeature[j].vecIntMinMax.resize(k+p+1);
	distFeat[j] = k+p ; // get distinct # of feature

      } // end if recursive discretization applies

      // set intData sets
      for ( i=0; i<numObs; ++i ) {
	obs = sortedObsIdx[i];
	intData[obs].X.resize(numAttrib);
	intData[obs].X[j]	= mapDblInt[origData[obs].X[j]] ;
      }

    }	// end for (j=0; j<numAttrib; ++j) for each attribute

    for (i=0; i<numObs ; ++i) {
      obs = sortedObsIdx[i];
      DEBUGPR(20, cout << "IntObs: " << obs << ": "
	       << intData[obs] << '\n');
    }
    DEBUGPR(0, cout << "distFeat: " << distFeat << "\n");

    maxL=0;
    for (j=0; j<numAttrib ; ++j) {
      numTotalCutPts += distFeat[j];
      if ( maxL-1 < distFeat[j] ) maxL = distFeat[j]+1;
    }

    W.resize(maxL);

  } // end integerizeData


  void RMA::integerizeFixedLengthData() {
    int i,j, glMaxL=-1;
    int sizeBin = fixedSizeBin();
    maxL=0;
    // fix X matrix
    for (i=0; i<numObs; ++i) {
      for (j=0; j<numAttrib; ++j) {
	if ( origData[i].X[j] < minX[j] )
	  minX[j] = origData[i].X[j] ;  // get minX[j]
	if ( origData[i].X[j] > maxX[j] )
	  maxX[j] = origData[i].X[j] ;  // get maxX[j]
      }
    }

    distFeat.resize(numAttrib);
    for (j=0; j<numAttrib; ++j) {
      maxL=-1;
      for (int i=0; i<numObs; ++i) {
	intData[i].X.resize(numAttrib);
        intData[i].X[j] = floor ( (origData[i].X[j]-minX[j])
				  / ((maxX[j]-minX[j])/(double)sizeBin) ) ;
	if (maxL<intData[i].X[j] ) maxL = intData[i].X[j];
      }

      distFeat[j] =  maxL;
      if ( glMaxL<maxL ) glMaxL=maxL;

      vecFeature[j].vecIntMinMax.resize(maxL);
      for (int i=0; i<maxL; ++i) {
	vecFeature[j].vecIntMinMax[0].minOrigVal
	  = (double) i * ((maxX[j]-minX[j])/(double)sizeBin) + minX[j];
	vecFeature[j].vecIntMinMax[0].maxOrigVal
	  = (double) (i+1) * ((maxX[j]-minX[j])/(double)sizeBin) + minX[j];
      }
    }

    W.resize(glMaxL);

  }


  // bucket sort each attribute
  void RMA::bucketSort(const int& j) {
    int tmpNumObs, bucketIndex, l=-1;
    int bucketSize = distFeat[j]+1;
    vector<vector<int> > buckets;
    buckets.resize(bucketSize);

    if (intData.size()==0) tmpNumObs = numObs;
    else tmpNumObs = numDistObs;

    for (int i=0; i<tmpNumObs; i++) {

      if (intData.size()==0)
	bucketIndex = intData[sortedObsIdx[i]].X[j];
      else
	bucketIndex = intData[sortedObsIdx[i]].X[j];

      if (bucketIndex<0)
	cerr<<"Error @ RMA: bucketIndex<0!";
      if (bucketIndex>=bucketSize)
	cerr<<"Error @ RMA: bucketIndex>=bucketSize!";

      buckets[bucketIndex].push_back(sortedObsIdx[i]);
    }

    // walk buckets to get sorted observation list on this attribute
    for (int v=0; v<bucketSize; ++v)
      for (int k=0; k<buckets[v].size(); ++k)
	sortedObsIdx[++l] = buckets[v][k];
  }


  void RMA::removeDuplicateObs() {

    if (sortedObsIdx.size()<=0) {
      DEBUGPR(0, cout << "intData is empty.\n");
      return;
    } else if (sortedObsIdx.size()==1) {
      DEBUGPR(15, cout << "There is only one observation" << "\n");
      return;
    }

    int obs1 = sortedObsIdx[0];
    int obs2 = sortedObsIdx[1];
    int k=0, l=-1;

    for (int i=1; i<sortedObsIdx.size(); ++i) { // for each sorted, covered observation

      for (int j=0; j<numAttrib; ++j) { // for each attribute

        if ( intData[obs1].X[j] == intData[obs2].X[j] ) {
          if (j==numAttrib-1) { // if it is in the same equivalent class
            intData[obs1].w += intData[obs2].w ;
            if (i!=intData.size()-1)  // if not the last observation
              obs2=sortedObsIdx[i+1];
          }

        } else {  // detected obs1 and obs2 are in different equivClass
          sortedObsIdx[++k] = sortedObsIdx[i];
          if (i!=sortedObsIdx.size()-1) { // if not the last observation
            obs1 = sortedObsIdx[i];
            obs2 = sortedObsIdx[i+1];
          }
          break; // as soon as we detect obs1 and obs2 are in different equivClass
                 // compare the next observation combinations
        }

      } // end for each attribute j
    } // end for each obs, i

    // remove observation with zero weight
    for (int i=0; i<k+1; ++i)
      if ( intData[sortedObsIdx[i]].w !=0 )
        sortedObsIdx[++l] = sortedObsIdx[i];

    // erase extra space
    sortedObsIdx.erase(sortedObsIdx.begin()+l+1, sortedObsIdx.end());

    DEBUGPR(1, ucout << "Size of sortedObsIdx: " << sortedObsIdx.size() << "\n");

    numDistObs = sortedObsIdx.size();

  }


  /************************** for greedy heurestic *****************************/
  //*
  solution* RMA::initialGuess() {
    if (!getInitialGuess()) return NULL;

#ifdef ACRO_HAVE_MPI
    if (uMPI::rank!=0) return NULL;
#endif //  ACRO_HAVE_MPI

    guess = new rmaSolution(this);

    startTime();
    double tmpMin=inf, tmpMax=-inf, minVal=inf, maxVal=-inf;
    double optLower, optUpper, maxObjValue;
    int optAttrib=-1, oldAttrib=-1, obs, i, j;
    bool fondNewBox;

    Lmin.clear(); Lmin.resize(numAttrib,0);
    Umin.resize(distFeat.size());
    copy(distFeat.begin(), distFeat.end(), Umin.begin());
    vecCoveredObs.resize(sortedObsIdx.size());
    copy(sortedObsIdx.begin(), sortedObsIdx.end(), vecCoveredObs.begin());

    DEBUGPR(10, cout << "Lmin " << Lmin );
    DEBUGPR(10, cout << "Umin " << Umin );

    ///////////////// Minimum Range ///////////////////////
    do {
      fondNewBox = false;
      for (j=0; j<numAttrib; ++j) { // for each feature
        if ( j != oldAttrib ) { // if this attribute is not restricted
          DEBUGPR(10, cout << "for featrue: " << j << "\n" );
          setObjVec(j);
          tmpMin = getMinRange(j);
          if (tmpMin<minVal) { // && ( tmpL!=0 || tmpU != distFeat[j] ) ) {
            minVal = tmpMin; fondNewBox = true;
            optAttrib = j; optLower = tmpL; optUpper = tmpU;
            DEBUGPR(10, cout << "optAttrib: (a,b): " << optAttrib << ": "
		     << optLower << ", " << optUpper << " min: " << minVal << "\n") ;
          } // end for each cut-value
        } // end if this attribute is not restricted
      } // end each feature
      if ( fondNewBox ) {
        dropObsNotCovered(optAttrib, optLower, optUpper);
        Lmin[optAttrib] = optLower;
        Umin[optAttrib] = optUpper;
        oldAttrib = optAttrib;
        DEBUGPR(10, cout << "vecCoveredObs: " << vecCoveredObs);
        DEBUGPR(10, cout << "final optAttrib: (a,b): " << optAttrib
		 << ": " << optLower << ", " << optUpper << " min: " << minVal << "\n" );
      }
    } while( fondNewBox );

    optAttrib=-1; oldAttrib=-1;
    Lmax.clear(); Lmax.resize(numAttrib);
    Umax.resize(distFeat.size());
    copy(distFeat.begin(), distFeat.end(), Umax.begin());
    vecCoveredObs.resize(sortedObsIdx.size());
    copy(sortedObsIdx.begin(), sortedObsIdx.end(), vecCoveredObs.begin());

    /////////////////////// Maximum Range ///////////////////////
    do {
      fondNewBox = false;
      for (j=0; j<numAttrib; ++j) { // for each feature
        if ( j != oldAttrib ) { // if this attribute is not restricted
          setObjVec(j);
          DEBUGPR(10, cout << "for attribute: " << j << "tempMax: "
		   << tmpMax << " tmpL: "<< tmpL << " tmpU: "<< tmpU <<"\n");
          tmpMax = getMaxRange(j);
          if (tmpMax>maxVal) {
            maxVal = tmpMax; fondNewBox = true;
            optAttrib = j; optLower = tmpL; optUpper = tmpU;
            DEBUGPR(10, cout << "optAttrib: (a,b): " << optAttrib << ": "
		     << optLower << ", " << optUpper << " max: " << maxVal << "\n" );
          } // end for each cut-value
        } // end if this attribute is not restricted
      } // end for each attribute
      if ( fondNewBox ) {
        dropObsNotCovered(optAttrib, optLower, optUpper);
        Lmax[optAttrib] = optLower;
        Umax[optAttrib] = optUpper;
        oldAttrib = optAttrib;
        DEBUGPR(10, cout << "vecCoveredObs: " << vecCoveredObs);
        DEBUGPR(10, cout << "final optAttrib: (a,b): " << optAttrib
		 << ": " << optLower << ", " << optUpper << " max: " << maxVal << "\n" );
      }
    } while( fondNewBox );

    /////////////////////// Final Optimal Range ///////////////////////
    if (maxVal>-minVal) {
      maxObjValue=maxVal;
      guess->a.resize(numAttrib);
      guess->b.resize(numAttrib);
      copy(Lmax.begin(), Lmax.end(), guess->a.begin());
      copy(Umax.begin(), Umax.end(), guess->b.begin());
      DEBUGPR(10, cout << "chose max\n");
    } else {
      maxObjValue=-minVal;
      guess->a.resize(numAttrib);
      guess->b.resize(numAttrib);
      copy(Lmin.begin(), Lmin.end(), guess->a.begin());
      copy(Umin.begin(), Umin.end(), guess->b.begin());
      DEBUGPR(10, cout << "chose min\n");
    }
    DEBUGPR(10, cout << "Lmin: " << Lmin << "Umin: " << Umin);
    DEBUGPR(10, cout << "Lmax: " << Lmax << "Umax: " << Umax);

    guess->value = maxObjValue;	// store the current incumbent!

    // Better incumbents may have been found along the way
    //this->rampUpIncumbentSync();

    cout << "initialGuess: " << guess->value << " ";
    endTime();
    DEBUGPR(10, cout << "initial a: " << guess->a
	     << "initial b: " << guess->b);
    DEBUGPR(10, cout << "Final vecCoveredObs: " << vecCoveredObs);
    return guess;

  } // end function solution* RMA::initialGuess()
  //*/

  // get Maximum range for the feature
  double RMA::getMaxRange(const int& j) {
    int s=Lmax[j]; tmpL=Lmax[j]; tmpU=Umax[j];
    double maxEndHere=0;
    double maxSoFar=-inf;    // min so far

    for (int i=Lmax[j]; i <= Umax[j]; ++i) {
      maxEndHere += W[i] ;
      if ( maxEndHere > maxSoFar ) { maxSoFar=maxEndHere; tmpL=s; tmpU=i; }
      if ( maxEndHere<0 ) { maxEndHere=0; s=i+1;}
    }

    DEBUGPR(10, cout << "Maximum contiguous sum is " << maxSoFar );
    DEBUGPR(10, cout << " feat (L,U): " << j << " ("
	     << tmpL << ", " << tmpU << ")\n" );
    return maxSoFar;
  }


  // get Miniumum range for the feature
  double RMA::getMinRange(const int &j) {
    int s=Lmin[j]; tmpL=Lmin[j]; tmpU=Umin[j];
    double minEndHere=0;
    double minSoFar=inf;

    for (int i=Lmin[j]; i <= Umin[j]; ++i) {
      minEndHere += W[i] ;
      if (minEndHere < minSoFar) {  minSoFar=minEndHere; tmpL=s; tmpU=i; }
      if ( minEndHere>0 ) { minEndHere=0; s=i+1; }
    }

    DEBUGPR(10, cout << "Minimum contiguous sum is " << minSoFar << " ");
    DEBUGPR(10, cout << " feat (L,U): " << j << " ("
	     << tmpL << ", " << tmpU << ")\n" );
    return minSoFar;
  }


  void RMA::setObjVec(const int &j) {
    int i, v, obs;
    for (i=0; i<=distFeat[j]; ++i) W[i] = 0;
    for (i=0; i<vecCoveredObs.size(); ++i) {
      obs = vecCoveredObs[i];
      v = intData[obs].X[j];
      W[v] += intData[obs].w;
    }
    DEBUGPR(10, cout << "W: " << W << endl);
  }


  void RMA::dropObsNotCovered(const int &j, const int &lower, const int &upper) {
    int obs, l=-1;
    DEBUGPR(10, cout << "before drop: " << vecCoveredObs;);
    for (int i=0 ; i<vecCoveredObs.size(); ++i) {
      obs = vecCoveredObs[i];
      if ( lower <= intData[obs].X[j] && intData[obs].X[j] <= upper )
        vecCoveredObs[++l] = obs;   // store covered observations
    }
    vecCoveredObs.resize(l+1);  // shrink the size of vecCoveredObs
    DEBUGPR(10, cout << "after drop: " << vecCoveredObs;);
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
      os << intData[i].w << ';';

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
      DEBUGPR(10, ucout << (*it).first << " => " << (*it).second << '\n');
    }

    DEBUGPR(10, ucout << "cut point (" << j << ", " << v << ") " );
    DEBUGPR(10, ucout << (isAlreadyInCache ? "is already in cache\n" : "is new\n"));

    // if not in the hash table, insert the cut point into the hash table.
    if (!isAlreadyInCache)
      mmapCachedCutPts.insert( make_pair(j, v) ) ;
      // mmapCachedCutPts.insert( make_pair<int, int>(j, v) ) ;

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
    bl<<al;
    au<<distFeat();
    bu<<au;
  };


  void RMASub::boundComputation(double* controlParam) {

    NumTiedSols=1;

    // if exceeded the # of maximum bouneded subproblem, stop
    if ( bGlobal()->subCount[2] >= global()->maxBoundedSP() ) {
      cout << "Not expand this branch since subproblems counted is : "
           << global()->maxBoundedSP() << "\n" ;
      setState(dead);
      return;
    }

    // sort each feature based on [au, bl]
    coveredObs = global()->sortedObsIdx;
    for (int j=0; j<numAttrib(); j++) bucketSortObs(j);
    setInitialEquivClass();	// set initial equivalence class, vecEquivClass

    // if there are enough discoverd cut points (storedCutPts) check only the list
    if ( global()->perCachedCutPts() < 1.0 &&  global()->binarySearchCutVal() )
      hybridBranching();
    else if ( global()->binarySearchCutVal() )
      binaryBranching();
    else if ( global()->perCachedCutPts() < 1.0 )
      setLiveCachedCutPts();
    else   //check all cut points
      strongBranching();

    DEBUGPR(10, printCurrentBounds());

    bound = _branchChoice.branch[0].roundedBound;	// look ahead bound
    setState(bounded);

    if (_branchChoice.branch[0].roundedBound < 0) {
      DEBUGPR(10, cout << "Bound < 0. \n");
      setState(dead);
      return;
    }

    if ( _branchChoice.branchVar > numAttrib() ) {
      DEBUGPR(10, ucout << "al: " << al << "au: " << au
	      << "bl: " << bl << "bu: " << bu );
      DEBUGPR(10, cout << "branchVar > numAttrib. \n");
      setState(dead);
      return;
    }

    // If (current objValue) >= (current bound), we found the solution.
    if ( workingSol()->value >= _branchChoice.branch[0].exactBound ) {
      DEBUGPR(5, workingSol()->printSolution());
      DEBUGPR(5, cout << "Bound: " << _branchChoice.branch[0].exactBound << "\n");
      foundSolution();
      setState(dead);
      return;
    }

    //////////////////////////////////// create listExclided list (start) ////////////////////////
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
      DEBUGPR(20, ucout << "ERROR: cutValue cannot be less than 0! (cutValue="
	      << _branchChoice.cutVal << ")\n" );
      exit (EXIT_FAILURE);
    }
    if (_branchChoice.cutVal >= bu[_branchChoice.branchVar]) {
      DEBUGPR(20, ucout << "ERROR: cutValue cannot be >= bu["
	      << _branchChoice.branchVar << "]! (cutValue="
	      << _branchChoice.cutVal << ")\n" );
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
    DEBUGPR(1, ucout << "catched cut-points: ");
    while (curr!=end) {
      j = curr->first;
      v = curr->second;
      DEBUGPR(1, ucout << j << ", " << v << "\n";);
      //if (j>numAttrib() || v<0) break;
      curr++;
      if (al[j]<=v && v<bu[j]) // if v in [al, bu)
        if ( !( au[j]<bl[j] && au[j]<=v && v<bl[j] ) )   // if not overlapping
          ++numLiveCachedCutPts;
    }
    DEBUGPR(1, ucout << "\n");
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

    globalPtr = parent->global();
    branchSubAsChildOf(parent);

    // set bound
    bound = parent->_branchChoice.branch[whichChild].roundedBound;
    whichChild = parent->_branchChoice.branch[whichChild].whichChild;

    DEBUGPR(10, ucout << "Bound: " << bound << "\n") ;

    int j = parent->_branchChoice.branchVar;
    int lowerBound, upperBound;

    if (j<0) {
      DEBUGPR(20, ucout << "ERROR: feature j cannot be < 0 (j=" << j << ")\n");
      cerr << "ERROR: feature j cannot be < 0 (j=" << j << ")\n";
      return;
    }
    else if (j>numAttrib()) {
      DEBUGPR(100, ucout << "ERROR: feature j cannot be > numAttrib (j=" << j << ")\n");
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
    DEBUGPR(10, ucout << "al: " << al << "au: " << au
	    << "bl: " << bl << "bu: " << bu );
  }


  // find a particular subproblem object to the problem description embodied in the object master
  void RMASub::RMASubFromRMA(RMA* master) {
    globalPtr = master;		// set a globalPtr
    //workingSol()->value = getObjectiveVal(); // set bound value as current solution
    DEBUGPR(20,ucout << "Created blank problem, out of rmaSub:::RMASubFromRMA" << "\n");
  };


  bool RMASub::candidateSolution() {
    DEBUGPR(100, ucout << "al: " << al << "au: " << au
	    << "bl: " << bl << "bu: " << bu );
    for (int j=0; j<numAttrib() ; ++j) {
      if ( al[j] != au[j] ) return false;
      if ( bl[j] != bu[j] ) return false;
    }
    workingSol()->a << al;
    workingSol()->b << bu;

    DEBUGPR(5, workingSol()->printSolution());
    return true;
  }


  // ******************* RMASub helper functions (start) *******************************

  void RMASub::branchingProcess(const int& j, const int& v) {

    vecBounds.resize(3);
    vecBounds[2]=-1;

    if ( al[j]<=v && v < min(au[j],bl[j]) ) {

      // Middle Child
      printSP(j, al[j], v, bl[j], bu[j]);
      if (global()->bruteForceEC()) {
        coveredObs1=coveredObs;
        setEquivClassBF(j, v, bl[j]);
      } else
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
      if (global()->bruteForceEC()) {
        coveredObs1=coveredObs;
        setEquivClassBF(j, v, v+1);
      } else
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
      if (global()->bruteForceEC()) {
        coveredObs1=coveredObs;
        setEquivClassBF(j, au[j], v+1);
      } else
        mergeEquivClass(j, al[j], au[j], v+1, bu[j]);
      vecBounds[1]=getBoundMerge();

    }

    branchChoice thisChoice(vecBounds[0], vecBounds[1], vecBounds[2], v, j);

    DEBUGPR(15, ucout << "Evaluating" << thisChoice << "\n");

    // select variable based on minimum of children
    // bounds given in lexicographically decreasing order
    thisChoice.sortBounds();

    for (int i=0; i<vecBounds.size(); i++)
      if (thisChoice.branch[i].exactBound <= workingSol()->value ) {
        thisChoice.branch[i].exactBound=-1;
        thisChoice.branch[i].roundedBound=-1;
      }

    DEBUGPR(10, ucout << "Sorted version is " << thisChoice << "\n");

    if (thisChoice < _branchChoice) {
      //cout << "branchBound: " << thisChoice.branch[0].exactBound << " "
      //     << _branchChoice.branch[0].exactBound;
      _branchChoice = thisChoice;
      DEBUGPR(50, ucout << "Improves best attribute: " << j << "\n");
      DEBUGPR(10, ucout << "Branch choice now: " << _branchChoice << "\n");
      NumTiedSols=1;
      //foundBound=true;
    } else if (thisChoice == _branchChoice) {
      //cout << "branchBound: " << thisChoice.branch[0].exactBound << " "
      //     << _branchChoice.branch[0].exactBound;
      if (global()->branchSelection()==0) {
        NumTiedSols++;
        srand (NumTiedSols*time(NULL)*100);
        double rand_num = (rand() % 10001 ) / 10000.0 ;
        //DEBUGPRX(0, global(), "rand: " << rand_num  << "\n");
        //DEBUGPRX(0, global(), "rand1: " << 1.0 /  NumTiedSols << "\n");
        if ( rand_num  <= 1.0 /  NumTiedSols ) {
          _branchChoice = thisChoice;
          DEBUGPR(50, ucout << "Improves best attribute: " << j << "\n");
          DEBUGPR(10, ucout << "Branch choice now is: " << _branchChoice << "\n");
        }
      } else if (global()->branchSelection()==2) {
        _branchChoice = thisChoice;
        DEBUGPR(50, ucout << "Improves best attribute: " << j << "\n");
        DEBUGPR(10, ucout << "Branch choice now is: " << _branchChoice << "\n");
      }
    }

  } // end function void RMASub::branchingProcess


  // strong branching
  void RMASub::strongBranching() {

    int numCutPtsInAttrib;
    checkIncumbent(numAttrib()-1);

    for (int j=0; j<numAttrib(); ++j ) {

      DEBUGPR(10, ucout << "original: ");
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
      if (!global()->bruteForceEC())
        (global()->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
      checkIncumbent(j);
    } // end for each feature

  } // end RMASub::strongBranching


  // branching using cut-point caching methods
  void RMASub::cachedBranching() {

    int k=0;

    sortCachedCutPtByAttrib();
    cachedCutPts = sortedCachedCutPts;
    checkIncumbent(numAttrib()-1);

    for (int j=0; j<numAttrib(); ++j) {
      while ( k<cachedCutPts.size() ) {
        if ( j == cachedCutPts[k].j ) {
          branchingProcess(cachedCutPts[k].j, cachedCutPts[k].v);
          ++k;
        } else break;
      }

      if (j==numAttrib()-1) break;
      if (!global()->bruteForceEC())
        (global()->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
      checkIncumbent(j);
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
    checkIncumbent(numAttrib()-1);

    for (int j=0; j<numAttrib(); ++j ) {  // for each attribute

      if (distFeat()[j]<30) {
        while ( k<cachedCutPts.size() ) {
          if ( j == cachedCutPts[k].j ) {
            branchingProcess(cachedCutPts[k].j, cachedCutPts[k].v);
            ++k;
          } else break;
        }
        if (j==numAttrib()-1) break;
        if (!global()->bruteForceEC())
	  (global()->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
        checkIncumbent(j);

      } else {  // binary search

        numCutValues = bu[j]-al[j];
        if (bl[j]>au[j]) numCutValues -= bl[j]-au[j];

        if (numCutValues==0)	{// if no cutValue in this feature,
          if (!global()->bruteForceEC())
	    (global()->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
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
              DEBUGPR(10, ucout << "cutValue>=bu[j] .\n");
              break;
            }
            if ( cutValue<al[j] ) {
              DEBUGPR(10, ucout << "cutValue<al[j] .\n");
              break;
            }
            if (vecCheckedCutVal[cutValue]) {
              DEBUGPR(10, ucout << "break since oldCutValue=cutValue.\n");
              break;
            }
            vecCheckedCutVal[cutValue] = true;

            printSP(j, al[j], au[j], bl[j], bu[j]);
            DEBUGPR(10, ucout << "j: " << j << " L: "<< L << " U: " << U
		    << " cutVal: " << cutValue << "\n");

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
          DEBUGPR(10, ucout << "sortedOBS1: " << coveredObs);
	  branchingProcess(j, cutValue);

          // compare objectives instead of bounds
          //if ( vecObjValue[0] > vecObjValue[1] ) L = cutValue;
          //else U = cutValue;

          // Compare bounds
          if ( vecBounds[0] < vecBounds[1] ) L = cutValue; else U = cutValue;

        }  // end while for each cut value of feature f

        if (j==numAttrib()-1) break;

        if (!global()->bruteForceEC())
          (global()->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
        checkIncumbent(j);

      } // end binary search
    }	// end for each attribute

  } // end RMASub::hybridBranching


  // split subproblems and choose cut value by binary search
  void RMASub::binaryBranching() {

    bool firstFewCutPts;
    int l, u, L, U, cutValue, numCutValues;
    vector<bool> vecCheckedCutVal;

    checkIncumbent(numAttrib()-1);

    for (int j=0; j<numAttrib(); ++j ) {  // for each attribute

      numCutValues = bu[j]-al[j];
      if (bl[j]>au[j]) numCutValues -= bl[j]-au[j];

      if (numCutValues==0)	{// if no cutValue in this feature,
        if (!global()->bruteForceEC())
          (global()->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
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
            DEBUGPR(10, ucout << "cutValue>=bu[j] .\n");
            break;
          }
          if ( cutValue<al[j] ) {
            DEBUGPR(10, ucout << "cutValue<al[j] .\n");
            break;
          }
          if (vecCheckedCutVal[cutValue]) {
            DEBUGPR(10, ucout << "break since oldCutValue=cutValue.\n");
            break;
          }
          vecCheckedCutVal[cutValue] = true;

          printSP(j, al[j], au[j], bl[j], bu[j]);
          DEBUGPR(10, ucout << "j: " << j << " L: "<< L << " U: " << U
		  << " cutVal: " << cutValue << "\n");

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
        DEBUGPR(10, cout << "sortedOBS1: " << coveredObs);
	branchingProcess(j, cutValue);

        // compare objectives instead of bounds
        //if ( vecObjValue[0] > vecObjValue[1] ) L = cutValue;
        //else U = cutValue;

        // Compare bounds
        if ( vecBounds[0] < vecBounds[1] ) U = cutValue; else L = cutValue;

      }  // end while for each cut value of feature f

      if (j==numAttrib()-1) break;

      if (!global()->bruteForceEC())
        (global()->countingSort()) ? countingSortEC(j) : bucketSortEC(j);
      checkIncumbent(j);
    }	// end for each attribute

  } // end RMASub::binarySplitSP


  void RMASub::setLiveCachedCutPts() {

    // numLiveCachedCutPts = (# of live cut points from the cache)
    int numLiveCachedCutPts = getNumLiveCachedCutPts();

    // if numCachedCutPts is less than the percentage, check all cut points
    if ( numLiveCachedCutPts
	 < global()->numTotalCutPts * global()->perCachedCutPts() )
      strongBranching();
    else { // if not, only check the storedCutPts
      // count number of subproblems only discovering cut-points from the chache
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
      cachedBranching();
    }

    // store cached cut-points
    globalPtr->setCachedCutPts(_branchChoice.branchVar, _branchChoice.cutVal);

  }


  void RMASub::bucketSortObs(const int& j) {
    int v, l=-1;
    int size = bu[j] - al[j] + 1 - max(0, bl[j]-au[j]);
    vector<vector<int> > buckets;
    buckets.resize(size);

    for (int i=0; i<coveredObs.size(); i++) {
      v = global()->intData[coveredObs[i]].X[j];
      if ( au[j] < bl[j] ) { 	// no overlapping
        if (v<au[j]) v -= al[j];
        else if ( au[j]<=v && v<=bl[j] ) v = au[j] - al[j];
        else if ( bl[j] < v ) v = v - (bl[j] - au[j]) - al[j];
      } else v -= al[j];	// overlapping

      if (v<0) {
	DEBUGPR(10, cout  << "below covered range \n"); continue;
      } else if (v>=size) {
	DEBUGPR(10, cout << "above covered range \n"); continue;
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
      v = global()->intData[obs].X[j];
      if ( au[j] < bl[j] ) { 	// no overlapping
        if (v<au[j]) v -= al[j];
        else if ( au[j]<=v && v<=bl[j] ) v = au[j] - al[j];
        else if ( bl[j] < v ) v = v - (bl[j] - au[j]) - al[j];
      } else v -= al[j];	// overlapping

      if (v<0) {
        DEBUGPR(10, cout  << "below covered range \n"); continue;
      } else if (v>=size) {
        DEBUGPR(10, cout << "above covered range \n"); continue;
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
    v = global()->intData[obs].X[j];
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
    v = global()->intData[obs].X[j];
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
      v = global()->intData[obs].X[j];
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
      v = global()->intData[obs].X[j];
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


  void RMASub::checkIncumbent(const int& j) {
    double tmpMin=inf, tmpMax=-inf;
    double minVal = - workingSol()->value;
    double maxVal =   workingSol()->value;
    double optMinLower, optMinUpper, optMaxLower, optMaxUpper;
    int optMinAttrib=-1, optMaxAttrib=-1;
    /*
      cout << "j: " << j << "; ";
      for (int i=0; i<coveredObs.size(); ++i)
      cout << global()->intData[coveredObs[i]].X[j] << " ";
      cout << "bound (" << al[j] << ", " << au[j] << ", "
      << bl[j] << ", " << bu[j] << ")\n ";
    */
    if (!global()->bruteForceIncumb()) curObs=0;
    tmpMin = getMinRange1(j);
    if (tmpMin<minVal) { // if better min incumbent was found
      minVal = tmpMin;
      optMinAttrib = j; optMinLower = aj; optMinUpper = bj;
      DEBUGPR(10, cout << "optAttrib: (a,b): " << optMinAttrib
              << ": " << optMinLower << ", " << optMinUpper
              << " min: " << minVal << "\n") ;
    }

    if (!global()->bruteForceIncumb()) curObs=0;
    tmpMax = getMaxRange1(j);
    if (tmpMax>maxVal) { // if better max incumbent was found
      maxVal = tmpMax;
      optMaxAttrib = j; optMaxLower = aj; optMaxUpper = bj;
      DEBUGPR(10, cout << "optAttrib: (a,b): " << optMaxAttrib
              << ": " << optMaxLower  << ", " << optMaxUpper
              << " max: " << maxVal << "\n" );
    }

    if ( max(maxVal, -minVal) > workingSol()->value ) {
      workingSol()->a << al;
      workingSol()->b << bu;
      if (maxVal>-minVal) {
        workingSol()->value =maxVal;
        workingSol()->a[optMaxAttrib]=optMaxLower;
        workingSol()->b[optMaxAttrib]=optMaxUpper;
        DEBUGPR(1, cout << "positive ");
      } else  {
        workingSol()->value =-minVal;
        workingSol()->a[optMinAttrib]=optMinLower;
        workingSol()->b[optMinAttrib]=optMinUpper;
        DEBUGPR(1, cout << "negative ");
      }
      foundSolution();
      DEBUGPR(1, cout << "new incumbent  " << workingSol()->value << '\n');
      DEBUGPR(10, workingSol()->printSolution());
    }

  }


  // get Maximum range for the feature
  double RMASub::getMaxRange1(const int& j) {
    double maxEndHere, maxSoFar, tmpObj;
    int v=al[j];

    tmpObj = getObjValue(j, v);
    aj=v; bj=v; maxEndHere=tmpObj; maxSoFar=tmpObj;
    if (!global()->bruteForceIncumb())
      if (au[j]<bl[j] && au[j]<=v && v<=bl[j]) {
        bj=bl[j]; v=bl[j];
      }
    ++v;

    for (; v <= bu[j]; ++v) { // for each value in this attribute
      if (!global()->bruteForceIncumb())
        if (au[j]<bl[j] && au[j]<v && v<=bl[j]) {
          v=bl[j]; continue;
        }
      tmpObj = getObjValue(j, v);
      if ( tmpObj > maxEndHere+tmpObj ) {
        maxEndHere = tmpObj;
        if (maxEndHere > maxSoFar) {
          maxSoFar=maxEndHere; aj=v; bj=v;
          if (!global()->bruteForceIncumb())
            if (au[j]<bl[j] && au[j]<v && v<=bl[j]) bj=bl[j];
        }
      } else {
        maxEndHere += tmpObj;
        if (maxEndHere > maxSoFar) {
          maxSoFar=maxEndHere; bj=v;
          if (!global()->bruteForceIncumb())
            if (au[j]<bl[j] && au[j]<v && v<=bl[j]) bj=bl[j];
        }
      }
    } // end for each value in this attribute

    DEBUGPR(10, cout << "Maximum contiguous sum is " << maxSoFar );
    DEBUGPR(10, cout << " attribute (L,U): " << j << " ("
	    << aj << ", " << bj << ")\n" );

    return maxSoFar;
  }


  // get Miniumum range for the feature
  double RMASub::getMinRange1(const int& j) {
    double minEndHere, minSoFar, tmpObj;
    int v=al[j];

    tmpObj = getObjValue(j, v);
    aj=v; bj=v; minEndHere=tmpObj; minSoFar=tmpObj;
    if (!global()->bruteForceIncumb())
      if (au[j]<bl[j] && au[j]<=v && v<=bl[j]) {
        bj=bl[j]; v=bl[j];
      }
    ++v;

    for (; v <= bu[j]; ++v) {
      if (!global()->bruteForceIncumb())
        if (au[j]<bl[j] && au[j]<v && v<=bl[j]) {
          v=bl[j]; continue;
        }
      tmpObj = getObjValue(j, v);

      if ( tmpObj < minEndHere+tmpObj ) {
        minEndHere = tmpObj;
        if (minEndHere < minSoFar) {
          minSoFar=minEndHere; aj=v; bj=v;
          if (!global()->bruteForceIncumb())
            if (au[j]<bl[j] && au[j]<v && v<=bl[j]) bj=bl[j];
        }
      } else {
        minEndHere += tmpObj;
        if (minEndHere < minSoFar) {
          minSoFar=minEndHere; bj=v;
          if (!global()->bruteForceIncumb())
            if (au[j]<bl[j] && au[j]<v && v<=bl[j]) bj=bl[j];
        }
      }
    }

    DEBUGPR(10, cout << "Minimum contiguous sum is " << minSoFar );
    DEBUGPR(10, cout << " attribute (L,U): " << j << " ("
	    <<  aj << ", " << bj << ")\n" );
    return minSoFar;
  }


  double RMASub::getObjValue(const int& j, const int& v) {
    int obs;
    double covgWt = 0.0;
    DEBUGPR(20, ucout << "j: " << j << " v:" << v << "\n");

    // for each covered obersvation
    if (global()->bruteForceIncumb())
      for (int i=0; i<coveredObs.size(); ++i) {
        obs = coveredObs[i];
        // if the observation's jth attribute value = cut-value
        if ( global()->intData[obs].X[j] == v )
          covgWt += global()->intData[obs].w;
      }
    else {

      for (int i=curObs; i<sortedECidx.size(); ++i) {

        obs = vecEquivClass[sortedECidx[i]].getObs();

        // if the observation's jth attribute value = cut-value
        if ( global()->intData[obs].X[j] == v )
          covgWt += vecEquivClass[sortedECidx[i]].getWt();

        else if (au[j]<bl[j]
		 && au[j]<=v && v<=bl[j]
		 && au[j]<=global()->intData[obs].X[j]
		 && global()->intData[obs].X[j]<=bl[j])
          covgWt += vecEquivClass[sortedECidx[i]].getWt();

        else if (global()->intData[obs].X[j]<v) {
          DEBUGPR(0, cout << "X[j] < v! ");
          //*
          DEBUGPR(20, cout << "curObs: " << curObs << " attribute: " << j << "; "
		  << global()->intData[obs].X[j]
		  << " < cutVal: " << v << "\n");
          for (int i=0; i<coveredObs.size(); ++i)
            DEBUGPR(20, cout << global()->intData[sortedECidx[i]].X[j]
		    << " bound (" << al[j] << ", " << au[j] << ", "
		    << bl[j] << ", " << bu[j] << ")\n )");
          //*/
        } else {
          curObs = i;
          break;
        }
        if (i==sortedECidx.size()-1) curObs = sortedECidx.size();
      }  // end for each covered observation

    } // end bruteForceIncumb way or not

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
      DEBUGPR(0, cout << "coveredObs is empty.\n");
      return;
    } else if (coveredObs.size()==1) {
      DEBUGPR(15, cout << "There is only one covered observation" << "\n");
      vecEquivClass.resize(1);
      vecEquivClass[0].addObsWt(coveredObs[0],
				global()->intData[coveredObs[0]].w);
      return;
    }

    vecEquivClass.resize(coveredObs.size());

    int obs1 = coveredObs[0];
    int obs2 = coveredObs[1];
    int k=0;
    vecEquivClass[0].addObsWt(obs1, global()->intData[obs1].w);

    for (int i=1; i<coveredObs.size(); ++i) { // for each sorted, covered observation

      for (int j=0; j<numAttrib(); ++j) { // for each attribute

        if ( isInSameClass(obs1, obs2, j, au[j], bl[j]) ) {
          if (j==numAttrib()-1) { // if it is in the same equivalent class
            vecEquivClass[k].addObsWt(obs2, global()->intData[obs2].w);
            if (i!=coveredObs.size()-1)  // if not the last observation
              obs2=coveredObs[i+1];
          }

        } else {  // detected obs1 and obs2 are in different equivClass
          vecEquivClass[++k].addObsWt(obs2, global()->intData[obs2].w);
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

    DEBUGPR(1,ucout << "Size of coveredObs: " << sortedECidx.size() << "\n");
    DEBUGPR(20,ucout << "Size of vecEquivClass: " << vecEquivClass.size() << "\n");
    for (int i=0; i<vecEquivClass.size(); ++i)
      DEBUGPR(30, cout << "EC: " << i << ": " << vecEquivClass[i] << "\n" );
  }


  void RMASub::mergeEquivClass(const int& j, const int& al_, const int& au_,
			       const int& bl_, const int& bu_) {

    if ( al_!=al[j] || bu_!=bu[j] )
      dropEquivClass(j, al_, bu_);
    else
      sortedECidx1 = sortedECidx;

    vecEquivClass1.clear();

    if (sortedECidx1.size()<=0){
      DEBUGPR(0, cout << "sortedECidx1 is empty. \n");
      return;
    }

    if (sortedECidx1.size()==1){
      DEBUGPR(15, cout << "There is only one equivalence class" << "\n");
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

    DEBUGPR(20, ucout << "Size of vecEquivClass1: " << vecEquivClass1.size() << "\n");
    DEBUGPR(25, ucout << "vecEquivClass1: \n");
    DEBUGPR(30, for (int i=0; i<vecEquivClass1.size(); ++i)
		  cout << "EC: " << i << ": " << vecEquivClass1[i] << "\n" );

  }  // end function RMASub::mergeEquivClass


  // drop some equivalence class from the initial equivalence class
  void RMASub::dropEquivClass(const int& j, const int&  al_, const int&  bu_) {

    int obs, idxEC, k=-1;	// observation number
    sortedECidx1.resize(sortedECidx.size());

    for (int i=0; i< sortedECidx.size(); i++) { // for each equivalence class
      idxEC = sortedECidx[i];
      obs = vecEquivClass[idxEC].getObs();
      // if covered, put the equiv class index to sortedECidx1
      if ( global()->intData[obs].X[j] >= al_
	   && global()->intData[obs].X[j] <= bu_ )
        sortedECidx1[++k] = idxEC;
    } // end each equivalence class

    sortedECidx1.resize(k+1); // erase extra space

  } // end function RMASub::dropEquivClass


  bool RMASub::isInSameClass(const int& obs1, const int& obs2,
			     const int& j, const int& au_, const int& bl_) {

    // if obs1.feat == obs2.feat
    if ( global()->intData[obs1].X[j]
	 == global()->intData[obs2].X[j] )
      return true;

    //  OR obs1.feat, obs2.feat in [au, bl]
    if ( au_<=bl_
	 && ( au_<=global()->intData[obs1].X[j]
	      && global()->intData[obs1].X[j]<=bl_)
	 && ( au_<=global()->intData[obs2].X[j]
	      && global()->intData[obs2].X[j]<=bl_)	)  {
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


  void RMASub::setEquivClassBF(const int& j_, const int& au_, const int& bl_) {

    int obs1, obs2, _au, _bl;
    int sizeEC = 1;
    bool foundEquivClass;
    vecEquivClass1.clear();
    vecEquivClass1.resize(coveredObs1.size());
    obs1 = coveredObs1[0];
    vecEquivClass1[0].addObsWt(obs1, global()->intData[obs1].w);

    for (int i=1; i<coveredObs1.size(); ++i) { // for each sorted, covered observation
      obs2 = coveredObs1[i];
      foundEquivClass=false;
      for (int k=0; k<sizeEC; ++k) {
        obs1 = vecEquivClass1[k].getObs();
	for (int j=0; j<numAttrib(); ++j) { // for each attribute

	  if (j==j_) { _au = au_; _bl =  bl_; } //DEBUGPR(50,ucout << j << _au << _bl << "\n");
	  else {_au = au[j]; _bl = bl[j]; }

	  if ( isInSameClass(obs1, obs2, j, _au, _bl) ) {
	    if (j==numAttrib()-1) { // if it is in the same equivalent class
	      vecEquivClass1[k].addObsWt(obs2, global()->intData[obs2].w);
              //cout << "added in the same equivClass\n";
              foundEquivClass=true;
	    }
	  } else break; // go to next equivClass

	} // end for each feat
        if (foundEquivClass) break; // go to next observation
      } // end for each equivClass
      if (!foundEquivClass) {
        ++sizeEC;
        vecEquivClass1[sizeEC-1].addObsWt
	  (obs2, global()->intData[obs2].w);
        //cout << "added in the end\n";
      }
    } // end for each obs

    vecEquivClass1.resize(sizeEC);

    DEBUGPR(20, ucout << "Size of sortedObs1: " << coveredObs1.size() << "\n");
    DEBUGPR(20, ucout << "Size of vecEquivClass1: " << vecEquivClass1.size() << "\n");
    DEBUGPR(25, ucout << "vecEquivClass1: \n" );
    DEBUGPR(30, for (int i=0; i<vecEquivClass1.size(); ++i)
		  cout << "EC: " << i << ": " << vecEquivClass1[i] << "\n" );
  }  // end function RMASub::setEquivClassBC


  void RMASub::printSP(const int& j, const int& al, const int& au,
		       const int& bl, const int& bu) const {
    DEBUGPR(10, cout << "j: " << j << " (al, au, bl, bu) = ("
	    << al << ", " << au << ", " << bl << ", " << bu << ")\n") ;
  }


  void RMASub::printCurrentBounds() {
    DEBUGPR(10, ucout << "Best local choice is " <<  _branchChoice << "\n");
    DEBUGPR(20, ucout << " optFeat=" << _branchChoice.branchVar
	     << " optCutValue=" << _branchChoice.cutVal
	     << " minBound=" << _branchChoice.branch[0].exactBound << endl);
  } // end junction RMASub::printCurrentBounds


  // ******************************************************************************
  // RMASolution methods

  rmaSolution::rmaSolution(RMA* global_) : solution(global_), global(global_) {
    DEBUGPRX(100, global, "Creating rmaSolution at " << (void*) this
	     << " with global=" << global << endl);
    // Only one solution representation in this application,
    // so typeId can just be 0.
    typeId = 0;
  }


  rmaSolution::rmaSolution(rmaSolution* toCopy) {
    DEBUGPRX(100,toCopy->global,"Copy constructing rmaSolution at "
	     << (void*) this << " from " << toCopy << endl);
    copy(toCopy);
    serial = ++(global->solSerialCounter);
  }


  void rmaSolution::copy(rmaSolution* toCopy) {
    solution::copy(toCopy);
    global = toCopy->global;
    a << toCopy->a;
    b << toCopy->b;
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

    if ( global->checkObjVal() ) checkObjValue();

  }  // end function printContents


  void const rmaSolution::printSolution() {
    cout << "printSolution: a: " << a << "printSolution: b: " << b << "\n";
  }


  void rmaSolution::checkObjValue() {
    int obs;
    double wt=0.0;

    for (int i=0; i<global->numDistObs; ++i) { // for each observation
      obs = global->sortedObsIdx[i];
      for (int j=0; j<global->numAttrib; ++j) { // for each attribute
      	if ( a[j] <= global->intData[obs].X[j]
	     && global->intData[obs].X[j] <= b[j] ) {
          // if this observation is covered by this solution
          if (j==global->numAttrib-1)
            wt+=global->intData[obs].w;
        } else break; // else go to the next observation
      }  // end for each attribute
    }  // end for each observation

    cout << "RMA ObjValue=" << wt << "\n";

  }  // end function rmaSolution::checkObjValue

  /*
    solution* rmaSolution::blankClone() {
    return new rmaSolution(this);
    }
  */

#ifdef ACRO_HAVE_MPI

  void rmaSolution::packContents(PackBuffer & outBuf) {
    outBuf << a << b;
  }

  void rmaSolution::unpackContents(UnPackBuffer &inBuf) {
    inBuf >> a >> b;
  }

  int rmaSolution::maxContentsBufSize() {
    return 2*(global->numAttrib+1) *sizeof(int) * 1.5 ;
  }

#endif

  double rmaSolution::sequenceData() {
    if (sequenceCursor < a.size())
      return a[sequenceCursor++];
    else
      return b[sequenceCursor++ - a.size()];
  }


} // ******************************* namespace pebbl (end) ***********************

ostream& operator<<(ostream& os, pebblRMA::branchChoice& bc)  {
  os << '(' << bc.branch[0].exactBound << ',' << bc.branch[1].exactBound << ','
     << bc.branch[2].exactBound << ")-(" << bc.branch[0].whichChild << ','
     << bc.branch[1].whichChild << ',' << bc.branch[2].whichChild
     << ")-<" << bc.branchVar << ">-(" << bc.cutVal << ')';
  return os;
}

// Operators to read and write RMA Objs to streams
ostream& operator<<(ostream& os, pebblRMA::Data& obj) {
  obj.write(os);
  return os;
}

istream& operator>>(istream& is, pebblRMA::Data& obj) {
  obj.read(is);
  return is;
}

ostream& operator<<(ostream& os, const deque<bool>& v)  {
  os << "(";
  for (deque<bool>::const_iterator i = v.begin(); i != v.end(); ++i)
    os << " " << *i;
  os << " )\n";
  return os;
}

ostream& operator<<(ostream& os, const vector<int>& v)  {
  os << "(";
  for (vector<int>::const_iterator i = v.begin(); i != v.end(); ++i)
    os << " " << *i;
  os << " )\n";
  return os;
}

ostream& operator<<(ostream& os, pebblRMA::EquivClass& obj)  {
  obj.write(os);
  return os;
}
