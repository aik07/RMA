/*
 *  File name:   Base.cpp
 *  Author:      Ai Kagawa
 *  Description: a source file for LPBase and Data classes
 */

#include "base.h"


namespace base {

  ////////////////////// Base class methods ////////////////////////

  // Standard serial read-in code.  Returns true if we can continue, false if
  // we have to bail out.
  bool Base::setup(int& argc, char**& argv) {

    if (!processParameters(argc,argv,min_num_required_args))
      return false;

    if (plist.size() == 0) {
        ucout << "Using default values for all solver options" << std::endl;
    } else {
      ucout << "User-specified solver options: " << std::endl;
      plist.write_parameters(ucout);
      ucout << std::endl;
    }

    set_parameters(plist,false);

    if ((argc > 0) && !checkParameters(argv[0]))
      return false;

    if (!setupProblem(argc,argv))
      return false;

    if (plist.unused() > 0) {
      ucout << "\nERROR: unused parameters: " << std::endl;
      plist.write_unused_parameters(ucout);
      ucout << utilib::Flush;
      return false;
    }

    return true;

  }


  bool Base::processParameters(int& argc, char**& argv,
  				  unsigned int min_num_required_args__) {

    if (argc > 0)
       solver_name = argv[0];
    else
       solver_name = "unknown";
    if (!parameters_registered) {
      register_parameters();
      parameters_registered=true;
    }

    plist.process_parameters(argc,argv,min_num_required_args__);

    // Set the name of the problem to be the last thing on the command
    // line. setName will extract the filename root. The setupProblem
    // method can overwrite this later.
    if ((argc > 1) && (argv[argc-1] != NULL))
      setName(argv[1]);

    return true;
  }


  bool Base::checkParameters(char const* progName) {

    if (help_parameter) {
      write_usage_info(progName,cout);
      return false;
    }

    if (debug_solver_params) {
      ucout << "---- LPBoost Parameters ----" << endl;
      write_parameter_values(ucout);
      ucout << endl << utilib::Flush;
    }

    return true;
  }


  void Base::write_usage_info(char const* progName,std::ostream& os) const {
    writeCommandUsage(progName,os);
    os << endl;
    plist.write_registered_parameters(os);
    os << endl;
  }


  void Base::writeCommandUsage(char const* progName,std::ostream& os) const {
    os << "\nUsage: " << progName << " { --parameter=value ... }";
    if (min_num_required_args == 1)
      os << " <problem data file>";
    else if (min_num_required_args == 1)
      os << " <" << min_num_required_args << " problem data files>";
    os << endl;
  }


  // This sets the official name of the problem by chewing up the
  // filename.  It can be overridden.  This version just finds the last
  // "/" or "\" in the name and removes it and everything before it.

  void Base::setName(const char* cname) {
  #if defined (TFLOPS)
    problemName = cname;
    int i=problemName.size();
    while (i >= 0) {
      if (cname[i] == '/') break;
      i--;
    }
    if (i >= 0)
       problemName.erase(0,i+1);
    // TODO: remove the .extension part for this case
  #else
    problemName = cname;
    size_type i = problemName.rfind("/");
    if (i == string::npos)
      i = problemName.rfind("\\");
    if (i != string::npos)
      problemName.erase(0,i+1);

    size_type n = problemName.length();

    if (n < 4)
      return;

    string endOfName(problemName,n-4,4);
    if ((endOfName == ".dat") || (endOfName == ".DAT"))
      problemName.erase(n-4,n);
    if ((endOfName == ".data") || (endOfName == ".DATA"))
        problemName.erase(n-5,n);
  #endif
  }


  ////////////////////// Data class methods ////////////////////////

  bool Data::readData(int argc, char** argv) {

    unsigned int i, j;
    double tmp;  string line;

    //setup(argc, argv); // setup function in Base class

    tc.startTime();

    numOrigObs=0;
    numAttrib=0;

    // read data from the data file
    if (argc <= 1) { cerr << "No filename specified\n"; return false;	}
    ifstream s(argv[1]); // open the data file
    // check whether or not the file is opened correctly
    if (!s) {	cerr << "Could not open file \"" << argv[1] << "\"\n"; return false; }

    // read how many columns and rows
    while (getline(s, line)) {
      if (numOrigObs==0) {
        istringstream streamCol(line);
        while ( streamCol >> tmp ) ++numAttrib;
      }
      ++numOrigObs;
    }
    --numAttrib; // last line is response value

  #ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
  #endif //  ACRO_HAVE_MPI
    cout << "(mxn): "<< numOrigObs << "\t" << numAttrib << "\n";
  #ifdef ACRO_HAVE_MPI
  	}
  #endif //  ACRO_HAVE_MPI

    s.clear();
    s.seekg(0, ios::beg);

    origData.resize(numOrigObs);
    for (i=0; i<numOrigObs; ++i) { // for each observation
      origData[i].X.resize(numAttrib);
      for (j=0; j<numAttrib; j++) // for each attribute
        s >> origData[i].X[j];
      s >> origData[i].y ;
    } // end while

    if (isLPBoost()) {
      for (i=0; i<numOrigObs; ++i)   // for each observation
        if (origData[i].y==0) origData[i].y=-1;
      for (i=0; i<numOrigObs; ++i)   // for each observation
        if (origData[i].y!=-1 && origData[i].y!=1)
          cerr << "y:" << origData[i].y
            << " Data contains not only -1(0) or +1!" << '\n';

    }

    s.close();  // close the data file

    // print out original obs info
    for (int i=0; i<numOrigObs; ++i)
      DEBUGPR(1, ucout << "obs: " << i << ": " << origData[i] << "\n" );

    if (readShuffledObs())
      readRandObs(argc, argv);

    DEBUGPRX(2, this, "setupProblem: \n");
    DEBUGPRX(2, this, tc.endCPUTime());
    DEBUGPRX(2, this, tc.endWallTime());

    setDataDimensions();

    return true;

  }


  // read shuffled observation from the data file
  bool Data::readRandObs(int argc, char** argv) {

  	ucout << "Use Shuffled Obs\n";

  	ifstream s(argv[2]); // open the data file

  	// check whether or not the file is opened correctly
  	if (!s) {	cerr << "Could not open file \"" << argv[2] << "\"\n"; return false; }

    vecRandObs.resize(numOrigObs);

  	// read data
  	for (int i=0; i < numOrigObs; ++i) s >> vecRandObs[i];

  	s.close();  // close the data file

    DEBUGPRX(2, this, "vecRandObs: " << vecRandObs);
  	return true;

  }


  void Data::setDataDimensions() {
    intData.resize(numOrigObs);
    distFeat.resize(numAttrib);
    vecFeature.resize(numAttrib);
    if (isREPR()) standData.resize(numOrigObs);
  }


  void Data::writeIntObs() {
    int i,j, obs;
    stringstream s;
    s << "int" << '.' ;
    ofstream os(s.str().c_str());
    for ( i=0; i<numTrainObs; ++i ) {
      for ( j=0; j<numAttrib; ++j ) {
        obs = vecTrainData[i];
        os << intData[obs].X[j] << " " ;
      }
      os << "\n";
    }
    os.close();
  }


  void Data::writeOrigObs() {
    int i,j, obs;
    stringstream s;
    s << "orig" << '.' ;
    ofstream os(s.str().c_str());
    for ( i=0; i<numTrainObs; ++i ) {
      for ( j=0; j<numAttrib; ++j ) {
        obs = vecTrainData[i];
        os << origData[obs].X[j] << " " ;
      }
      os << "\n";
    }
    os.close();
  }


  void Data::setXStat() {

    int i, j, obs;
    avgX.resize(numAttrib);
    sdX.resize(numAttrib);

    for (j=0; j<numAttrib; ++j) { avgX[j]=0, sdX[j]=0;}

    for (i=0; i<numTrainObs; ++i) {
      obs = vecTrainData[i];
      for (j=0; j<numAttrib; ++j) {
        avgX[j] += origData[obs].X[j];
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    // get std dev of X in each attribute
    for (j=0; j<numAttrib; ++j) {
      avgX[j] /= numTrainObs;		// get avgX  for each attribute
      for (i=0; i<numTrainObs; ++i) {
        obs = vecTrainData[i];
        sdX[j] += pow(origData[obs].X[j]-avgX[j], 2);
      }
      sdX[j] /= numTrainObs;
      sdX[j] = sqrt(sdX[j]);
    }

  }


  void Data::integerizeData() {

    bool isSplit, flag;
    int i, j, k, l, r, p, q, o, obs;
    double tmpL, tmpU, tmpL1, tmpU1, tmp1U;

    double interval;  // confidence interval range
    double eps, eps0; // episilon, aggregation level

    vector<double> vecTemp(numOrigObs);

    set<double> setDistVal;         // a set continas all distinct values for each attribute
    set<double>::iterator it, itp;  // iterator for the set

    map<double, int> mapDblInt;     // a container maps from an original value to an integeried value
    map<double, int>::iterator itm; // iterator for the map

    vector<IntMinMax> copyIntMinMax;  // a vector contins min and max for each integerized value

    tc.startTime();

    if (isLPBoost()) setXStat();

    for (j=0; j<numAttrib; ++j) {

      DEBUGPRX(2, this, "feat: " << j << "\n");
      setDistVal.clear();
      for (i=0; i<numTrainObs; ++i) {
        obs = vecTrainData[i];
        setDistVal.insert(origData[obs].X[j]);
      }

      DEBUGPRX(2, this, "setDistVal: " ;
        for (it=setDistVal.begin(); it!=setDistVal.end(); ++it)
          cout << *it << " ";
        cout << '\n');

      // get 95% confidence interval range
      interval = min(4.0*sdX[j], *setDistVal.rbegin() - *setDistVal.begin()) ;

      // episiolon, aggregation level, for integerization
      eps = min(getDelta(), getLimitInterval()) * interval ;


      eps0 = eps;
      DEBUGPRX(2, this, "delta: " << getDelta() << "\n"
        << "max: " << *setDistVal.rbegin()
        << " min: " << *setDistVal.begin() << "\n"
        << "eps: " << eps << "\n"
        << "limitInterval: " << getLimitInterval()*interval << endl);

      /************ assign integer without recursive integerization ************/
      k=0;
      mapDblInt.clear();
      vecFeature[j].vecIntMinMax.resize(setDistVal.size());

      // the min value is equal to the maximum value for each integer assigned
      itp = setDistVal.begin();
      vecFeature[j].vecIntMinMax[0].minOrigVal = *itp;
      vecFeature[j].vecIntMinMax[0].maxOrigVal = *itp;

      // walk thorugh the set of distincet values
      // some value can be aggregated by the level of the episilon
      for (it=setDistVal.begin(); it!=setDistVal.end(); ++it) {
        DEBUGPRX(2, this, "tmpL: " << *itp << " tmpU: " << *it
          << " diff: " << (*it-*itp) << endl);
        if ( (*it-*itp)>eps ) { // aggregating some value
          vecFeature[j].vecIntMinMax[++k-1].maxOrigVal = *(--it);
          vecFeature[j].vecIntMinMax[k].minOrigVal     = *(++it);
        }
        itp = it;
        mapDblInt[*it] = k;
      }
      vecFeature[j].vecIntMinMax[k].maxOrigVal = *(--it);

      DEBUGPRX(2, this, "mapDblInt contains:";
        for (itm = mapDblInt.begin(); itm != mapDblInt.end(); ++itm)
          cout << " [" << itm->first << ':' << itm->second << ']';
        cout << '\n');
      vecFeature[j].vecIntMinMax.resize(k+1);
      distFeat[j] = k ; // get distinct # of feature

      /************************ recursive integerization ************************/

      // if there is interval limit
      // and the size of distince value is not same as the number of integers assigned
      if (getLimitInterval()!=inf || k!=setDistVal.size()-1) {

        copyIntMinMax.resize(k+1);
        for (i=0; i<=k; ++i) {
          copyIntMinMax[i].minOrigVal = vecFeature[j].vecIntMinMax[i].minOrigVal;
          copyIntMinMax[i].maxOrigVal = vecFeature[j].vecIntMinMax[i].maxOrigVal;
        }

        DEBUGPRX(2, this, "\nvecIntMin ";
          for (i=0; i<=k; ++i)
            cout << copyIntMinMax[i].minOrigVal << ' ';
          cout << "\nvecIntMax ";
          for (i=0; i<=k; ++i)
            cout << copyIntMinMax[i].maxOrigVal << ' ';
          cout << '\n');

        p=0;
        for (i=0; i<=k; ++i) {

          isSplit=true;
          eps=eps0;
          r=0;

          tmpL = copyIntMinMax[i].minOrigVal;
          tmpU = copyIntMinMax[i].maxOrigVal;

          while ( (tmpU-tmpL) > getLimitInterval()*interval
                  && isSplit && eps>.0001) {

            isSplit=false;
            eps *= shrinkDelta() ;

            DEBUGPRX(2, this, "new eps: " << eps << '\n');

            for (q=0; q<=r; ++q) {
              l=0;
              tmpL1 = vecFeature[j].vecIntMinMax[i+p+q].minOrigVal;
              tmpU1 = vecFeature[j].vecIntMinMax[i+p+q].maxOrigVal;

              DEBUGPRX(2, this, " q: " << q
                 << " tmpL2: " << tmpL1 << " tmpU2: " << tmpU1
                 << " diff: " << tmpU1 - tmpL1 << endl);

              if ( ( tmpU1 - tmpL1 ) < 0 ) {

                DEBUGPRX(2, this, "Something Wrong!!!!!!!!!!!!!!!!!!!!!\n");

                DEBUGPRX(2, this, endl << "vecIntMin2 ";
                  for (o=0; o<=k+p; ++o)
                    cout << vecFeature[j].vecIntMinMax[o].minOrigVal << ' ';
                  cout << "\nvecIntMax2 ";
                  for (o=0; o<=k+p; ++o)
                    cout << vecFeature[j].vecIntMinMax[o].maxOrigVal << ' ';
                  cout << '\n');

              } else if ( ( tmpU1 - tmpL1 ) > getLimitInterval()*interval ) {

                isSplit=true;

                for (it=setDistVal.find(tmpL1); ; ++it) {
                  tmp1U = *it;
                  DEBUGPRX(2, this,
                    "tmpL1: " << tmpL1 << " tmpU1: " << tmp1U
                    << " diff: " << tmp1U-tmpL1 << endl);

                  if ( ( tmp1U-tmpL1 ) > eps ) {
                    ++l; ++r; flag=true;
                    vecFeature[j].vecIntMinMax[i+p+l+q-1].maxOrigVal = tmpL1;
                    vecFeature[j].vecIntMinMax[i+p+l+q].minOrigVal   = tmp1U;
                    DEBUGPRX(2, this, " idx: " << i+p+l+q-1
                         << " tmpL4: " << vecFeature[j].vecIntMinMax[i+p+l+q-1].maxOrigVal
                         << " tmpU4: " << vecFeature[j].vecIntMinMax[i+p+l+q].minOrigVal
                         << endl);
                    DEBUGPRX(2, this, " i: " << i << "p: " << p << " r: " << r
                             << " l: " << l << " q: " << q    << endl);
                  }

                  tmpL1 = tmp1U ;
                  vecFeature[j].vecIntMinMax[i+p+l+q].maxOrigVal = tmpU;

                  if ( tmp1U==tmpU1 ) break;

                } // end for each inner sub interval
              }  // end if each interval is less than the threthold
            } // end for (p=0; p<=r; ++p)
          } // end while ( (tmpU-tmpL) > getLimitInterval()*interval && isSplit)

          p+=r;

          if ( (tmpU-tmpL) <= getLimitInterval()*interval && p>0 ) {
            vecFeature[j].vecIntMinMax[i+p].minOrigVal = copyIntMinMax[i].minOrigVal;
            vecFeature[j].vecIntMinMax[i+p].maxOrigVal = copyIntMinMax[i].maxOrigVal;
          }

        } // end for (i=0; i<=k; ++i), each original interval

        DEBUGPRX(2, this, "\nvecIntMin1 ";
          for (i=0; i<=k+p; ++i)
            cout << vecFeature[j].vecIntMinMax[i].minOrigVal << ' ';
          cout << "\nvecIntMax1 ";
          for (i=0; i<=k+p; ++i)
            cout << vecFeature[j].vecIntMinMax[i].maxOrigVal << ' ';
          cout << '\n');

        o=0;
        for (it = setDistVal.begin(); it != setDistVal.end(); ++it) {
          if ( *it > vecFeature[j].vecIntMinMax[o].maxOrigVal ) ++o;
          mapDblInt[*it] = o;
        }

        DEBUGPRX(2, this, "mapDblInt1 contains:";
          for (itm = mapDblInt.begin(); itm != mapDblInt.end(); ++itm)
            cout << " [" << itm->first << ':' << itm->second << ']';
          cout << '\n');

        vecFeature[j].vecIntMinMax.resize(k+p+1);
        distFeat[j] = k+p ; // get distinct # of feature

      } // end if recursive discretization applies

      // set intData sets
      for ( i=0; i<numTrainObs; ++i ) {
        obs = vecTrainData[i];
        intData[obs].X.resize(numAttrib);
        intData[obs].X[j]	= mapDblInt[origData[obs].X[j]] ;
      }

    }	// end for each attribute

    /*
    for (i=0; i<numTrainObs ; ++i) {
      obs = vecTrainData[i];
      DEBUGPRX(20, this, "IntObs: " << obs << ": "
        << intData[obs] << '\n');
    }*/

    DEBUGPRX(1, this, "distFeat: " << distFeat << "\n");

  #ifdef ACRO_HAVE_MPI
    if (uMPI::rank==0) {
  #endif //  ACRO_HAVE_MPI
    if (writePred()) {
      writeIntObs();
      writeOrigObs();
    }
  #ifdef ACRO_HAVE_MPI
    }
  #endif //  ACRO_HAVE_MPI

    maxL=0;
    for (j=0; j<numAttrib ; ++j) {
      numTotalCutPts += distFeat[j];
      if ( maxL-1 < distFeat[j] )
        maxL = distFeat[j]+1;
    }

    ////////////////////////////////////////////////////////////////////////////
    for (int i=0; i<numTrainObs; ++i) {
      obs = vecTrainData[i];
      DEBUGPRX(1, this, "obs: " << obs << ": "
        << intData[obs] << "\n" );
    }

    DEBUGPRX(1, this, "integerizeProblem: \t" );
    DEBUGPRX(1, this, tc.endCPUTime());
    DEBUGPRX(2, this, tc.endWallTime());

  } // end integerizeData


  void Data::setStandData(){

    int i, j, obs;
    avgY=0,
    sdY=0;
    avgX.resize(numAttrib);
    sdX.resize(numAttrib);
    maxX.resize(numAttrib);
    minX.resize(numAttrib);

    for (j=0; j<numAttrib; ++j) {
      avgX[j]=0;
      sdX[j]=0;
      minX[j] = inf ;
      maxX[j] = -inf;
    }

    ////////////////////////////////////////////////////////////////////////////
    for (i=0; i<numTrainObs; ++i) {
      obs = vecTrainData[i];
      avgY += origData[obs].y ;
      for (j=0; j<numAttrib; ++j) {
        avgX[j] += origData[obs].X[j];
      }
    }

    avgY /= numTrainObs; // get average response value

    // get std dev of y in each attribute
    for (i=0; i<numTrainObs; ++i) {
      obs = vecTrainData[i];
      sdY += pow(origData[obs].y-avgY, 2);
    }
    sdY /= numTrainObs;  sdY = sqrt(sdY);

    // standardize response value
    for (i=0; i<numTrainObs; ++i) {
      obs = vecTrainData[i];
      standData[obs].y = (origData[obs].y - avgY) / sdY ;
      standData[obs].X.resize(numAttrib);
    }

    ////////////////////////////////////////////////////////////////////////////
    // get std dev of X in each attribute
    for (j=0; j<numAttrib; ++j) {
      avgX[j] /= numTrainObs;
      for (i=0; i<numTrainObs; ++i) {
        obs = vecTrainData[i];
        sdX[j] += pow(origData[obs].X[j]-avgX[j], 2);
      }
      sdX[j] /= numTrainObs; sdX[j] = sqrt(sdX[j]);
    }

    // standardize X in each attribute
    for (j=0; j<numAttrib; ++j)
      for (i=0; i<numTrainObs; ++i) {
        obs = vecTrainData[i];
        standData[obs].X[j]
          = (origData[obs].X[j] - avgX[j]) / sdX[j] ;
      }

    ////////////////////////////////////////////////////////////////////////////
    for (int i=0; i<numTrainObs; ++i) {
      obs = vecTrainData[i];
      DEBUGPRX(1, this, "obs: " << obs << ": "
        << standData[obs] << "\n" );
    }

  }

  void Data::setPosNegObs() {

    numPosTrainObs=0;
    numNegTrainObs=0;
    for (int i=0; i<numTrainObs; ++i) {
      if (origData[vecTrainData[i]].y==1) ++numPosTrainObs;
      else                                ++numNegTrainObs;
    }
    ucout << "m^+ m^-: " << numPosTrainObs << "\t" << numNegTrainObs << "\n";

  }

  //  integerize into fixed bin
  void Data::integerizeFixedLengthData() {

    int i,j, obs, glMaxL=-1;
    int sizeBin = fixedSizeBin();
    maxL=0;

    // fix X matrix
    for (i=0; i<numTrainObs; ++i) {
      obs = vecTrainData[i];
      for (j=0; j<numAttrib; ++j) {
        if ( standData[obs].X[j] < minX[j] )
          minX[j] = standData[obs].X[j] ;  // get minX[j]
        if ( standData[obs].X[j] > maxX[j] )
          maxX[j] = standData[obs].X[j] ;  // get maxX[j]
      }
    }

    distFeat.resize(numAttrib);
    for (j=0; j<numAttrib; ++j) {
      maxL=-1;
      for (int i=0; i<numTrainObs; ++i) {
        obs = vecTrainData[i];
        intData[obs].X.resize(numAttrib);
        intData[obs].X[j] = floor ( (origData[obs].X[j]-minX[j])
                                / ((maxX[j]-minX[j])/(double)sizeBin) ) ;
        if (maxL<intData[obs].X[j] ) maxL = intData[obs].X[j];
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

  }


} // namespace boosting


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

ostream& operator<<(ostream& os, const vector<double>& v)  {
	os << "(";
	for (vector<double>::const_iterator i = v.begin(); i != v.end(); ++i)
		os << " " << *i;
	os << " )\n";
	return os;
}

ostream& operator<<(ostream& os, const vector<vector<int> >& v)  {
	for (int i=0; i<v.size(); ++i) {
	  os << "(";
	  for (int j=0; j<v[i].size(); ++j) os <<  v[i][j] << " ";
	  os << " )\n";
	}
	return os;
}

ostream& operator<<(ostream& os, const vector<vector<double> >& v)  {
	for (int i=0; i<v.size(); ++i) {
	  os << "(";
	  for (int j=0; j<v[i].size(); ++j) os <<  v[i][j] << " ";
	  os << " )\n";
	}
	return os;
}

// Operators to read and write RMA Objs to streams
ostream& operator<<(ostream& os, base::DataXy& obj) {
	obj.write(os);
	return os;
}

istream& operator>>(istream& is, base::DataXy& obj) {
	obj.read(is);
	return is;
}

// Operators to read and write RMA Objs to streams
ostream& operator<<(ostream& os, base::DataXw& obj) {
	obj.write(os);
	return os;
}

istream& operator>>(istream& is, base::DataXw& obj) {
	obj.read(is);
	return is;
}
