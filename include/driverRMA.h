/**********************************************************
* File name:   driverRMA.h
* Author:      Ai Kagawa
* Description: a header file for RMA driver class
***********************************************************/

#ifndef RMA_h
#define RMA_h

#include <iostream>

#include <pebbl_config.h>
#include <pebbl/utilib/ParameterList.h>
#include <pebbl/utilib/memdebug.h>
#include <pebbl/utilib/seconds.h>
#include <pebbl/utilib/CommonIO.h>
#include <pebbl/bb/pebblParams.h>
#include <pebbl/pbb/parPebblParams.h>

#include "dataRMA.h"
#include "serRMA.h"
#include "greedyRMA.h"

#ifdef ACRO_HAVE_MPI
  #include <pebbl/pbb/parBranching.h>
  #include "parRMA.h"
  #define outstream ucout
  #define IO(action) if (uMPI::iDoIO) { CommonIO: end_tagging(); action; }
#else
  typedef void parRMA;
  #define outstream cout;
  #define IO(action) action;
#endif


namespace pebblRMA {

//static double inf = numeric_limits<double>::infinity();
//static int intInf = numeric_limits<int>::max();

using namespace utilib;
using namespace argRMA;
using namespace greedyRMA;


class DriverRMA : public ArgRMA,
                  virtual public pebbl::pebblParams,
                  virtual public pebbl::parallelPebblParams {

public:

  DriverRMA(int argc, char** arg);

  ~DriverRMA() {
    #ifdef ACRO_HAVE_MPI
      if (parallel) { CommonIO::end(); uMPI::done(); }
    #endif // ACRO_HAVE_MPI
  }

  void solveRMA();

  bool setup(int& argc, char**& argv);

  // Parameter-related methods

	void write_usage_info(char const* progName,std::ostream& os) const;

  void writeCommandUsage(char const* progName,std::ostream& os) const;

  bool processParameters(int& argc, char**& argv,
         unsigned int min_num_required_args);

  /// Register the parameters into a ParameterList object
	void register_parameters() { plist.register_parameters(*this); }
  //void register_parameters() { plist.register_parameters(args); }

  /// Check parameters for setup problems and perform debugging I/O
  bool checkParameters(char const* progName = "");

	bool setupProblem(int argc, char** argv) { true; }

	virtual void setName(const char* cname);

//private:
  bool          parallel;

  //ArgRMA*       args;
  Data*         data;
  GreedyRMA*    grma;
  RMA*          rma ;
  parRMA*       prma;

  ParameterList plist;
	bool          parameters_registered;
	string        problemName;
	string        solver_name;
	unsigned int  min_num_required_args;

	Time          tc;
	double        wallTime;
  double        cpuTime;

};

} // end namespace pebblRMA

#endif // RMA_h
