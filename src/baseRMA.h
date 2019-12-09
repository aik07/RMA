/**********************************************************
 * File name:   argRMA.h
 * Author:      Ai Kagawa
 * Description: a header file for RMA argument class
**********************************************************/

#ifndef BASE_RMA_h
#define BASE_RMA_h

#include "argRMA.h"
//*
#include <limits>
#include <string>

#include <pebbl_config.h>
#include <pebbl/utilib/ParameterSet.h>
#include <pebbl/utilib/ParameterList.h>
#include <pebbl/utilib/CommonIO.h>
#include <pebbl/utilib/memdebug.h>
#include <pebbl/utilib/seconds.h>
#include <pebbl/bb/branching.h>
#include <pebbl/bb/pebblParams.h>
#include <pebbl/pbb/parPebblParams.h>
//*/

using namespace std;
using namespace utilib;
using namespace pebbl;
using namespace arg;

namespace base {

class BaseRMA : public ArgRMA,
                virtual public pebblParams,
                virtual public parallelPebblParams {

public:

  BaseRMA(): parameters_registered(false), min_num_required_args(0) { }

  virtual ~BaseRMA() {}

  bool   setup(int& argc, char**& argv);

  // Parameter-related ethods
  void   write_usage_info(char const* progName, ostream& os) const;
  void   writeCommandUsage(char const* progName, ostream& os) const;
  bool   processParameters(int& argc, char**& argv,
                           unsigned int min_num_required_args_=0);

  // Register the parameters into a ParameterList object
  void   register_parameters() { plist.register_parameters(*this); }

  /// Check parameters for setup problems and perform debugging I/O
  bool   checkParameters(char const* progName = "");

  bool   setupProblem(int argc, char** argv) { true; }

  virtual void setName(const char* cname);

 //////////////////////////////////////////////////////////////////
 
private:
  ParameterList plist;
  bool          parameters_registered;
  string        problemName;
  string        solver_name;
  unsigned int  min_num_required_args;

};

} // namespace base

#endif // BASE_RMA_h
