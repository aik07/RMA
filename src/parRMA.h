/**********************************************************
*  Author:      Ai Kagawa
*  Description: a header file for the parallel RMA solver using PEBBL
**********************************************************/


#ifndef pebbl_paralleRMA_h
#define pebbl_paralleRMA_h

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <utility>
#include <map>
#include <stack>
#include <cmath>
#include <algorithm>    // std::in
#include <mpi.h>

#include <pebbl_config.h>

#ifdef ACRO_HAVE_MPI

#include <pebbl/utilib/logEvent.h>
#include <pebbl/utilib/stl_auxiliary.h>
#include <pebbl/utilib/exception_mngr.h>
#include <pebbl/utilib/comments.h>
#include <pebbl/utilib/mpiUtil.h>
#include <pebbl/utilib/std_headers.h>
#include <pebbl/utilib/PackBuf.h>
#include <pebbl/utilib/BitArray.h>
#include <pebbl/comm/mpiComm.h>
#include <pebbl/comm/coTree.h>
#include <pebbl/comm/outBufferQ.h>
#include <pebbl/sched/ThreadObj.h>
#include <pebbl/sched/SelfAdjustThread.h>

#include <pebbl/pbb/parBranching.h>
#include <pebbl/pbb/packedSolution.h>
#include <pebbl/pbb/parPebblBase.h>

#include "serRMA.h"
#include "utility.h"


using namespace std;
using namespace pebbl;
using namespace utilib;


namespace pebblRMA {


  class RMA;
  class RMASub;
  class CutPtThd;


  //**************************************************************************
  //  The parallel branching class...
  class parRMA : virtual public parallelBranching, virtual public RMA {

  public:

    parRMA(MPI_Comm comm_ = MPI_COMM_WORLD);
    ~parRMA();

    parallelBranchSub * blankParallelSub();
    //loadBalDebug = data->loadBalDebug;

    // Need this to make sure the extra thread is set up
    void placeTasks();

    void pack(PackBuffer &outBuf);
    void unpack(UnPackBuffer &inBuf);
    int spPackSize();

    /*
      virtual bool continueRampUp() {
      return (spCount() <= rampUpFeatureFac * data->numAttrib)
      && parallelBranching::continueRampUp();
      }
    */

    /// Note: use VB flag?
    void reset(bool VBflag=true);

    // In parallel, restrict writing to verification log to processor
    // 0 when ramping up.
    bool verifyLog() {
      return _verifyLog && (!rampingUp() || (uMPI::rank == 0));
    };

    ostream* openVerifyLogFile();

    void setCachedCutPts(const unsigned int& j, const unsigned int& v) ;

    CutPtThd* cutPtCaster;		    // Thread to broadcast cut point data
    MessageID cutPtBroadcastTag;	// Message tag

  protected:
    double rampUpFeatureFac;

  }; // end reset function


  //**************************************************************************
  //  The parallel branchSub class...
  class parRMASub : virtual public parallelBranchSub, virtual public RMASub {

  public:

    parRMASub() {} //RMASub()
    virtual ~parRMASub() {}

    // Return a pointer to the global branching object
    parRMA* global() const { return globalPtr; }

    // Return a pointer to the parallel global base class object
    parallelBranching* pGlobal() const { return global(); }

    void setGlobalInfo(parRMA* global_) {
      globalPtr = global_;
      RMASub::setGlobalInfo(global_);	// set serial layer pointer etc.
    };

    virtual parallelBranchSub* makeParallelChild(int whichChild);

    void pack(utilib::PackBuffer &outBuffer);
    void unpack(utilib::UnPackBuffer & inBuffer);

    void boundComputation(double* controlParam);
    void parStrongBranching(const unsigned int& firstIdx, const unsigned int& lastIdx);
    void setLiveCachedCutPts();
    void parCachedBranching(unsigned int firstIdx, unsigned int lastIdx);

    void setNumLiveCutPts();

  protected:
    parRMA* globalPtr;  // A pointer to the global parallel branching object

  private:
    int numLiveCutPts;
    bool isCachedCutPts;

  }; // end parRMASub function


  // **********************************************************
  // CutPtThd
  class CutPtThd : public broadcastPBThread {
  public:
    CutPtThd(parRMA* global_, MessageID msgID);

    // virtual functions
    bool unloadBuffer();
    void initialLoadBuffer(PackBuffer* buf) { relayLoadBuffer(buf); };
    void relayLoadBuffer(PackBuffer* buf);

    void setCutPtThd(const unsigned int& f, const unsigned int& v);
    void preBroadcastMessage(const int& owningProc);

    int j, v;
    parRMA* ptrParRMA;
  }; // end CutPtThd class

} // namespace pebblRMA

#endif // ACRO_HAVE_MPI

#endif // pebbl_paralleRMA_h
