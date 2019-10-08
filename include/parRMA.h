/**
 * file parRMA.h
 * author Ai Kagawa
 *
 *  Example class to use object-oriented branching framework.
 *  Solves Rectangle Maximum Agreement.
 */

#ifndef pebbl_paralleRMA_h
#define pebbl_paralleRMA_h

#include <iostream>
#include <pebbl_config.h>

#ifdef ACRO_HAVE_MPI
#include <pebbl/pbb/parBranching.h>

#include <vector>
#include "serRMA.h"

using namespace utilib;
using namespace std;
using namespace pebbl;


namespace pebblRMA {

class RMA;
class RMASub;

#ifdef ACRO_HAVE_MPI
	class CutPtThd;
#endif

	//**************************************************************************
	//  The parallel branching class...
	class parRMA : public parallelBranching, public RMA {

	public:

		parRMA();
		~parRMA();

		parallelBranchSub * blankParallelSub();

		bool setup(int& argc,char**& argv) {
		  return parallelBranching::setup(argc,argv);
		}

		void setParameter(Data* data, const int& deb_int) {
			debug = deb_int;
			//////////////////////////////////////////
			loadBalDebug = data->loadBalDebug;
		}

		virtual void printSolutionTime() const {
			ucout << "ERMA Solution: " << incumbentValue
						<< "\tCPU time: " << totalCPU << "\n";
		}

		// Need this to make sure the extra thread is set up
		void placeTasks();

		void pack(PackBuffer &outBuf);
		void unpack(UnPackBuffer &inBuf);
		int spPackSize();

	  virtual bool continueRampUp() {
		  return (spCount() <= rampUpFeatureFac*numAttrib)
							&& parallelBranching::continueRampUp();
	  }

		/// Note: use VB flag?
		void reset(const bool& VBflag=true);

	  // In parallel, restrict writing to verification log to processor
	  // 0 when ramping up.
	  bool verifyLog() {
	  	return _verifyLog && (!rampingUp() || (uMPI::rank == 0));
	  };

	  ostream* openVerifyLogFile();

	  void setCachedCutPts(const int& j, const int& v) ;

	  CutPtThd* cutPtCaster;		    // Thread to broadcast cut point data
	  MessageID cutPtBroadcastTag;	// Message tag

	protected:
		double rampUpFeatureFac;

	};//************************************************************************


	//**************************************************************************
	//  The parallel branchSub class...
	class parRMASub : public parallelBranchSub, public RMASub {

	public:

		parRMASub() : RMASub() {}
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
		void parStrongBranching(const int& firstIdx, const int& lastIdx);
		void setLiveCachedCutPts();
		void parCachedBranching(int firstIdx, int lastIdx);

		void setNumLiveCutPts();

	protected:
		parRMA* globalPtr;  // A pointer to the global parallel branching object

	private:
		int numLiveCutPts;
		bool isCachedCutPts;

	};// **********************************************************


	// **********************************************************
	// CutPtThd
	class CutPtThd : public broadcastPBThread {
	public:
		CutPtThd(parRMA* global_, MessageID msgID);

		// virtual functions
		bool unloadBuffer();
		void initialLoadBuffer(PackBuffer* buf) { relayLoadBuffer(buf); };
		void relayLoadBuffer(PackBuffer* buf);

		void setCutPtThd(const int& f, const int& v);
		void preBroadcastMessage(const int& owningProc);

		int j, v;
		parRMA* ptrParRMA;
	}; // **********************************************************

} // namespace lpboost

#endif // ACRO_HAVE_MPI
#endif // pebbl_paralleRMA_h
