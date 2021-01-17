//**********************************************************
// Author:      Ai Kagawa
// Description: a serial or parallel driver
// ***********************************************************


#include "solveRMA.h"
using namespace rma;


int main(int argc, char** argv) {

  SolveRMA rma_sol;
  rma_sol.setupSolveRMA(argc, argv);
  rma_sol.solveRMA();

  return 0;

}
