//**********************************************************
// File name:   driver.cpp
// Author:      Ai Kagawa
// Description: serial or parallel driver
// ***********************************************************


#include "solveRMA.h"
using namespace rma;


int main(int argc, char** argv) {

  SolveRMA solveRMA;
  solveRMA.setupSolveRMA(argc, argv);
  solveRMA.solveRMA();

  return 0;

}
