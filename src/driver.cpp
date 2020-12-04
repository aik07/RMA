//**********************************************************
// File name:   driver.cpp
// Author:      Ai Kagawa
// Description: serial or parallel driver
// ***********************************************************

#include "argRMA.h"
#include "dataRMA.h"
#include "solveRMA.h"
using namespace arg;
using namespace data;
using namespace rma;


int main(int argc, char** argv) {

  SolveRMA solveRMA;
  solveRMA.setupSolveRMA(argc, argv);
  solveRMA.solveRMA();

  return 0;

}
