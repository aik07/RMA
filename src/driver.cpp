//**********************************************************
// File name:   driver.cpp
// Author:      Ai Kagawa
// Description: serial or parallel driver
// ***********************************************************

#include "argRMA.h"
#include "dataRMA.h"
#include "driverRMA.h"
using namespace arg;
using namespace data;
using namespace rma;


int main(int argc, char** argv) {

  DriverRMA driverRMA;
  driverRMA.setupDriverRMA(argc, argv);
  driverRMA.solveRMA();

  // Using the following code, I get an error
  // DriverRMA* driverRMA = new DriverRMA();
  // driverRMA->setupDriverRMA(argc, argv);
  // driverRMA->solveRMA();

  return 0;

}
