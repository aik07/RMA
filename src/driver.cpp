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

  DriverRMA *driver = new DriverRMA(argc, argv);
  driver->solveRMA();
  return 0;

}
