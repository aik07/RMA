/**********************************************************
* File name:   driver.cpp
* Author:      Ai Kagawa
* Description: serial or parallel driver
***********************************************************/

#include "driverRMA.h"
using namespace pebblRMA;


int main(int argc, char** argv) {
  Arguments *args   = new Arguments(argc, argv);
  Data      *data   = new Data(argc, argv, args);
  DriverRMA *driver = new DriverRMA(args, data);
  //driver->data->readData(argc, argv);
  //driver->setup(argc, argv);
  driver->solveRMA();
  return 0;
}
