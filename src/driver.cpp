/**********************************************************
* File name:   driver.cpp
* Author:      Ai Kagawa
* Description: serial or parallel driver
***********************************************************/

#include "driverRMA.h"
using namespace pebblRMA;


int main(int argc, char** argv) {
  DriverRMA *driver = new DriverRMA(argc, argv);
  //driver->data->readData(argc, argv);
  //driver->setup(argc, argv);
  driver->solveRMA();
  return 0;
}
