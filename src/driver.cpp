/**********************************************************
* File name:   driver.cpp
* Author:      Ai Kagawa
* Description: serial or parallel driver
***********************************************************/
#include "argRMA.h"
#include "dataRMA.h"
#include "driverRMA.h"
using namespace arg;
using namespace data;
using namespace rma;

int main(int argc, char** argv) {
  // Arguments *args   = new Arguments(argc, argv);
  // Data      *data   = new Data(argc, argv, args);

  Arguments args(argc, argv);
  Data      data(argc, argv, &args);
  DriverRMA *driver = new DriverRMA(&args, &data);
  //driver->data->readData(argc, argv);
  //driver->setup(argc, argv);

  driver->solveRMA();
  return 0;
}
