/*
 * driver.cpp: parallel rma driver
 *
 * Author: Ai Kagawa
 */

#include "rma.h"
using namespace pebblRMA;

int main(int argc, char** argv) {
  DriverRMA rma;
  rma.solveRMA();
  return 0;
}
