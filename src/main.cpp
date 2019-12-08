/*
 * main.cpp: parallel rma driver
 *
 * Author: Ai Kagawa
 */

#include <pebbl_config.h>
#include "serRMA.h"

#ifdef ACRO_HAVE_MPI
#include "parRMA.h"
#define outstream ucout
#define IO(action) if (uMPI::iDoIO) { CommonIO::end_tagging(); action; }
#else
#include "serRMA.h"
typedef void parRMA;
#define outstream cout
#define IO(action) action;
#endif

using namespace pebblRMA;

int main(int argc, char** argv) {
  return driver<RMA,parRMA>(argc,argv);
}
