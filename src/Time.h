/*
 *  Author:     Ai Kagawa
 * Description: timer related functiosn
 */


#ifndef TIME_h
#define TIME_h

#include <limits>
#include <ctime>        // std::time
#include <sys/time.h>
#include <iostream>
#include <string>
#include <pebbl/utilib/CommonIO.h>

using namespace utilib;


class Time {

public:

  // set the start time
  void startTime() {
    timeStartWall = get_wall_time();
    timeStartCPU  = clock();
  }

  // Wall time utility function
  double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
  }

  // get run time in CPU Time
  double getCPUTime() {
    timeEndCPU = clock();
    clockTicksTaken = timeEndCPU - timeStartCPU;
    return clockTicksTaken / (double) CLOCKS_PER_SEC ;
  }

  // get run time in Wall time
  double getWallTime() {
  	timeEndWall = get_wall_time();
    return  timeEndWall - timeStartWall ;
  }

  // print out the CPU time
  void printCPUTime() {
    ucout <<  "CPU Time: " << clockTicksTaken / (double) CLOCKS_PER_SEC <<"\n";
  }

  //print out Wall time
  void printWallTime() {
    ucout <<  "Wall Time: " << timeEndWall - timeStartWall <<"\n";
  }

private:
  double timeStartWall, timeEndWall;
  clock_t timeStartCPU, timeEndCPU, clockTicksTaken;

}; // class Time


#endif  // TIME_h
