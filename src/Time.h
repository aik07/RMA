/*
 *  File name: Time.h
 *  Author:    Ai Kagawa
 */


#ifndef TIME_h
#define TIME_h

#include <limits>
#include <ctime>        // std::time
#include <sys/time.h>
#include <iostream>
#include <pebbl/utilib/CommonIO.h>


using namespace utilib;

class Time {

public:

  double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
  }

  void startTime() {
    timeStartWall = get_wall_time();
    timeStartCPU = clock();
  }

  double getCPUTime() {
    timeEndCPU = clock();
    clockTicksTaken = timeEndCPU - timeStartCPU;
    return clockTicksTaken / (double) CLOCKS_PER_SEC ;
  }

  double getWallTime() {
  	timeEndWall = get_wall_time();
    return  timeEndWall - timeStartWall ;
  }

  void printCPUTime() {
    ucout <<  "CPU Time: " << clockTicksTaken / (double) CLOCKS_PER_SEC <<"\n";
  }

  void printWallTime() {
    ucout <<  "Wall Time: " << timeEndWall - timeStartWall <<"\n";
  }

private:
  double timeStartWall, timeEndWall;
  clock_t timeStartCPU, timeEndCPU, clockTicksTaken;

};


#endif  // TIME_h
