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

  double endCPUTime() {
    timeEndCPU = clock();
    clockTicksTaken = timeEndCPU - timeStartCPU;
    ucout <<  "CPU Time: " << clockTicksTaken / (double) CLOCKS_PER_SEC <<"\n";
    return clockTicksTaken / (double) CLOCKS_PER_SEC ;
  }

  double endWallTime() {
  	timeEndWall = get_wall_time();
  	ucout <<  "Wall Time: " << timeEndWall - timeStartWall <<"\n";
    return  timeEndWall - timeStartWall ;
  }

  double timeStartWall, timeEndWall;
  clock_t timeStartCPU, timeEndCPU, clockTicksTaken;

};


#endif  // TIME_h
