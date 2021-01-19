/*
 *  Author:     Ai Kagawa
 * Description: a header file of utility functions for RMA
 */


#ifndef UTILITY_h
#define UTILITY_h

#include <limits>
#include <vector>
#include <deque>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

double getInf();    // get the double infinity value
int    getIntInf(); // get the integer inifinitiy value

// for 1D deque ooutput
template<class T>
ostream& operator<<(ostream& os, const deque<T>& v);

// for 1D set output
template<class T>
ostream& operator<<(ostream& os, const set<T>& v);

// for map output
ostream& operator<<(ostream& os, const map<double, int>& v);

// for 1D vector output
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v);

// for 1D vector output
// template<class T>
// ostream& operator<<(ostream& os, vector<T>& v);

// for 2D vector output
template<class T>
ostream& operator<<(ostream& os, const vector<vector<T> >& v);

// for 2D deque output
template<class T>
ostream& operator<<(ostream& os, const deque<deque<T> >& v);

#endif  // UTILITY_h
