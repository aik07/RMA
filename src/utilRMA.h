#ifndef UTILITY_h
#define UTILITY_h

#include <limits>
#include <vector>
#include <deque>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

double getInf();    // get the double infinity value
int    getIntInf(); // get the integer inifinitiy value

template<class T>
ostream& operator<<(ostream& os, const deque<T>& v);

template<class T>
ostream& operator<<(ostream& os, const vector<T>& v);

template<class T>
ostream& operator<<(ostream& os, const vector<vector<T> >& v);

// ostream& operator<<(ostream& os, const deque<bool>& v);
// ostream& operator<<(ostream& os, const vector<int>& v);
// ostream& operator<<(ostream& os, const vector<double>& v);
// ostream& operator<<(ostream& os, const vector<vector<int> >& v);
// ostream& operator<<(ostream& os, const vector<vector<double> >& v);

#endif  // UTILITY_h
