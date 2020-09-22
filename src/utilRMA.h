#ifndef UTILITY_h
#define UTILITY_h

#include <limits>
#include <vector>
#include <deque>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

//int intInf = numeric_limits<int>::max();

double getInf();

ostream& operator<<(ostream& os, const deque<bool>& v);
ostream& operator<<(ostream& os, const vector<int>& v);
ostream& operator<<(ostream& os, const vector<double>& v);
ostream& operator<<(ostream& os, const vector<vector<int> >& v);
ostream& operator<<(ostream& os, const vector<vector<double> >& v);

#endif  // UTILITY_h
