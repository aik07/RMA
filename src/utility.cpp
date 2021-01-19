/*
 *  Author:     Ai Kagawa
 * Description: a source file of utility functions for RMA
 */


#include "utility.h"


// get the double infinity value
double getInf() {
 return std::numeric_limits<double>::infinity();
}


// get the integer inifinitiy value
int getIntInf() {
 return std::numeric_limits<int>::max();
}


// for 1D deque output
template<class T>
ostream &operator<<(ostream &os, const deque<T> &v) {
  os << "(";
  for (typename deque<T>::const_iterator i = v.begin(); i != v.end(); ++i)
    os << " " << *i;
  os << " )\n";
  return os;
}
template ostream &operator<< <bool>(ostream &os, const deque<bool> &v);


// for 1D set output
template<class T>
ostream &operator<<(ostream &os, const set<T> &v) {
  os << "(";
  for (typename set<T>::const_iterator i = v.begin(); i != v.end(); ++i)
    os << " " << *i;
  os << " )\n";
  return os;
}


// for map output
ostream &operator<<(ostream &os, const map<double, int> &m) {
  for (typename map<double, int>::const_iterator i =m.begin(); i !=m.end(); ++i)
    os << " [" << i->first << ':' << i->second << ']';
  os << " \n";
  return os;
}


// for 1D const vector output
template<class T>
ostream &operator<<(ostream &os, const vector<T> &v) {
  os << "(";
  for (typename vector<T>::const_iterator i = v.begin(); i != v.end(); ++i)
    os << " " << *i;
  os << " )\n";
  return os;
}
template ostream &operator<< <bool>(ostream &os, const vector<bool> &v);
template ostream &operator<< <double>(ostream &os, const vector<double> &v);
template ostream &operator<< <int>(ostream &os, const vector<int> &v);
template ostream &operator<< <unsigned int>(ostream &os,
                                            const vector<unsigned int> &v);


// // for 1D vector output
// template<class T>
// ostream &operator<<(ostream &os, vector<T> &v) {
// os << "(";
// for (typename vector<T>::const_iterator i = v.begin(); i != v.end(); ++i)
//   os << " " << *i;
// os << " )\n";
// return os;
// }
// template ostream &operator<< <bool>(ostream &os, vector<bool> &v);
// template ostream &operator<< <double>(ostream &os, vector<double> &v);
// template ostream &operator<< <int>(ostream &os, vector<int> &v);
// template ostream &operator<< <unsigned int>(ostream &os,
//                                             vector<unsigned int> &v);


// for 2D vector output
template<class T>
ostream &operator<<(ostream &os, const vector<vector<T> > &v) {
  for (unsigned int i = 0; i < v.size(); ++i) {
    os << "(";
    for (unsigned int j = 0; j < v[i].size(); ++j)
      os << v[i][j] << " ";
    os << " )\n";
  }
  return os;
}


template ostream &operator<< <unsigned int>(ostream &os,
                 const vector<vector<unsigned int> > &v);

template ostream &operator<< <double>(ostream &os,
                 const vector<vector<double> > &v);

template ostream &operator<< <bool>(ostream &os,
                 const vector<vector<bool> > &v);


// for 2D deque output
template<class T>
ostream &operator<<(ostream &os, const deque<deque<T> > &v) {
for (unsigned int i = 0; i < v.size(); ++i) {
  os << "(";
  for (unsigned int j = 0; j < v[i].size(); ++j)
    os << v[i][j] << " ";
  os << " )\n";
}
return os;
}

template ostream &operator<< <bool>(ostream &os,
               const deque<deque<bool> > &v);
