#include "utilRMA.h"

double getInf() {
 return std::numeric_limits<double>::infinity();
}

ostream &operator<<(ostream &os, const deque<bool> &v) {
  os << "(";
  for (deque<bool>::const_iterator i = v.begin(); i != v.end(); ++i)
    os << " " << *i;
  os << " )\n";
  return os;
}

ostream &operator<<(ostream &os, const vector<int> &v) {
  os << "(";
  for (vector<int>::const_iterator i = v.begin(); i != v.end(); ++i)
    os << " " << *i;
  os << " )\n";
  return os;
}

ostream &operator<<(ostream &os, const vector<double> &v) {
  os << "(";
  for (vector<double>::const_iterator i = v.begin(); i != v.end(); ++i)
    os << " " << *i;
  os << " )\n";
  return os;
}

ostream &operator<<(ostream &os, const vector<vector<int>> &v) {
  for (unsigned int i = 0; i < v.size(); ++i) {
    os << "(";
    for (unsigned int j = 0; j < v[i].size(); ++j)
      os << v[i][j] << " ";
    os << " )\n";
  }
  return os;
}

ostream &operator<<(ostream &os, const vector<vector<double>> &v) {
  for (unsigned int i = 0; i < v.size(); ++i) {
    os << "(";
    for (unsigned int j = 0; j < v[i].size(); ++j)
      os << v[i][j] << " ";
    os << " )\n";
  }
  return os;
}
