#include <cstdio>
#include <cmath>
#include "est/sparse.hpp"
#include "bicgstab.hpp"

vector<double> sparse__bicgstab(sparse::matrix<double>& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;

  sparse__BiCGSTAB(A, &x[0], &b[0], max_iter, tol);

  return x;
}
