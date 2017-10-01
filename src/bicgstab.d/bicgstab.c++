#include <cstdio>
#include <cmath>
#include "est/sparse.hpp"
#include "bicgstab.hpp"

vector<double> sparse__bicgstab(sparse::matrix<double>& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;
  A[0][0] = 1.0;
  sparse__BiCGSTAB(A, x, b, max_iter, tol);
  fprintf(stderr,"sparse__BiCGSTAB is successed\n");
  return x;
}
