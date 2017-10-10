#include <cstdio>
#include <cmath>
#include "est/sparse.hpp"
#include "cl_bicgstab.hpp"

vector<double> cl_bicgstab(sparse::matrix<double>& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;

  vector<double> dinv(A.size());

  for ( int k=1, n=A.size(); k<n; k++){
    if (A[k][k] != 0.0) dinv[k] = 1.0/A[k][k];
    else { dinv[1] = 0.0; break; }
  }
  /* for TDMA */
  int i=1, n= A.size();
  vector<double> a_(n), b_(n), c_(n);

  if ( dinv[1] != 0.0 ) {
    for ( i=1; i<n; i++)   a_[i] =  A[i][i];
    for ( i=1; i<n-1; i++) b_[i] = -A[i][i+1];
    for ( i=2; i<n;   i++) c_[i] = -A[i][i-1];
  }

  sparse__BiCGSTAB(A, &x[0], &b[0], max_iter, tol, &dinv[0],
		   &a_[0], &b[0], &c_[0]);

  return x;
}
