#include <cstdio>
#include <cmath>
#include "est/sparse.hpp"
#include "cl_bicgstab.hpp"


int isTridiagonal(sparse::matrix<double>&A);
vector<double> TDMA(sparse::matrix<double>&A, vector<double>&b);


vector<double> cl_bicgstab(sparse::matrix<double>& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;

  vector<double> dinv(A.size());

  for ( int k=1, n=A.size(); k<n; k++){
    if (A[k][k] != 0.0) dinv[k] = 1.0/A[k][k];
    else { dinv[1] = 0.0; break; }
  }

  if ( dinv[1] != 0.0 ) if ( isTridiagonal(A) ) return TDMA(A,b);

  sparse__BiCGSTAB(A, &x[0], &b[0], max_iter, tol, &dinv[0]);
	
  return x;
}
