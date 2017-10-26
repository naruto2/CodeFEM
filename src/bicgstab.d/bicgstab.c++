#include <cstdio>
#include <cmath>
#include "est/sparse.hpp"
#include "cl_bicgstab.hpp"


int isTridiagonal(sparse::matrix<double>&A);
vector<double> TDMA(sparse::matrix<double>&A, vector<double>&b);


vector<double> cl_bicgstab(sparse::matrix<double>& A, vector<double>& b){
  A[0][0] = 1.0; b[0] = 0.0;
  static vector<double>x(b.size());
  int max_iter = A.size();
  double tol = 0.0000001;

  vector<double> dinv(A.size());

  for ( int k=1, n=A.size(); k<n; k++){
    if (A[k][k] != 0.0) dinv[k] = 1.0/A[k][k];
    else { dinv[1] = 0.0; break; }
  }

  if ( dinv[1] != 0.0 ) if ( isTridiagonal(A) ) return TDMA(A,b);



  printf("b[0] = %f\n",b[0]);
  printf("b[1] = %f\n",b[1]);
  printf("b[b.size()-3] = %f\n",b[b.size()-3]);
  printf("b[b.size()-2] = %f\n",b[b.size()-2]);
  printf("b[b.size()-1] = %f\n",b[b.size()-1]);


  sparse__BiCGSTAB(A, &x[0], &b[0], max_iter, tol, &dinv[0]);

  printf("x[0] = %f\n",x[0]);
  printf("x[1] = %f\n",x[1]);
  printf("x[x.size()-3] = %f\n",x[x.size()-3]);
  printf("x[x.size()-2] = %f\n",x[x.size()-2]);
  printf("x[x.size()-1] = %f\n",x[x.size()-1]);
  
  return x;
}
