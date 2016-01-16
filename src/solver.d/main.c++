#include <cstdlib>
#include "solver.hpp"
#include "cg.h"
#include "bicgstab.h"
#include "cgs.h"
#include "bicg.h"
#include "qmr.h"
#include "jacobi.h"
#include "incholesky.h"
#include "Preconditioner.hpp"
#include "psc98.hpp"

int main(int argc, char **argv){
  matrix A;
  vector<double> b;
  getprob(A,b);
  int n = A.size();
  vector<double> x(n);
  Preconditioner M, M2;
  int max_iter = 100000;
  double tol = 0.000001;
  
#if 0
  A.resize(3);
  A[0][0] = 1.0; A[0][1] = 0.5; A[0][2] = 0.0;
  A[1][0] = 0.5; A[1][1] = 1.0; A[1][2] = 0.5;
  A[2][0] = 0.0; A[2][1] = 0.5; A[2][2] = 1.0;

  x.resize(3);
  b.resize(3);
  b[0] = 3.0;
  b[1] = 4.0;
  b[2] = 3.0;
#endif
#if 1
  M.ic(A);
#endif  
  
  //A.jacobi(A);

  A.sync();
  CG(A, x, b, M, max_iter, tol);
  //Jacobi(A, x, b, M, max_iter, tol);
  A.T(); QMR(A, x, b, M, M2, max_iter, tol);
  check(x);

  for ( int i = 0; i<3; i++ ) printf("%f\n",x[i]);
  return 0;
}
