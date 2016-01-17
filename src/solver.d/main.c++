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
#include "est/matrix.hpp"


void blockmatrix(matrix &A, matrix &B);


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
  A.resize(4);
  A[0][0] = 1.80; A[0][1] = 2.88; A[0][2] = 2.05; A[0][3] =-0.89;
  A[1][0] = 5.25; A[1][1] =-2.95; A[1][2] =-0.95; A[1][3] =-3.80;
  A[2][0] = 1.58; A[2][1] =-2.69; A[2][2] =-2.90; A[2][3] =-1.04;
  A[3][0] =-1.11; A[3][1] =-0.66; A[3][2] =-0.59; A[3][3] = 0.80;

  x.resize(4);
  b.resize(4);
  b[0] = 3.0;
  b[1] = 4.0;
  b[2] = 3.0;
  A.sync();
#endif
#if 1
  M.ic(A);
#endif  
  matrix B;
  blockmatrix(A,B);
  A.sync();
  CG(A, x, b, M, max_iter, tol);
  //Jacobi(A, x, b, M, max_iter, tol);
  //A.T(); QMR(A, x, b, M, M2, max_iter, tol);
  check(x);
  return 0;
}
