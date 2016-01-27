#include <cstdlib>
#include "est/solver.hpp"
#include "cg.h"
#include "bicgstab.h"
#include "cgs.h"
#include "bicg.h"
#include "qmr.h"
#include "gmres.h"
#include "ir.h"
#include "cheby.h"
#include "jacobi.h"
#include "incholesky.h"
#include "lu.h"
#include "Preconditioner.hpp"
#include "est/psc98.hpp"
#include "est/matrix.hpp"

void blockmatrix(matrix &A, matrix &B);

int main(int argc, char **argv){
  
#if 0
  matrix B(3);
  B[0][0] = 1.00; B[0][1] = 0.50; B[0][2] = 0.00; 
  B[1][0] = 0.50; B[1][1] = 1.00; B[1][2] = 0.50; 
  B[2][0] = 0.00; B[2][1] = 0.50; B[2][2] = 1.00;
  LUDecomp(B,B.size());
  cout << B << endl;
  return 0;
#endif
  matrix A;
  vector<double> b;
  getprob(A,b);
  int n = A.size();
  vector<double> x(n);
  Preconditioner M, M2;
  int max_iter = 100000;
  double tol = 0.000001;

#if 1
  matrix H;
  H.resize(n);
  int m = 1000;
  M.ic(A);
  //GMRES(A, x, b, M, H, m, max_iter, tol);
  CHEBY(A,x,b,M,max_iter, tol, 19.738218,129105.765540);
  //CG(A, x, b, M, max_iter, tol);
  //M2.ic(A); A.T(); QMR(A, x, b, M, M2, max_iter, tol);
#endif
  check(x);
  return 0;
}
