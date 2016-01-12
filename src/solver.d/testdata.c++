#include "solver.hpp"

int main()
{
  Matrix A(2);
  Vector x(2), b(2);
  Preconditioner M,M1,M2;
  int max_iter = 100000;
  double double_tol = 0.000001;
  Real tol;

  tol = double_tol;

  A[0][0] =  2.0; A[0][1] =  1.0;
  A[1][0] =  1.0; A[1][1] =  2.0; 

  x[0] = 0.0;
  x[1] = 0.0;
  
  b[0] = 1.0;
  b[1] = 5.0;

  
  CG(A, x, b, M, max_iter, tol);

  printf("%f\n",x[0]);
  printf("%f\n",x[1]);

  return 0;
}
