#include "est/solver.hpp"
#include "est/esolver.hpp"
#include "est/psc98.hpp"

int main(int argc, char ** argv) {
  matrix A, A1, A2; vector<double> b, b1, b2, x;
  getprob(A,b);
#if 0
  A.resize(3);
  b.resize(A.size());
  A[0][0] = 1.00; A[0][1] = 0.50; A[0][2] = 0.00;
  A[1][0] = 0.50; A[1][1] = 1.00; A[1][2] = 0.50;
  A[2][0] = 0.00; A[2][1] = 0.50; A[2][2] = 1.00;

  b[0] = 2.0;
  b[1] = 3.0;
  b[2] = 2.0;
#endif
  Preconditioner M;
  A1 = A; b1 = b;
  A2 = A; b2 = b;
  double mine = minesolver(A1,b1);
  double maxe = maxesolver(A2,b2);
  printf("mine = %f maxe = %f\n", mine, maxe);
  x = cheby(M,A,b,mine,maxe);

  //for ( int i=0; i< x.size(); i++) printf("%f\n",x[i]); return 0;
  check(x);
  return 0;
}

