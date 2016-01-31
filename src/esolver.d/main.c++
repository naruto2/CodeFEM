#include "est/solver.hpp"
#include "est/esolver.hpp"
#include "est/psc98.hpp"

//#include "../solver.d/cg.h"
//#include "../solver.d/incholesky.h"
//#include "../solver.d/Preconditioner.hpp"

int main(int argc, char ** argv) {
  matrix A;
  vector<double> b,x;
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
  //Preconditioner M;
  int max_iter = 100000;
  double tol = 0.000001;

  x.resize(b.size());  
  printf("%f\n",minesolver(A,b)); return 0;
  //LUDecomp(A); LUSolver(A,b,x);
  //CG(A, x, b, M, max_iter, tol);
  //for ( int i=0; i< x.size(); i++) printf("%f\n",x[i]);
  check(x);
  //printf("%f\n",maxesolver(A,x));
  return 0;
}

