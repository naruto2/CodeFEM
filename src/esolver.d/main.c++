#include "est/solver.hpp"
#include "est/psc98.hpp"


double maxesolver(Matrix & B, Vector & x);

int main(int argc, char ** argv) {
  matrix A;
  vector<double> b,x;
  getprob(A,b);

  //plotmatrix(A);
  
  printf("%f\n",maxesolver(A,x));
  
  return 0;
}

