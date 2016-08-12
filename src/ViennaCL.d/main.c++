#include <stdio.h>
#include "est/matrix.hpp"
#include "est/psc98.hpp"
#include "ViennaCLinf.hpp"

int main(int argc, char **argv){

  matrix A; vector<double> x, b;
  getprob(A,b);

  x = gpucgilut(A,b);
  
  check(x);

  return 0;
}
