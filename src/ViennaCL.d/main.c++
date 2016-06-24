#include <stdio.h>
#include "est/matrix.hpp"
#include "est/psc98.hpp"
#include "ViennaCL.hpp"

int main(int argc, char **argv){

  matrix A; vector<double> x, b;
  getprob(A,b);

  x = gpubicgstab(A,b);
  
  check(x);
  return 0;
}
