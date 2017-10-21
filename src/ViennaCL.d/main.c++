#include <stdio.h>
#include "est/sparse.hpp"
#include "est/psc98.hpp"
#include "est/ViennaCL.hpp"
#include "est/op.hpp"

int main(int argc, char **argv){
  initop(argc,argv);
  
  sparse::matrix<double> A; vector<double> x, b;
  psc98_init(A,b);
  x = vcl_bicgstab(A,b);
  psc98_check(x);
  return 0;
}
