#include <vector>
#include "est/sparse.hpp"
#include "est/bicgstab.hpp"
#include "est/psc98.hpp"


int main(int argc, char **argv){
  cl_bicgstab_init(argc,argv);

  sparse::matrix<double> A; vector<double> x, b;
  psc98_init(A,b);
  x = cl_bicgstab(A,b);
  psc98_check(x);
  return 0;
}
