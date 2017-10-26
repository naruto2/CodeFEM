#include <vector>
#include "est/op.hpp"
#include "est/sparse.hpp"
#include "est/bicgstab.hpp"
#include "est/Eigen.hpp"
#include "est/psc98.hpp"
#include "est/glirulus.hpp"
#include "est/ViennaCL.hpp"


int main(int argc, char **argv){
  //cl_bicgstab_init(argc,argv);
  glirulus_init(argc,argv);
  sparse::matrix<double> A; vector<double> x, b;
  psc98_init(A,b);
  x = cl_bicgstab(A,b);
  psc98_check(x);
  return 0;
}
