#include <vector>
#include "est/op.hpp"
#include "est/sparse.hpp"
#include "est/bicgstab.hpp"
#include "est/Eigen.hpp"
#include "est/psc98.hpp"
#include "est/glirulus.hpp"
#include "est/ViennaCL.hpp"
#include "est/perfectpivot.hpp"


int main(int argc, char **argv){
  glirulus_init(argc,argv);
  sparse::matrix<double> A; vector<double> x, b;
  psc98_init(A,b);
  x = glirulus(A,b);
  psc98_check(x);
  return 0;
}
