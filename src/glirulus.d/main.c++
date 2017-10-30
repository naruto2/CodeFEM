#include <cstdio>
#include "est/sparse.hpp"
#include <vector>
#include "est/glirulus.hpp"
#include "est/ViennaCL.hpp"
#include "est/bicgstab.hpp"

int main(int argc, char **argv)
{
  double res;
  glirulus_init(argc,argv);
    
  sparse::matrix<double> A;  vector<double> b,x;

  int ret=glirulus_mm(A,b); if(ret) return ret;

  x = glirulus(A,b);

  res = glirulus_check(A,x,b);
  
  return 0; 
}
