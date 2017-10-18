#include <cstdio>
#include <est/sparse.hpp>
#include <vector>
#include "glirulus.hpp"

int main(int argc, char **argv)
{
  glirulus_init(argc,argv);
    
  sparse::matrix<double> A;  vector<double> b,x;

  int ret=glirulus_mm(A,b); if(ret) return ret;

  x = glirulus(A,b);

  return glirulus_check(A,x,b);
}
