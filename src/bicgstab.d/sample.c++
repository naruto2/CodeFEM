#include <cstdio>
#include <vector>
#include "est/sparse.hpp"
#include "est/bicgstab.hpp"

int main(int argc, char **argv)
{
  cl_bicgstab_init(argc,argv);

  sparse::matrix<double> A(3);
  vector<double> b(3), x;

  A[1][1] = 2.0; A[1][2] = 1.0; b[1] = 3.0;
  A[2][2] = 1.0; b[2] = 1.0;

  x = cl_bicgstab(A,b);
  printf("x[1]=%f\n",x[1]);
  printf("x[2]=%f\n",x[2]);
  return 0;
}
