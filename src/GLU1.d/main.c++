#include <cstdio>
#include "est/sparse.hpp"

int GLU1(sparse::matrix<double>&A);

int main()
{
  sparse::matrix<double> A(5);

  A[1][1] = 1.0;
  A[2][2] = 1.0;
  A[3][4] = 1.0;
  A[4][3] = 1.0;
  
  printf("GLU1(A) = %d\n",GLU1(A));
  return 0;
}
