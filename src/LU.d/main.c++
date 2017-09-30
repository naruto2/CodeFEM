#include <cstdio>
#include "est/sparse.hpp"

int LU(sparse::matrix<double>&A);

int main()
{
  sparse::matrix<double>A(2);
  A[0][0]= 1.0;
  A[1][1] = 1.0;
  printf("LU=%d\n",LU(A));
  
  return 0;
}
