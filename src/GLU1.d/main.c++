#include <cstdio>
#include "est/sparse.hpp"

int GLU1(sparse::matrix<double>&A);
int GSLV1(sparse::matrix<double>&A, vector<double>&b);

int main()
{
  sparse::matrix<double> A(5);
  vector<double> b(5);
  A[1][1] = 1.0; b[1] = 1.0;
  A[2][2] = 1.0; b[2] = 2.0;
  A[3][4] = 1.0; b[3] = 3.0;
  A[4][3] = 1.0; b[4] = 4.0;
  
  printf("GLU1(A) = %d\n",GLU1(A));
  GSLV1(A,b);
  for ( int i=1; i<b.size(); i++)  printf("%f ",b[i]);
  printf("\n");
  return 0;
}
