#include <cstdio>


#include "matrix.hpp"



void plotmatrix( matrix A );

void mulmatrix( vector<double> &x, matrix &A, vector<double> &b)
{
  int i, j;
  
  for ( i = 0; i < dim1(A); i++) {
    x[i] = 0.0;
    for ( j = 0; j < dim1(A); j++) {
      x[i]+=A[i][j]*b[j];
    }
  }
}

int main() {
  matrix A(2);
  vector<double> b(2), x(2);
  
  A[0][0] =  2.0;   A[0][1] = -1.0;
  A[1][0] = -1.0;   A[1][1] =  2.0;

  b[0] = 2.0;
  b[1] = 1.0;

  mulmatrix(x,A,b);
  
  printf("%f %f\n",x[0],x[1]);

  cout << A ;

  plotmatrix(A);
  return 0;
}

