#include <cstdio>

#include "matrix.hpp"

void plotmatrix( matrix A );

int main() {
  int n = 10;
  matrix A(n+1);
  matrix B(2);
  matrix C(2);

  A[0][0] = 1.1;
  A[1][0] = 2.1;
  C[0][1] = 1.2;  C[0][1] = 1.2;
  B[0][1] = 1.2;
  A[0][1] = 1.2;
  A[1][1] = 2.2;

  B[0][0] = 1.1;
  B[1][0] = 2.1;

  B[1][1] = 2.2;

  C[0][0] = 1.1;
  C[1][0] = 2.1;

  C[1][1] = 2.2;
  
  A[2][2] = 3.3;
  
  for(int i=0; i<2; i++ ) for (int j=0; j<2; j++ ) {
      A[i][j] += 10.0;
      C[i][j] += 30.0;
      B[i][j] += 20.0;
    }

  int k;
  for ( k = 1; k<A.size(); k++) {
    A[k][k] = 1.0;
  }
  cout << A << B << C;

  plotmatrix(A);
  
  return 0;
}

