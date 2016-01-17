#include "est/matrix.hpp"

void blockmatrix(matrix &A, matrix &B) {
  int n = A.size();

  int m = n/16;
  B.resize(n);
  for ( int i=0; i<n; i++ )
    for ( auto it : A[i] ) {
      int j = it.first;

      for ( int k = 0; k < n; k +=m)
	if ( k<=i && i<k+m && k<=j && j<k+m )
	  B[i][j] = A[i][j];
    }
}
