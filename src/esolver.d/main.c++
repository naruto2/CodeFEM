#include "est/solver.hpp"
#include "est/psc98.hpp"


double maxesolver(Matrix & B, Vector & x);

static int halfbw(Matrix& a)
{
  int n, i, j, w ;

  n = a.size();
  for( i=0; i<n; i++ )
    for( auto it : a[i] ) {
      j = it.first;
      if ( a[i][j] != 0.0 ) if( w < abs(i-j) ) w = abs(i-j);
    }
  return w;
}

#include "lu.h"

int main(int argc, char ** argv) {
  matrix A;
  vector<double> b,x;
  getprob(A,b);

  printf("halfbw = %d\n",halfbw(A));

  LUDecomp(A,A.size());
  plotmatrix(A);
  
  //printf("%f\n",maxesolver(A,x));
  
  return 0;
}

