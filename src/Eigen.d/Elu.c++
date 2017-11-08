#include <stdio.h>
#include "est/sparse.hpp"
#include "est/psc98.hpp"
#include "est/op.hpp"
#include "est/TDMA.hpp"
#include "solvers.h"

vector<double> jacobi(sparse::matrix<double>&A,
		      vector<double>&x,vector<double>&b);

int enough(sparse::matrix<double>&A, vector<double>&x, vector<double>&b);

vector<double> Elu(sparse::matrix<double>&A, vector<double>&b)
{
  int diag = 1;

  for ( int k=1, n=A.size(); k<n; k++)
    if (A[k][k] == 0.0) { diag = 0; break; }

  if ( diag && isTridiagonal(A) ) return TDMA(A,b);


  if ( A.size() > 16000 ) {
    fprintf(stderr,"Warning: Elu() can't solve n < 16000 matrix\n");
    vector<double> x(A.size());
    for(int k=0;k<A.size();k++) x[k] = b[k];
    return x;
  }


  for ( int i=1; i<A.size(); i++)
    for ( auto it: A[i]) { int j = it.first;
      Tri(i,j,A[i][j]);
    }
  Tri(0,0,1.0); b[0] = 0.0;
  Smatrix Aa = MapSmatrix(A.size(),A.size());

  Vector bb(A.size());
  for(int k=0;k<A.size();k++) bb[k] = b[k];

  Vector xx = Elu(Aa,bb);
  
  vector<double> x(A.size());
  for(int k=0;k<A.size();k++) x[k] = xx[k];

  if(!enough(A,x,b)) jacobi(A,x,b);
  return x;
}

#if 0
int main(int argc, char **argv){
  initop(argc,argv);

  sparse::matrix<double> A; vector<double> x,b;
  psc98_init(A,b);
  x = Elu(A,b);
  psc98_check(x);
  return 0;
}
#endif
