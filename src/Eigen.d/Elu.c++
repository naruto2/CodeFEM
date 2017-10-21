#include <stdio.h>
#include "est/sparse.hpp"
#include "est/psc98.hpp"
#include "est/op.hpp"
#include "solvers.h"


vector<double> Elu(sparse::matrix<double>&A, vector<double>&b)
{
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