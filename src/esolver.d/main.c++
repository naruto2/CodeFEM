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


int LUSolver(matrix &A, vector<double> &b, vector<double> &x, int n)
{
  if(n <= 0) return 0;
  int j;
  // 前進代入(forward substitution)
  //  LY=bからYを計算
  for(int i = 0; i < n; ++i){
    double bly = b[i];
    //for(int j = 0; j < i; ++j){
    for ( auto it : A[i] ) {
      j = it.first;
      if ( j < i )
	bly -= A[i][j]*x[j];
    }
    x[i] = bly/A[i][i];
  }

  // 後退代入(back substitution)
  //  UX=YからXを計算
  for(int i = n-1; i >= 0; --i){
    double yux = x[i];
    //for(int j = i+1; j < n; ++j){
    for ( auto it : A[i] ) {
      j = it.first;
      if ( i < j && j < n )
	yux -= A[i][j]*x[j];
    }
    x[i] = yux;
  }

  return 1;
}

#include "lu.h"
#include "../solver.d/cg.h"
#include "../solver.d/incholesky.h"
#include "../solver.d/Preconditioner.hpp"

int main(int argc, char ** argv) {
  matrix A;
  vector<double> b,x;
  getprob(A,b);

#if 0
  A.resize(3);
  b.resize(A.size());
  A[0][0] = 1.00; A[0][1] = 0.50; A[0][2] = 0.00;
  A[1][0] = 0.50; A[1][1] = 1.00; A[1][2] = 0.50;
  A[2][0] = 0.00; A[2][1] = 0.50; A[2][2] = 1.00;

  b[0] = 2.0;
  b[1] = 3.0;
  b[2] = 2.0;
#endif
  Preconditioner M;
  int max_iter = 100000;
  double tol = 0.000001;

  x.resize(b.size());  
  LUDecomp(A,A.size()); LUSolver(A,b,x,x.size());
  //CG(A, x, b, M, max_iter, tol);
  //for ( int i=0; i< x.size(); i++) printf("%f\n",x[i]);
  check(x);
  //printf("%f\n",maxesolver(A,x));
  return 0;
}

