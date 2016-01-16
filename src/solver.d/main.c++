#include <cstdlib>
#include "solver.hpp"
#include "cg.h"
#include "bicgstab.h"
#include "cgs.h"
#include "bicg.h"
#include "qmr.h"
#include "jacobi.h"

#include "incholesky.h"

void backward(matrix&L, vector<double>&y, vector<double>& b, vector<double>& d, vector<double>& x) {
  int n = L.size();
  
  for(int i = 0; i < n; ++i){
    double rly = b[i];
    for ( auto it : L[i] ) { int j = it.first;
      if ( j<i ) rly -= L[i][j]*y[j];
    }
    y[i] = rly/L[i][i];
  }

  for(int i = n-1; i >= 0; --i){
    double lu = 0.0;
    //for(int j = i+1; j < n; ++j){
    for ( auto it : L[i] ) { int j = it.first;
      if ( j>i ) lu += L[j][i]*x[j];
    }
    x[i] = y[i]-d[i]*lu;
  }

#if 0
  y[0] = b[0];
  int n = L.size();
  for ( int i = 1; i < n; i++ ){
    y[i] = b[i];
    internal_map Li = L[i];
    for ( internal_map::iterator k = Li.begin(); k != Li.end(); k++) y[i] -= L[i][k->first]*y[k->first];
    //for (int  k = 0; k < i; k++ ) y[i] -= L[i][k]*y[k];
  }

  for ( int i = 0; i<n; i++ ) y[i] /= d[i];

  matrix LT(L.size());

  static int init = 1;
  if (init = 1 ) { // LT trans L
    init = 0;
    int n = L.size();
    for(int i = 0; i<n; i++ ){
      internal_map Li = L[i];
      for ( internal_map::iterator  j = Li.begin(); j != Li.end() ; j++ ) LT[j->first][i] = L[i][j->first];
    }
    LT.sync();
  }

  x[n-1] = y[n-1];
  
  for ( int i = n-2; i >= 0; i-- ) {
    x[i] = y[i];
    internal_map LTi = LT[i];
    for ( internal_map::iterator j = LTi.begin(); j != LTi.end(); j++) if(j->first != i) x[i] -= LT[i][j->first]*x[j->first];
    //for ( int j = i+1; j<n; j++) x[i] -= LT[i][j]*x[j];
  }
#endif     

}

#include "Preconditioner.hpp"
extern "C" {
  void genmat(int,int*,double*,double*);
  void chkval(FILE*,int,double*);
}

int main(int argc, char **argv){
  vector<double> AA(10);
  double B;
  int    JA[10], i, j, n, w;

  genmat(-1,&JA[0],&AA[0],&B);
  n = JA[0]; w = JA[2];

  matrix A(n); vector<double> x(n), b(n);

  for ( i=1; i<=n; i++) {
    for (j=0;j<=w-1;j++) {
      JA[j] =  -1;
      AA[j] = 0.0;
    }
    genmat(i,&JA[0],&AA[0],&B);
    for ( j=0; j<w; j++) if (JA[j] != -1){
	if( JA[j] <=0 || n < JA[j] ) {
	  ;
	} else {
	  A[i-1][JA[j]-1] = AA[j];
	  b[i-1] = B;
	}
      }
  }

  Preconditioner M,M2;
  int max_iter = 100000;
  double double_tol = 0.000001;
  Real tol;

  tol = double_tol;

#if 0
  A.resize(3);
  A[0][0] = 1.0; A[0][1] = 0.5; A[0][2] = 0.0;
  A[1][0] = 0.5; A[1][1] = 1.0; A[1][2] = 0.5;
  A[2][0] = 0.0; A[2][1] = 0.5; A[2][2] = 1.0;

  x.resize(3);
  b.resize(3);
  b[0] = 3.0;
  b[1] = 4.0;
  b[2] = 3.0;
#endif
  A.sync();



#if 1
  matrix L(A.size());
  vector<double> d(A.size());
  IncompleteCholeskyDecomp2(A, L, d, A.size());
  M.setic(L,d);
#endif  
  //A.T();
  //A.jacobi(A);
  //vector<double> y(A.size());
  //backward(L,y,b,d,x);
  //cout <<A<< L;
  //Jacobi(A, x, b, M, max_iter, tol);

  CG(A, x, b, M, max_iter, tol);
  for (int i = 0; i < 3; i++) printf("%f\n",x[i]);
  return 0;
  //QMR(A, x, b, M, M2, max_iter, tol);
  
  chkval(stdout,n,&x[0]);
  return 0;
}
