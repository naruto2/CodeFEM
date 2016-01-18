#include <cstdlib>
#include "solver.hpp"
#include "cg.h"
#include "bicgstab.h"
#include "cgs.h"
#include "bicg.h"
#include "qmr.h"
#include "jacobi.h"
#include "incholesky.h"
#include "Preconditioner.hpp"
#include "psc98.hpp"
#include "est/matrix.hpp"

void blockmatrix(matrix &A, matrix &B);


long dim2(matrix& a){
  return a.size()-1;		    
}

long halfbw(matrix& a){
  return a.size()-1;
}

#define ary1(p,n) static vector<double> p(n)
#define forall(a,x,b) for(x=a; x<=b; x++)

double absv(double v) {
  return fabs(v);
}

long lessv(long a, long b) {
  if ( a < b ) return a;
  return b;
}





vector<double>& LU(matrix& a)
{
  double nrm,aik,apkk;
  internal_map apk, ai;
  static long i,j,k,w,n,pk,pklimit,klimit;
  n=dim2(a);w=halfbw(a);

  ary1(p,n+1);

  for(k=1;k<=n;k++){
    pk=k; nrm=absv(a[k][k]);

    if(nrm == 0.0) return p;
    p[k]=pk;
    klimit=lessv(k+w,n); apkk= -a[pk][k];
    forall(k+1,i,klimit){
      a[i][k] /= apkk;
      pklimit=lessv(pk+w,n); aik=a[i][k]; apk=a[pk]; ai=a[i];
      forall(k+1,j,pklimit) ai[j] += aik*apk[j];
    }
  }

  return p;
}


static long gauss_c(matrix& a, vector<double> b)
{
  double s,bpi;
  long i,j,n,w,pi,pilimit;
  internal_map api;

  vector<double> p;
   
  n = dim2(a);
  w = halfbw(a);

  p = LU(a);

  for(i=0;i<=n;i++){
    pilimit=lessv(i+w,n); bpi=b[p[i]];
    //forall(i+1,j,pilimit) b[j] += a[j][i]*bpi;
    forall(i,j,pilimit) b[j] += a[j][i]*bpi;
  }
  for(i=n;i>=0;i--){
    s = b[(pi=p[i])]; api=a[pi]; pilimit=lessv(pi+w,n);
    //forall(i+1,j,pilimit) s -= api[j]*b[j];
    forall(i,j,pilimit) s -= api[j]*b[j];
    b[pi] = s/api[i];
  }

  return 1;
}




int main(int argc, char **argv){
  matrix A;
  vector<double> b;
  getprob(A,b);
  int n = A.size();
  vector<double> x(n);
  Preconditioner M, M2;
  int max_iter = 100000;
  double tol = 0.000001;
  
#if 0
  A.resize(4);
  A[0][0] = 1.00; A[0][1] = 0.00; A[0][2] = 0.00; A[0][3] = 0.00;
  A[1][0] = 0.00; A[1][1] = 1.00; A[1][2] = 0.00; A[1][3] = 0.00;
  A[2][0] = 0.00; A[2][1] = 0.00; A[2][2] = 1.00; A[2][3] = 0.00;
  A[3][0] = 0.00; A[3][1] = 0.00; A[3][2] = 0.00; A[3][3] = 1.00;

  x.resize(4);
  b.resize(4);
  b[0] = 1.0;
  b[1] = 1.0;
  b[2] = 1.0;
  b[3] = 2.0;
  A.sync();
#endif
#if 1
  A.resize(3);
  A[0][0] = 1.00; A[0][1] = 0.00; A[0][2] = 0.00; 
  A[1][0] = 0.00; A[1][1] = 1.00; A[1][2] = 1.00; 
  A[2][0] = 0.00; A[2][1] = 0.00; A[2][2] = 1.00; 

  x.resize(3);
  b.resize(3);
  b[0] = 1.0;
  b[1] = 1.0;
  b[2] = 2.0;
  A.sync();
#endif
#if 1
  M.ic(A);
#endif  
  gauss_c(A,b);
  cout << A;
  for ( int i=0; i<2; i++) printf("%f\n",b[i]);
  return 0;
  //matrix B;
  //blockmatrix(A,B);
  A.sync();
  CG(A, x, b, M, max_iter, tol);
  //Jacobi(A, x, b, M, max_iter, tol);
  //A.T(); QMR(A, x, b, M, M2, max_iter, tol);
  check(x);
  return 0;
}
