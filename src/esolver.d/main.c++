#include "est/solver.hpp"
#include "est/psc98.hpp"

double maxesolver(Matrix & B, Vector & x);

#include "lu.h"
#include "../solver.d/cg.h"
#include "../solver.d/incholesky.h"
#include "../solver.d/Preconditioner.hpp"

/* function invert */
static double invert(Matrix& A,Vector& x,double tol,int n,int *itr)
{
  //static Matrix L;
  //static Matrix U;
  static Vector y(n);
  static Vector s(n);
  double c, d, xmin=0.0,lambda=10000000.0;
  int m, i;

  //ary2(L,n+1,n+1);
  //ary2(U,n+1,n+1);
  //ary1(y,n+1);
  //ary1(s,n+1);
  
  
  //lu(A,L,U,tol,n);           // AをLU分解する。結果はLとUに入れる
  LUDecomp(A);
  for (m=0; m<itr[0]; m++){
    for (i=0; i<n; i++)
      s[i]=x[i];
    //subs(L,U,s,y,tol,n); //　LUy = s の解yを求める処理。
    LUSolver(A,s,y);
    c=0; d=0;
    for (i=0; i<n; i++){
      c+=y[i]*x[i];
      d+=y[i]*y[i];
    }
    lambda =c/d;
    itr[1]=m;
    if (fabs(xmin-lambda)<tol){
      return lambda;  /* 収束した */
    }
    xmin=lambda;
    for (i=0; i<n; i++)
      x[i]=y[i]/sqrt(d);  // xにyを正規化して代入
  }
  itr[1]=m;
  return lambda;  //Itr回では収束しなかった
}


double minesolver(matrix &A, vector<double> &x)
{
  int n=x.size(), itr[2]={1000000, 0};
  return invert(A,x,0.0000001,n,itr);
}




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
  printf("%f\n",minesolver(A,b)); return 0;
  //LUDecomp(A); LUSolver(A,b,x);
  //CG(A, x, b, M, max_iter, tol);
  //for ( int i=0; i< x.size(); i++) printf("%f\n",x[i]);
  check(x);
  //printf("%f\n",maxesolver(A,x));
  return 0;
}

