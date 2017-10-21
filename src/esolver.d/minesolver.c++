#include "est/sparse.hpp"
#include "est/solver.hpp"
#include "lu.h"

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
