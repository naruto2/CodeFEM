#include <cmath>
#include <vector>
#include "est/sparse.hpp"
#include "est/solver.hpp"


static const double tol=1e-6;  // 収束を判定する小さい数
static const int Itr=100;     //反復の最大値


static void mvsub(Matrix& A,Vector& b,Vector& c,int n)
{
  int i, j;
  // c = A・bを計算する部分
  for (i=0; i<n; i++){
    c[i]=0.0;
    for (j=0; j<n; j++){
      c[i]+=A[i][j]*b[j];
    }
  }
}

//大きさを１にする変換
static void normarize(Vector& x, int n)
{
  double c=0.0;
  int i;

  for (i=0; i<n; i++)
    c+=x[i]*x[i];
  c = sqrt(c);
  for (i=0; i<n; i++)
    x[i]/=c;
}

/* function power */
static double power(Matrix& A,Vector& x,int n,int *itr)
{
  int i=0, m=0;
  double c=0.0, s=0.0, lambda=0.0, xmax=0.0;
  static Vector y(n);

  normarize(x,n); //xを正規化する。
  for (m=0; m<itr[0]; m++){
    printf("m = %d\n",m);
    y = A * x;
    s=0.0; c=0.0;
    for (i=0; i<n; i++){
      s+=y[i]*y[i];
      c+=y[i]*x[i];
    }
    lambda=s/c;
    itr[1]=m;
    if (fabs(xmax-lambda)<tol)
      return lambda;
    xmax=lambda;
    //  cout << m << " " << lambda << "\n" ;
    for (i=0; i<n; i++)
      x[i]=y[i]/sqrt(s);  //x←y ,更に正規化する。
  }
  return lambda;
}


double maxesolver(matrix &A, Vector& x)
{
  int itr[2]={Itr, 0};

  x.resize(A.size());

  for ( int i=0; i< x.size(); i++ )
    x[i] = 1.0;
  return power(A,x,x.size(),itr);
}

