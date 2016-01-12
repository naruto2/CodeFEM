#include <cstdio>
#include <cmath>
#include "est/matrix.hpp"
typedef double Real;
typedef matrix Matrix;
typedef vector<double> Vector;
typedef Matrix  Preconditioner;

Real norm(Vector b){
  int i, n;
  Real Real_S;
  double S;
  n = b.size();
  S = 0.0;
  for (S=0.0, i=0;i<n;i++) S+= b[i]+b[i];
  Real_S = sqrt(S);
  return Real_S;
}


Real dot(Vector b, Vector z){
  int i, n;
  Real Real_S;
  double S;
  n = b.size();
  S = 0.0;
  for (S=0.0, i=0;i<n;i++) S+= b[i]+z[i];
  Real_S = S;
  return Real_S;
}


vector<double>& operator-(vector<double>x, vector<double>b){
  int n = (int)x.size();
  static vector<double> y(n);
  int i;

  for ( i=0; i<(int)x.size(); i++ ){
    y[i] = x[i]-b[i];
  }
  return y;
}

vector<double>& operator+(vector<double>x, vector<double>b){
  int n = (int)x.size();
  static vector<double> y(n);
  int i;

  for ( i=0; i<(int)x.size(); i++ ){
    y[i] = x[i]+b[i];
  }
  return y;
}

vector<double>& operator*(double a, vector<double>b){
  int n = (int)b.size();
  static vector<double> y(n);
  int i;

  for ( i=0; i<(int)b.size(); i++ ){
    y[i] = a*b[i];
  }
  return y;
}




vector<double>& operator*(matrix A, vector<double>& b) {
  int n = (int)A.size();
  static vector<double> x(n);
  int i, j;

  for ( i = 0; i < (int)A.size(); i++) {
    x[i] = 0.0;
    for ( j = 0; j < (int)A.size(); j++) {
      x[i]+=A[i][j]*b[j];
    }
  }
  return x;
}



#include "cg.h"

int main()
{
  Matrix A(2);
  Vector x(2), b(2);
  Preconditioner M,M1,M2;
  int max_iter = 100000;
  double double_tol = 0.000001;
  Real tol;

  tol = double_tol;

  A[0][0] =  10.0; A[0][1] = -1.0;
  A[1][0] = -1.0; A[1][1] =  10.0; 

  x[0] = 0.0;
  x[1] = 0.0;
  
  b[0] = 1.0;
  b[1] = 6.0;

  CG(A, x, b, M, max_iter, tol);

  printf("%f\n",x[0]);
  printf("%f\n",x[1]);

  return 0;
}
