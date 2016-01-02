#include <cstdio>
#include <est/All>

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

#if 0
Real operator/(double x, Real y){
  Real z;
  z = x / y.x;
  return z;
}

double operator+(double x, Real y){
  return x + y.x;
}
#endif

int main()
{
  Matrix A(2);
  Vector x(2), b(2);
  Preconditioner M,M1,M2;
  int max_iter = 100000;
  double double_tol = 0.000001;
  Real tol;

  tol = double_tol;

  A[0][0] = 2.0; A[0][1] = 0.0;
  A[1][0] = 0.0; A[1][1] = 2.0; 

  x[0] = 0.0;
  x[1] = 0.0;
  
  b[0] = 1.0;
  b[1] = 6.0;

  CG(A, x, b, M, max_iter, tol);
  //BiCGSTAB(A, x, b, M, max_iter, tol);
  //CGS(A, x, b, M, max_iter, tol);
  //BiCG(A, x, b, M, max_iter, tol);
  //QMR(A,x,b,M1,M2,max_iter,tol);

  printf("%f\n",x[0]);
  printf("%f\n",x[1]);

  return 0;
}
