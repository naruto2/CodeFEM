#include <cstdio>
#include <cmath>
#include "est/sparse.hpp"

int enough(sparse::matrix<double>&A, vector<double>&x, const vector<double>&b);
  

double dot(const vector<double>&x, const vector<double>&y){
  double S = 0.0;
  int n = x.size();
  for ( int i=1; i<n; i++ ) S += x[i]*y[i];
  return S;
}

double norm(const vector<double>&x){
  return sqrt(dot(x,x));
}

vector<double>& operator-(const vector<double>&x, vector<double>&y){
  int n = x.size();
  static vector<double> z(n);
  for (int i=1; i<n; i++ ) z[i] = x[i]-y[i];
  return z;
}

vector<double>& operator+(vector<double>&x, vector<double>&y){
  int n = x.size();
  static vector<double> z(n);
  for (int i=1; i<n; i++ ) z[i] = x[i]+y[i];
  return z;
}

vector<double>& operator*(double a, vector<double>&x){
  int n = x.size();
  static vector<double> z(n);
  for (int i=1; i<n; i++ ) z[i] = a*x[i];
  return z;
}

vector<double>& operator*(const sparse::matrix<double>&A, vector<double>&b) {
  int n = A.size();
  static vector<double> x(n);

  for ( int i=1; i<n; i++ ) {
    x[i] = 0.0;
    for(auto  Ai : A[i]) {
      int j = Ai.first;
      x[i] += Ai.second * b[j];
    }
  }
  return x;
}

//*****************************************************************
// Iterative template routine -- Jacobi
//
// Jacobi solves the symmetric positive definite linear
// system Ax=b using the Jacobi method.
//
// Jacobi follows the algorithm described on p. 12 in the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************

typedef sparse::matrix<double> Matrix;
typedef vector<double> Vector;
typedef double Real;
int 
Jacobi(Matrix &A, Vector &x, const Vector &b,
   int &max_iter, Real &tol)
{
  Real resid;
  Vector p, z, q;
  Vector alpha(1), beta(1), rho(1), rho_1(1);

  Real normb = norm(b);


  Vector r = b - A*x;
  printf("enter\n");
  
  if (normb == 0.0) 
    normb = 1;

  if ((resid = norm(r) / normb) <= tol && enough(A,x,b) ) {
    tol = resid;
    printf("resid=%f,normb=%f\n",resid,normb);
    max_iter = 0;
    return 0;
  }

  printf("world\n");

  int n = A.size();
  
  for (int k = 1; k <= max_iter; k++) {
    printf("loop\n");

    x = A*r;
    for (int i=1; i<n; i++) x[i] = (b[i]-x[i])/A[i][i];
#if 0
    for ( int i=0; i<n; i++ ) {
      internal_map Ai = A[i];
      internal_map::iterator j = Ai.begin();
      for ( x[i]=0.0; j != Ai.end(); j++ ) if ( i != j->first ) {
	  x[i] += j->second * r[j->first];
	}
      x[i] = (b[i]-x[i])/A[i][i];
    }
#endif
    
    r = x;    

    if ((resid = norm(b - A*x) / normb) <= tol && enough(A,x,b) ) {
      tol = resid;
      max_iter = k;
      return 0;     
    }
    printf("k = %d    resid = %f\n",k,resid);
  }
  tol = resid;
  return 1;
}



vector<double> jacobi(sparse::matrix<double>& A, vector<double>&x,
		      vector<double>& b){
  int max_iter = 1000000;
  double tol = 0.00000001;
  printf("hello jacobi\n");
  Jacobi(A, x, b, max_iter, tol);
  return x;
}

