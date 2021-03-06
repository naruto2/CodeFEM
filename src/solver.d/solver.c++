#include <cstdio>
#include <cmath>
#include "est/matrix.hpp"


double dot(const vector<double>&x, const vector<double>&y){
  double S = 0.0;
  int n = x.size();
  for ( int i=0; i<n; i++ ) S += x[i]*y[i];
  return S;
}

double norm(const vector<double>&x){
  return sqrt(dot(x,x));
}

vector<double>& operator-(const vector<double>&x, vector<double>&y){
  int n = x.size();
  static vector<double> z(n);
  for (int i=0; i<n; i++ ) z[i] = x[i]-y[i];
  return z;
}

vector<double>& operator+(vector<double>&x, vector<double>&y){
  int n = x.size();
  static vector<double> z(n);
  for (int i=0; i<n; i++ ) z[i] = x[i]+y[i];
  return z;
}

vector<double>& operator*(double a, vector<double>&x){
  int n = x.size();
  static vector<double> z(n);
  for (int i=0; i<n; i++ ) z[i] = a*x[i];
  return z;
}

vector<double>& operator*(const matrix&A, vector<double>&b) {
  int n = A.size();
  static vector<double> x(n);

  for ( int i=0; i<n; i++ ) {
    internal_map Ai = A[i];
    internal_map::iterator j = Ai.begin();
    for ( x[i]=0.0; j != Ai.end(); j++ )
      x[i] += j->second * b[j->first];
  }
  return x;
}

typedef double Real;
typedef matrix Matrix;
typedef vector<double> Vector;

#include "incholesky.h"
#include "Preconditioner.hpp"
#include "cg.h"

vector<double> cg(Preconditioner& M, matrix& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;
  CG(A, x, b, M, max_iter, tol);
  return x;
}


#include "ir.h"

vector<double> ir(Preconditioner& M, matrix& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;
  IR(A, x, b, M, max_iter, tol);
  return x;
}

#include "cgs.h"

vector<double> cgs(Preconditioner& M, matrix& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;
  CGS(A, x, b, M, max_iter, tol);
  return x;
}


#include "bicgstab.h"

vector<double> bicgstab(Preconditioner& M, matrix& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;
  BiCGSTAB(A, x, b, M, max_iter, tol);
  return x;
}



#include "bicg.h"

vector<double> bicg(Preconditioner& M, matrix& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;
  A.T();
  BiCG(A, x, b, M, max_iter, tol);
  return x;
}


#include "qmr.h"

vector<double> qmr(Preconditioner& M, Preconditioner& M2, matrix& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;
  A.T();
  QMR(A, x, b, M, M2, max_iter, tol);
  return x;
}

#include "gmres.h"

vector<double> gmres(Preconditioner& M, matrix& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;
  int m = 1000;
  matrix H;
  H.resize(A.size());
  GMRES(A, x, b, M, H, m, max_iter, tol);
  return x;
}


#include "cheby.h"

vector<double> cheby(Preconditioner& M, matrix& A, vector<double>& b, double mine, double maxe){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;
  CHEBY(A, x, b, M, max_iter, tol, mine, maxe);
  return x;
}



#include "jacobi.h"

vector<double> jacobi(Preconditioner& M, matrix& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;
  Jacobi(A, x, b, M, max_iter, tol);
  return x;
}

