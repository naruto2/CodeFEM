#include <cstdio>
#include <cmath>
#include "est/matrix.hpp"
#include "est/sparse.hpp"

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

vector<double>& operator*(const sparse::matrix<double>&A, vector<double>&b) {
  int n = A.size();
  static vector<double> x(n);

  for ( int i=0; i<n; i++ ) {
    auto j = A[i].begin();
    for ( x[i]=0.0; j != A[i].end(); j++ )
      x[i] += j->second * b[j->first];
  }
  return x;
}

#include "Preconditioner.hpp"
#include "bicgstab.hpp"

vector<double> sparse__bicgstab(Preconditioner& M, sparse::matrix<double>& A, vector<double>& b){
  static vector<double>x(b.size());
  int max_iter = 1000000;
  double tol = 0.0000001;
  sparse__BiCGSTAB(A, x, b, M, max_iter, tol);
  fprintf(stderr,"sparse__BiCGSTAB is successed\n");
  return x;
}
