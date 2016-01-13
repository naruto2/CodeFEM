#include <cstdio>
#include <cmath>
#include "est/matrix.hpp"
#include "Preconditioner.hpp"

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
