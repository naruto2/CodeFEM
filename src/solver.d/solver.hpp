#include <cstdio>
#include <cmath>
#include "est/matrix.hpp"
typedef double Real;
typedef matrix Matrix;
typedef vector<double> Vector;


Real norm(Vector b){
  int i, n;
  Real Real_S;
  double S;
  n = b.size();
  S = 0.0;
  for (S=0.0, i=0;i<n;i++) S+= b[i]*b[i];
  Real_S = sqrt(S);
  return Real_S;
}


Real dot(Vector b, Vector z){
  int i, n;
  Real Real_S;
  double S;
  n = b.size();
  S = 0.0;
  for (S=0.0, i=0;i<n;i++) S+= b[i]*z[i];
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




vector<double>& operator*(matrix A, vector<double> b) {
  int n = (int)A.size();
  static vector<double> x(n);
  int i, j;


  for ( i = 0; i < (int)A.size(); i++) {
    x[i] = 0.0;
    
#if 1
    internal_map Ai;
    Ai = A[i];    
    internal_map::iterator it;
    it = Ai.begin();
    while ( it != Ai.end() ) {
      x[i]+=it->second*b[it->first];
      it++;
    }

#endif

#if 0
      for ( j = 0; j < (int)A.size(); j++) if (A[i][j]!=0.0){
	  printf ("%d %.1f ",j, A[i][j]);
	x[i]+=A[i][j]*b[j];
      }    
    printf("\n");
#endif
    
#if 0    
    for ( j = 0; j < (int)A.size(); j++) 
      x[i]+=A[i][j]*b[j];
#endif 
  }
  return x;
}



#include "cg.h"

