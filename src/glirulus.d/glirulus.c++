#include <cstdio>
#include <cmath>
#include "mmio.h"
#include "est/sparse.hpp"
#include "est/bicgstab.hpp"
#include "est/ViennaCL.hpp"
#include "est/Eigen.hpp"
#include "est/op.hpp"
#include <vector>
#include "est/GLU1.hpp"

double L_inf(vector<double>&x)
{
  int N = x.size();
  double max=0.0;
  for (int k=1;k<=N;k++) if ( max < abs(x[k])) max=abs(x[k]);
  return max;
}

int glirulus_init(int argc, char **argv)
{
  initop(argc,argv);
  cl_bicgstab_init(argc,argv);
  return 0;
}


int glirulus_mm(sparse::matrix<double>&A,vector<double>&b)
{
  if ( !defop("-A") ) {
    fprintf(stderr,"Usage: glirulus -A file.mtx\n");
    return 1;
  }

  FILE *fp;
  fp = fopen(getop("-A").c_str(),"r");
  if ( fp == NULL ) {
    fprintf(stderr,"Error: Can't open file %s\n",getop("-A").c_str());
    return 2;
  }
  fclose(fp); fp = NULL;

  static double *val; static int *I, *J;
  int M, N, nz, ret;

  ret = mm_read_unsymmetric_sparse(getop("-A").c_str(),&M,&N,&nz,&val,&I,&J);

  if ( ret ){
    fprintf(stderr,"Error: mm_read_unsymmetric_sparse()=%d file %s\n",ret,
	    getop("-A").c_str());
    return ret;
  }
  if ( M != N ) {
    fprintf(stderr,"Error: %s's matrix is not square.\n",getop("-A").c_str());
    return 3;
  }


  A.resize(N+1);
  int i, j, k;

  for (int k=0;k<nz;k++) {
    if ( 0<=I[k]&&I[k]<N+1&&0<=J[k]&&J[k]<N+1)
			   A[I[k]+1][J[k]+1] = val[k];
  }

  b.resize(N+1);
  for (int i=1; i<=N; i++){
    b[i] = 0.0;
    for (auto it: A[i]) { int j = it.first;
      b[i] += A[i][j];
    }
  }

  return 0;
}

vector<double> residual(sparse::matrix<double>&A,vector<double>&x,vector<double>&b)
{
  int n = A.size();
  vector<double> r(n);

  for (int i=1; i<n; i++){
    r[i] = 0.0;
    for (auto it: A[i]) { int j = it.first;
      r[i] += A[i][j]*x[j];
    }
  }

  for (int k=1; k<n; k++) r[k] -= b[k];
  return r;
}



double glirulus_check(sparse::matrix<double>&A,vector<double>&x,vector<double>&b)
{
  vector<double> r = residual(A,x,b);
  return L_inf(r);
}


int isReallyNaN(double x) {
  return x != x;    // xがNaNであればtrue, それ以外ではfalse
}

vector<double> glirulus(sparse::matrix<double>&A,vector<double>&b)
{
  double res;
  vector<double> x;
  x = vcl_bicgstab(A,b);

  for ( int i=0; i<x.size(); i++)
    if ( isReallyNaN(x[i])) {
      x = cl_bicgstab(A,b);
      break;
    }
  res = glirulus_check(A,x,b);
  if(defop("-v")) fprintf(stderr,"res= %f\n",res);
  for (int k=0; k<8; k++ ) {
    if ( res < 0.0000004 ) break;
    if (k<4) x = cl_bicgstab(A,b);
    else x = Elu(A,b);
    res = glirulus_check(A,x,b);
    if(defop("-v")) fprintf(stderr,"res= %f\n",res);
  }
  return x;
}
