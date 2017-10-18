#include <cstdio>
#include <cmath>
#include "mmio.h"
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>
#include <est/op.hpp>
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



int glirulus_check(sparse::matrix<double>&A,vector<double>&x,vector<double>&b)
{
  vector<double> r;
  r = residual(A,x,b);
  printf("%s L_inf%f\n",getop("-A").c_str(),L_inf(r));
  return 0;
}



vector<double> glirulus(sparse::matrix<double>&A,vector<double>&b)
{
  vector<double> x, r;
  x = cl_bicgstab(A,b);

  r = residual(A,x,b);

  if ( L_inf(r) > 0.000999 ){
    if (A.size()<10000) x = GSLV1(A,b);
  }
  return x;
}
