#include <cstdio>
#include <cmath>
#include "mmio.h"
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>
#include <est/op.hpp>

int main(int argc, char **argv)
{
  if ( argc < 2 ) return 0;
  initop(argc,argv);
  cl_bicgstab_init(argc,argv);

  static double *val; static int *I, *J;
  int M, N, nz, ret;

  ret = mm_read_unsymmetric_sparse(getop("-f").c_str(),&M,&N,&nz,&val,&I,&J);

  if ( ret != 0 ) printf("mm_read_unsymmetric_sparse()=%d %s\n",ret,
			 getop("-f").c_str());
  if ( M != N ) return 0;

  sparse::matrix<double> A(N+1);
  int i, j, k;

  for (int k=0;k<nz;k++) {
    //printf("%d %d %f\n",I[k],J[k],val[k]);
    if ( 0<=I[k]&&I[k]<N+1&&0<=J[k]&&J[k]<N+1)
			   A[I[k]+1][J[k]+1] = val[k];

  }

  vector<double> x,b(N+1);

  for (int i=1; i<=N; i++){
    b[i] = 0.0;
    for (auto it: A[i]) { int j = it.first;
      b[i] += A[i][j];
    }
  }
  x = cl_bicgstab(A,b);
  
  //for (k=1;k<=N;k++) printf("%d %f\n",k,x[k]);
  for (k=1;k<=N;k++) x[k] -= 1.0;

  double max=0.0;
  for ( k=1;k<=N;k++) if ( max < abs(x[k])) max=abs(x[k]);
  printf("%s L_inf%f\n",getop("-f").c_str(),max);
  return 0;
}
