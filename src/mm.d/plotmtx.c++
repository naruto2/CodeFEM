#include <cstdio>
#include <cmath>
#include "mmio.h"
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>
#include <est/op.hpp>


void matrix2png(sparse::matrix<double>&A)
{
  FILE *pp;
  long i, j, n=A.size();

  pp = popen("gnuplot","w");
  fprintf(pp,"unset border\n");
  fprintf(pp,"unset xtics\n");
  fprintf(pp,"set xrange [0:%d]\n",n);
  fprintf(pp,"set yrange [%d:0]\n",n);
  fprintf(pp,"set size square\n");
  
  for(i=1;i<n;i++){
    for( auto it : A[i]){
      j = it.first;
      if ( A[i][j] != 0 ) fprintf(pp,"set label \".\" at %d, %d;\n",j,i);
    }
  }
  fprintf(pp,"plot '-' with lines title \"\"\n");
  fprintf(pp,"1 1\n");
  fprintf(pp,"%d %d\n",n-1,n-1);
  fprintf(pp,"e\n\n");
  fflush(pp);
  sleep(60);
  pclose(pp);
};

int main(int argc, char **argv)
{
  if ( argc < 2 ) return 0;
  initop(argc,argv);
  cl_bicgstab_init(argc,argv);

  static double *val; static int *I, *J;
  int M, N, nz;

  mm_read_unsymmetric_sparse(getop("-f").c_str(),&M,&N,&nz,&val,&I,&J);

  if ( M != N ) return 0;

  sparse::matrix<double> A(N+1);
  int i, j, k;

  for (int k=0;k<nz;k++) {
    //printf("%d %d %f\n",I[k],J[k],val[k]);
    if ( 0<=I[k]&&I[k]<N+1&&0<=J[k]&&J[k]<N+1)
			   A[I[k]+1][J[k]+1] = val[k];

  }
  matrix2png(A);
  return 0;
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
