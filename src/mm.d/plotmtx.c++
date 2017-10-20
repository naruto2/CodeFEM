#include <cstdio>
#include <cmath>
#include "mmio.h"
#include "est/sparse.hpp"
#include "est/bicgstab.hpp"
#include "est/op.hpp"


void matrix2png(sparse::matrix<double>&A, const char *fname)
{
  FILE *pp;
  long i, j, n=A.size();

  pp = popen("gnuplot","w");
  fprintf(pp,"unset border\n");
  fprintf(pp,"unset xtics\n");
  fprintf(pp,"set xrange [0:%d]\n",n);
  fprintf(pp,"set yrange [%d:0]\n",n);
  fprintf(pp,"set size square\n");
  fprintf(pp,"set terminal png\n");
  fprintf(pp,"set output \"%s\"\n",fname);
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
  pclose(pp);
};

int main(int argc, char **argv)
{
  if ( argc < 2 ) return 0;
  initop(argc,argv);


  static double *val; static int *I, *J;
  int M, N, nz;

  mm_read_unsymmetric_sparse(getop("-f").c_str(),&M,&N,&nz,&val,&I,&J);

  if ( M != N ) return 0;

  sparse::matrix<double> A(N+1);
  int i, j, k;

  for (int k=0;k<nz;k++) {
    if ( 0<=I[k]&&I[k]<N+1&&0<=J[k]&&J[k]<N+1)
			   A[I[k]+1][J[k]+1] = val[k];

  }
  if ( defop("-o")) matrix2png(A,getop("-o").c_str());
  else plotmatrix(A);
  return 0;
}
