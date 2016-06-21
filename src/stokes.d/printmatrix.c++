#include <stdio.h>
#include "est/matrix.hpp"

void printmatrix(matrix &A, const char *name)
{
  FILE *fp = fopen(name,"w");
  
  fprintf(fp,"size = %ld\n",A.size()-1);

  unsigned long i, j, n;
  n = A.size()-1;

  for ( i = 1; i <= n; i++ ){
    for ( j = 1; j <=n; j++ ){
      fprintf(fp,"%f ",A[i][j]);
    }
    fprintf(fp,"\n");
  }
}



void printvector(vector<double> &b, const char *name)
{
  FILE *fp = fopen(name,"w");

  fprintf(fp,"size = %ld\n",b.size()-1);

  unsigned long i, n;
  n = b.size()-1;

  for ( i = 1; i<=n; i++ ) {
    fprintf(fp,"%f\n",b[i]);
  }
  fclose(fp);
}
