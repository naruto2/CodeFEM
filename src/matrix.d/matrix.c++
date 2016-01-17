#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>


#include "est/matrix.hpp"


void plotmatrix(matrix A) {
  int i, j;


#if 0
  for(i=0; i<dim1(A); i++) {
    for(j=0; j<dim1(A); j++)
      if (A[i][j] != 0.0 )
	printf(".");
      else
	printf(" ");
    printf("\n");
  }
#endif

  FILE *pp;

  pp = popen("gnuplot","w");
  fprintf(pp,"unset border\n");
  fprintf(pp,"unset xtics\n");
 
  fprintf(pp,"set xrange [-1:%d]\n",dim1(A));
  fprintf(pp,"set yrange [%d:-1]\n",dim1(A));

  for(i=0; i<dim1(A); i++)
    for ( auto it : A[i] ) {
      j = it.first;
      if( A[i][j] != 0.0 )
	fprintf(pp,"set label \".\" at %d, %d;\n",j,i);
    }
  fprintf(pp,"plot '-' with lines title \"\"\n");
  fprintf(pp,"0 0\n");
  fprintf(pp,"%d %d\n",dim1(A)-1,dim1(A)-1);
  fprintf(pp,"e\n\n");
  fflush(pp);
  sleep(60*3);
  pclose(pp);
  

}
