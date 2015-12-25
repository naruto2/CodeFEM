#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>


#include "est/matrix.hpp"


void plotmatrix(matrix A) {
  int i, j;

  for(i=0; i<A.size(); i++) {
    for(j=0; j<A.size(); j++)
      if (A[i][j] != 0.0 )
	printf(".");
      else
	printf(" ");
    printf("\n");
  }


  FILE *pp;

  pp = popen("gnuplot","w");
  fprintf(pp,"unset border\n");
  fprintf(pp,"unset xtics\n");
 
  fprintf(pp,"set xrange [-1:%ld]\n",A.size());
  fprintf(pp,"set yrange [%ld:-1]\n",A.size());

  for(i=0; i<A.size(); i++)
    for (j=0; j<A.size(); j++)
      if( A[i][j] != 0.0 )
	fprintf(pp,"set label \".\" at %ld, %ld;\n",j,i);

  fprintf(pp,"plot '-' with lines title \"\"\n");
  fprintf(pp,"0 0\n");
  fprintf(pp,"%d %d\n",A.size()-1,A.size()-1);
  fprintf(pp,"e\n\n");
  fflush(pp);
  sleep(60*3);
  pclose(pp);
  

}
