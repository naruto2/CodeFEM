#include <stdio.h>
#include "ary.h"
#include "mx.h"
#include "solver.h"


int main2(void){
  int i, n = 512;
  static MX *A; static double *x, *b;

  initmx(A,n+1,8); ary1(x, n+1); ary1(b, n+1);

  for ( i =1; i<=n; i++){
    mx(A,i,i) = 2.0; b[i] = 1.0;
  }
  for (i=1;i<n;i++) mx(A,i,i+1) = -1.0;
  for (i=1;i<n;i++) mx(A,i+1,i) = -1.0;

  estiva_pcgssolver2(A,x,b);

  for (i=1;i<=n;i++) printf("%f\n", x[i] );
}
