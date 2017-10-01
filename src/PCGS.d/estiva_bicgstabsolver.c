#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "estiva/ary.h"
#include "estiva/mx.h"
#include "estiva/op.h"
#include "estiva/solver.h"
#include "estiva/std.h"
#include "estiva/vec.h"
#include "estiva/eblas.h"
int estiva_bicgstabsolver(void *Apointer, double *xk, double *b);
  
int estiva_bicgstab(long *Ai, long *Aj, double *Aa, double *X, double *B)
{
  long i, k, n;

  for (i=1;Ai[i] != 0; i++ );
  n = Ai[i-1];
  
  static MX *A;
  initmx(A,n+1,30);

  for(k=0; Ai[k]!=0; k++) mx(A,Ai[k],Aj[k]) = Aa[k];

  static double *x, *b;
  ary1(x,n+1); ary1(b,n+1);
  
  for (i=1; i<=n; i++) b[i] = B[i];
  estiva_bicgstabsolver(A, x, b);

  for(int i=1;i<=n;i++) X[i] = x[i];
  return 0;
}
