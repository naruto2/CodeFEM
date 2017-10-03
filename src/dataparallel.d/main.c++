#include <stdlib.h>
#include <CL/cl.h>
#include <stdio.h>


double cl_norm(int n, double *x);
void   cl_copy(int n, double *y, double *x);
double  cl_dot(int n, double *y, double *x);
void   cl_init(int argc, char **argv);
void cl_finalize(void);

int main(int argc, char **argv)
{
  cl_init(argc,argv);

  int k, n = 65537;
  double *x = (double*)malloc(sizeof(double)*n);
  double *y = (double*)malloc(sizeof(double)*n);

  x[0] = 0.0; for(k=1;k<n;k++) x[k] = 1.0;

  printf("norm(n,x) = %f\n",cl_norm(n,x));

  cl_copy(n,y,x);
  printf("dot(n,y,x) = %f\n",cl_dot(n,y,x));

  cl_finalize();
  return 0;
}
