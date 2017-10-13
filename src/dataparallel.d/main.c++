#include <stdlib.h>
#include <CL/cl.h>
#include <stdio.h>
#include <vector>
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>

double cl_norm(int n, double *x);
void   cl_copy(int n, double *y, double *x);
double  cl_dot(int n, double *y, double *x);
void    cl_finalize(void);

double gp_norm(int n, double *x);
double gp_dot(int n, double *x, double *y);
double gp_copy(int n, double *x, double *y);

int main(int argc, char **argv)
{
  cl_bicgstab_init(argc,argv);

  int k, n = 65537;
  double *x = (double*)malloc(sizeof(double)*n);
  double *y = (double*)malloc(sizeof(double)*n);

  x[0] = 0.0; for(k=1;k<n;k++) x[k] = 1.0;

  printf("gp_norm(n,x) = %f\n",gp_norm(n,x));

  gp_copy(n,y,x);
  
  printf("gp_dot(n,y,x) = %f\n",gp_dot(n,y,x));

  cl_finalize();
  return 0;
}
