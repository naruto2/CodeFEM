#include "kernel.h"
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

/*--------------------
  kernel.cl
  --------------------*/

__kernel void vec_add(int n, __global double *w, __global double* z,
		      __global double* x, __global double* y)
{
  int i = n/16;
  int g = get_global_id(0);
  int gid;
  w[g] = 0.0;
  for (gid = get_global_id(0); gid<i; gid +=128) {
    w[g] += x[gid] * y[gid]
      + x[gid+1] * y[gid+1]
      + x[gid+2] * y[gid+2]
      + x[gid+3] * y[gid+3]
      + x[gid+4] * y[gid+4]
      + x[gid+5] * y[gid+5]
      + x[gid+6] * y[gid+6]
      + x[gid+7] * y[gid+7]
      + x[gid+8] * y[gid+8]
      + x[gid+9] * y[gid+9]
      + x[gid+10] * y[gid+10]
      + x[gid+11] * y[gid+11]
      + x[gid+12] * y[gid+12]
      + x[gid+13] * y[gid+13]
      + x[gid+14] * y[gid+14]
      + x[gid+15] * y[gid+15]
      ;
  }
}




__kernel void vec_dot(int n, __global double* z,
		      __global double* x, __global double* y)
{
  __global double16 *z16 = (__global double16*)z;
  __global double16 *x16 = (__global double16*)x;
  __global double16 *y16 = (__global double16*)y;

  int range=get_global_id(0)*2048;
    do {
      if ( range <= n-16 )             z16[range] = x16[range] * y16[range];
      else for (  ; range<n; range++ )   z[range] =   x[range] *   y[range];
      range +=2048;
    } while (range < n-16);
    z[2046] = 3.14;
}



