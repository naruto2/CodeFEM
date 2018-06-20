/*--------------------
  kernel.cl
    --------------------*/

__kernel void vec_add(int n, __global float* z,
		      __global float* x, __global float* y)
{
  int gid = get_global_id(0);
  
  if(gid < n){
    z[gid] = x[gid] + y[gid];
  }
}
