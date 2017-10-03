#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void cl_norm(int n,__global float *x)
{
  int k, rank = get_global_id(0);
  if(!rank){
    x[0] = 0.0;
    for(k=1;k<n;k++) x[0] += x[k]*x[k];
    x[0] = sqrt(x[0]);
  }
}

__kernel void cl_copy(int n,__global float *y, __global float *x)
{
  int j, k;
  int np   = get_global_size(0);
  int size = n/np;
  int rank = get_global_id(0);

  for(j = rank, k=j*size;k<(j+1)*size;k++) y[k] = x[k];
  if ( rank == 0 ) for(k=np*size;k<n;k++)  y[k] = x[k];
}

__kernel void cl_dot(int n,__global float *y, __global float *x)
{
  int k, rank = get_global_id(0);
  if(!rank){
    x[0] = 0.0;
    for(k=1;k<n;k++) x[0] += y[k]*x[k];
  }
}
