#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define THREADSPERBLOCK 256

__kernel void dotProduct(__global const float* a,
             __global const float* b,
             __global float *partial_c,
             int n,
             int threadsPerBlock,
             int blocksPerGrid)
{
  int i;
  local float cache[THREADSPERBLOCK];
  int tid = get_global_id(0);
  int cacheIndex = get_local_id(0);
  int gid = get_group_id(0);
  float temp = 0;

  while(tid < n){
    temp += a[tid] * b[tid];
    tid += threadsPerBlock * blocksPerGrid;
  }

  cache[cacheIndex] = temp;
  barrier(CLK_LOCAL_MEM_FENCE);

  i = threadsPerBlock / 2;

  while( i != 0 ){
    if(cacheIndex < i)
      cache[cacheIndex] += cache[cacheIndex + i];
    barrier(CLK_LOCAL_MEM_FENCE);
    i /= 2;
  }

  if(cacheIndex == 0){
    partial_c[gid] = cache[0];
  }
}
