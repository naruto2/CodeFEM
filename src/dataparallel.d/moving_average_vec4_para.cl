__kernel void moving_average_vec4_para(__global int4  *values,
                                       __global int *out,
                                       int length,
                                       int name_num,
                                       int width,
				       __global float *x,
				       __global float *y,
				       __global float *z)
{
  out[0] = get_global_size(0);
}


__kernel void mynorm(int n,__global float *x)
{
  int k, rank = get_global_id(0);
  if(!rank){
    x[0] = 0.0;
    for(k=1;k<n;k++) x[0] += x[k]*x[k];
    x[0] = sqrt(x[0]);
  }
}
