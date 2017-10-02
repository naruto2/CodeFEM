__kernel void moving_average_vec4_para(__global int4  *values,
                                       __global int *out,
                                       int length,
                                       int name_num,
                                       int width,
				       __global float *x,
				       __global float *y,
				       __global float *z)
{
  int i, j, k;
  int np = out[0]=get_global_size(0);
  int n = 8192;
  int size = n/np;
  int rank = get_global_id(0);
  float tmp;
  
  for ( j=0;j<np;j++ )if( j==rank ) {
      for(k=j*size;k<(j+1)*size;k++) z[k] = 0.0;
    }
#if 0
  for ( j=0;j<np;j++ )if( j==rank ) {
      for(k=j*size;k<(j+1)*size;k++) z[k] = x[k]+y[k];
    }
#endif

  for ( j=0;j<np;j++ )if( j==rank ) {
      z[j*size] = 0.0;
      for(k=j*size;k<(j+1)*size;k++) z[j*size] += x[k]*y[k];
    }
  barrier(CLK_GLOBAL_MEM_FENCE);
  if ( rank == 0 ) {
    tmp = 0.0;
    for ( k=0;k<np;k++) tmp += z[k*size];
  }  
  barrier(CLK_GLOBAL_MEM_FENCE);
  if ( rank == 0 ) z[0] = tmp;
}
