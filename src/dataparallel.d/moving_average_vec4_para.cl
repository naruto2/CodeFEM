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
  int np   = out[0]=get_global_size(0);
  int n    = 8192;
  int size = n/np;
  int rank = get_global_id(0);
  float tmp;

  for(j = rank, k=j*size;k<(j+1)*size;k++) z[k] = 0.0;
  if ( rank == 0 ) for(k=np*size;k<n;k++)  z[k] = 0.0;

  for(j = rank, k=j*size;k<(j+1)*size;k++) z[k] = x[k]+y[k];
  if ( rank == 0 ) for(k=np*size;k<n;k++)  z[k] = x[k]+y[k];

  if ( rank == 0 ){
    tmp = 0.0;
    for(k=0;k<n;k++) tmp += x[k]*y[k];
    out[1] = tmp;
    out[2] = tmp;
    out[3] = tmp;
    out[4] = tmp;
  }

  if ( rank == 1){
    out[1] = out[2]+out[3];
  }
}


__kernel void mynorm(int n,__global float *x)
{
  int np   = get_global_size(0);
  int rank = get_global_id(0);
  int k;

  if ( rank == 0 ) {
    x[0] = 0.0;
    for (k=1;k<n;k++) x[0] += x[k]*x[k];
  }

}
