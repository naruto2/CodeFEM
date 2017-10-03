#pragma OPENCL EXTENSION cl_khr_fp64 : enable


__kernel void cl_norm(int n,__global double *x,__global double *npa)
{
  int np   = get_global_size(0);
  int size = n/np;
  int j, k, rank = get_global_id(0);

  npa[rank] = 0.0;
  for(j = rank, k=j*size;k<(j+1)*size;k++) npa[j] += x[k]*x[k];
  if ( rank == 0 ) for(k=np*size;k<n;k++)  npa[0] += x[k]*x[k];
}


__kernel void cl_copy(int n,__global double *y, __global double *x)
{
  int j, k;
  int np   = get_global_size(0);
  int size = n/np;
  int rank = get_global_id(0);

  for(j = rank, k=j*size;k<(j+1)*size;k++) y[k] = x[k];
  if ( rank == 0 ) for(k=np*size;k<n;k++)  y[k] = x[k];
}


__kernel void cl_dot(int n,__global double *y, __global double *x,
		     __global double *npa)
{
  int np   = get_global_size(0);
  int size = n/np;
  int j, k, rank = get_global_id(0);

  npa[rank] = 0.0;
  for(j = rank, k=j*size;k<(j+1)*size;k++) npa[j] += y[k]*x[k];
  if ( rank == 0 ) for(k=np*size;k<n;k++)  npa[0] += y[k]*x[k];
}


__kernel void cl_phase1(int n,__global double *p, __global double *r,
			__global double *v, double beta, double omega)
{
  int j, k;
  int np   = get_global_size(0);
  int size = n/np;
  int rank = get_global_id(0);

  for(j = rank, k=j*size;k<(j+1)*size;k++)
     p[k] = r[k] + beta * (p[k] - omega *v[k]);
  if ( rank == 0 ) for(k=np*size;k<n;k++)
		     p[k] = r[k] + beta * (p[k] - omega *v[k]);
}


__kernel void cl_phase3(int n,__global double *s, __global double *r,
			__global double *v, double alpha)
{
  int j, k;
  int np   = get_global_size(0);
  int size = n/np;
  int rank = get_global_id(0);

  for(j = rank, k=j*size;k<(j+1)*size;k++)
    s[k] = r[k] - alpha * v[k];
  if ( rank == 0 ) for(k=np*size;k<n;k++)
		     s[k] = r[k] - alpha * v[k];
}


__kernel void cl_phase4(int n,__global double *x, __global double *phat,
			double alpha)
{
  int j, k;
  int np   = get_global_size(0);
  int size = n/np;
  int rank = get_global_id(0);

  for(j = rank, k=j*size;k<(j+1)*size;k++)
    x[k] = x[k] + alpha*phat[k];
  
  if ( rank == 0 ) for(k=np*size;k<n;k++)
		     x[k] = x[k] + alpha*phat[k];
}
