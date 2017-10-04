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


__kernel void cl_phase0(int n, __global double *r,
		   __global double *Aa, __global int *col_ind,
		   __global int *row_ptr, __global double *x,
		   __global double *rtilde, __global double *b,
		   __global double *npa)
{
  int i,j=0, k;
  int np   = get_global_size(0);
  int size = n/np;
  int rank = get_global_id(0);
  i = rank;

  npa[i]=0.0;
  for(j = rank, k=j*size;k<(j+1)*size;k++) if(k!=0){
    r[k] = 0.0;
    for (int j=row_ptr[k];j<row_ptr[k+1];j++) r[k] += Aa[j] * x[col_ind[j]];
    rtilde[k] = r[k] = b[k] - r[k];
    npa[i] += r[k]*r[k];
  }
  if ( rank == 0 )
    for(k=np*size;k<n;k++)
    {
    r[k] = 0.0;
    for (int j=row_ptr[k];j<row_ptr[k+1];j++) r[k] += Aa[j] * x[col_ind[j]];
    rtilde[k] = r[k] = b[k] - r[k];
    npa[i] += r[k]*r[k];
  }

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


__kernel void cl_phase6(int n,__global double *x, __global double *s,
			__global double *r, __global double *t,
			__global double *phat, __global double *shat,
			double alpha, double omega)
{
  int j, k;
  int np   = get_global_size(0);
  int size = n/np;
  int rank = get_global_id(0);

  for(j = rank, k=j*size;k<(j+1)*size;k++){
    x[k] = x[k] + alpha*phat[k] + omega*shat[k];
    r[k] = s[k] - omega * t[k];
  }
  if ( rank == 0 ) for(k=np*size;k<n;k++){
      x[k] = x[k] + alpha*phat[k] + omega*shat[k];
      r[k] = s[k] - omega * t[k];
    }
}
