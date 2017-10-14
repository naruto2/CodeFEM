#if 0
#define __kernel
#define __global
#define __local
#define CLK_LOCAL_MEM_FENCE
#include <math.h>
int get_local_size(int);
int get_local_id(int);
void barrier();
#endif

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define NP 128
#define LOOP1 int k=i*size, end=(i+1)*size; k<end; k++
#define LOOP2 int j=row_ptr[k]; j<row_ptr[k+1]; j++
#define LOOP3 int k=np*size; k<n; k++


__kernel void _norm(int n,__global double *x,__global double *npa)
{
  int   np = get_local_size(0);
  int    i = get_local_id(0);
  int size = n/np;

  npa[i] = 0.0;
  for (LOOP1)  if( k)  npa[i] += x[k]*x[k];
  if (!i) for (LOOP3)  npa[i] += x[k]*x[k]; 
}


__kernel void _dot(int n,__global double *y, __global double *x,
       	       __global double *npa)
{
  int   np = get_local_size(0);
  int    i = get_local_id(0);
  int size = n/np;

  npa[i] = 0.0;
  for (LOOP1) if( k) npa[i] += y[k]*x[k];
  if(!i) for (LOOP3) npa[i] += y[k]*x[k];  
}


static void _copy(int n,__global double *y, __global double *x)
{
  int   np = get_local_size(0);
  int    i = get_local_id(0);
  int size = n/np;

  for (LOOP1) if( k) y[k] = x[k];
  if(!i) for (LOOP3) y[k] = x[k];
}


__kernel void gp_copy(int n,__global double *y, __global double *x)
{
  _copy(n,y,x);
}


static void _presolve_pointjacobi(int n,__global double *x,
	 	     __global double *dinv,   __global double *d)
{
  int   np = get_local_size(0);
  int    i = get_local_id(0);
  int size = n/np;

  for (LOOP1) if( k) x[k] = dinv[k]*d[k];
  if(!i) for (LOOP3) x[k] = dinv[k]*d[k];
}


__kernel void gp_presolve_pointjacobi(int n,__global double *x,
	 	     __global double *dinv,   __global double *d)
{
	  _presolve_pointjacobi(n,x,dinv,d);
}


static void _presolve(int n,__global double *x,
	 	     __global double *dinv,   __global double *d)
{
	if ( dinv[1] == 0.0 ) {
	   _copy(n,x,d);
    	}
	else {
	   _presolve_pointjacobi(n,x,dinv,d);
	}
}


__kernel void gp_presolve(int n,__global double *x,
	 	     __global double *dinv,   __global double *d)
{
	  _presolve(n,x,dinv,d);
}


static double _phase0(int n, __global double *r,
		   __global double *Aa, __global int *col_ind,
		   __global int *row_ptr, __global double *x,
		   __global double *rtilde, __global double *b)
{
  __local double npa[NP];	
  int   np = get_local_size(0);
  int    i = get_local_id(0);
  int size = n/np;
  double tmpa;

  npa[i] = 0.0;
  for (LOOP1) if( k) {
    r[k] = 0.0;
    for (LOOP2) r[k] += Aa[j] * x[col_ind[j]];
    rtilde[k] = r[k] = b[k] - r[k];
    npa[i] += r[k]*r[k];
  }
  if(!i) for (LOOP3) {
    r[k] = 0.0;
    for (LOOP2) r[k] += Aa[j] * x[col_ind[j]];
    rtilde[k] = r[k] = b[k] - r[k];
    npa[i] += r[k]*r[k];
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if (!i) {
        tmpa = 0.0;
	for (int k=0;k<np;k++) tmpa += npa[k]; 
	npa[0] = sqrt(tmpa);	
      }
  return npa[0];
}

__kernel void gp_phase0(int n, __global double *r,
		   __global double *Aa, __global int *col_ind,
		   __global int *row_ptr, __global double *x,
		   __global double *rtilde, __global double *b,
		   __global double *result)
{
	double tmpa = _phase0(n,r,Aa,col_ind,row_ptr,x,rtilde,b);
	if (!get_local_id(0) ) result[0] = tmpa;
}


static void _phase1(int n,__global double *p, __global double *r,
       		   __global double *v, double beta, double omega)
{
  int   np = get_local_size(0);
  int    i = get_local_id(0);
  int size = n/np;

  for (LOOP1) if( k) p[k] = r[k] + beta * (p[k] - omega *v[k]);
  if(!i) for (LOOP3) p[k] = r[k] + beta * (p[k] - omega *v[k]);
}


__kernel void gp_phase1(int n,__global double *p, __global double *r,
			__global double *v, double beta, double omega)
{
   _phase1(n,p,r,v,beta,omega);
}


static double _phase2(int n, __global double *v,
			__global double *Aa, __global int *col_ind,
			__global int *row_ptr, __global	double *phat,
			__global double *rtilde)
{
  __local double npa[NP];	
  int   np = get_local_size(0);
  int    i = get_local_id(0);
  int size = n/np;
  double tmpa;

  npa[i] = 0.0;
  for (LOOP1) if( k) {
     v[k] = 0.0;
     for (LOOP2) v[k] += Aa[j] * phat[col_ind[j]];
     npa[i] += rtilde[k]*v[k];
  }
  if(!i) for (LOOP3) {
     v[k] = 0.0;
     for (LOOP2) v[k] += Aa[j] * phat[col_ind[j]];
     npa[i] += rtilde[k]*v[k];
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if (!i) {
        tmpa = 0.0;
	for (int k=0;k<np;k++) tmpa += npa[k]; 
        npa[0] = tmpa;
    }
  return npa[0];
}


__kernel void gp_phase2(int n, __global double *v,
			__global double *Aa, __global int *col_ind,
			__global int *row_ptr, __global	double *phat,
			__global double *rtilde, __global double *result)
{

  double tmpa =  _phase2(n,v,Aa,col_ind,row_ptr,phat,rtilde);
  if (!get_local_id(0) ) result[0] = tmpa;
}


static double _phase3(int n,__global double *s, __global double *r,
		__global double *v, double alpha)
{
  __local double npa[NP];	
  int   np = get_local_size(0);
  int    i = get_local_id(0);
  int size = n/np;
  double tmpa;

  npa[i] = 0.0;
  for (LOOP1) if( k) {
    s[k] = r[k] - alpha * v[k];
    npa[i] += s[k]*s[k];
  }
  if(!i) for (LOOP3) {
    s[k] = r[k] - alpha * v[k];
    npa[i] += s[k]*s[k];
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if (!i) {
        tmpa = 0.0;
	for (int k=0;k<np;k++) tmpa += npa[k]; 
	npa[0] = sqrt(tmpa);
    }
  return npa[0];
}


__kernel void gp_phase3(int n,__global double *s, __global double *r,
		__global double *v, double alpha, __global double *result)
{
	double tmpa = _phase3(n,s,r,v,alpha);
	if ( !get_local_id(0) ) result[0] = tmpa;
}


static void _phase4(int n,__global double *x, __global double *phat,
			double alpha)
{
  int   np = get_local_size(0);
  int    i = get_local_id(0);
  int size = n/np;

  for (LOOP1) if( k) x[k] = x[k] + alpha*phat[k];
  if(!i) for (LOOP3) x[k] = x[k] + alpha*phat[k];
}


__kernel void gp_phase4(int n,__global double *x, __global double *phat,
			double alpha)
{
     _phase4(n,x,phat,alpha);
}


static double _phase5(int n,__global double *t,__global double *Aa,
	 	 __global int *col_ind,__global int *row_ptr,
		 __global double *shat,__global double *s)
{
  __local double npa[NP];
  __local double npb[NP];	
  int   np = get_local_size(0);
  int    i = get_local_id(0);
  int size = n/np;

  npa[i] = 0.0;
  npb[i] = 0.0;
  for (LOOP1) if( k) {
    t[k] = 0.0;
    for (LOOP2) t[k] += Aa[j] * shat[col_ind[j]];
    npa[i] += t[k]*s[k];
    npb[i] += t[k]*t[k];
  }
  if(!i) for (LOOP3) {
    t[k] = 0.0;
    for (LOOP2) t[k] += Aa[j] * shat[col_ind[j]];
    npa[i] += t[k]*s[k];
    npb[i] += t[k]*t[k];
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if (!i) {
     	double tmpa = 0.0, tmpb = 0.0;
	for (int k=0;k<np;k++) { tmpa += npa[k]; tmpb += npb[k]; }
	npa[0] = tmpa/tmpb;
    }
  return npa[0];
}


__kernel void gp_phase5(int n,__global double *t,__global double *Aa,
	 	 __global int *col_ind,__global int *row_ptr,
		 __global double *shat,__global double *s,
		 __global double *result,__global double *resuls)
{
    double tmpa = _phase5(n,t,Aa,col_ind,row_ptr,shat,s);
    if (!get_local_id(0) ) result[0] = tmpa;
}


static double _phase6(int n,__global double *x, __global double *s,
			__global double *r, __global double *t,
			__global double *phat, __global double *shat,
			double alpha, double omega)
{
  __local double npa[NP];
  int   np = get_local_size(0);
  int    i = get_local_id(0); 
  int size = n/np;

  npa[i] = 0.0;
  for (LOOP1) if( k) {
    x[k] = x[k] + alpha*phat[k] + omega*shat[k];
    r[k] = s[k] - omega * t[k];
    npa[i] += r[k]*r[k];
  }
  if(!i) for (LOOP3) {
    x[k] = x[k] + alpha*phat[k] + omega*shat[k];
    r[k] = s[k] - omega * t[k];
    npa[i] += r[k]*r[k];
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if (!i) {
        double tmpa = 0.0;
	for (int k=0;k<np;k++) tmpa += npa[k]; 
	npa[0] = sqrt(tmpa);
    }
  return npa[0];
}


__kernel void gp_phase6(int n,__global double *x, __global double *s,
			__global double *r, __global double *t,
			__global double *phat, __global double *shat,
			double alpha, double omega, __global double *result)
{
    double tmpa = _phase6(n,x,s,r,t,phat,shat,alpha,omega);
    if (!get_local_id(0) ) result[0] = tmpa;
}
