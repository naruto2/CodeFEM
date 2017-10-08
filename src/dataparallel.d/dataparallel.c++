#include <stdlib.h>
#include <CL/cl.h>
#include <stdio.h>
#include <cmath>

#define MAX_SOURCE_SIZE (0x100000)

static cl_context context = NULL;
static cl_command_queue command_queue = NULL;
static cl_uint work_dim = 1;
static size_t global_item_size[3];
static size_t local_item_size[3];
static size_t np;
static cl_program program = NULL;
static char *kernel_src_str;

static void check_np(int n)
{
  if ( global_item_size[0] < 0 )
    global_item_size[0] = 1024;
  if ( 1024 < global_item_size[0])
    global_item_size[0] = 1024;
  if ( n<=global_item_size[0])
    global_item_size[0] = n-1;
}

 
static int cl_load(cl_kernel &kernel, const char *name)
{
  int ret;
  kernel = clCreateKernel(program, name, &ret);
  if (ret) { fprintf(stderr,"%s=%d\n",name,ret); abort(); }
  return ret;
}

static int cl_mem_r(int size, cl_mem &mem)
{
  int ret;
  if (!mem)  mem = clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret);
  if (ret) { fprintf(stderr,"clCreateBuffer()=%d\n",ret); abort();}
  return ret;
}

static int cl_mem_w(int size, cl_mem &mem)
{
  int ret;
  if (!mem) mem = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size, NULL, &ret);
  if (ret) { fprintf(stderr,"clCreateBuffer()=%d\n",ret); abort();}
  return ret;
}

static int cl_mem_rw(int size, cl_mem &mem)
{
  int ret;
  if(!mem) mem = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
  if (ret) { fprintf(stderr,"clCreateBuffer()=%d\n",ret); abort();}
  return ret;
}


static int cl_run(cl_kernel kernel)
{
  int ret;
  ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
				  global_item_size, local_item_size,
				  0, NULL, NULL);
  if (ret) { fprintf(stderr,"clEnqueueNDRangeKernel()=%d\n",ret); abort();}
  return ret;
}


static int cl_send(size_t size, cl_mem &mem_x, void *x)
{
  int ret;
  ret = clEnqueueWriteBuffer(command_queue, mem_x, CL_TRUE, 0,
			     size, x, 0, NULL, NULL);
  if (ret) { fprintf(stderr,"clEnqueWriteBuffer()=%d\n",ret); abort();}
  return ret;
}


static int cl_get(size_t size, cl_mem &mem_x, void *x)
{
  int ret;
  ret = clEnqueueReadBuffer(command_queue, mem_x, CL_TRUE, 0,
			    size, x,0, NULL, NULL);
  if (ret) { fprintf(stderr,"clEnqueReadBuffer()=%d\n",ret); abort(); }
  return ret;
}


double cl_norm(int n, double *x)
{
  check_np(n);

  static cl_kernel kernel;
  cl_load(kernel,"cl_norm");

  static double *npa;
  if (!npa) npa = (double*)malloc(np*sizeof(double));

  static cl_mem mem_x;
  static cl_mem mem_npa;

  cl_mem_r(n*sizeof(double), mem_x);
  cl_mem_w(np*sizeof(double), mem_npa);
  
  cl_send(n*sizeof(double), mem_x, x);

  clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_x);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_npa);

  cl_run(kernel);

  cl_get(np*sizeof(double), mem_npa, npa);

  for(int k=1;k<np;k++) npa[0] += npa[k];
  
  return sqrt(npa[0]);
}


void cl_copy(int n, double *y, double *x)
{
  check_np(n);
  
  static cl_kernel  kernel;
  cl_load(kernel,"cl_copy");

  static cl_mem mem_y;
  static cl_mem mem_x;
  cl_mem_w(n*sizeof(double), mem_y);
  cl_mem_r(n*sizeof(double), mem_x);

  cl_send(n*sizeof(double), mem_x, x);

  clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_y);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_x);

  cl_run(kernel);
          
  cl_get(n*sizeof(double), mem_y, y);
}


double cl_dot(int n, double *y, double *x)
{
  check_np(n);
  
  static cl_kernel kernel;
  cl_load(kernel,"cl_dot");

  static double *npa;
  if (!npa) npa = (double*)malloc(np*sizeof(double));

  static cl_mem mem_y;
  static cl_mem mem_x;
  static cl_mem mem_npa;
  cl_mem_r(n*sizeof(double), mem_y);
  cl_mem_r(n*sizeof(double), mem_x);
  cl_mem_w(np*sizeof(double), mem_npa);
  
  cl_send(n*sizeof(double), mem_y, y);
  cl_send(n*sizeof(double), mem_x, x);

  clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_y);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_x);
  clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_npa);

  cl_run(kernel);

  cl_get(np*sizeof(double), mem_npa, npa);            

  for(int k=1;k<np;k++) npa[0] += npa[k];
    
  return npa[0];
}


void cl_matrixvector(int n, double *r, double *Aa, int *col_ind,
		     int *row_ptr, double *b, int w)
{
  static int ww;
  int ret;
  check_np(n);

  static cl_kernel kernel;
  ret = cl_load(kernel,"cl_matrixvector");

  static cl_mem mem_r;
  static cl_mem mem_Aa;
  static cl_mem mem_col_ind;
  static cl_mem mem_row_ptr;
  static cl_mem mem_b;

  cl_mem_w(n*sizeof(double), mem_r);

  if ( ww < w){
    if(mem_Aa)
      { ret = clReleaseMemObject(mem_Aa); mem_Aa = NULL; }
    ret = cl_mem_r(w*sizeof(double), mem_Aa);
    if(mem_col_ind)
      { ret = clReleaseMemObject(mem_col_ind); mem_col_ind = NULL; }
    ret = cl_mem_r(w*sizeof(double), mem_col_ind);
    ww = w;
  }

  cl_mem_r((n+1)*sizeof(double), mem_row_ptr);
  cl_mem_r(n*sizeof(double), mem_b);

  ret = clEnqueueWriteBuffer(command_queue, mem_Aa, CL_TRUE, 0,
			       w*sizeof(double),
                               Aa, 0, NULL, NULL);
  ret = clEnqueueWriteBuffer(command_queue, mem_col_ind, CL_TRUE, 0,
			       w*sizeof(int),
                               col_ind, 0, NULL, NULL);
  ret = clEnqueueWriteBuffer(command_queue, mem_row_ptr, CL_TRUE, 0,
			       (n+1)*sizeof(int),
                               row_ptr, 0, NULL, NULL);

  ret = clEnqueueWriteBuffer(command_queue, mem_b, CL_TRUE, 0,
			       n*sizeof(double),
                               b, 0, NULL, NULL);

  ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
  ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_r);
  ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_Aa);
  ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_col_ind);
  ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&mem_row_ptr);
  ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&mem_b);
    
  ret = cl_run(kernel);

  ret = clEnqueueReadBuffer(command_queue, mem_r, CL_TRUE, 0,
                              n * sizeof(double),
                              r,0, NULL, NULL);
}


static int ww = 0;
static cl_mem mem_Aa= NULL;
static cl_mem mem_col_ind = NULL;
static cl_mem mem_row_ptr = NULL;

static int cl_send_A(int n,int w, double *Aa, int *col_ind, int *row_ptr)
{
  if ( ww < w ) {
    if(mem_Aa)
      { clReleaseMemObject(mem_Aa); mem_Aa = NULL; }
    cl_mem_r(w*sizeof(double), mem_Aa);
    if(mem_col_ind)
      { clReleaseMemObject(mem_col_ind); mem_col_ind = NULL; }
    cl_mem_r(w*sizeof(double), mem_col_ind);

    cl_send(w*sizeof(double), mem_Aa, Aa);
    cl_send(w*sizeof(int), mem_col_ind, col_ind);
    
    cl_mem_r((n+1)*sizeof(int), mem_row_ptr);
    cl_send((n+1)*sizeof(int), mem_row_ptr, row_ptr);

    ww = w;
  }
  return 0;
}


double cl_phase0(int n, double *r, double *Aa, int *col_ind,
	       int *row_ptr, double *x, double *rtilde,
	       double *b, int w)
{
  check_np(n);

  static cl_kernel kernel;
  cl_load(kernel,"cl_phase0");

  cl_send_A(n,w,Aa,col_ind,row_ptr);

  static double *npa;
  if (!npa) npa = (double*)malloc(np*sizeof(double));

  static cl_mem mem_r;
  static cl_mem mem_x;
  static cl_mem mem_rtilde;
  static cl_mem mem_b;
  static cl_mem mem_npa;

  cl_mem_rw(n*sizeof(double), mem_r);
  cl_mem_r(n*sizeof(double), mem_x);
  cl_mem_rw(n*sizeof(double), mem_rtilde);
  cl_mem_r(n*sizeof(double), mem_b);
  cl_mem_w(np*sizeof(double), mem_npa);
  
  cl_send(n*sizeof(double), mem_r, r);
  cl_send(n*sizeof(double), mem_x, x);
  cl_send(n*sizeof(double), mem_rtilde, rtilde);
  cl_send(n*sizeof(double), mem_b, b);

  clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_r);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_Aa);
  clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_col_ind);
  clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&mem_row_ptr);
  clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&mem_x);
  clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&mem_rtilde);
  clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&mem_b);
  clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&mem_npa);
    
  cl_run(kernel);

  cl_get(n*sizeof(double), mem_r, r);
  cl_get(n*sizeof(double), mem_rtilde, rtilde);
  cl_get(np*sizeof(double), mem_npa, npa);

  for(int k=1;k<np;k++) npa[0] += npa[k];

  return sqrt(npa[0]);
}


void cl_phase1(int n, double *p, double *r, double *v,
	       double beta, double omega)
{
  check_np(n);

  static cl_kernel kernel;
  cl_load(kernel,"cl_phase1");

  static cl_mem mem_p;
  static cl_mem mem_r;
  static cl_mem mem_v;

  cl_mem_rw(n*sizeof(double), mem_p);
  cl_mem_r(n*sizeof(double), mem_r);
  cl_mem_r(n*sizeof(double), mem_v);

  cl_send(n*sizeof(double), mem_p, p);
  cl_send(n*sizeof(double), mem_r, r);
  cl_send(n*sizeof(double), mem_v, v);

  clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_p);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_r);
  clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_v);
  clSetKernelArg(kernel, 4, sizeof(double), (void *)&beta);
  clSetKernelArg(kernel, 5, sizeof(double), (void *)&omega);

  cl_run(kernel);

  cl_get(n*sizeof(double), mem_p, p);
}


double cl_phase2(int n, double *v,
		 double *Aa, int *col_ind, int *row_ptr,
		 double *phat, double *rtilde, int w)
{
  int k, ret;
  static cl_kernel kernel = NULL;  
  check_np(n);
  ret = cl_load(kernel,"cl_phase2");

  static double *npa = NULL;
  static cl_mem mem_v = NULL;
  static cl_mem mem_phat = NULL;
  static cl_mem mem_rtilde = NULL;
  static cl_mem mem_npa = NULL;

  if(!npa) npa = (double*)malloc(sizeof(double)*global_item_size[0]);

  ret = cl_mem_w(n*sizeof(double), mem_v);

    if ( ww < w) {
    if(mem_Aa)
      { ret = clReleaseMemObject(mem_Aa); mem_Aa = NULL; }
    ret = cl_mem_r(w*sizeof(double), mem_Aa);
    if(mem_col_ind)
      { ret = clReleaseMemObject(mem_col_ind); mem_col_ind = NULL; }
    ret = cl_mem_r(w*sizeof(double), mem_col_ind);

    ret = clEnqueueWriteBuffer(command_queue, mem_Aa, CL_TRUE, 0,
			       w*sizeof(double),
                               Aa, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_col_ind, CL_TRUE, 0,
			       w*sizeof(int),
                               col_ind, 0, NULL, NULL);

    ret = cl_mem_r((n+1)*sizeof(int), mem_row_ptr);

    ret = clEnqueueWriteBuffer(command_queue, mem_row_ptr, CL_TRUE, 0,
				 (n+1)*sizeof(int),
				 row_ptr, 0, NULL, NULL);
    ww = w;
  }

    ret = cl_mem_r(n*sizeof(double), mem_phat);
    ret = cl_mem_r(n*sizeof(double), mem_rtilde);
    ret = cl_mem_r(global_item_size[0]*sizeof(double), mem_npa);
  
    ret = clEnqueueWriteBuffer(command_queue, mem_phat, CL_TRUE, 0,
			       n*sizeof(double),
                               phat, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_rtilde, CL_TRUE, 0,
			       n*sizeof(double),
                               rtilde, 0, NULL, NULL);

    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_v);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_Aa);
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_col_ind);
    ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&mem_row_ptr);
    ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&mem_phat);
    ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&mem_rtilde);
    ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&mem_npa);
    
    ret = cl_run(kernel);

    ret = clEnqueueReadBuffer(command_queue, mem_v, CL_TRUE, 0,
                              n * sizeof(double),
                              v,0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_npa, CL_TRUE, 0,
			      global_item_size[0] * sizeof(double),
			      npa,0, NULL, NULL);

    for(k=1;k<global_item_size[0];k++) npa[0] += npa[k];

    return npa[0];
}


double cl_phase3(int n, double *s, double *r, double *v,
	       double alpha)
{
  int k, ret;
  static cl_kernel kernel = NULL;
  check_np(n);
  ret = cl_load(kernel,"cl_phase3");

  static double *npa = NULL;
  static cl_mem mem_s = NULL;
  static cl_mem mem_r = NULL;
  static cl_mem mem_v = NULL;
  static cl_mem mem_npa=NULL;

  if(!npa) npa = (double*)malloc(sizeof(double)*global_item_size[0]);

  ret = cl_mem_rw(n*sizeof(double), mem_s);
  ret = cl_mem_r(n*sizeof(double), mem_r);
  ret = cl_mem_r(n*sizeof(double), mem_v);
  ret = cl_mem_w(global_item_size[0]*sizeof(double), mem_npa);
  
    ret = clEnqueueWriteBuffer(command_queue, mem_s, CL_TRUE, 0,
			       n*sizeof(double),
                               s, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_r, CL_TRUE, 0,
			       n*sizeof(double),
                               r, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_v, CL_TRUE, 0,
			       n*sizeof(double),
                               v, 0, NULL, NULL);
  
    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_s);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_r);
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_v);
    ret = clSetKernelArg(kernel, 4, sizeof(double), (void *)&alpha);
    ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&mem_npa);
    
    
    ret = cl_run(kernel);

    ret = clEnqueueReadBuffer(command_queue, mem_s, CL_TRUE, 0,
                              n * sizeof(double),
                              s,0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_npa, CL_TRUE, 0,
                              global_item_size[0] * sizeof(double),
                              npa,0, NULL, NULL);

    for(k=1;k<global_item_size[0];k++) npa[0] += npa[k];
    return sqrt(npa[0]);
}


void cl_phase4(int n, double *s, double *phat, double alpha)
{
  int ret;
  static cl_kernel kernel = NULL;
  check_np(n);
  ret = cl_load(kernel,"cl_phase4");
  static cl_mem mem_s = NULL;
  static cl_mem mem_phat = NULL;

  ret = cl_mem_rw(n*sizeof(double), mem_s);
  ret = cl_mem_r(n*sizeof(double), mem_phat);
  
    ret = clEnqueueWriteBuffer(command_queue, mem_s, CL_TRUE, 0,
			       n*sizeof(double),
                               s, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_phat, CL_TRUE, 0,
			       n*sizeof(double),
                               phat, 0, NULL, NULL);

    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_s);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_phat);
    ret = clSetKernelArg(kernel, 3, sizeof(double), (void *)&alpha);

    ret = cl_run(kernel);

    ret = clEnqueueReadBuffer(command_queue, mem_s, CL_TRUE, 0,
                              n * sizeof(double),
                              s,0, NULL, NULL);
}

double cl_phase5(int n, double *t, double *Aa, int *col_ind,
	       int *row_ptr, double *shat, double *s, int w)
{
  int k, ret;
  static cl_kernel kernel = NULL;
  check_np(n);
  ret = cl_load(kernel,"cl_phase5");
  
  static double *npa = NULL;
  static double *npb = NULL;
  
  static cl_mem mem_t = NULL;
  static cl_mem mem_r = NULL;
  static cl_mem mem_shat = NULL;
  static cl_mem mem_s = NULL;
  static cl_mem mem_npa = NULL;
  static cl_mem mem_npb = NULL;
  if(!npa) npa = (double*)malloc(sizeof(double)*global_item_size[0]);
  if(!npb) npb = (double*)malloc(sizeof(double)*global_item_size[0]);

  ret = cl_mem_rw(n*sizeof(double), mem_t);

  if ( ww < w) {
    if(mem_Aa)
      { ret = clReleaseMemObject(mem_Aa); mem_Aa = NULL; }
    ret = cl_mem_r(w*sizeof(double), mem_Aa);
    if(mem_col_ind)
      { ret = clReleaseMemObject(mem_col_ind); mem_col_ind = NULL; }
    ret = cl_mem_r(w*sizeof(double), mem_col_ind);

    ret = clEnqueueWriteBuffer(command_queue, mem_Aa, CL_TRUE, 0,
			       w*sizeof(double),
                               Aa, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_col_ind, CL_TRUE, 0,
			       w*sizeof(int),
                               col_ind, 0, NULL, NULL);

    ret = cl_mem_r((n+1)*sizeof(int), mem_row_ptr);

    ret = clEnqueueWriteBuffer(command_queue, mem_row_ptr, CL_TRUE, 0,
				 (n+1)*sizeof(int),
				 row_ptr, 0, NULL, NULL);
    ww = w;
  }
  ret = cl_mem_r(n*sizeof(double), mem_shat);
  ret = cl_mem_r(n*sizeof(double), mem_s);
  ret = cl_mem_w(global_item_size[0]*sizeof(double), mem_npa);
  ret = cl_mem_w(global_item_size[0]*sizeof(double), mem_npb);
  
    ret = clEnqueueWriteBuffer(command_queue, mem_r, CL_TRUE, 0,
			       n*sizeof(double),
                               t, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_shat, CL_TRUE, 0,
			       n*sizeof(double),
                               shat, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_s, CL_TRUE, 0,
			       n*sizeof(double),
                               s, 0, NULL, NULL);

    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_t);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_Aa);
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_col_ind);
    ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&mem_row_ptr);
    ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&mem_shat);
    ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&mem_s);
    ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&mem_npa);
    ret = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&mem_npb);

    ret = cl_run(kernel);

    ret = clEnqueueReadBuffer(command_queue, mem_t, CL_TRUE, 0,
                              n * sizeof(double),
                              t,0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_npa, CL_TRUE, 0,
			      global_item_size[0] * sizeof(double),
			      npa,0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_npb, CL_TRUE, 0,
			      global_item_size[0] * sizeof(double),
			      npb,0, NULL, NULL);

    for(k=1;k<global_item_size[0];k++) {
      npa[0] += npa[k];
      npb[0] += npb[k];
    }
      
    return npa[0]/npb[0];
}


double cl_phase6(int n, double *x, double *s, double *r, double *t,
		 double *phat, double *shat, double alpha, double omega)
{
  int k, ret;
  static cl_kernel kernel = NULL;
  check_np(n);
  ret = cl_load(kernel,"cl_phase6");
  static double *npa = NULL;

  static cl_mem mem_x = NULL;
  static cl_mem mem_s = NULL;
  static cl_mem mem_r = NULL;
  static cl_mem mem_t = NULL;
  static cl_mem mem_phat = NULL;
  static cl_mem mem_shat = NULL;
  static cl_mem mem_npa  = NULL;
  

  if(!npa) npa = (double*)malloc(sizeof(double)*global_item_size[0]);
  
  ret = cl_mem_w(global_item_size[0]*sizeof(double), mem_npa);
  ret = cl_mem_rw(n*sizeof(double), mem_x);
  ret = cl_mem_r(n*sizeof(double), mem_s);
  ret = cl_mem_w(n*sizeof(double), mem_r);
  ret = cl_mem_r(n*sizeof(double), mem_t);
  ret = cl_mem_r(n*sizeof(double), mem_phat);
  ret = cl_mem_r(n*sizeof(double), mem_shat);

    ret = clEnqueueWriteBuffer(command_queue, mem_x, CL_TRUE, 0,
			       n*sizeof(double),
                               x, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_s, CL_TRUE, 0,
			       n*sizeof(double),
                               s, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_t, CL_TRUE, 0,
			       n*sizeof(double),
                               t, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_phat, CL_TRUE, 0,
			       n*sizeof(double),
                               phat, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_shat, CL_TRUE, 0,
			       n*sizeof(double),
                               shat, 0, NULL, NULL);

    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_x);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_s);
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_r);
    ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&mem_t);
    ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&mem_phat);
    ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&mem_shat);
    ret = clSetKernelArg(kernel, 7, sizeof(double), (void *)&alpha);
    ret = clSetKernelArg(kernel, 8, sizeof(double), (void *)&omega);
    ret = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *)&mem_npa);
    
    ret = cl_run(kernel);

    ret = clEnqueueReadBuffer(command_queue, mem_x, CL_TRUE, 0,
                              n * sizeof(double),
                              x,0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_r, CL_TRUE, 0,
                              n * sizeof(double),
                              r,0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_npa, CL_TRUE, 0,
			      global_item_size[0] * sizeof(double),
			      npa,0, NULL, NULL);

    for(k=1;k<global_item_size[0];k++) npa[0] += npa[k];
    return sqrt(npa[0]);
}


void cl_bicgstab_init(int argc, char **argv)
{

  if (kernel_src_str != NULL ) return;
  cl_platform_id platform_id = NULL;
  cl_uint ret_num_platforms;
  cl_device_id device_id = NULL;
  cl_uint ret_num_devices;
  size_t kernel_code_size;
  cl_int ret;
  FILE *fp;
    
  /* Set parameters for data parallel processing (work item) */

  global_item_size[0] = 64; /* np Global number of work items */
  local_item_size[0] = 1;  /* Number of work items per work group */

  if ( argc > 1 ) global_item_size[0] = atoi(argv[1]);    
  np = global_item_size[0];

  /* --> global_item_size[0] / local_item_size[0] becomes 2, which indirectly sets the number of workgroups to 2*/


  /* Allocate space to read in kernel code */
  kernel_src_str = (char *)malloc(MAX_SOURCE_SIZE);
    
  /* Get Platform*/
  ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);

  /* Get Device */
  ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id,
		       &ret_num_devices);

  /* Create Context */
  context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);


#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  /* Create Command Queue */  
  command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
#pragma GCC diagnostic warning "-Wdeprecated-declarations"
  
  
  /* Read kernel source code */     
  fp = fopen("/usr/include/est/cl_bicgstab_kernel.cl", "r");
  kernel_code_size = fread(kernel_src_str, 1, MAX_SOURCE_SIZE, fp);
  fclose(fp);
  
  /* Create Program Object */
  program = clCreateProgramWithSource(context, 1, (const char **)&kernel_src_str,
				      (const size_t *)&kernel_code_size, &ret);

  /* Compile kernel */         
  ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
}


void cl_finalize(void)
{
    /* OpenCL Object Finalization */
    int ret;
    ret = clReleaseProgram(program);
    ret = clReleaseCommandQueue(command_queue);
    ret = clReleaseContext(context);

    /* Deallocate memory on the host */
    free(kernel_src_str);
}

