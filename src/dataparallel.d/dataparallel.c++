#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <CL/cl.h>

#define MAX_SOURCE_SIZE (0x100000)
#define NP 128

static cl_context context = NULL;
static cl_command_queue command_queue = NULL;
static cl_uint work_dim = 1;
static size_t global_item_size[3];
static size_t local_item_size[3];
static size_t np;
static cl_program program = NULL;
static char *kernel_src_str;
static int ret;
#define cl_load(kernel, name) if(!kernel) kernel=clCreateKernel(program, name, &ret)
#define cl_mem_r(size,mem) if(!mem) mem=clCreateBuffer(context, CL_MEM_READ_ONLY, size, NULL, &ret)
#define cl_mem_w(size,mem) if(!mem) mem=clCreateBuffer(context, CL_MEM_WRITE_ONLY, size, NULL, &ret)
#define cl_mem_rw(size,mem) if(!mem) mem=clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret)


static void check_np(int n)
{
  if ( global_item_size[0] < 0 )
    global_item_size[0] = NP;
  if ( NP < global_item_size[0])
    global_item_size[0] = NP;
  if ( n<=global_item_size[0])
    np = global_item_size[0] = n-1;
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

static int ww = 0;
static cl_mem mem_Aa= NULL;
static cl_mem mem_col_ind = NULL;
static cl_mem mem_row_ptr = NULL;

int cl_send_A(int n,int w, double *Aa, int *col_ind, int *row_ptr)
{
  if ( ww < w ) {
    if(mem_Aa)
      { clReleaseMemObject(mem_Aa); mem_Aa = NULL; }
    cl_mem_r(w*sizeof(double), mem_Aa);
    if(mem_col_ind)
      { clReleaseMemObject(mem_col_ind); mem_col_ind = NULL; }
    cl_mem_r(w*sizeof(double), mem_col_ind);
    cl_mem_r((n+1)*sizeof(int), mem_row_ptr);
    ww = w;
  }
  cl_send(w*sizeof(double), mem_Aa, Aa);
  cl_send(w*sizeof(int), mem_col_ind, col_ind);
  cl_send((n+1)*sizeof(int), mem_row_ptr, row_ptr);
  return 0;
}


void cl_matrixvector(int n, double *r, double *Aa, int *col_ind,
		     int *row_ptr, double *b, int w)
{
  check_np(n);
  
  static cl_kernel kernel;
  cl_load(kernel,"cl_matrixvector");

  static cl_mem mem_r;
  static cl_mem mem_b;

  cl_mem_w(n*sizeof(double), mem_r);
  cl_mem_r(n*sizeof(double), mem_b);

  cl_send(n*sizeof(double), mem_b, b);

  clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_r);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_Aa);
  clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_col_ind);
  clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&mem_row_ptr);
  clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&mem_b);
    
  cl_run(kernel);

  cl_get(n*sizeof(double), mem_r, r);
}


double cl_phase0(int n, double *r, double *Aa, int *col_ind,
		 int *row_ptr, double *x, double *rtilde,
		 double *b, int w)
{
  check_np(n);

  static cl_kernel kernel;
  cl_load(kernel,"cl_phase0");

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
  check_np(n);

  static cl_kernel kernel;
  cl_load(kernel,"cl_phase2");

  static double *npa = NULL;
  if (!npa) npa = (double*)malloc(np*sizeof(double));

  static cl_mem mem_v = NULL;
  static cl_mem mem_phat = NULL;
  static cl_mem mem_rtilde = NULL;
  static cl_mem mem_npa = NULL;

  cl_mem_w(n*sizeof(double), mem_v);
  cl_mem_r(n*sizeof(double), mem_phat);
  cl_mem_r(n*sizeof(double), mem_rtilde);
  cl_mem_r(np*sizeof(double), mem_npa);
  
  cl_send(n*sizeof(double), mem_phat, phat);
  cl_send(n*sizeof(double), mem_rtilde, rtilde);

  clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_v);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_Aa);
  clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_col_ind);
  clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&mem_row_ptr);
  clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&mem_phat);
  clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&mem_rtilde);
  clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&mem_npa);
  
  cl_run(kernel);

  cl_get(n*sizeof(double), mem_v, v);
  cl_get(np*sizeof(double), mem_npa, npa);

  for(int k=1;k<np;k++) npa[0] += npa[k];

  return npa[0];
}


double cl_phase3(int n, double *s, double *r, double *v, double alpha)
{
  check_np(n);
  
  static cl_kernel kernel;
  cl_load(kernel,"cl_phase3");

  static double *npa;
  if (!npa) npa = (double*)malloc(np*sizeof(double));

  static cl_mem mem_s = NULL;
  static cl_mem mem_r = NULL;
  static cl_mem mem_v = NULL;
  static cl_mem mem_npa=NULL;

  cl_mem_rw(n*sizeof(double), mem_s);
  cl_mem_r(n*sizeof(double), mem_r);
  cl_mem_r(n*sizeof(double), mem_v);
  cl_mem_w(np*sizeof(double), mem_npa);
  
  cl_send(n*sizeof(double), mem_s, s);
  cl_send(n*sizeof(double), mem_r, r);
  cl_send(n*sizeof(double), mem_v, v);
  
  clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_s);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_r);
  clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_v);
  clSetKernelArg(kernel, 4, sizeof(double), (void *)&alpha);
  clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&mem_npa);
    
  cl_run(kernel);

  cl_get(n*sizeof(double), mem_s, s);
  cl_get(np*sizeof(double), mem_npa, npa);

  for(int k=1;k<np;k++) npa[0] += npa[k];
  return sqrt(npa[0]);
}


void cl_phase4(int n, double *s, double *phat, double alpha)
{
  check_np(n);
  
  static cl_kernel kernel;
  cl_load(kernel,"cl_phase4");

  static cl_mem mem_s;
  static cl_mem mem_phat;

  cl_mem_rw(n*sizeof(double), mem_s);
  cl_mem_r(n*sizeof(double), mem_phat);
  
  cl_send(n*sizeof(double), mem_s, s);
  cl_send(n*sizeof(double), mem_phat, phat);
  
  clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_s);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_phat);
  clSetKernelArg(kernel, 3, sizeof(double), (void *)&alpha);

  cl_run(kernel);

  cl_get(n*sizeof(double), mem_s, s);
}


double cl_phase5(int n, double *t, double *Aa, int *col_ind,
	       int *row_ptr, double *shat, double *s, int w)
{
  check_np(n);

  static cl_kernel kernel;
  cl_load(kernel,"cl_phase5");

  static double *npa;
  static double *npb;

  if (!npa) npa = (double*)malloc(np*sizeof(double));
  if (!npb) npb = (double*)malloc(np*sizeof(double));
  
  static cl_mem mem_t;
  static cl_mem mem_shat;
  static cl_mem mem_s;
  static cl_mem mem_npa;
  static cl_mem mem_npb;

  cl_mem_rw(n*sizeof(double), mem_t);
  cl_mem_r(n*sizeof(double), mem_shat);
  cl_mem_r(n*sizeof(double), mem_s);
  cl_mem_w(np*sizeof(double), mem_npa);
  cl_mem_w(np*sizeof(double), mem_npb);
  
  cl_send(n*sizeof(double), mem_t, t);
  cl_send(n*sizeof(double), mem_shat, shat);
  cl_send(n*sizeof(double), mem_s, s);

  clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_t);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_Aa);
  clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_col_ind);
  clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&mem_row_ptr);
  clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&mem_shat);
  clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&mem_s);
  clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&mem_npa);
  clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&mem_npb);

  cl_run(kernel);

  cl_get(n*sizeof(double), mem_t, t);
  cl_get(np*sizeof(double), mem_npa, npa);
  cl_get(np*sizeof(double), mem_npb, npb);

  for(int k=1;k<np;k++) {
    npa[0] += npa[k];
    npb[0] += npb[k];
  }
  return npa[0]/npb[0];
}


double cl_phase6(int n, double *x, double *s, double *r, double *t,
		 double *phat, double *shat, double alpha, double omega)
{
  check_np(n);

  static cl_kernel kernel;
  cl_load(kernel,"cl_phase6");

  static double *npa;
  if (!npa) npa = (double*)malloc(np*sizeof(double));

  static cl_mem mem_x;
  static cl_mem mem_s;
  static cl_mem mem_r;
  static cl_mem mem_t;
  static cl_mem mem_phat;
  static cl_mem mem_shat;
  static cl_mem mem_npa; 
  
  cl_mem_rw(n*sizeof(double), mem_x);
  cl_mem_r(n*sizeof(double), mem_s);
  cl_mem_w(n*sizeof(double), mem_r);
  cl_mem_r(n*sizeof(double), mem_t);
  cl_mem_r(n*sizeof(double), mem_phat);
  cl_mem_r(n*sizeof(double), mem_shat);
  cl_mem_w(np*sizeof(double), mem_npa);
  
  cl_send(n*sizeof(double), mem_x, x);
  cl_send(n*sizeof(double), mem_s, s);
  cl_send(n*sizeof(double), mem_t, t);
  cl_send(n*sizeof(double), mem_phat, phat);
  cl_send(n*sizeof(double), mem_shat, shat);

  clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_x);
  clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_s);
  clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_r);
  clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&mem_t);
  clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&mem_phat);
  clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&mem_shat);
  clSetKernelArg(kernel, 7, sizeof(double), (void *)&alpha);
  clSetKernelArg(kernel, 8, sizeof(double), (void *)&omega);
  clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *)&mem_npa);
    
  cl_run(kernel);

  cl_get(n*sizeof(double), mem_x, x);
  cl_get(n*sizeof(double), mem_r, r);
  cl_get(np*sizeof(double), mem_npa, npa);

  for(int k=1;k<np;k++) npa[0] += npa[k];
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

  global_item_size[0] = NP; /* np Global number of work items */
  local_item_size[0] =  NP;

  if ( argc > 1 ) global_item_size[0] = atoi(argv[1]);    
  np = local_item_size[0] = global_item_size[0];

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

  /* clGetDeviceInfoで設定可能なインデックス空間の情報を得る */
  /* 株式会社フィックスターズ: OpenCL入門, インプレスジャパン p187, 2010. */
  cl_uint work_item_dim;
  size_t work_item_sizes[3];
  size_t work_group_size;
  cl_uint compute_unit = 0;
  cl_ulong local_mem_size = 0;
  
  ret = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
			sizeof(cl_uint), &work_item_dim, NULL);
  if (ret) abort();
  ret = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_SIZES,
			sizeof(work_item_sizes), work_item_sizes, NULL);
  if (ret) abort();
  ret = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE,
			sizeof(size_t), &work_group_size, NULL);
  if (ret) abort();
  ret = clGetDeviceInfo(device_id, CL_DEVICE_MAX_COMPUTE_UNITS,
			sizeof(cl_uint), &compute_unit, NULL);
  if (ret) abort();
  ret = clGetDeviceInfo(device_id, CL_DEVICE_LOCAL_MEM_SIZE,
			sizeof(cl_ulong), &local_mem_size, NULL);
  if (ret) abort();

  if ( 0 ) {
    printf("work_item_dim: %d ",work_item_dim);
    printf("work_item_sizes: %d %d %d ",
	   work_item_sizes[0],work_item_sizes[1],work_item_sizes[2]);
    printf("work_group_size: %d ",work_group_size);
    printf("compute_unit: %d\n",compute_unit);
    printf("local_mem_size: %d\n",local_mem_size);
  }
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

double gp_norm(int n, double *x)
{
  check_np(n);

  static cl_kernel kernel;
  cl_load(kernel,"gp_norm");

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
  return npa[0];
}


double gp_dot(int n, double *y, double *x)
{
  check_np(n);
  
  static cl_kernel kernel;
  cl_load(kernel,"gp_dot");

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
    
  return npa[0];
}

