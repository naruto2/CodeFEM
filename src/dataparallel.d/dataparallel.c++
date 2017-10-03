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
static cl_program program = NULL;
static char *kernel_src_str;

double cl_norm(int n, double *x)
{
  static cl_kernel kernel = NULL;
  static cl_mem mem_x = NULL;
  static cl_mem mem_npa = NULL;
  static double *npa;
  
  int k, ret;
  if(!npa)
    npa = (double*)malloc(sizeof(double)*global_item_size[0]);

  if(!kernel) {
    kernel = clCreateKernel(program, "cl_norm", &ret);
    printf("cl_norm=%d\n",ret);
  }
  if(!mem_x)
    mem_x  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			    n * sizeof(double), NULL, &ret);
  if(!mem_npa)
    mem_npa  = clCreateBuffer(context, CL_MEM_READ_WRITE,
			      global_item_size[0] * sizeof(double),
			      NULL, &ret);

  
    ret = clEnqueueWriteBuffer(command_queue, mem_x, CL_TRUE, 0,
			       n*sizeof(double),
                               x, 0, NULL, NULL);

    ret = clEnqueueWriteBuffer(command_queue, mem_npa, CL_TRUE, 0,
			       global_item_size[0]*sizeof(double),
                               npa, 0, NULL, NULL);

    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_x);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_npa);

    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);

    ret = clEnqueueReadBuffer(command_queue, mem_npa, CL_TRUE, 0,
                              global_item_size[0] * sizeof(double),
                              npa,0, NULL, NULL);

    for(k=1;k<global_item_size[0];k++) npa[0] += npa[k];
    npa[0] = sqrt(npa[0]);

    return npa[0];
}


void cl_copy(int n, double *y, double *x)
{
  static cl_mem mem_y = NULL;
  static cl_mem mem_x = NULL;
  static cl_kernel  kernel = NULL;
  int k, ret;

  if(!kernel) {
    kernel = clCreateKernel(program, "cl_copy", &ret);
    printf("cl_copy=%d\n",ret);
  }
  if(!mem_y)
    mem_y  = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
				 n * sizeof(double), NULL, &ret);
  if(!mem_x)
    mem_x  = clCreateBuffer(context, CL_MEM_READ_ONLY,
				 n * sizeof(double), NULL, &ret);

    ret = clEnqueueWriteBuffer(command_queue, mem_x, CL_TRUE, 0,
			       n*sizeof(double),
                               x, 0, NULL, NULL);
    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_y);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_x);

    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_y, CL_TRUE, 0,
                              n * sizeof(double),
                              y,0, NULL, NULL);
}


double cl_dot(int n, double *y, double *x)
{
  static cl_mem mem_y = NULL;
  static cl_mem mem_x = NULL;
  static cl_mem mem_npa = NULL;
  static cl_kernel  kernel = NULL;
  static double *npa = NULL;
  int k, ret;

  if(!npa)
    npa = (double*)malloc(sizeof(double)*global_item_size[0]);

  if(!kernel) {
    kernel = clCreateKernel(program, "cl_dot", &ret);
    printf("cl_dot=%d\n",ret);
  }
  if(!mem_y)
    mem_y  = clCreateBuffer(context, CL_MEM_READ_ONLY,
				 n * sizeof(double), NULL, &ret);
  if(!mem_x)
    mem_x  = clCreateBuffer(context, CL_MEM_READ_ONLY,
				 n * sizeof(double), NULL, &ret);
  if(!mem_npa)
    mem_npa  = clCreateBuffer(context, CL_MEM_READ_WRITE,
			      global_item_size[0] * sizeof(double),
			      NULL, &ret);

  
    ret = clEnqueueWriteBuffer(command_queue, mem_y, CL_TRUE, 0,
			       n*sizeof(double),
                               y, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_x, CL_TRUE, 0,
			       n*sizeof(double),
                               x, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_npa, CL_TRUE, 0,
			       global_item_size[0]*sizeof(double),
                               npa, 0, NULL, NULL);


    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_y);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_x);
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_npa);

    
    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_npa, CL_TRUE, 0,
                              global_item_size[0] * sizeof(double),
                              npa,0, NULL, NULL);

    for(k=1;k<global_item_size[0];k++) npa[0] += npa[k];
    
    return npa[0];
}


double cl_phase0(int n, double *r, double *Aa, int *col_ind,
	       int *row_ptr, double *x, double *rtilde,
	       double *b, int w)
{
  static int ww = 0;
  static double *npa = NULL;
  static cl_mem mem_r = NULL;
  static cl_mem mem_Aa= NULL;
  static cl_mem mem_col_ind = NULL;
  static cl_mem mem_row_ptr = NULL;
  static cl_mem mem_x = NULL;
  static cl_mem mem_rtilde = NULL;
  static cl_mem mem_b = NULL;
  static cl_mem mem_npa = NULL;

  static cl_kernel kernel = NULL;
  int k, ret;

  if(!npa)
    npa = (double*)malloc(sizeof(double)*global_item_size[0]);

  if(!kernel) {
    kernel = clCreateKernel(program, "cl_phase0", &ret);
    printf("cl_phase0=%d\n",ret);
 }
  if(!mem_r)
    mem_r  = clCreateBuffer(context, CL_MEM_READ_WRITE,
			    n * sizeof(double), NULL, &ret);

  if ( ww < w){
    if(mem_Aa) ret = clReleaseMemObject(mem_Aa);
    mem_Aa  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			     w * sizeof(double), NULL, &ret);
    if(mem_col_ind) ret = clReleaseMemObject(mem_col_ind);
    mem_col_ind  = clCreateBuffer(context, CL_MEM_READ_ONLY,
				  w * sizeof(int), NULL, &ret);
    ww = w;
  }

    if(!mem_row_ptr)
    mem_row_ptr  = clCreateBuffer(context, CL_MEM_READ_ONLY,
				  (n+1) * sizeof(double), NULL, &ret);
  if(!mem_x)
    mem_x  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			    n * sizeof(double), NULL, &ret);
  if(!mem_rtilde)
    mem_rtilde  = clCreateBuffer(context, CL_MEM_READ_WRITE,
			    n * sizeof(double), NULL, &ret);
  if(!mem_b)
    mem_b  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			    n * sizeof(double), NULL, &ret);
  if(!mem_npa)
    mem_npa  = clCreateBuffer(context, CL_MEM_READ_WRITE,
			      global_item_size[0] * sizeof(double),
			      NULL, &ret);
  
    ret = clEnqueueWriteBuffer(command_queue, mem_r, CL_TRUE, 0,
			       n*sizeof(double),
                               r, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_Aa, CL_TRUE, 0,
			       w*sizeof(double),
                               Aa, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_col_ind, CL_TRUE, 0,
			       w*sizeof(int),
                               col_ind, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_row_ptr, CL_TRUE, 0,
			       (n+1)*sizeof(int),
                               row_ptr, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_x, CL_TRUE, 0,
			       n*sizeof(double),
                               x, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_rtilde, CL_TRUE, 0,
			       n*sizeof(double),
                               rtilde, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_b, CL_TRUE, 0,
			       n*sizeof(double),
                               b, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_npa, CL_TRUE, 0,
			       global_item_size[0]*sizeof(double),
			       npa, 0, NULL, NULL);

    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_r);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_Aa);
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_col_ind);
    ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&mem_row_ptr);
    ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&mem_x);
    ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&mem_rtilde);
    ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&mem_b);
    ret = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&mem_npa);

    
    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);

    ret = clEnqueueReadBuffer(command_queue, mem_r, CL_TRUE, 0,
                              n * sizeof(double),
                              r,0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_rtilde, CL_TRUE, 0,
                              n * sizeof(double),
                              rtilde,0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_npa, CL_TRUE, 0,
			      global_item_size[0] * sizeof(double),
			      npa,0, NULL, NULL);

    for(k=1;k<global_item_size[0];k++) npa[0] += npa[k];

    return sqrt(npa[0]);
}


void cl_phase1(int n, double *p, double *r, double *v,
	       double beta, double omega)
{
  static cl_mem mem_p = NULL;
  static cl_mem mem_r = NULL;
  static cl_mem mem_v = NULL;
  static cl_kernel kernel = NULL;
  int k, ret;


  if(!kernel) {
    kernel = clCreateKernel(program, "cl_phase1", &ret);
    printf("cl_phase1=%d\n",ret);
  }
  if(!mem_p)
    mem_p  = clCreateBuffer(context, CL_MEM_READ_WRITE,
			    n * sizeof(double), NULL, &ret);
  if(!mem_r)
    mem_r  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			    n * sizeof(double), NULL, &ret);
  if(!mem_v)
    mem_v  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			    n * sizeof(double), NULL, &ret);

    ret = clEnqueueWriteBuffer(command_queue, mem_p, CL_TRUE, 0,
			       n*sizeof(double),
                               p, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_r, CL_TRUE, 0,
			       n*sizeof(double),
                               r, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_v, CL_TRUE, 0,
			       n*sizeof(double),
                               v, 0, NULL, NULL);


    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_p);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_r);
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&mem_v);
    ret = clSetKernelArg(kernel, 4, sizeof(double), (void *)&beta);
    ret = clSetKernelArg(kernel, 5, sizeof(double), (void *)&omega);
    
    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);

    ret = clEnqueueReadBuffer(command_queue, mem_p, CL_TRUE, 0,
                              n * sizeof(double),
                              p,0, NULL, NULL);
}


void cl_phase3(int n, double *s, double *r, double *v,
	       double alpha)
{
  static cl_mem mem_s = NULL;
  static cl_mem mem_r = NULL;
  static cl_mem mem_v = NULL;
  static cl_kernel kernel = NULL;
  int ret;


  if(!kernel) {
    kernel = clCreateKernel(program, "cl_phase3", &ret);
    printf("cl_phase3=%d\n",ret);
  }
  if(!mem_s)
    mem_s  = clCreateBuffer(context, CL_MEM_READ_WRITE,
			    n * sizeof(double), NULL, &ret);
  if(!mem_r)
    mem_r  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			    n * sizeof(double), NULL, &ret);
  if(!mem_v)
    mem_v  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			    n * sizeof(double), NULL, &ret);

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
    
    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);

    ret = clEnqueueReadBuffer(command_queue, mem_s, CL_TRUE, 0,
                              n * sizeof(double),
                              s,0, NULL, NULL);
}


void cl_phase4(int n, double *s, double *phat, double alpha)
{
  static cl_mem mem_s = NULL;
  static cl_mem mem_phat = NULL;
  static cl_kernel kernel = NULL;
  int ret;


  if(!kernel) {
    kernel = clCreateKernel(program, "cl_phase4", &ret);
    printf("cl_phase4=%d\n",ret);
  }
  if(!mem_s)
    mem_s  = clCreateBuffer(context, CL_MEM_READ_WRITE,
			    n * sizeof(double), NULL, &ret);
  if(!mem_phat)
    mem_phat  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			    n * sizeof(double), NULL, &ret);

  
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
    
    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);

    ret = clEnqueueReadBuffer(command_queue, mem_s, CL_TRUE, 0,
                              n * sizeof(double),
                              s,0, NULL, NULL);
}


void cl_phase6(int n, double *x, double *s, double *r, double *t, double *phat,
	       double *shat, double alpha, double omega)
{
  static cl_mem mem_x = NULL;
  static cl_mem mem_s = NULL;
  static cl_mem mem_r = NULL;
  static cl_mem mem_t = NULL;
  static cl_mem mem_phat = NULL;
  static cl_mem mem_shat = NULL;
  static cl_kernel kernel = NULL;
  int ret;


  if(!kernel) {
    kernel = clCreateKernel(program, "cl_phase6", &ret);
    printf("cl_phase6=%d\n",ret);
  }
  if(!mem_x)
    mem_x  = clCreateBuffer(context, CL_MEM_READ_WRITE,
			    n * sizeof(double), NULL, &ret);
  if(!mem_s)
    mem_s  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			    n * sizeof(double), NULL, &ret);
  if(!mem_r)
    mem_r  = clCreateBuffer(context, CL_MEM_READ_WRITE,
			    n * sizeof(double), NULL, &ret);
  if(!mem_t)
    mem_t  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			    n * sizeof(double), NULL, &ret);
  if(!mem_phat)
    mem_phat  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			       n * sizeof(double), NULL, &ret);
  if(!mem_shat)
    mem_shat  = clCreateBuffer(context, CL_MEM_READ_ONLY,
			       n * sizeof(double), NULL, &ret);

    ret = clEnqueueWriteBuffer(command_queue, mem_x, CL_TRUE, 0,
			       n*sizeof(double),
                               x, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_s, CL_TRUE, 0,
			       n*sizeof(double),
                               s, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_r, CL_TRUE, 0,
			       n*sizeof(double),
                               r, 0, NULL, NULL);
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
    
    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);

    ret = clEnqueueReadBuffer(command_queue, mem_x, CL_TRUE, 0,
                              n * sizeof(double),
                              x,0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_r, CL_TRUE, 0,
                              n * sizeof(double),
                              r,0, NULL, NULL);
}


void cl_init(int argc, char **argv)
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
  fp = fopen("moving_average_vec4_para.cl", "r");
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

