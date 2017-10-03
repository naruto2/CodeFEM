#include <stdlib.h>
#include <CL/cl.h>
#include <stdio.h>

#define MAX_SOURCE_SIZE (0x100000)

static cl_context context = NULL;
static cl_command_queue command_queue = NULL;
static cl_uint work_dim = 1;
static size_t global_item_size[3];
static size_t local_item_size[3];
static cl_program program = NULL;
static char *kernel_src_str;

static float cl_norm(int n, float *x)
{
  static cl_mem mem_x = NULL;
  static cl_kernel kernel_norm = NULL;
  int k, ret;

  if(!kernel_norm) {
    kernel_norm = clCreateKernel(program, "cl_norm", &ret);
    printf("cl_norm=%d\n",ret);
  }
  if(!mem_x)
    mem_x  = clCreateBuffer(context, CL_MEM_READ_WRITE,
				 n * sizeof(float), NULL, &ret);

    ret = clEnqueueWriteBuffer(command_queue, mem_x, CL_TRUE, 0,
			       n*sizeof(float),
                               x, 0, NULL, NULL);
    ret = clSetKernelArg(kernel_norm, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel_norm, 1, sizeof(cl_mem), (void *)&mem_x);
    ret = clEnqueueNDRangeKernel(command_queue, kernel_norm, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_x, CL_TRUE, 0,
                              1 * sizeof(float),
                              x,0, NULL, NULL);
    return x[0];
}


static void cl_copy(int n, float *y, float *x)
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
				 n * sizeof(float), NULL, &ret);
  if(!mem_x)
    mem_x  = clCreateBuffer(context, CL_MEM_READ_ONLY,
				 n * sizeof(float), NULL, &ret);

    ret = clEnqueueWriteBuffer(command_queue, mem_x, CL_TRUE, 0,
			       n*sizeof(float),
                               x, 0, NULL, NULL);
    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_y);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_x);

    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_y, CL_TRUE, 0,
                              n * sizeof(float),
                              y,0, NULL, NULL);
}


static float cl_dot(int n, float *y, float *x)
{
  static cl_mem mem_y = NULL;
  static cl_mem mem_x = NULL;
  static cl_kernel  kernel = NULL;
  int k, ret;

  if(!kernel) {
    kernel = clCreateKernel(program, "cl_dot", &ret);
    printf("cl_dot=%d\n",ret);
  }
  if(!mem_y)
    mem_y  = clCreateBuffer(context, CL_MEM_READ_WRITE,
				 n * sizeof(float), NULL, &ret);
  if(!mem_x)
    mem_x  = clCreateBuffer(context, CL_MEM_READ_WRITE,
				 n * sizeof(float), NULL, &ret);

    ret = clEnqueueWriteBuffer(command_queue, mem_y, CL_TRUE, 0,
			       n*sizeof(float),
                               y, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_x, CL_TRUE, 0,
			       n*sizeof(float),
                               x, 0, NULL, NULL);
    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&n);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_y);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_x);

    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_x, CL_TRUE, 0,
                              1 * sizeof(float),
                              x,0, NULL, NULL);
    return x[0];
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

  global_item_size[0] = 4; /* np Global number of work items */
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

  /* Create Command Queue */  
  command_queue = clCreateCommandQueue(context, device_id, 0, &ret);

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


int main(int argc, char **argv)
{
  cl_init(argc,argv);

  int k, n = 65537;
  float *x = (float*)malloc(sizeof(float)*n);
  float *y = (float*)malloc(sizeof(float)*n);

  x[0] = 0.0; for(k=1;k<n;k++) x[k] = 1.0;

  printf("norm(n,x) = %f\n",cl_norm(n,x));

  cl_copy(n,y,x);
  printf("dot(n,y,x) = %f\n",cl_dot(n,y,x));

  cl_finalize();
  return 0;
}
