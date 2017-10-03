#include <stdlib.h>
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
#include <stdio.h>

#define NAME_NUM (8)  /* Number of stocks */
#define DATA_NUM (21) /* Number of data to process for each stock */

/* Read Stock data */
int stock_array_many[NAME_NUM*DATA_NUM];



/* Moving average width */
#define WINDOW_SIZE (13)

#define MAX_SOURCE_SIZE (0x100000)

static cl_context context = NULL;
static cl_command_queue command_queue = NULL;
static cl_uint work_dim = 1;
static size_t global_item_size[3];
static size_t local_item_size[3];
static cl_program program = NULL;


static float cl_norm(int nn, float *x_norm)
{
  static cl_mem mem_x_norm = NULL;
  static cl_kernel kernel_norm = NULL;
  int k, ret;

  if(!kernel_norm) {
    kernel_norm = clCreateKernel(program, "mynorm", &ret);
    printf("kernel_norm=%d\n",ret);
  }
  if(!mem_x_norm)
    mem_x_norm  = clCreateBuffer(context, CL_MEM_READ_WRITE,
				 nn * sizeof(float), NULL, &ret);

    ret = clEnqueueWriteBuffer(command_queue, mem_x_norm, CL_TRUE, 0,
			       nn*sizeof(float),
                               x_norm, 0, NULL, NULL);
    ret = clSetKernelArg(kernel_norm, 0, sizeof(int), (void *)&nn);
    ret = clSetKernelArg(kernel_norm, 1, sizeof(cl_mem), (void *)&mem_x_norm);
    ret = clEnqueueNDRangeKernel(command_queue, kernel_norm, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_x_norm, CL_TRUE, 0,
                              1 * sizeof(float),
                              x_norm,0, NULL, NULL);
    return x_norm[0];
}


static void cl_copy(int nn, float *y_copy, float *x_norm)
{
  static cl_mem mem_y_copy = NULL;
  static cl_mem mem_x_norm = NULL;
  static cl_kernel  kernel = NULL;
  int k, ret;

  if(!kernel) {
    kernel = clCreateKernel(program, "cl_copy", &ret);
    printf("cl_copy=%d\n",ret);
  }
  if(!mem_y_copy)
    mem_y_copy  = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
				 nn * sizeof(float), NULL, &ret);
  if(!mem_x_norm)
    mem_x_norm  = clCreateBuffer(context, CL_MEM_READ_ONLY,
				 nn * sizeof(float), NULL, &ret);

    ret = clEnqueueWriteBuffer(command_queue, mem_x_norm, CL_TRUE, 0,
			       nn*sizeof(float),
                               x_norm, 0, NULL, NULL);
    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&nn);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_y_copy);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_x_norm);

    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_y_copy, CL_TRUE, 0,
                              nn * sizeof(float),
                              y_copy,0, NULL, NULL);
}


static float cl_dot(int nn, float *y_copy, float *x_norm)
{
  static cl_mem mem_y_copy = NULL;
  static cl_mem mem_x_norm = NULL;
  static cl_kernel  kernel = NULL;
  int k, ret;

  if(!kernel) {
    kernel = clCreateKernel(program, "cl_dot", &ret);
    printf("cl_dot=%d\n",ret);
  }
  if(!mem_y_copy)
    mem_y_copy  = clCreateBuffer(context, CL_MEM_READ_WRITE,
				 nn * sizeof(float), NULL, &ret);
  if(!mem_x_norm)
    mem_x_norm  = clCreateBuffer(context, CL_MEM_READ_WRITE,
				 nn * sizeof(float), NULL, &ret);

    ret = clEnqueueWriteBuffer(command_queue, mem_y_copy, CL_TRUE, 0,
			       nn*sizeof(float),
                               y_copy, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, mem_x_norm, CL_TRUE, 0,
			       nn*sizeof(float),
                               x_norm, 0, NULL, NULL);
    ret = clSetKernelArg(kernel, 0, sizeof(int), (void *)&nn);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mem_y_copy);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&mem_x_norm);

    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, mem_x_norm, CL_TRUE, 0,
                              1 * sizeof(float),
                              x_norm,0, NULL, NULL);
    return x_norm[0];
}



int main(int argc, char **argv)
{
    cl_platform_id platform_id = NULL;
    cl_uint ret_num_platforms;
    cl_device_id device_id = NULL;
    cl_uint ret_num_devices;

    cl_mem memobjIn = NULL;
    cl_mem memobjIn1 = NULL;
    cl_mem memobjIn2 = NULL;
    cl_mem memobjOut = NULL;
    cl_mem memobjOut1 = NULL;

    cl_kernel kernel = NULL;
    size_t kernel_code_size;
    char *kernel_src_str;
    int *result;
    float *z, *x, *y;
    cl_int ret;
    FILE *fp;
    int k;

    int window_num = (int)WINDOW_SIZE;
    int point_num = NAME_NUM * DATA_NUM;
    int data_num = (int)DATA_NUM;
    int name_num = (int)NAME_NUM;

    int i,j;

    /* Allocate space to read in kernel code */
    kernel_src_str = (char *)malloc(MAX_SOURCE_SIZE);

    /* Allocate space for the result on the host side */
    result = (int *)malloc(point_num*sizeof(int));
    z = (float *)malloc(8192*sizeof(float));
    x = (float *)malloc(8192*sizeof(float));
    y = (float *)malloc(8192*sizeof(float));

    for(k=0;k<8192;k++) x[k] = 1.0;
    for(k=0;k<8192;k++) y[k] = 1.0;

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

    /* Create kernel */
    kernel = clCreateKernel(program, "moving_average_vec4_para", &ret);

    /* Create buffer for the input data on the device */
    memobjIn  = clCreateBuffer(context, CL_MEM_READ_WRITE,
                               point_num * sizeof(int), NULL, &ret);

    /* Create buffer for the result on the device */    
    memobjOut = clCreateBuffer(context, CL_MEM_READ_WRITE,
                               point_num * sizeof(float), NULL, &ret);

    /* Copy input data to the global memory on the device*/
    ret = clEnqueueWriteBuffer(command_queue, memobjIn, CL_TRUE, 0,
                               point_num * sizeof(int),
                               stock_array_many, 0, NULL, NULL);

    memobjIn1  = clCreateBuffer(context, CL_MEM_READ_WRITE,
				8192 * sizeof(float), NULL, &ret);
    memobjIn2  = clCreateBuffer(context, CL_MEM_READ_WRITE,
				8192 * sizeof(float), NULL, &ret);
    memobjOut1 = clCreateBuffer(context, CL_MEM_READ_WRITE,
				8192 * sizeof(float), NULL, &ret);

    ret = clEnqueueWriteBuffer(command_queue, memobjIn1, CL_TRUE, 0,
                               8192 * sizeof(float),
                               x, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(command_queue, memobjIn2, CL_TRUE, 0,
                               8192 * sizeof(float),
                               y, 0, NULL, NULL);

    
    /* Set Kernel Arguments */
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&memobjIn);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&memobjOut);
    ret = clSetKernelArg(kernel, 2, sizeof(int),    (void *)&data_num);
    ret = clSetKernelArg(kernel, 3, sizeof(int),    (void *)&name_num);
    ret = clSetKernelArg(kernel, 4, sizeof(int),    (void *)&window_num);
    ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&memobjIn1);
    ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&memobjIn2);
    ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&memobjOut1);

    /* Set parameters for data parallel processing (work item) */

    global_item_size[0] = 4; /* np Global number of work items */
    local_item_size[0] = 1;  /* Number of work items per work group */

    if ( argc > 1 ) global_item_size[0] = atoi(argv[1]);    
    
    /* --> global_item_size[0] / local_item_size[0] becomes 2, which indirectly sets the number of workgroups to 2*/

    /* Execute Data Parallel Kernel */ 
    ret = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL,
                                 global_item_size, local_item_size,
                                 0, NULL, NULL);

    /* Copy result from device to host */
    ret = clEnqueueReadBuffer(command_queue, memobjOut, CL_TRUE, 0,
                              point_num * sizeof(int),
                              result, 0, NULL, NULL);

    ret = clEnqueueReadBuffer(command_queue, memobjOut1, CL_TRUE, 0,
                              8192 * sizeof(float),
                              z, 0, NULL, NULL);
    ret = clReleaseKernel(kernel);
    ret = clReleaseMemObject(memobjIn);
    ret = clReleaseMemObject(memobjOut);
    ret = clReleaseMemObject(memobjOut1);

#ifndef mynorm // mynorm用 ################################################## 
    int nn = 1024;
    float *x_norm = (float*)calloc(sizeof(float),nn);
    x_norm[0] = 0.0;
    for(k=1;k<nn;k++) x_norm[k] = 1.0;
    printf("norm(n,x) = %f\n",cl_norm(nn,x_norm));

    float *y_copy = (float*)calloc(sizeof(float),nn);
    cl_copy(nn,y_copy,x_norm);
    printf("dot(n,y,x) = %f\n",cl_dot(nn,y_copy,x_norm));
#endif // #####################################################################
    /* OpenCL Object Finalization */

    ret = clReleaseProgram(program);
    ret = clReleaseCommandQueue(command_queue);
    ret = clReleaseContext(context);

    printf("np=%d\n",result[0]);    
    /* Deallocate memory on the host */
    free(result);
    free(kernel_src_str);
    return 0;
}

