#include "CL.h"
#include <chrono>

int main(int argc, char** argv)
{
  std::size_t const n = 32;

  cl_platform_id platform_id = NULL;
  cl_device_id device_id = NULL;
  cl_context context = NULL;
  cl_command_queue command_queue = NULL;

  cl_kernel kernel = NULL;
  cl_uint ret_num_devices;
  cl_uint ret_num_platforms;
  cl_int ret;

  cl_mem x_dev, y_dev, z_dev, w_dev;//for device memory

  std::chrono::system_clock::time_point  start, end;

  double* x = (double*)malloc(sizeof(double) * n);
  double* y = (double*)malloc(sizeof(double) * n);
  double* z = (double*)malloc(sizeof(double) * n);
  double* w = (double*)malloc(sizeof(double) * n);

  for (int i = 0; i < n; i++){
    x[i] = 1.0;
    y[i] = 1.0;
  }

  //Get the platforms
  ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);

  //Get the device
  ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1,
		       &device_id, &ret_num_devices);

  //Create context
  context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);

  //Create Command Queue
  cl_command_queue_properties prop[8];
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
#pragma GCC diagnostic warning "-Wdeprecated-declarations"
  cl_program program = cc(device_id,context);



  //Create Kernel object
  kernel = clCreateKernel(program, "vec_add", &ret);
  if (ret) printf("error: clCreateKernel vec_add\n");

  //Create memory object
  x_dev = clCreateBuffer(context, CL_MEM_READ_WRITE,
			 n * sizeof(double), NULL, &ret);

  y_dev = clCreateBuffer(context, CL_MEM_READ_WRITE,
			 n * sizeof(double), NULL, &ret);

  z_dev = clCreateBuffer(context, CL_MEM_READ_WRITE,
			 n * sizeof(double), NULL, &ret);

  w_dev = clCreateBuffer(context, CL_MEM_READ_WRITE,
			 128 * sizeof(double), NULL, &ret);



  
  //memory copy host to device
  ret = clEnqueueWriteBuffer(command_queue, x_dev, CL_TRUE,
			     0, n * sizeof(double), x,
			     0, NULL, NULL);
  ret = clEnqueueWriteBuffer(command_queue, y_dev, CL_TRUE,
			     0, n * sizeof(double), y,
			     0, NULL, NULL);
  ret = clEnqueueWriteBuffer(command_queue, z_dev, CL_TRUE,
			     0, n * sizeof(double), z,
			     0, NULL, NULL);
  ret = clEnqueueWriteBuffer(command_queue, w_dev, CL_TRUE,
			     0, 128 * sizeof(double), w,
			     0, NULL, NULL);
  //Set args for kernel object
  ret = clSetKernelArg(kernel, 0, sizeof(int), &n);
  ret += clSetKernelArg(kernel, 1, sizeof(cl_mem), &w_dev);
  ret += clSetKernelArg(kernel, 2, sizeof(cl_mem), &z_dev);
  ret += clSetKernelArg(kernel, 3, sizeof(cl_mem), &x_dev);
  ret += clSetKernelArg(kernel, 4, sizeof(cl_mem), &y_dev);

  size_t global_work_size[3] = {128, 1, 1};
  size_t local_work_size[3] = {32, 1, 1};
  start = std::chrono::system_clock::now();
  // 計測したい処理

  //Execute Kernel Program
  ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL,
			       global_work_size, local_work_size,
			       0, NULL, NULL);

  //memory copy device to host
  ret = clEnqueueReadBuffer(command_queue, w_dev, CL_TRUE, 0,
			    128 * sizeof(double), w, 0,
			    NULL, NULL);
  end = std::chrono::system_clock::now();
  double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

  printf("%.3f\n",elapsed/1000);
  
  
  for(std::size_t i=0; i<n; ++i){
    //std::cout << x[i] << "*" << y[i] << " = " << w[i]<< std::endl;
  }
#if 0

  double* zd = (double*)malloc(sizeof(double) * n);
  double* xd = (double*)malloc(sizeof(double) * n);
  double* yd = (double*)malloc(sizeof(double) * n);

  for (int i = 0; i < n; i++){
    zd[i] = 0.123;
    xd[i] = 1.0;
    yd[i] = (double)i;
  }
  cl_kernel kernel_vec_dot = clCreateKernel(program, "vec_dot", &ret);
  if (ret) printf("error: clCreateKernel vec_dot\n");
  
  //Create memory object
  z_dev = clCreateBuffer(context, CL_MEM_READ_WRITE,
			 n * sizeof(double), NULL, &ret);
  if (ret) printf("error: clCreateBuffer\n");
  x_dev = clCreateBuffer(context, CL_MEM_READ_WRITE,
			 n * sizeof(double), NULL, &ret);
  if (ret) printf("error: clCreateBuffer\n");
  y_dev = clCreateBuffer(context, CL_MEM_READ_WRITE,
			 n * sizeof(double), NULL, &ret);
  if (ret) printf("error: clCreateBuffer\n");

  //memory copy host to device
  ret = clEnqueueWriteBuffer(command_queue, z_dev, CL_TRUE,
			     0, n * sizeof(double), z,
			     0, NULL, NULL);
  if (ret) printf("error: clEnqueueWriteBuffer\n");
  ret = clEnqueueWriteBuffer(command_queue, x_dev, CL_TRUE,
			     0, n * sizeof(double), x,
			     0, NULL, NULL);
  if (ret) printf("error: clEnqueueWriteBuffer\n");
  ret = clEnqueueWriteBuffer(command_queue, y_dev, CL_TRUE,
			     0, n * sizeof(double), y,
			     0, NULL, NULL);
  if (ret) printf("error: clEnqueueWriteBuffer\n");
  //Set args for kernel object
  ret  = clSetKernelArg(kernel_vec_dot, 0, sizeof(int), &n);
  if (ret) {
    printf("error: clSetKernelArg\n");
    if ( ret == CL_INVALID_KERNEL ) printf(" CL_INVALID_KERNEL\n");
    if ( ret == CL_INVALID_ARG_INDEX )      printf(" CL_INVALID_ARG_INDEX\n");
    if ( ret == CL_INVALID_ARG_VALUE )      printf(" CL_ARG_VALUE\n");
    if ( ret == CL_INVALID_MEM_OBJECT )     printf(" CL_MEM_OBJECT\n");
    if ( ret == CL_INVALID_SAMPLER )        printf(" CL_SAMPLER\n");
    if ( ret == CL_INVALID_ARG_SIZE )       printf(" CL_ARG_SIZE\n");

  }
  ret = clSetKernelArg(kernel_vec_dot, 1, sizeof(cl_mem), &z_dev);
  if (ret) printf("error: clSetKernelArg\n");
  ret = clSetKernelArg(kernel_vec_dot, 2, sizeof(cl_mem), &x_dev);
  if (ret) printf("error: clSetKernelArg\n");
  ret = clSetKernelArg(kernel_vec_dot, 3, sizeof(cl_mem), &y_dev);
  if (ret) printf("error: clSetKernelArg\n");

  // size_t global_work_size[3] = {n, 0, 0};
  // size_t local_work_size[3] = {n, 0, 0};

  //Execute Kernel Program
  ret = clEnqueueNDRangeKernel(command_queue, kernel_vec_dot, 1, NULL,
			       global_work_size, local_work_size,
			       0, NULL, NULL);
  if ( ret ) printf("error: clEnqueueNDRangeKernel\n");
  
  //memory copy device to host
  ret = clEnqueueReadBuffer(command_queue, z_dev, CL_TRUE, 0,
			    n * sizeof(double), z, 0,
			    NULL, NULL);

  if ( ret ) printf("error: clEnqueueReadBuffer\n");

  for(std::size_t i=0; i<n; ++i){
    //std::cout << xd[i] << "*" << yd[i] << " = " << zd[i]<< std::endl;
  }
#endif

  ret = clFlush(command_queue);
  if ( ret ) printf("error: clFlush\n");
  ret = clFinish(command_queue);
  if ( ret ) printf("error: clFinish\n");

  
  //memory free
  ret = clReleaseKernel(kernel);
  ret = clReleaseProgram(program);
  ret = clReleaseCommandQueue(command_queue);
  ret = clReleaseContext(context);

  ret = clReleaseMemObject(x_dev);
  ret = clReleaseMemObject(y_dev);
  ret = clReleaseMemObject(z_dev);
  free(x);
  free(y);
  free(z);

  return 0;
}
