#include <iostream>
#include <stdio.h>
#include <malloc.h>
#include <CL/cl.h>
#include "opencl.h"


cl_uint ret_num_platforms;
cl_uint ret_num_devices;
cl_platform_id platform_id;
cl_device_id device_id;
cl_context context;
cl_command_queue command_queue;
cl_program program;


void OpenCL::init(void)
{
  //Read a Kernel code
  FILE* fp = fopen("kernel.cl", "r");
  if(!fp){
    std::cerr << "Failed to load kernel" << std::endl;
    exit(-1);
  }
#define MAX_SOURCE_SIZE (0x100000)
  char* source_str = (char*)malloc(MAX_SOURCE_SIZE);
  std::size_t source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
  fclose(fp);

  //Get the platforms
  cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);

  //Get the device
  ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1,
		       &device_id, &ret_num_devices);

  //Create context
  context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);

  //Create Command Queue
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
#pragma GCC diagnostic warning "-Wdeprecated-declarations"

  //Create Program Object
  program = clCreateProgramWithSource(context, 1, (const char**)&source_str,
				      (const size_t*)&source_size, &ret);

  //Compile of kernel's source code
  ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
  if ( ret != CL_SUCCESS ) {
    if ( ret == CL_INVALID_PROGRAM )
      printf("clBuildProgram error: CL_INVALID_PROGRAM\n");
    if ( ret == CL_INVALID_VALUE )
      printf("clBuildProgram error: CL_INVALID_VALUE\n");
    if ( ret == CL_INVALID_DEVICE )
      printf("clBuildProgram error: CL_INVALID_DEVICE\n");
    if ( ret == CL_INVALID_BINARY )
      printf("clBuildProgram error: CL_INVALID_BINARY\n");
    if ( ret == CL_INVALID_BUILD_OPTIONS )
      printf("clBuildProgram error: CL_INVALID_BUILD_OPTIONS\n");
    if ( ret == CL_INVALID_OPERATION )
      printf("clBuildProgram error: CL_INVALID_OPERATION\n");
    if ( ret == CL_COMPILER_NOT_AVAILABLE )
      printf("clBuildProgram error: CL_COMPILER_NOT_AVAILABLE\n");
    if ( ret == CL_BUILD_PROGRAM_FAILURE )
      printf("clBuildProgram error: CL_BUILD_PROGRAM_FAILURE2\n");
    if ( ret == CL_OUT_OF_HOST_MEMORY )
      printf("clBuildProgram error: CL_OUT_OF_HOST_MEMORY\n");

    char errorlog[1000];
    size_t errorlogsize;
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG,
			  1000,errorlog,&errorlogsize);
    puts(errorlog);
    exit(0);
  }
  free(source_str);
}


cl_kernel IsKernelFunction(std::string function)
{
  cl_int ret;
  return clCreateKernel(program, function.c_str(), &ret);
}


cl_mem vec[32];


cl_mem *array_function(size_t size, void *a, int n)
{
  static int k;
  cl_int ret;

  vec[k] = clCreateBuffer(context, CL_MEM_READ_WRITE,
			  n * size, NULL, &ret);
  ret = clEnqueueWriteBuffer(command_queue, vec[k], CL_TRUE,
			     0, n * size, a,
			     0, NULL, NULL);
  k++;
  return &vec[k-1];
}
