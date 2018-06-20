#ifndef OPENCL_H
#define OPENCL_H
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <sys/time.h>
#include <CL/cl.h>

#define checkError(openclFunction)	        \
	if (cl_int err = (openclFunction))	    \
	{                                       \
		printf("error : %d\n", err);        \
	}

namespace OpenCL
{
  void init(void);
	void initialize(int platformIdx, int deviceIdx);
	void detectLine(unsigned char* result, unsigned char* origin);
	void release();
	cl_program compileProgram(char* fileName);
	cl_kernel createKernel(cl_program program, char* kernelName);

}

extern cl_uint ret_num_devices;
extern cl_platform_id platform_id;
extern cl_device_id device_id;
extern cl_context context;
extern cl_command_queue command_queue;
extern cl_program program;
extern cl_kernel IsKernelFunction(std::string function);
extern cl_mem *array_function(size_t size, void *a, int n);
extern cl_mem vec[32];


#define array(type, a, n) type* a = (type*)malloc(sizeof(type)*n)
#define scalar(n) sizeof(int),&n
#define vector(a,n) sizeof(cl_mem), array_function(sizeof(*a), a, n)
#define arg(kernel,k,vectorf) ret = clSetKernelArg(kernel,k,vectorf)

#endif
