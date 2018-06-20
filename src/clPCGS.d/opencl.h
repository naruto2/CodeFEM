#ifndef OPENCL_H
#define OPENCL_H

#include <stdio.h>
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

#endif
