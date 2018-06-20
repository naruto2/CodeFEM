#ifndef _CL_H_
#define _CL_H_
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <sys/time.h>
#include <CL/cl.h>
#define MAX_SOURCE_SIZE (0x100000)


cl_program cc(cl_device_id device_id, cl_context context)
{
  cl_program program = NULL;
  const char filename[] = "./kernel.cl";

  //Read a Kernel code
  FILE* fp = fopen(filename, "r");
  if(!fp){
    std::cerr << "Failed to load kernel" << std::endl;
    exit(0);
  }
  char* source_str = (char*)malloc(MAX_SOURCE_SIZE);
  std::size_t source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
  fclose(fp);

  int ret;
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
  return program;
}

#endif
