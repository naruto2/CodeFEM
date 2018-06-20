#include "opencl.h"


int main(int argc, char** argv){

  OpenCL::init();

  int n = 170;
  cl_int ret;

  array(float, x, n); array(float, y, n); array(float, z, n);
  
  for(int i=0;i<n;i++) { x[i] = 1.0;  y[i] = i; }
  
  cl_kernel vec_add = IsKernelFunction("vec_add");
  arg(vec_add,0,scalar(n));
  arg(vec_add,1,vector(z,n));
  arg(vec_add,2,vector(x,n));
  arg(vec_add,3,vector(y,n));
  
  size_t global_work_size[3] = {(size_t)n, (size_t)1, 1};
  size_t local_work_size[3]  = {(size_t)n, (size_t)1, 1};

  ret = clEnqueueNDRangeKernel(command_queue, vec_add, 1, NULL,
			       global_work_size, local_work_size,
			       0, NULL, NULL);

  ret = clEnqueueReadBuffer(command_queue, vec[0], CL_TRUE, 0,
			    n * sizeof(float), z, 0,
			    NULL, NULL);
  
  for(std::size_t i=0; i<n; ++i){
    std::cout << x[i] << "+" << y[i] << " = " << z[i]<< std::endl;
  }

  ret = clFlush(command_queue);
  ret = clFinish(command_queue);
  ret = clReleaseKernel(vec_add);
  ret = clReleaseProgram(program);
  ret = clReleaseCommandQueue(command_queue);
  ret = clReleaseMemObject(vec[0]);
  ret = clReleaseMemObject(vec[1]);
  ret = clReleaseMemObject(vec[2]);

  free(x);
  free(y);
  free(z);

  return 0;
}
