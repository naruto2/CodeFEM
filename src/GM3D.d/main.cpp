#include <iostream>
#include <CL/opencl.h>
#include "myCLUtils.h"
#include "myCLManager2.h"

#define N (33 * 1024)
#define THREADSPERBLOCK 256
#define imin(a, b) (a < b ? a : b)
#define expectedSum(x) (x * (x + 1)*(2 * x + 1) / 3)

int main(int argc, char **argv)
{
  cl_int ret;
  int i, n = N;
  const int queueNum = 1, debug = 0;
  const int threadsPerBlock = THREADSPERBLOCK;
  const int blocksPerGrid 
    = imin(32, (N + THREADSPERBLOCK - 1) / THREADSPERBLOCK);
  const size_t globalWorkSize[1] = {N};
  const size_t localWorkSize[1] = {threadsPerBlock};

  float *a, *b, *partial_c;
  cl_mem dev_a, dev_b, dev_partial_c;
  myCLManager2 clm(queueNum, debug);
  clm.register_program("dotSource", "dotProduct.cl");
  cl_program program = clm.get_program("dotSource");
  cl_kernel dotKernel = clCreateKernel(program, "dotProduct", &ret);
  myCLUtils::report_error(__FUNCTION__, "clCreateKernel", ret);
  a = new float[N];
  b = new float[N];
  partial_c = new float[blocksPerGrid];
  for(i = 0; i < N; i++){
    a[i] = i;
    b[i] = i * 2;
  }
  dev_a = clCreateBuffer(clm.context,
             CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
             N * sizeof(float), a, &ret);
  myCLUtils::report_error(__FUNCTION__, "clCreateBuffer", ret);
  dev_b = clCreateBuffer(clm.context,
             CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
             N * sizeof(float), b, &ret);
  myCLUtils::report_error(__FUNCTION__, "clCreateBuffer", ret);
  dev_partial_c = clCreateBuffer(clm.context,
                 CL_MEM_WRITE_ONLY,
                 blocksPerGrid * sizeof(float), NULL, &ret);
  myCLUtils::report_error(__FUNCTION__, "clCreateBuffer", ret);
  i = 0;
  ret = clSetKernelArg(dotKernel, i++, sizeof(dev_a), &dev_a);
  myCLUtils::report_error(__FUNCTION__, "clSetKernelArg", ret);
  ret = clSetKernelArg(dotKernel, i++, sizeof(dev_b), &dev_b);
  myCLUtils::report_error(__FUNCTION__, "clSetKernelArg", ret);
  ret = clSetKernelArg(dotKernel, i++, sizeof(dev_partial_c), &dev_partial_c);
  myCLUtils::report_error(__FUNCTION__, "clSetKernelArg", ret);
  ret = clSetKernelArg(dotKernel, i++, sizeof(n), &n);
  myCLUtils::report_error(__FUNCTION__, "clSetKernelArg", ret);
  ret = clSetKernelArg(dotKernel, i++, 
               sizeof(threadsPerBlock), &threadsPerBlock);
  myCLUtils::report_error(__FUNCTION__, "clSetKernelArg", ret);
  ret = clSetKernelArg(dotKernel, i++, 
               sizeof(blocksPerGrid), &blocksPerGrid);
  myCLUtils::report_error(__FUNCTION__, "clSetKernelArg", ret);
  ret = clEnqueueNDRangeKernel(clm.queues[0],
                   dotKernel,
                   1,
                   NULL,
                   globalWorkSize,
                   localWorkSize,
                   0, NULL, NULL);
  myCLUtils::report_error(__FUNCTION__, "clEnqueueNDRangeKernel", ret);
  ret = clEnqueueReadBuffer(clm.queues[0],
                dev_partial_c,
                CL_TRUE,
                0,
                blocksPerGrid * sizeof(float),
                partial_c,
                0, NULL, NULL);
  myCLUtils::report_error(__FUNCTION__, "clEnqueueReadBuffer", ret);
  float c = 0;
  for(int i = 0; i < blocksPerGrid; i++){
    c += partial_c[i];
  }
  delete a;
  delete b;
  delete partial_c;
  ret = clReleaseMemObject(dev_a);
  myCLUtils::report_error(__FUNCTION__, "clReleaseMemObject", ret);
  clReleaseMemObject(dev_b);
  myCLUtils::report_error(__FUNCTION__, "clReleaseMemObject", ret);
  clReleaseMemObject(dev_partial_c);
  myCLUtils::report_error(__FUNCTION__, "clReleaseMemObject", ret);
  cout << "Blocks per Grid = " << blocksPerGrid << "\n";
  cout << "Threads per Block = " << threadsPerBlock << "\n";
  cout << "GPU sum  = " << c << "\n";
  cout << "Expected = " << expectedSum((float)(N - 1)) << endl;
  return 0;
}
