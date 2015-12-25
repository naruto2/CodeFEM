#define __CL_ENABLE_EXCEPTIONS
#pragma warning(push)
#pragma warning(disable:4290)
#include <CL/cl.hpp>
#pragma warning(pop)
#include <iostream>
#include <cstdio>
//#include <conio.h>
#define nullptr  NULL

#pragma comment(lib, "opencl.lib")

const char* CLHelloSrcCodeStr = "__kernel void hello(void) { }";

void ExecuteMySimpleOpenCLProgram(cl::Context& context)
{
  cl_int err = CL_SUCCESS;
  const std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();

  cl::Program::Sources source(1, std::make_pair(CLHelloSrcCodeStr, strlen(CLHelloSrcCodeStr)));
  cl::Program program = cl::Program(context, source);
  program.build(devices);

  cl::Kernel kernel(program, "hello", &err);

  cl::Event event;
  cl::CommandQueue queue(context, devices[0], 0, &err);
  puts("Enqueue.");
  queue.enqueueNDRangeKernel(
			     kernel, cl::NullRange, cl::NDRange(4, 4), cl::NullRange, nullptr, &event);
  puts("Waiting...");
  event.wait();
  puts("Finished.");
}

int main()
{
  cl_int err = CL_SUCCESS;
    try
      {
	std::vector<cl::Platform> platforms;
	cl::Platform::get(&platforms);
	if (platforms.empty())
	  {
	    std::cerr << "NO OpenCL Platform." << std::endl;
	    return -1;
	  }
	std::cout << "OpenCL Platform Count = " << platforms.size() << std::endl;

	for (std::vector<cl::Platform>::iterator it = platforms.begin(); it != platforms.end(); ++it)
	  {
	    // 何番目のプラットフォームがどのベンダーの GPU / CPU になるかは不定。
	    // SDK やドライバーのインストール順序などにも左右される。
	    std::string platformName;
	    err = it->getInfo(CL_PLATFORM_NAME, &platformName);
	    std::cout << "PlatformName = " << platformName << std::endl;
	    // とりあえずすべての OpenCL デバイスを列挙させてみる。
	    std::vector<cl::Device> devs;
	          try
		    {
		      it->getDevices(CL_DEVICE_TYPE_ALL, &devs);
		    }
		  catch (const cl::Error&)
		    {
		      puts("NO OpenCL device in the platform!!");
		      continue;
		    }
		  for (std::vector<cl::Device>::iterator jt = devs.begin(); jt != devs.end(); ++jt)
		    {
		      std::string devName;
		      err = jt->getInfo(CL_DEVICE_NAME, &devName);
		      std::cout << "DeviceName = " << devName << std::endl;
		    }

		  // GPU デバイスの取得を試みて、取得できれば OpenCL プログラムを作成して実行させる。
		        try
			  {
			    it->getDevices(CL_DEVICE_TYPE_GPU, &devs); // デバイスの取得に失敗したら例外が発生する。
			  }
			catch (const cl::Error&)
			  {
			    puts("NO GPU device in the platform!!");
			    devs.clear();
			  }
			// 誤って余計な例外を捕捉しないよう、try-catch ブロックは分ける。
			if (!devs.empty())
			  {
			            cl_context_properties properties[] =
				      {
					CL_CONTEXT_PLATFORM, reinterpret_cast<cl_context_properties>((*it)()),
					          0
				      };
				    cl::Context context(CL_DEVICE_TYPE_GPU, properties);
				    ExecuteMySimpleOpenCLProgram(context);
			  }

			// CPU デバイスの取得を試みて、取得できれば OpenCL プログラムを作成して実行させる。
			      try
				{
				  it->getDevices(CL_DEVICE_TYPE_CPU, &devs); // デバイスの取得に失敗したら例外が発生する。
				}
			      catch (const cl::Error&)
				{
				  puts("NO CPU device in the platform!!");
				  devs.clear();
				}
			      // 誤って余計な例外を捕捉しないよう、try-catch ブロックは分ける。
			      if (!devs.empty())
				{
				          cl_context_properties properties[] =
					    {
					      CL_CONTEXT_PLATFORM, reinterpret_cast<cl_context_properties>((*it)()),
					                0
					    };
					  cl::Context context(CL_DEVICE_TYPE_CPU, properties);
					  ExecuteMySimpleOpenCLProgram(context);
				}
	  }
      }
    catch (const cl::Error& errObj)
      {
	std::cerr << "ERROR: " << errObj.what() << "(" << errObj.err() << ")" << std::endl;
      }

    puts("Press any...");
    //_getch();
    return 0;
}
