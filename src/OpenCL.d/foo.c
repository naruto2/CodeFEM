#include <CL/cl.h>
#include <stdio.h>


int main()
{
  char buffer[2048];
  cl_platform_id platform_id = NULL;//プラットフォームID
  cl_device_id device_id = NULL;
  size_t size_ret;

  cl_uint num_platforms;
  cl_platform_id clSelectedPlatformID = NULL;
  clGetPlatformIDs (1, &platform_id, &num_platforms);//プラットフォームIDの取得
  printf("プラットフォーム数 : %d\n",(int) num_platforms);

  clGetPlatformInfo (platform_id, CL_PLATFORM_NAME, sizeof(buffer), buffer, NULL);
  printf("CL_PLATFORM_NAME:    %s\n",buffer);
  clGetPlatformInfo (platform_id, CL_PLATFORM_PROFILE, sizeof(buffer), buffer, NULL);
  printf("CL_PLATFORM_PROFILE: %s\n",buffer);
  clGetPlatformInfo (platform_id, CL_PLATFORM_VERSION, sizeof(buffer), buffer, NULL);
  printf("CL_PLATFORM_VERSION: %s\n",buffer);
  clGetPlatformInfo (platform_id, CL_PLATFORM_VENDOR, sizeof(buffer), buffer, NULL);
  printf("CL_PLATFORM_VENDOR:  %s\n",buffer);
  clGetPlatformInfo (platform_id, CL_PLATFORM_EXTENSIONS, sizeof(buffer), buffer, NULL);
  printf("CL_PLATFORM_EXTENSIONS: %s\n",buffer);


  cl_uint devnum;
  clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_ALL, 1, &device_id, &devnum);
  printf("deviceは %d\n",devnum);

  cl_uint ret_uint;

  // CPUのデバイスIDを取る
  clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_CPU, 1, &device_id, &devnum);
  printf("\n***** CPUは %dこ *****\n",devnum);

  clGetDeviceInfo( device_id, CL_DEVICE_NAME, sizeof(buffer), &buffer, &size_ret);
  printf("CL_DEVICE_NAME: %s\n",buffer);
  clGetDeviceInfo( device_id, CL_DEVICE_VENDOR, sizeof(buffer), &buffer, (size_t*)buffer);
  printf("CL_DEVICE_VENDOR: %s\n",buffer);
  clGetDeviceInfo( device_id, CL_DEVICE_OPENCL_C_VERSION, sizeof(buffer), &buffer, &size_ret);
  printf("CL_DEVICE_OPENCL_C_VERSION: %s\n",buffer);
  clGetDeviceInfo( device_id, CL_DEVICE_PROFILE , sizeof(buffer), &buffer, &size_ret);
  printf("CL_DEVICE_PROFILE: %s\n",buffer);
  clGetDeviceInfo( device_id, CL_DEVICE_VERSION, sizeof(buffer), &buffer, &size_ret);
  printf("CL_DEVICE_VERSION: %s\n",buffer);
  clGetDeviceInfo( device_id, CL_DRIVER_VERSION, sizeof(buffer), &buffer, &size_ret);
  printf("CL_DRIVER_VERSION: %s\n",buffer);
  clGetDeviceInfo( device_id, CL_DEVICE_EXTENSIONS , sizeof(buffer), &buffer, &size_ret);
  printf("CL_DEVICE_EXTENSIONS: %s\n",buffer);



  // GPUのデバイスIDを取る
  clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &devnum);
  printf("\n***** GPUは %dこ *****\n",devnum);

  clGetDeviceInfo( device_id, CL_DEVICE_NAME, sizeof(buffer), &buffer, &size_ret);
  printf("CL_DEVICE_NAME: %s\n",buffer);
  clGetDeviceInfo( device_id, CL_DEVICE_VENDOR, sizeof(buffer), &buffer, &size_ret);
  printf("CL_DEVICE_VENDOR: %s\n",buffer);
  clGetDeviceInfo( device_id, CL_DEVICE_OPENCL_C_VERSION, sizeof(buffer), &buffer, &size_ret);
  printf("CL_DEVICE_OPENCL_C_VERSION: %s\n",buffer);
  clGetDeviceInfo( device_id, CL_DEVICE_PROFILE , sizeof(buffer), &buffer, &size_ret);
  printf("CL_DEVICE_PROFILE: %s\n",buffer);
  clGetDeviceInfo( device_id, CL_DEVICE_VERSION, sizeof(buffer), &buffer, &size_ret);
  printf("CL_DEVICE_VERSION: %s\n",buffer);
  clGetDeviceInfo( device_id, CL_DRIVER_VERSION, sizeof(buffer), &buffer, &size_ret);
  printf("CL_DRIVER_VERSION: %s\n",buffer);
  clGetDeviceInfo( device_id, CL_DEVICE_EXTENSIONS , sizeof(buffer), &buffer, &size_ret);
  printf("CL_DEVICE_EXTENSIONS: %s\n",buffer);



  puts("***");
  clGetDeviceInfo( device_id, CL_DEVICE_ADDRESS_BITS, sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_ADDRESS_BITS : %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE : %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_MAX_CLOCK_FREQUENCY  , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_MAX_CLOCK_FREQUENCY: %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_MEM_BASE_ADDR_ALIGN , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_MEM_BASE_ADDR_ALIGN: %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE  , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE: %i\n",ret_uint);

  clGetDeviceInfo( device_id, CL_DEVICE_NATIVE_VECTOR_WIDTH_CHAR  , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_NATIVE_VECTOR_WIDTH_CHAR: %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_NATIVE_VECTOR_WIDTH_SHORT , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_NATIVE_VECTOR_WIDTH_SHORT: %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_NATIVE_VECTOR_WIDTH_INT , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_NATIVE_VECTOR_WIDTH_INT: %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_NATIVE_VECTOR_WIDTH_LONG: %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT: %i\n",ret_uint);
  // cl_khr_fp64がサポートされてないと0を返す
  clGetDeviceInfo( device_id, CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE: %i\n",ret_uint);
  // cl_khr_fp16がサポートされてないと0を返す
  clGetDeviceInfo( device_id, CL_DEVICE_NATIVE_VECTOR_WIDTH_HALF , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_NATIVE_VECTOR_WIDTH_HALF: %i\n",ret_uint);

  clGetDeviceInfo( device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR: %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT: %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT: %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG: %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT: %i\n",ret_uint);
  // cl_khr_fp64がサポートされてないと0を返す
  clGetDeviceInfo( device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE: %i\n",ret_uint);
  // cl_khr_fp16がサポートされてないと0を返す
  clGetDeviceInfo( device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF: %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_VENDOR_ID , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_VENDOR_ID: %i\n",ret_uint);
  clGetDeviceInfo( device_id, CL_DEVICE_MAX_COMPUTE_UNITS , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_MAX_COMPUTE_UNITS: %i\n",ret_uint);

  clGetDeviceInfo( device_id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS , sizeof(ret_uint), &ret_uint, &size_ret);
  printf("CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: %i\n",ret_uint);

}
