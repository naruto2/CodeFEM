/**--------------------
   deviceQuery.cpp
   --------------------*/

#include <iostream>
#include <cstring>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

//for platform information
struct platform_info{
  char info_name[1024];
  cl_platform_info info_type;
};

//for device information
struct device_info{
  char info_name[1024];
  cl_device_info info_type;
};

struct platform_info pinfo_list[] = {
  {"CL_PLATFORM_PROFILE", CL_PLATFORM_PROFILE},
  {"CL_PLATFORM_VERSION", CL_PLATFORM_VERSION},
  {"CL_PLATFORM_NAME", CL_PLATFORM_NAME},
  {"CL_PLATFORM_VENDOR", CL_PLATFORM_VENDOR},
  {"CL_PLATFORM_EXTENSIONS", CL_PLATFORM_EXTENSIONS}
};

struct device_info dinfo_list[] = {
  {"CL_DEVICE_TYPE", CL_DEVICE_TYPE},
  {"CL_DEVICE_VENDOR_ID", CL_DEVICE_VENDOR_ID},
  {"CL_DEVICE_MAX_COMPUTE_UNITS", CL_DEVICE_MAX_COMPUTE_UNITS},
  {"CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS", CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS},
  {"CL_DEVICE_MAX_WORK_ITEM_SIZES", CL_DEVICE_MAX_WORK_ITEM_SIZES},
  {"CL_DEVICE_MAX_WORK_GROUP_SIZE", CL_DEVICE_MAX_WORK_GROUP_SIZE},
  {"CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR", CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR},
  {"CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT", CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT},
  {"CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT", CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT},
  {"CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG", CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG},
  {"CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT", CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT},
  {"CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE", CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE},
  {"CL_DEVICE_MAX_CLOCK_FREQUENCY", CL_DEVICE_MAX_CLOCK_FREQUENCY},
  {"CL_DEVICE_ADDRESS_BITS", CL_DEVICE_ADDRESS_BITS},
  {"CL_DEVICE_MAX_MEM_ALLOC_SIZE", CL_DEVICE_MAX_MEM_ALLOC_SIZE},
  {"CL_DEVICE_IMAGE_SUPPORT", CL_DEVICE_IMAGE_SUPPORT},
  {"CL_DEVICE_MAX_READ_IMAGE_ARGS", CL_DEVICE_MAX_READ_IMAGE_ARGS},
  {"CL_DEVICE_MAX_WRITE_IMAGE_ARGS", CL_DEVICE_MAX_WRITE_IMAGE_ARGS},
  {"CL_DEVICE_IMAGE2D_MAX_WIDTH", CL_DEVICE_IMAGE2D_MAX_WIDTH},
  {"CL_DEVICE_IMAGE2D_MAX_HEIGHT", CL_DEVICE_IMAGE2D_MAX_HEIGHT},
  {"CL_DEVICE_IMAGE3D_MAX_WIDTH", CL_DEVICE_IMAGE3D_MAX_WIDTH},
  {"CL_DEVICE_IMAGE3D_MAX_HEIGHT", CL_DEVICE_IMAGE3D_MAX_HEIGHT},
  {"CL_DEVICE_IMAGE3D_MAX_DEPTH", CL_DEVICE_IMAGE3D_MAX_DEPTH},
  {"CL_DEVICE_MAX_SAMPLERS", CL_DEVICE_MAX_SAMPLERS},
  {"CL_DEVICE_MAX_PARAMETER_SIZE", CL_DEVICE_MAX_PARAMETER_SIZE},
  {"CL_DEVICE_MEM_BASE_ADDR_ALIGN", CL_DEVICE_MEM_BASE_ADDR_ALIGN},
  {"CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE", CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE},
  {"CL_DEVICE_SINGLE_FP_CONFIG", CL_DEVICE_SINGLE_FP_CONFIG},
  {"CL_DEVICE_GLOBAL_MEM_CACHE_TYPE", CL_DEVICE_GLOBAL_MEM_CACHE_TYPE},
  {"CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE", CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE},
  {"CL_DEVICE_GLOBAL_MEM_CACHE_SIZE", CL_DEVICE_GLOBAL_MEM_CACHE_SIZE},
  {"CL_DEVICE_GLOBAL_MEM_SIZE", CL_DEVICE_GLOBAL_MEM_SIZE},
  {"CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE", CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE},
  {"CL_DEVICE_MAX_CONSTANT_ARGS", CL_DEVICE_MAX_CONSTANT_ARGS},
  {"CL_DEVICE_LOCAL_MEM_TYPE", CL_DEVICE_LOCAL_MEM_TYPE},
  {"CL_DEVICE_LOCAL_MEM_SIZE", CL_DEVICE_LOCAL_MEM_SIZE},
  {"CL_DEVICE_ERROR_CORRECTION_SUPPORT", CL_DEVICE_ERROR_CORRECTION_SUPPORT},
  {"CL_DEVICE_PROFILING_TIMER_RESOLUTION", CL_DEVICE_PROFILING_TIMER_RESOLUTION},
  {"CL_DEVICE_ENDIAN_LITTLE", CL_DEVICE_ENDIAN_LITTLE},
  {"CL_DEVICE_AVAILABLE", CL_DEVICE_AVAILABLE},
  {"CL_DEVICE_EXECUTION_CAPABILITIES", CL_DEVICE_EXECUTION_CAPABILITIES},
  {"CL_DEVICE_QUEUE_PROPERTIES", CL_DEVICE_QUEUE_PROPERTIES},
  {"CL_DEVICE_PLATFORM", CL_DEVICE_PLATFORM},
  {"CL_DEVICE_NAME", CL_DEVICE_NAME},
  {"CL_DEVICE_VENDOR", CL_DEVICE_VENDOR},
  {"CL_DEVICE_VERSION", CL_DEVICE_VERSION},
  {"CL_DEVICE_PROFILE", CL_DEVICE_PROFILE},
  {"CL_DEVICE_VERSION", CL_DEVICE_VERSION},
  {"CL_DEVICE_EXTENSIONS", CL_DEVICE_EXTENSIONS}
};

std::size_t const pinfo_size = sizeof(pinfo_list) / sizeof(pinfo_list[0]);
std::size_t const dinfo_size = sizeof(dinfo_list) / sizeof(dinfo_list[0]);

struct platform_data{
  cl_uint size;
  cl_platform_id* platform;
};

struct device_data{
  cl_uint size;
  cl_device_id* device;
};

platform_data* allocate_platform(void){
  cl_int ret;
  cl_uint size;
  struct platform_data* pdata = NULL;

  do{
    //Get: Number of Platforms
    ret = clGetPlatformIDs(0, NULL, &size);
    if(ret != CL_SUCCESS){
      std::cerr << "Error - clGetPlatformIDs(" << ret << ")" << std::endl;
      break;
    }

    //Allocate memory
    pdata = (platform_data*)malloc(sizeof(platform_data));
    if(pdata == NULL){
      std::cerr << "Error - Can not allocate memory" << std::endl;
      break;
    }
    pdata->size = size;
    pdata->platform = (cl_platform_id*)malloc(sizeof(cl_platform_id) * size);
    if(pdata->platform == NULL){
      std::cerr << "Error - Can not allocate memory" << std::endl;
      break;
    }

    //Get: Platform's information
    ret = clGetPlatformIDs(size, pdata->platform, &pdata->size);
    if(ret != CL_SUCCESS){
      std::cerr << "Error - clGetPlatformIDs(" << ret << ")" << std::endl;
      break;
    }
  } while(0);

  return pdata;
}

void free_platform(platform_data* pdata){
  if(pdata != NULL){
    if(pdata->platform != NULL){
      free(pdata->platform);
      pdata->platform = NULL;
    }
    free(pdata);
    pdata = NULL;
  }
}

device_data* allocate_device(cl_platform_id* platform_id){
  cl_int ret;
  cl_uint size;
  struct device_data* ddata = NULL;

  do{
    if(platform_id == NULL){
      break;
    }

    //Get: Number of Devices
    ret = clGetDeviceIDs(platform_id[0], CL_DEVICE_TYPE_ALL,
			 0, NULL, &size);
    if(ret != CL_SUCCESS){
      break;
    }

    //Allocate memory
    ddata = (device_data*)malloc(sizeof(device_data) * size);
    if(ddata == NULL){
      break;
    }
    ddata->size = size;
    ddata->device = (cl_device_id*)malloc(sizeof(cl_device_id) * size);
    if(ddata->device == NULL){
      break;
    }

    //Get: Device's information
    ret = clGetDeviceIDs(platform_id[0], CL_DEVICE_TYPE_ALL,
			 size, ddata->device, &ddata->size);
    if(ret != CL_SUCCESS){
      break;
    }
  } while(0);

  return ddata;
}

void free_device(device_data* ddata){
  if(ddata != NULL){
    if(ddata->device != NULL){
      free(ddata->device);
      ddata->device = NULL;
    }
    free(ddata);
    ddata = NULL;
  }
}

void show_platform_id(cl_platform_id* platform_id){
  std::size_t const buf_size = 1025;
  char info[buf_size];

  cl_int ret;

  for(std::size_t i=0; i<pinfo_size; ++i){
    memset(info, 0x00, sizeof(info));

    ret = clGetPlatformInfo(platform_id[0], pinfo_list[i].info_type,
			    sizeof(info), info, NULL);
    if(ret != CL_SUCCESS){
      break;
    }

    std::cout << "  " << pinfo_list[i].info_name << " = " << info << std::endl;
  }
}

void show_device_id(cl_device_id* device_id){
  std::size_t const buf_size = 1025;
  char info[buf_size];

  cl_int ret;

  for(std::size_t i=0; i<dinfo_size; ++i){
    memset(info, 0x00, sizeof(info));

    ret = clGetDeviceInfo(device_id[0], dinfo_list[i].info_type,
			  sizeof(info), info, NULL);
    if(ret != CL_SUCCESS){
      break;
    }

    std::cout << "  " << dinfo_list[i].info_name << " = " << info << std::endl;
  }
}

int main(void){
  platform_data* pdata;
  device_data* ddata;

  do{
    pdata = allocate_platform();
    if(pdata == NULL){
      break;
    }
    std::cout << "Number of Platforms = " << pdata->size << std::endl;

    for(cl_uint i=0; i<pdata->size; ++i){
      std::cout << "Platform :" << i << std::endl;
      show_platform_id(&pdata->platform[i]);

      ddata = allocate_device(&pdata->platform[i]);
      if(ddata == NULL){
	break;
      }
      std::cout << "  Number of Devices = " << ddata->size << std::endl;
      for(cl_uint j=0; j<ddata->size; ++j){
	std::cout << "------------------------------" << std::endl;
	std::cout << "  Device: " << j << std::endl;
	show_device_id(&ddata->device[j]);
      }
      free_device(ddata);
    }
  } while(0);

  free_platform(pdata);

  return 0;
}
