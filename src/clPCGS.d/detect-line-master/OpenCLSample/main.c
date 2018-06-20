#include <stdio.h>
#include <CL/cl.h>

int main(int argc, char **argv)
{

	// プラットフォーム取得
	cl_uint platformNumber = 0;
	cl_platform_id platformIds[8];
	clGetPlatformIDs(8, platformIds, &platformNumber);

	char string[256];
	cl_device_type type;
	cl_uint value;
	size_t sizes[3];
	cl_ulong ulvalue;
	for (int i = 0; i < platformNumber; i++)
	{
		printf("platform idx : %d\n", i);
		cl_platform_id platform = platformIds[i];
		clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, 256, string, NULL);
		printf("platform vendor : %s\n", string);
		clGetPlatformInfo(platform, CL_PLATFORM_NAME, 256, string, NULL);
		printf("platform name : %s\n", string);
		clGetPlatformInfo(platform, CL_PLATFORM_VERSION, 256, string, NULL);
		printf("platform version : %s\n", string);

		// デバイス取得
		cl_uint deviceNumber = 0;
		cl_device_id deviceIds[8];
		clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 8, deviceIds, &deviceNumber);
		for (int j = 0; j < deviceNumber; j++)
		{
			printf("	device idx : %d\n", j);
			cl_device_id device = deviceIds[j];
			clGetDeviceInfo(device, CL_DEVICE_NAME, 256, string, NULL);
			printf("	device name : %s\n", string);
			clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof(cl_device_type), &type, NULL);
			if (type == CL_DEVICE_TYPE_CPU) printf("	device type : CPU\n");
			if (type == CL_DEVICE_TYPE_GPU) printf("	device type : GPU\n");
			if (type == CL_DEVICE_TYPE_ACCELERATOR) printf("	device type : ACCELERATOR\n");
			clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &value, NULL);
			printf("	device max compute units : %d\n", value);
			clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t) * 3, sizes, NULL);
			printf("	device max work item sizes : [%d][%d][%d]\n", sizes[0], sizes[1], sizes[2]);
			clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(cl_uint), &value, NULL);
			printf("	device max work group size : %d\n", value);
			clGetDeviceInfo(device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), &ulvalue, NULL);
			printf("	device max mem alloc size : %d\n", ulvalue);
			clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(cl_ulong), &ulvalue, NULL);
			printf("	device max constant buffer size : %d\n", ulvalue);
		}
	}

}