#include <stdio.h>
#include <malloc.h>
#include <CL/cl.h>
#include "opencl.h"

cl_device_id gDevice;
cl_command_queue gCommandQueue;
cl_context gContext;

cl_program gProgram;
cl_kernel gKernel;

cl_mem gResult;
cl_mem gOrigin;

/*
OpenCLを使用するために必要なものを作ります。
*/
void OpenCL::initialize(int platformIdx, int deviceIdx)
{
#if 0
  // プラットフォーム取得
	cl_uint platformNumber = 0;
	cl_platform_id platformIds[8];
	checkError(clGetPlatformIDs(8, platformIds, &platformNumber));
	cl_platform_id platform = platformIds[platformIdx];

	// デバイス取得
	cl_uint deviceNumber = 0;
	cl_device_id deviceIds[8];
	checkError(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 8, deviceIds, &deviceNumber));
	gDevice = deviceIds[deviceIdx];

	// コンテキスト（メモリ確保などに使用）とコマンドキュー（カーネルの実行などに使用）を作成
	gContext = clCreateContext(NULL, 1, &gDevice, NULL, NULL, NULL);

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
	gCommandQueue = clCreateCommandQueue(gContext, gDevice, 0, NULL);
#pragma GCC diagnostic warning "-Wdeprecated-declarations"
	
	// カーネルプログラムのコンパイル
	gProgram = compileProgram((char*)"kernel.cl");
	gKernel = createKernel(gProgram, (char*)"detectLine");

	// デバイスで使用するメモリオブジェクトを確保
	gResult = clCreateBuffer(gContext, CL_MEM_READ_WRITE, sizeof(char) * 640 * 480, NULL, NULL);
	gOrigin = clCreateBuffer(gContext, CL_MEM_READ_WRITE, sizeof(char) * 640 * 480 * 3, NULL, NULL);
}

/*
OpenCLカーネルを呼び出して計算を実行します。
*/
void OpenCL::detectLine(unsigned char* result, unsigned char* origin)
{
	cl_int2 size = { 640, 480 };
	// ホストからデバイスにメモリ転送
	checkError(clEnqueueWriteBuffer(gCommandQueue, gOrigin, CL_TRUE, 0, sizeof(char) * 640 * 480 * 3, origin, 0, NULL, NULL));
	// メモリオブジェクトをカーネル関数の引数にセット
	checkError(clSetKernelArg(gKernel, 0, sizeof(cl_mem), &gResult));
	checkError(clSetKernelArg(gKernel, 1, sizeof(cl_mem), &gOrigin));
	checkError(clSetKernelArg(gKernel, 2, sizeof(cl_int2), &size));
	// カーネルの並列実行数を設定
	size_t workSize[2] = { 640, 480 };
	// カーネルの呼び出し
	checkError(clEnqueueNDRangeKernel(gCommandQueue, gKernel, 2, NULL, workSize, NULL, 0, NULL, NULL));
	// デバイスからホストにメモリ転送
	checkError(clEnqueueReadBuffer(gCommandQueue, gResult, CL_TRUE, 0, sizeof(char) * 640 * 480, result, 0, NULL, NULL));
#endif
}

void OpenCL::release()
{
	clFinish(gCommandQueue);
	clReleaseMemObject(gResult);
	clReleaseMemObject(gOrigin);
	clReleaseCommandQueue(gCommandQueue);
	clReleaseContext(gContext);
}

/*
OpenCLのカーネルプログラムをコンパイルして、生成されたプログラムオブジェクトを返します。
*/
cl_program OpenCL::compileProgram(char* fileName)
{
	// プログラムの読み込み
	FILE* fp;
	fp = fopen(fileName, "r");
	if (fp == NULL)
	{
		printf("%s load failed\n", fileName);
		return NULL;
	}

	fseek(fp, 0, SEEK_END);
	const int filesize = ftell(fp);

	fseek(fp, 0, 0);
	char* sourceString = (char*)malloc(filesize);
	size_t sourceSize = fread(sourceString, sizeof(char), filesize, fp);
	fclose(fp);

	// プログラムのコンパイル
	cl_program program = clCreateProgramWithSource(gContext, 1, (const char**)&sourceString, (const size_t*)&sourceSize, NULL);
	cl_int err = clBuildProgram(program, 1, &gDevice, NULL, NULL, NULL);
	// コンパイルに失敗した場合はエラー内容を表示
	if (err != CL_SUCCESS)
	{
		size_t logSize;
		clGetProgramBuildInfo(program, gDevice, CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
		char* buildLog = (char*)malloc((logSize + 1));
		clGetProgramBuildInfo(program, gDevice, CL_PROGRAM_BUILD_LOG, logSize, buildLog, NULL);
		printf("%s", buildLog);
		free(buildLog);
	}
	free(sourceString);
	return program;
}

/*
プログラムオブジェクトからカーネルオブジェクトを作ります。
*/
cl_kernel OpenCL::createKernel(cl_program program, char* kernelName)
{
	return clCreateKernel(program, kernelName, NULL);
}
