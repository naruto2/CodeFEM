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
OpenCL���g�p���邽�߂ɕK�v�Ȃ��̂����܂��B
*/
void OpenCL::initialize(int platformIdx, int deviceIdx)
{
	// �v���b�g�t�H�[���擾
	cl_uint platformNumber = 0;
	cl_platform_id platformIds[8];
	checkError(clGetPlatformIDs(8, platformIds, &platformNumber));
	cl_platform_id platform = platformIds[platformIdx];

	// �f�o�C�X�擾
	cl_uint deviceNumber = 0;
	cl_device_id deviceIds[8];
	checkError(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 8, deviceIds, &deviceNumber));
	gDevice = deviceIds[deviceIdx];

	// �R���e�L�X�g�i�������m�ۂȂǂɎg�p�j�ƃR�}���h�L���[�i�J�[�l���̎��s�ȂǂɎg�p�j���쐬
	gContext = clCreateContext(NULL, 1, &gDevice, NULL, NULL, NULL);
	gCommandQueue = clCreateCommandQueue(gContext, gDevice, 0, NULL);

	// �J�[�l���v���O�����̃R���p�C��
	gProgram = compileProgram("kernel.cl");
	gKernel = createKernel(gProgram, "detectLine");

	// �f�o�C�X�Ŏg�p���郁�����I�u�W�F�N�g���m��
	gResult = clCreateBuffer(gContext, CL_MEM_READ_WRITE, sizeof(char) * 640 * 480, NULL, NULL);
	gOrigin = clCreateBuffer(gContext, CL_MEM_READ_WRITE, sizeof(char) * 640 * 480 * 3, NULL, NULL);
}

/*
OpenCL�J�[�l�����Ăяo���Čv�Z�����s���܂��B
*/
void OpenCL::detectLine(unsigned char* result, unsigned char* origin)
{
	cl_int2 size = { 640, 480 };
	// �z�X�g����f�o�C�X�Ƀ������]��
	checkError(clEnqueueWriteBuffer(gCommandQueue, gOrigin, CL_TRUE, 0, sizeof(char) * 640 * 480 * 3, origin, 0, NULL, NULL));
	// �������I�u�W�F�N�g���J�[�l���֐��̈����ɃZ�b�g
	checkError(clSetKernelArg(gKernel, 0, sizeof(cl_mem), &gResult));
	checkError(clSetKernelArg(gKernel, 1, sizeof(cl_mem), &gOrigin));
	checkError(clSetKernelArg(gKernel, 2, sizeof(cl_int2), &size));
	// �J�[�l���̕�����s����ݒ�
	size_t workSize[2] = { 640, 480 };
	// �J�[�l���̌Ăяo��
	checkError(clEnqueueNDRangeKernel(gCommandQueue, gKernel, 2, NULL, workSize, NULL, 0, NULL, NULL));
	// �f�o�C�X����z�X�g�Ƀ������]��
	checkError(clEnqueueReadBuffer(gCommandQueue, gResult, CL_TRUE, 0, sizeof(char) * 640 * 480, result, 0, NULL, NULL));
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
OpenCL�̃J�[�l���v���O�������R���p�C�����āA�������ꂽ�v���O�����I�u�W�F�N�g��Ԃ��܂��B
*/
cl_program OpenCL::compileProgram(char* fileName)
{
	// �v���O�����̓ǂݍ���
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

	// �v���O�����̃R���p�C��
	cl_program program = clCreateProgramWithSource(gContext, 1, (const char**)&sourceString, (const size_t*)&sourceSize, NULL);
	cl_int err = clBuildProgram(program, 1, &gDevice, NULL, NULL, NULL);
	// �R���p�C���Ɏ��s�����ꍇ�̓G���[���e��\��
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
�v���O�����I�u�W�F�N�g����J�[�l���I�u�W�F�N�g�����܂��B
*/
cl_kernel OpenCL::createKernel(cl_program program, char* kernelName)
{
	return clCreateKernel(program, kernelName, NULL);
}