#include <CL/cl.h>

int main(int argc, char **argv)
{

  int errcode;
  uint platformCount;

  //OpenCLが実行可能なプラットフォームが存在するかどうかを検査します。
  errcode = clGetPlatformIDs(0, null, out platformCount);
  if (errcode != CL_SUCCESS)
    throw new Exception("Error at clGetPlatformIDs : " + errcode);

  //プラットフォームのIDを取得します。
  IntPtr[] platforms = new IntPtr[platformCount];
  clGetPlatformIDs(platformCount, platforms, out platformCount);

  return 0;
}
