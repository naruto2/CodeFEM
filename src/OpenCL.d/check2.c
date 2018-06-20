#include <stdio.h>
// Appleの場合
//#include <OpenCL/opencl.h>
// その他の場合
#include <CL/cl.h>

int main (){
  // プラットフォームに関するID・情報の取得
  // プラットフォームIDを取得し，格納する変数の宣言
  cl_platform_id platforms[10];
  // プラットフォームの数を格納する変数の宣言
  cl_uint num_platforms;

  // プラットフォームのリストを取得する
  // clGetPlatformIDs(cl_uint num_entries,         : 変数platformsに出力するcl_platform_idエントリの数を指定
  //                  cl_platform_id *platforms,   : 返り値(OpenCLプラットフォームのリスト)の出力先となるポインタ
  //                  cl_uint *num_platforms);     : 返り値(利用可能なOpenCLプラットフォーム数)の出力先となるポインタ
  clGetPlatformIDs(10, platforms, &num_platforms);
  printf("OpenCLプラットフォームの数: %d\n", num_platforms);

  char message[1024];
  // プラットフォームの情報を取得する
  // clGetPlatformInfo(cl_platform_id   platform,             : clGetPlatformIDsで取得したプラットフォームIDを指定
  //                   cl_platform_info param_name,           : 取得する情報を下表から指定する
  //                   size_t           param_value_size,     : param_valueが指定するメモリのサイズをバイトで指定
  //                   void            *parma_value,          : param_nameで指定した情報についてのポインタを指定，NULLなら無視
  //                   size_t          *param_value_size_ret  : param_valueにコピーされるデータサイズをバイトで返すNULLなら無視
  //
  // ・CL_PLATFORM_PROFILE: 実装がサポートするOpenCLプロファイルの名称(FULL_PROFILE or EMBEDDED_PROFILE)を返す．
  // ・CL_PLATFORM_VERSION: 実装がサポートするOpenCLのバージョンを文字列で返す．
  // ・CL_PLATFORM_NAME: プラットフォーム名を文字列で返す
  // ・CL_PLATFORM_VENDOR: プラットフォームベンダ名を文字列で返す．
  // ・CL_PLATFORM_EXTENSIONS: プラットフォームがサポートする拡張機能の名称を，リスト形式の文字列で返す．
  clGetPlatformInfo(platforms[0], CL_PLATFORM_NAME, sizeof(message), message, NULL);
  printf("名前: %s\n", message);

  // デバイスに関するID・情報の取得
  // デバイスIDを取得し，格納する変数の宣言
  cl_device_id device_id[10];
  // デバイスの数を格納する変数の宣言
  cl_uint num_devices;

  // プラットフォーム上で利用可能なデバイスのリストを取得する
  // clGetDeviceIDs( cl_platform_id  platform,      : OpenCLプラットフォームを指定
  //                 cl_device_type  device_type,   : 取得するデバイスの種類を指定
  //                                                   ・CL_DEVICE_TYPE_CPU
  //                                                   ・CL_DEVICE_TYPE_GPU
  //                                                   ・CL_DEVICE_TYPE_ACCELERATOR
  //                                                   ・CL_DEVICE_TYPE_DEFAULT
  //                                                   ・CL_DEVICE_TYPE_ALL
  //                 cl_uint         num_entries,   : devicesに追加可能なcl_device_idエントリの数を指定
  //                 cl_device_id   *devices,       : OpenCLデバイスを指定するリストの出力先のポインタ．NULLなら無視
  //                 cl_uint        *num_devices);  : device_typeとして利用可能なデバイス数の出力先のポインタ
  clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, 10, device_id, &num_devices);
  printf("デバイスの数: %d\n", num_devices);
  for(int i = 0; i < num_devices; i++){
    // clGetDeviceInfo(cl_device_id  device,                : clGetDeviceIDsで取得したデバイスを指定します
    //                 cl_device_ifo param_name,            : 取得する情報を指定
    //                 size_t        param_value_size,      : param_valueが指すメモリのサイズをバイトで指定
    //                 void         *param_value,           : 取得する情報についての値が返されるポインタ NULLならば無視
    //                 size_t       *param_value_size_ret); : param_valueにコピーされるデータのサイズをバイトで返すNULL無視

    // CL_DEVICE_NAME: デバイス名を文字列で返す
    clGetDeviceInfo(device_id[i], CL_DEVICE_NAME, sizeof(message), message, NULL);
    printf("デバイス名: %s\n", message);

    // CL_DEVICE_MAX_COMPUTE_UNITS: OpenCLデバイス上の並列演算コア数を返す
    cl_uint units;
    clGetDeviceInfo(device_id[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(units), &units, NULL);
    printf("最大演算ユニット数: %u\n", units);

    // CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: データ並列実行モデルので用いるグローバル・ローカルワークアイテムIDの次元数の最大値
    cl_uint max_dimensions;
    clGetDeviceInfo(device_id[i], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(max_dimensions), &max_dimensions, NULL);
    printf("最大次元数: %d\n", max_dimensions);

    // CL_DEVICE_MAX_WORK_ITEM_SIZE:
    size_t max_cores[3];
    clGetDeviceInfo(device_id[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(max_cores), max_cores, NULL);
    printf("最大並列動作数:");

    for (int dimension = 0; dimension < max_dimensions; dimension++) printf("%ld ", max_cores[dimension]);
    printf("\n");

    cl_uint max_workitemsize[3];
    clGetDeviceInfo(device_id[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(max_workitemsize), max_workitemsize, NULL);
    printf("最大work group size:");
    printf("%u \n", max_workitemsize[0]);
#if 0
    for (int dimension = 0; dimension < max_dimensions; dimension++)
      printf("%u ", max_workitemsize[dimension]);
    printf("\n");
#endif
    cl_uint compute_unit = 0;
    clGetDeviceInfo(device_id[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint),&compute_unit, NULL);
    printf("演算ユニット数:");
    printf("%u \n", compute_unit);
    int ret;
    ret = clGetDeviceInfo(device_id[i], CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, sizeof(cl_uint),&compute_unit, NULL);
    if ( ret != CL_SUCCESS ) printf("error\n");
    printf("int型ベクタの望ましい長さ:");
    printf("%u \n", compute_unit);

  }
  return 0;
}
