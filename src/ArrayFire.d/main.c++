#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

double gettimeofday_sec()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}


#define PRINT_MAT(X) cout << #X << ":\n" << X << endl << endl

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
  size_t num = 6;
  if(argc >= 2) {
    num = atoi(argv[1]);
  }

  MatrixXf A = MatrixXf::Zero(num, num);
  for(int j = 0; j < num; ++j) {
    for(int i = 0; i < num; ++i) {
      A(j,i) = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    }
  }

  volatile double t1 = gettimeofday_sec();

  //FullPivLU< MatrixXf > flu(A);

  volatile double t2 = gettimeofday_sec();

  PartialPivLU< MatrixXf > plu(A);

  volatile double t3 = gettimeofday_sec();

  std::cout << "FullPiv time: " << t2 - t1 << " sec." << std::endl;
  std::cout << "PartPiv time: " << t3 - t2 << " sec." << std::endl;

  #ifdef ENABLE_PRINT
  PRINT_MAT(A);
  #endif

  return 0;
}
