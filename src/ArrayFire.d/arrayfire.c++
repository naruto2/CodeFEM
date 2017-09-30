#include <iostream>
#include <cstdio>
#include <arrayfire.h>
#include <vector>
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


int main(int argc, char** argv)
{

  size_t num = 6;
  if(argc >= 2) {
    num = atoi(argv[1]);
  }

  std::vector<float> vec(num * num);
  for(int i = 0; i < num; ++i) {
    vec[i] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
  }

  af::array A(num, num, &vec[0]);
  af::array lower(num,num);
  af::array upper(num,num);
  af::array out(num, num);
  af::array pivot(num);
  const bool is_lapack_piv=true;
  
  volatile double t1 = gettimeofday_sec();

  af::lu(lower,upper,pivot,A);

  volatile double t2 = gettimeofday_sec();

  #ifdef ENABLE_PRINT
  af::print(A);
  //      af::print(lu);
  #endif

  std::cout << "time: " << t2 - t1 << " sec." << std::endl;

  return 0;
}
