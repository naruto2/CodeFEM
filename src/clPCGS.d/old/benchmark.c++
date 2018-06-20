#include <cstdio>
#include <chrono>

#define n (65535*512)

int main( int argc, char ** argv)
{
  static double z[n], x[n], y[n];

  for (int i = 0; i < n; i++){
    z[i] = 0.123;
    x[i] = 1.0;
    y[i] = (double)i;
  }

  std::chrono::system_clock::time_point  start, end;
  start = std::chrono::system_clock::now();
  
  for (int i = 0; i < n; i++){
    z[i] = x[i] * y[i];
  }


  end = std::chrono::system_clock::now();
  double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

  printf("%.3f\n",elapsed/1000);
  
  return 0;
}
