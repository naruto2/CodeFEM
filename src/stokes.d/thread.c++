#include <iostream>
#include <thread>
#include <mutex>
#include <unistd.h>
#include <chrono>
#include <cstdio>
#include <unistd.h>
#include <time.h>

void baton2(void)
{
  static clock_t wall=0;
  static long i=-1;
  static char* s[4];
  if ( i == -1 ) {
    setbuf(stderr,NULL);
    fprintf(stderr," ");
    s[0] = (char*)"|";
    s[1] = (char*)"/";
    s[2] = (char*)"-";
    s[3] = (char*)"\\";
    i = 0;
  }
  i++;
  i %= 4;
  fprintf(stderr,"\b%s",s[ i ]);
  usleep(100000);
}

using namespace std;

int main(int argc, char **argv) {

  thread t1([]() { for(int i=0; i<=10000000; i++) {baton2();}});

  t1.detach();   // do not wait t1 to finish
  sleep(10);
  t1.~thread();
  return 0;
}
