#include <iostream>
#include <thread>
#include <mutex>
#include <unistd.h>
#include <chrono>
#include <cstdio>
#include <unistd.h>
#include <time.h>

extern void start_baton(void);
extern int end_baton(void);
#define batonth() for(start_baton();end_baton();)

using namespace std;

static int exit_flag=0;
static int flag = -1;

void baton2(void)
{
  static char* s[4];
  setbuf(stdout,NULL);
  s[0] = (char*)"|";
  s[1] = (char*)"/";
  s[2] = (char*)"-";
  s[3] = (char*)"\\";


  for ( int i=0; ; i++, i %= 4){
    printf("\b%s",s[ i ]);
    usleep(80000);
    if ( exit_flag==0 ) continue;
    break;
  }
}


void start_baton(void)
{
  exit_flag=0;
  flag = -1;
  thread t1([]() {baton2();});
  t1.detach();
}

int end_baton(void)
{
  if (flag==0){ exit_flag=1; return 0;}
  flag =  0;
  return 1;
}

#if 0
int main(int argc, char **argv) {

  batonth() sleep(10);
  printf("hello\n");

  batonth() sleep(10);
  printf("world\n");

  batonth() sleep(10);
  printf("hello\n");

  batonth() sleep(10);
  printf("world\n");

  batonth() sleep(10);


  return 0;
}
#endif
