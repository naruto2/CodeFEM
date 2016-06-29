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

void baton2(long i)
{
  static char* s[4];
  if ( i == 0 ) {
    setbuf(stdout,NULL);
    s[0] = (char*)"|";
    s[1] = (char*)"/";
    s[2] = (char*)"-";
    s[3] = (char*)"\\";
  }
  i %= 4;
  printf("\b%s",s[ i ]);
  usleep(8000000);
}

static thread t;
static int flag = 1;
void start_baton(void)
{
  flag = 1;
  thread t1([]() { for(int i=0; i<=10000000; i++) {baton2(i);}});
  t1.detach();
  t = move(t1);
}

int end_baton(void)
{
  t.~thread();
  return flag--;
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
