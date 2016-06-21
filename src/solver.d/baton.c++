#include <cstdio>
#include <unistd.h>
#include <time.h>

void baton(void)
{
  static clock_t wall=0;
  if ( wall < clock()) {
    wall = clock()+CLOCKS_PER_SEC/8;
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
  }
}
