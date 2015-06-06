#include <cstdio>
#include "xyc_nde.h"
#include <vector>

#define     fp2xyc(fp,Z)       estiva_fp2xyc(fp,&(Z))
#define     fprintmesh(fp,Z,N) estiva_fprintmesh(fp,Z,N)

extern "C" {
void estiva_fp2xyc(void*, xyc**);
void estiva_fprintmesh(void*, xyc*,nde*);
}

extern void delaunay1(xyc**, nde**);

int main(int argc, char **argv)
{ 
  FILE *fp; static xyc*Z; static nde*N; 
  
  fp = fopen(argv[1],"r");
  fp2xyc(fp,Z); 
  fclose(fp);
  
  delaunay1(&Z,&N);
  
  fp = stdout;
  fprintmesh(fp,Z,N);
  fclose(fp);
  return 0;
}
