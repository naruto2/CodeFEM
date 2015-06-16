#include <cstdio>
#include "xyc_nde.h"
#include <vector>
using namespace std;
#include "ary.h"

#define     fp2xyc(fp,Z)       estiva_fp2xyc(fp,&(Z))
#define     fprintmesh(fp,Z,N) estiva_fprintmesh(fp,Z,N)

extern "C" {
  void estiva_fp2xyc(void*, xyc**);
  void estiva_fprintmesh(void*, xyc*,nde*);
}

extern void vectoary_xyc(vector<xyc> &, xyc *);
extern void vectoary_nde(vector<nde> &, nde *);
extern void arytovec_xyc(xyc *,vector<xyc> &);
extern void arytovec_nde(nde *,vector<nde> &);
extern void delaunay1(vector<xyc>&, vector<nde>&);
extern void fprintmesh1(FILE *, vector<xyc>&, vector<nde>&);

int main(int argc, char **argv)
{ 
  FILE *fp; static xyc*Z;
  
  fp = fopen(argv[1],"r");
  fp2xyc(fp,Z); 
  fclose(fp);

  vector<xyc> Zov;
  vector<nde> Nov;
  arytovec_xyc(Z,Zov);
  delaunay1(Zov, Nov);

  fp = stdout;
  fprintmesh1(fp,Zov,Nov);
  fclose(fp);
  return 0;
}
