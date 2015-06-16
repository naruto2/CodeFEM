#include <cstdio>
#include <vector>
using namespace std;
#include "xyc_nde.h"

extern void fp2xyc(FILE *,vector<xyc>&);
extern void delaunay(vector<xyc>&, vector<nde>&);
extern void fprintmesh(FILE *, vector<xyc>&, vector<nde>&);

int main(int argc, char **argv)
{ 
  FILE *fp;
  vector<xyc> Z;
  vector<nde> N;
  
  if ( NULL == (fp = fopen(argv[1],"r"))) return 0;
  fp2xyc(fp,Z); 
  fclose(fp);

  delaunay(Z, N);

  fprintmesh(stdout,Z,N);

  return 0;
}
