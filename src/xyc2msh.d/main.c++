#include <cstdio>
#include <vector>
using namespace std;
#include "xyc_nde.h"

extern void fp2xyc(FILE *,vector<xyc>&);
extern void delaunay(vector<xyc>&, vector<nde>&);
extern void fprintmesh1(FILE *, vector<xyc>&, vector<nde>&);

int main(int argc, char **argv)
{ 
  FILE *fp;
  vector<xyc> Z;
  vector<nde> N;
  
  fp = fopen(argv[1],"r");
  fp2xyc(fp,Z); 
  fclose(fp);

  delaunay(Z, N);

  fp = stdout;
  fprintmesh1(fp,Z,N);
  fclose(fp);

  return 0;
}
