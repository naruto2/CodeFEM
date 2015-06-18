#include <cstdio>
#include <vector>
#include <fstream>
using namespace std;
#include "xyc_nde.h"

extern void ifs2xyc(ifstream &, vector<xyc>&);
extern void delaunay(vector<xyc>&, vector<nde>&);
extern void fprintmesh(FILE *, vector<xyc>&, vector<nde>&);

int main(int argc, char **argv)
{ 
  FILE *fp;
  vector<xyc> Z;
  vector<nde> N;
  ifstream ifs(argv[1]);
  
  if ( NULL == (fp = fopen(argv[1],"r"))) return 0;
  ifs2xyc(ifs,Z); 
  fclose(fp);

  delaunay(Z, N);

  fprintmesh(stdout,Z,N);

  return 0;
}
