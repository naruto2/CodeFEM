#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include "est/op.hpp"
#include "est/xmesh.hpp"
using namespace std;

extern vector<xyc> ncpolynomial1(vector<xyc> Z, vector<nde> &N );


int main(int argc, char ** argv)
{
  initop(argc,argv);
  vector<xyc> Z;
  vector<nde> N;
  ifstream ifs(argv[1]);

  if (!ifs) {
    cerr << "Error: 入力ストリームを開けませんでした(xmesh)" << endl;
    return 0;
  }

  in2xyc(ifs,Z);
  ifs.close();

  delaunay(Z, N);
  
  vector<xyc> Mid = ncpolynomial1(Z,N);

  for ( unsigned long i=1; i< Mid.size(); i++ ) {
    printf("%d %f %f %s\n",i, Mid[i].x, Mid[i].y, Mid[i].label );
  }
  printf("-------------------------------------\n");

  for ( unsigned long e=1; e<N.size()-1; e++ ) {
    printf("%d %d %d %d %d %d %d\n",e,N[e].a,N[e].b,N[e].c,N[e].A,N[e].B,N[e].C);
  }    
  return 0;
}
