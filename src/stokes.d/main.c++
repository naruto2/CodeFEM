#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include "est/op.hpp"
#include "est/xmesh.hpp"
using namespace std;

extern vector<xyc> ncpolynomial1(vector<xyc> Z, vector<nde> &N );


static vector<double> S_(vector<xyc> Z, vector<nde> N)
{

  long i, e;
  double xi,xj,xk,yi,yj,yk;

  e = N.size()-1;
  vector<double> S;
  S.resize(e+1);
  
  for(i=1;i<=e;i++){
    xi = Z[N[i].a].x;
    xj = Z[N[i].b].x;
    xk = Z[N[i].c].x;
    yi = Z[N[i].a].y;
    yj = Z[N[i].b].y;
    yk = Z[N[i].c].y;
    S[i] = xi*yj+xj*yk+xk*yi-yi*xj-yj*xk-yk*xi;
    S[i] /= 2.0;
  }
  return S;
}


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
  printf("N.size()-1 = %ld\n",N.size()-1);
  
  vector<xyc> Mid = ncpolynomial1(Z,N);

  for ( unsigned long i=1; i< Mid.size(); i++ ) {
    printf("%d %f %f %s\n",i, Mid[i].x, Mid[i].y, Mid[i].label );
  }
  printf("-------------------------------------\n");

  for ( unsigned long e=1; e<N.size(); e++ ) {
    printf("%d %d %d %d %d %d %d\n",e,N[e].a,N[e].b,N[e].c,N[e].A,N[e].B,N[e].C);
  }    
  printf("-------------------------------------\n");
  
  vector<double> S = S_(Z,N);
  
  for ( unsigned long e=1; e<N.size(); e++ ) {
    printf("%ld %f\n",e,S[e]);
  }

  unsigned long m = Mid.size()-1;
  unsigned long n = S.size()-1;
  unsigned long NUM = 2*m+n;
  printf("m = %ld  n = %ld\n",m,n);
  return 0;
}
