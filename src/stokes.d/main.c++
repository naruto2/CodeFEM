#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include "est/op.hpp"
#include "est/xmesh.hpp"
#include "est/matrix.hpp"
using namespace std;

extern vector<xyc> ncpolynomial1(vector<xyc> Z, vector<nde> &N );


#define length(a,b) \
  ((Z[b].x-Z[a].x)*(Z[b].x-Z[a].x)+(Z[b].y-Z[a].y)*(Z[b].y-Z[a].y))

#define inner(a,b,c) \
  ((Z[c].x-Z[a].x)*(Z[b].x-Z[c].x)+(Z[c].y-Z[a].y)*(Z[b].y-Z[c].y))


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


static matrix M__(vector<xyc> Mid, vector<nde> N, vector<double> S)
{
  matrix M;
  long  e, a, b, c, A, B, C, m, n;

  m = Mid.size()-1;
  n = N.size()-1;
  M.resize(m+1);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;

    M[A][A] += S[e]/3.0;
    M[B][B] += S[e]/3.0;
    M[C][C] += S[e]/3.0;
  }
  return M;
}

static matrix K__(vector<xyc> Mid, vector<xyc> Z, vector<nde> N, vector<double> S)
{
  matrix K;
  long e, a, b, c, A, B, C, m, n;
  double s;

  m = Mid.size()-1;
  n = N.size()-1;
  K.resize(m+1);
  
  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;
    s = S[e];

    K[A][A] += length(b,c)/s; K[A][B] +=inner(a,b,c)/s; K[A][C] +=inner(a,c,b)/s;
    K[B][A] +=inner(b,a,c)/s; K[B][B] += length(c,a)/s; K[B][C] +=inner(b,c,a)/s;
    K[C][A] +=inner(c,a,b)/s; K[C][B] +=inner(c,b,a)/s; K[C][C] += length(a,b)/s;

  }
  return K;
}

static matrix Hx__(vector<xyc> Mid, vector<xyc> Z, vector<nde> N)
{
  matrix Hx;
  long e, a, b, c, A, B, C, m, n;

  m = Mid.size()-1;
  n = N.size()-1;
  Hx.resize(m+1);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;

    Hx[A][e] += -(Z[c].y-Z[b].y);
    Hx[B][e] += -(Z[a].y-Z[c].y);
    Hx[C][e] += -(Z[b].y-Z[a].y);
  }
  return Hx;
}


static matrix Hy__(vector<xyc> Mid, vector<xyc> Z, vector<nde> N)
{
  matrix Hy;
  long e, a, b, c, A, B, C, m, n;

  m = Mid.size()-1;
  n = N.size()-1;
  Hy.resize(m+1);

  for(e=1;e<=n;e++){
    a = N[e].a, b = N[e].b, c = N[e].c;
    A = N[e].A, B = N[e].B, C = N[e].C;

    Hy[A][e] += (Z[c].x-Z[b].x);
    Hy[B][e] += (Z[a].x-Z[c].x);
    Hy[C][e] += (Z[b].x-Z[a].x);

  }
  return Hy;
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



  matrix M = M__(Mid, N, S);
  matrix K = K__(Mid, Z, N, S);
  matrix Hx= Hx__(Mid, Z, N);
  matrix Hy= Hy__(Mid, Z, N);
  printf("m = %ld  n = %ld\n",m,n);
  plotmatrix(Hy);
  return 0;
}
