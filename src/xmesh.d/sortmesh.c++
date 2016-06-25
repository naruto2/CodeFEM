#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <functional>
#include "est/xmesh.hpp"

using namespace std;

typedef struct { int i, a, b, c, A, B, C; double x, y; } ende;
typedef struct { int i; double x, y; char *label; } exyc;

bool operator<(const exyc& left, const exyc& right)
{
  if ( left.y == right.y ) return left.x < right.x;
  return left.y < right.y;
}

bool operator>(const exyc& left, const exyc& right)
{
  if ( left.y == right.y ) return left.x > right.x;
  return left.y > right.y;
}


bool operator<(const ende& left, const ende& right)
{
  if ( left.y == right.y ) return left.x < right.x;
  return left.y < right.y;
}

bool operator>(const ende& left, const ende& right)
{
  if ( left.y == right.y ) return left.x > right.x;
  return left.y > right.y;
}

extern void sortmesh(vector<xyc> &Z0, vector<nde> &N0)
{
  unsigned long i;

  vector<exyc> Z(Z0.size());
  vector<ende> N(N0.size());
  
  // 重心を求めるフェーズ
  for ( i=0; i<Z.size(); i++) {
    Z[i].i = i;
    Z[i].x = Z0[i].x;
    Z[i].y = Z0[i].y;
    Z[i].label = Z0[i].label;
  }
  for ( i=0; i<N.size(); i++) {
    N[i].i = i;
    N[i].a = N0[i].a;
    N[i].b = N0[i].b;
    N[i].c = N0[i].c;
    N[i].A = N0[i].A;
    N[i].B = N0[i].B;
    N[i].C = N0[i].C;
    N[i].x = (Z[N0[i].a].x + Z[N0[i].b].x + Z[N0[i].c].x ) / 3.0;
    N[i].y = (Z[N0[i].a].y + Z[N0[i].b].y + Z[N0[i].c].y ) / 3.0;
  }

  // ソートするフェーズ
  sort(&Z[1],&Z[Z.size()]);
  sort(&N[1],&N[N.size()]);


  // 逆関数を構成するフェーズ
  vector<int> f(Z.size());
  for ( i = 0; i<Z.size(); i++) { f[Z[i].i] =i; }

  vector<int> g(N.size());
  for ( i = 0; i<N.size(); i++) { g[N[i].i] =i; }
  

  // 要素-節点対応表を変換するフェーズ
  for ( i = 0; i<N.size(); i++ )
    {
      N[i].a = f[N[i].a];
      N[i].b = f[N[i].b];
      N[i].c = f[N[i].c];
    }
  for ( i = 0; i<N.size(); i++ )
    {
      N[i].A = g[N[i].A];
      N[i].B = g[N[i].B];
      N[i].C = g[N[i].C];
    }
  
  // 書き戻すフェーズ
  for ( i=0; i<Z.size(); i++) {
    Z0[i].x = Z[i].x;
    Z0[i].y = Z[i].y;
    Z0[i].label = Z[i].label;
  }
  for ( i=0; i<N.size(); i++) {
    N0[i].a = N[i].a;
    N0[i].b = N[i].b;
    N0[i].c = N[i].c;
    N0[i].A = N[i].A;
    N0[i].B = N[i].B;
    N0[i].C = N[i].C;
  }

#if 0 // データデバッグプリントフェーズ  
  printf("<xyc>\n");
  for ( i=0; i<Z.size(); i++)
    printf("%ld %f %f %s\n",i,Z[i].x,Z[i].y,Z[i].label);

  printf("<nde>\n");
  for ( i=0; i<N.size(); i++)
    printf("%ld %d %d %d %d %d %d\n",i,
	   N[i].a,N[i].b,N[i].c,N[i].A,N[i].B,N[i].C);
  exit(0);
#endif 
}

