#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <vector>
#include <estiva/ary.h>
#include <estiva/foreach.h>
#include <estiva/mesh.h>

using namespace std;

static void makeMid(vector<nde> &N)
{
  int i, e, m, n, u;
  e = N.size(); 

  for(i=1;i<e;i++)
    foreach(m,&N[i].A,&N[i].B,&N[i].C)
      m = -m;

  n=1;
  for(i=1;i<e;i++) 
    foreach(m,&N[i].A,&N[i].B,&N[i].C)
      if(m<=0){
	foreach(u,&N[-m].A,&N[-m].B,&N[-m].C) if(u== -i) u=n;
	m=n++;
      }
}

#define max(x,y) (x>y?x:y)

static char *bound(char *s1, char *s2)
{

  if(s1 == NULL || s2 == NULL) return strdup("");
  return strcmp(s1,s2)<0?s1:s2;
}

static vector<xyc> makeMV(vector<xyc> Z, vector<nde> &N)
{
  
  int i, e, m, k;
  // tsukud-y e = N.size()-2;
  e = N.size()-1;
  
  m = 1;
  for(i=1;i<=e;i++)
    foreach(k,&N[i].A,&N[i].B,&N[i].C)
      m = max(m,k);


  vector<xyc> MV;
  MV.resize(m+1);

  for(i=1;i<=e;i++){
    MV[N[i].A].label = bound(Z[N[i].b].label,Z[N[i].c].label);
    MV[N[i].B].label = bound(Z[N[i].c].label,Z[N[i].a].label);
    MV[N[i].C].label = bound(Z[N[i].a].label,Z[N[i].b].label);

    MV[N[i].A].x = (Z[N[i].b].x + Z[N[i].c].x)/2.0;
    MV[N[i].B].x = (Z[N[i].c].x + Z[N[i].a].x)/2.0;
    MV[N[i].C].x = (Z[N[i].a].x + Z[N[i].b].x)/2.0;

    MV[N[i].A].y = (Z[N[i].b].y + Z[N[i].c].y)/2.0;
    MV[N[i].B].y = (Z[N[i].c].y + Z[N[i].a].y)/2.0;
    MV[N[i].C].y = (Z[N[i].a].y + Z[N[i].b].y)/2.0;
  }

  return MV;
}

static vector<xyc> makeM(vector<xyc> Z, vector<nde> &N)
{
  makeMid(N);
  return makeMV(Z,N);
}

vector<xyc> ncpolynomial1(vector<xyc> Z, vector<nde> &N)
{
  vector<xyc> Mid;
  Mid = makeM(Z,N);

  return Mid;
}
