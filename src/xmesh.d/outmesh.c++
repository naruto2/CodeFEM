#include <cstdio>
#include <iostream>
#include <vector>
using namespace std;
#include "est/xmesh.hpp"

static char *label(vector<xyc>&Z,int i)
{ static char *normal = (char*)"inner";
  return Z[i].label==NULL? normal:Z[i].label;
}

void outmesh(ostream &os, vector<xyc>&Zv,vector<nde>&Nv)
{
  char ss[999];
  int i;

  os << "<xyc>" << endl;

  for(i=1;i<(int)Zv.size();i++) { 
    snprintf(ss,990,"%4d %.9f %.9f %s",i,Zv[i].x,Zv[i].y,label(Zv,i));
    os << ss << endl;
  }

  os << "<nde>" << endl;

  for(i=1;i<(int)Nv.size();i++){
    snprintf(ss,990,"%4d  %4d %4d %4d  %4d %4d %4d",
            i, Nv[i].a,Nv[i].b,Nv[i].c,Nv[i].A,Nv[i].B,Nv[i].C);
    os << ss << endl;
  }
}
