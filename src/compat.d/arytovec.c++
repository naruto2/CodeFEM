#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;
#include "xyc_nde.h"
#include "ary.h"

void arytovec_nde(nde *N, vector<nde> &Nv)
{
  int i;
  Nv.resize(dim1(N)+1);
  for(i =0 ;i <= dim1(N); i++) Nv[i] = N[i];
}

void vectoary_nde(vector<nde> &Nv, nde *N)
{
  int i;
  for(i=0; i<=dim1(N); i++ ) N[i] = Nv[i];
}

void arytovec_xyc(xyc *Z, vector<xyc> &Zv)
{
  int i;
  Zv.resize(dim1(Z)+1);
  for(i =0 ;i <= dim1(Z); i++) Zv[i] = Z[i];
}

void vectoary_xyc(vector<xyc> &Zv, xyc *Z)
{
  int i;
  for(i=0; i<=dim1(Z); i++ ) Z[i] = Zv[i];
}


