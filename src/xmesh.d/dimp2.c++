#include <iostream>
#include <cstdio>
#include <vector>
#include "est/xmesh.hpp"

using namespace std;

long dimp2(vector<nde>&N)
{
  long e, E;
  int n=1;
  E = N.size();
  for(e=1; e<E; e++) {
    n = max(n,N[e].A);
    n = max(n,N[e].B);
    n = max(n,N[e].C);
  }
  return n;
}
