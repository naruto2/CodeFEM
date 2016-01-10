#include <vector>
using namespace std;
#include "est/xmesh.hpp"

int count_nodes(vector<nde>&N)
{
  int n=0, i;
  for (i=1; i<(int)N.size(); i++)
    if (N[i].a!=0) n++;
  return n;
}
