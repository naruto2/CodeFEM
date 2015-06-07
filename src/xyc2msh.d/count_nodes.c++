#include <cstdlib>
#include <vector>
using namespace std;
#include "xyc_nde.h"

int count_nodes(vector<nde>&N)
{
  int n=0, i;
  for (i=1; i<=(int)(N.size()-1); i++)
    if (N[i].a!=0) n++;
  return n;
}
