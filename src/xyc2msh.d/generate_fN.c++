#include <cstdlib>
#include <vector>
using namespace std;
#include "xyc_nde.h"
#include "ary.h"

int generate_fN(vector<nde>&N,vector<int>&fNv)
{
  int i, n=0;
  for (i=1; i<(int)fNv.size(); i++) {
    while(N[n].a==0) n++;
    fNv[i] = n++;
  }
  return n;
}
