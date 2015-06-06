#include <cstdlib>
#include <vector>
using namespace std;
#include "xyc_nde.h"
#include "ary.h"


#define T double

int generate_fN(vector<nde>&N,int*fN)
{
  int i, n=0;
  for (i=1; i<=dim1(fN); i++) {
    while(N[n].a==0) n++;
    fN[i] = n++;
  }
  return n;
}
