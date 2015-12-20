#include <vector>
using namespace std;
#include "est/xyc2msh.hpp"

int generate_fN(vector<nde>&N,vector<int>&fN)
{
  int i, n=0;
  for (i=1; i<(int)fN.size(); i++) {
    while(N[n].a==0) n++;
    fN[i] = n++;
  }
  return n;
}
