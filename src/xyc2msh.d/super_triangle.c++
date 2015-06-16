#include <cstdlib>
#include <vector>
using namespace std;
#include "xyc_nde.h"
#include "max2.hpp"
#include "min2.hpp"

#define T double

void super_triangle(vector<xyc>&Z)
{
  T xmin,ymin,xmax,ymax,length;
  int z, i;
  static char super_node[] = "super_node";
  z=Z.size(); z -=4;
  xmin = Z[1].x; ymin = Z[1].y;
  xmax = Z[1].x; ymax = Z[1].y;
  for (i=2; i<=z; i++){
    xmin=min2(xmin,Z[i].x); xmax=max2(xmax,Z[i].x);
    ymin=min2(ymin,Z[i].y); ymax=max2(ymax,Z[i].y);
  }

  length = max2(xmax-xmin, ymax-ymin);
  Z[z+1].x = xmin - length/2.0;
  Z[z+1].y = ymin - length/2.0;
  Z[z+2].x = xmin + length*3.0;
  Z[z+2].y = ymin - length/2.0;
  Z[z+3].x = xmin - length/2.0;
  Z[z+3].y = ymin + length*3.0;
  Z[z+1].label = super_node;
  Z[z+2].label = super_node;
  Z[z+3].label = super_node;
}   
