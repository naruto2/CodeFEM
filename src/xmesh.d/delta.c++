#include "est/xmesh.hpp"

double delta(int i, vector<xyc>&Z, vector<nde>&N)
{
  double xi, xj, xk, yi, yj, yk;

  xi = Z[N[i].a].x;
  xj = Z[N[i].b].x;
  xk = Z[N[i].c].x;
  yi = Z[N[i].a].y;
  yj = Z[N[i].b].y;
  yk = Z[N[i].c].y;
  return (xi*yj+xj*yk+xk*yi-yi*xj-yj*xk-yk*xi)/2.0;
}
