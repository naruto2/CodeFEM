#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "est/xmesh.hpp"
#include "estiva/ary.h"
#include "estiva/std.h"

#define n static_bind(long,x)

void estiva_forgammap1(long *x)
{
  n = 1;
}

int estiva_forgammap1_loop(long *x, const char *name, vector<xyc> &Z)
{
  while ( n < (long )Z.size() ) {
    if (Z[n].label && !strcmp(Z[n].label, name) ) {
      *x = n;
      n++;
      return 1;
    }
    else {
      n++;
    }
  }
  static_free(x);
  return 0;
}

