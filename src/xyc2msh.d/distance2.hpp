#ifndef DISTANCE2_HPP_
#define DISTANCE2_HPP_

template<class Real>
Real distance2(Real x0,Real y0,Real x1,Real y1)
{
  return ((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
}
#endif
