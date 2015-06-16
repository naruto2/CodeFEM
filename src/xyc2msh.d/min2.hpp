#ifndef MIN2_HPP_
#define MIN2_HPP_

template<class Real>
Real min2(Real a, Real b)
{
  if ( a < b ) return a;
  return b;
}
#endif
