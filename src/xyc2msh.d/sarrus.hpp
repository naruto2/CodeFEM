#ifndef SARRUS_HPP_
#define SARRUS_HPP_

template<class Real>
Real sarrus(Real a11,Real a12,Real a13,
	    Real a21,Real a22,Real a23,
	    Real a31,Real a32,Real a33)
{
  return 
 ((a11)*(a22)*(a33)+(a21)*(a32)*(a13)+(a31)*(a12)*(a23)\
  -(a13)*(a22)*(a31)-(a23)*(a32)*(a11)-(a33)*(a12)*(a21));
}
#endif
