#ifndef CRAMER3_HPP_
#define CRAMER3_HPP_
#include "sarrus.hpp"

template<class Real>
void cramer3(Real *px,Real *py,Real *pz, 
             Real a11,Real a12,Real a13,
	     Real a21,Real a22,Real a23, 
	     Real a31,Real a32,Real a33,
	     Real  b1,Real b2, Real b3 )
{ Real det;

  *px = sarrus(b1 ,a12,a13, b2 ,a22,a23, b3 ,a32,a33);
  *py = sarrus(a11,b1 ,a13, a21,b2 ,a23, a31,b3 ,a33);
  *pz = sarrus(a11,a12,b1 , a21,a22,b2 , a31,a32,b3 );
  det = sarrus(a11,a12,a13, a21,a22,a23, a31,a32,a33);
  if(det != 0.0){ *px/=det;*py/=det;*pz/=det;}
}
#endif
