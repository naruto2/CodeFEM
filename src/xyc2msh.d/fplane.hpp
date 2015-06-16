#ifndef FPLANE_HPP_
#define FPLANE_HPP_

template<class Real>
Real fplane(Real x, Real y, Real z,
	    Real x0,Real y0,Real z0,
	    Real x1,Real y1,Real z1,
	    Real x2,Real y2,Real z2)
{ Real xn,yn,zn;

  xn = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);
  yn = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);
  zn = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);

  return xn*x+yn*y+zn*z-(xn*x0+yn*y0+zn*z0);
}
#endif
