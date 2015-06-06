#define T double

T fplane(T x, T y, T z,
	 T x0,T y0,T z0,
	 T x1,T y1,T z1,
	 T x2,T y2,T z2)
{ T xn,yn,zn;

  xn = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);
  yn = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);
  zn = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);

  return xn*x+yn*y+zn*z-(xn*x0+yn*y0+zn*z0);
}
