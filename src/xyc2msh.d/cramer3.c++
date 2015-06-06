#define T double

T sarrus(T,T,T,T,T,T,T,T,T);

void cramer3(T *px,T *py,T *pz, 
             T a11,T a12,T a13,
	     T a21,T a22,T a23, 
	     T a31,T a32,T a33,
	     T  b1,T b2, T b3 )
{ T det;

  *px = sarrus(b1 ,a12,a13, b2 ,a22,a23, b3 ,a32,a33);
  *py = sarrus(a11,b1 ,a13, a21,b2 ,a23, a31,b3 ,a33);
  *pz = sarrus(a11,a12,b1 , a21,a22,b2 , a31,a32,b3 );
  det = sarrus(a11,a12,a13, a21,a22,a23, a31,a32,a33);
  if(det != 0.0){ *px/=det;*py/=det;*pz/=det;}
}

