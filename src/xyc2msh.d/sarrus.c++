#define T double

T sarrus(T a11,T a12,T a13,T a21,T a22,T a23, T a31,T a32,T a33)
{
  return 
 ((a11)*(a22)*(a33)+(a21)*(a32)*(a13)+(a31)*(a12)*(a23)\
  -(a13)*(a22)*(a31)-(a23)*(a32)*(a11)-(a33)*(a12)*(a21));
}
