static double axij(long i, long j, double *u, double B1, double B2, double B3)
{
  double Au=0.0, Bu=0.0, Cu=0.0, Du=0.0, Eu=0.0, Fu=0.0,
         alphaBj, betaBj, gammaBj, axij;

  for (long k = 1; k <= 6; k++) {
    Au += a(k)*u[k];
    Bu += b(k)*u[k];
    Cu += c(k)*u[k];
    Du += d(k)*u[k];
    Eu += e(k)*u[k];
    Fu += f(k)*u[k];
  }
  alphaBj =     b(j)*B1 +     c(j)*B2;
  betaBj  = 2.0*e(j)*B1 +     d(j)*B2;
  gammaBj =     d(j)*B1 + 2.0*f(j)*B2;

  axij =
  + 420.0 * (a(i)*Au                     ) * (3.0*alphaBj +     betaBj +     gammaBj)
  + 105.0 * (a(i)*Bu + b(i)*Au           ) * (4.0*alphaBj + 2.0*betaBj +     gammaBj)
  + 105.0 * (a(i)*Cu + c(i)*Au           ) * (4.0*alphaBj +     betaBj + 2.0*gammaBj)
  +  21.0 * (a(i)*Du + d(i)*Au           ) * (5.0*alphaBj + 2.0*betaBj + 2.0*gammaBj)
  +  21.0 * (b(i)*Cu + c(i)*Bu           ) * (5.0*alphaBj + 2.0*betaBj + 2.0*gammaBj)
  +  42.0 * (a(i)*Eu + e(i)*Au + b(i)*Bu ) * (5.0*alphaBj + 3.0*betaBj +     gammaBj)
  +  42.0 * (a(i)*Fu + f(i)*Au + c(i)*Cu ) * (5.0*alphaBj +     betaBj + 3.0*gammaBj)
  +   7.0 * (b(i)*Du + d(i)*Bu           ) * (6.0*alphaBj + 3.0*betaBj + 2.0*gammaBj)
  +   7.0 * (c(i)*Eu + e(i)*Cu           ) * (6.0*alphaBj + 3.0*betaBj + 2.0*gammaBj)
  +   7.0 * (b(i)*Fu + f(i)*Bu           ) * (6.0*alphaBj + 2.0*betaBj + 3.0*gammaBj)
  +   7.0 * (c(i)*Du + d(i)*Cu           ) * (6.0*alphaBj + 2.0*betaBj + 3.0*gammaBj)
  +  21.0 * (b(i)*Eu + e(i)*Bu           ) * (6.0*alphaBj + 4.0*betaBj +     gammaBj)
  +  21.0 * (c(i)*Fu + f(i)*Cu           ) * (6.0*alphaBj +     betaBj + 4.0*gammaBj)
  +   3.0 * (d(i)*Eu + e(i)*Du           ) * (7.0*alphaBj + 4.0*betaBj + 2.0*gammaBj)
  +   3.0 * (d(i)*Fu + f(i)*Du           ) * (7.0*alphaBj + 2.0*betaBj + 4.0*gammaBj)
  +   2.0 * (e(i)*Fu + f(i)*Eu + d(i)*Du ) * (7.0*alphaBj + 3.0*betaBj + 3.0*gammaBj)
  +  12.0 * (e(i)*Eu                     ) * (7.0*alphaBj + 5.0*betaBj +     gammaBj)
  +  12.0 * (f(i)*Fu                     ) * (7.0*alphaBj +     betaBj + 5.0*gammaBj)
  ;
  return axij/1260.0;
}
