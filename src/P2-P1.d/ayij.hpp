static double ayij(long i, long j, double *v, double C1, double C2, double C3)
{
  double Av=0.0, Bv=0.0, Cv=0.0, Dv=0.0, Ev=0.0, Fv=0.0,
    alphaCj, betaCj, gammaCj, ayij;

  for (long k = 1; k <= 6; k++) {
    Av += a(k)*v[k];
    Bv += b(k)*v[k];
    Cv += c(k)*v[k];
    Dv += d(k)*v[k];
    Ev += e(k)*v[k];
    Fv += f(k)*v[k];
  }
  alphaCj =     b(j)*C1 +     c(j)*C2;
  betaCj  = 2.0*e(j)*C1 +     d(j)*C2;
  gammaCj =     d(j)*C1 + 2.0*f(j)*C2;

  ayij =
    420.0 * (a(i)*Av                    ) * (3.0*alphaCj +     betaCj +     gammaCj)
  + 105.0 * (a(i)*Bv + b(i)*Av          ) * (4.0*alphaCj + 2.0*betaCj +     gammaCj)
  + 105.0 * (a(i)*Cv + c(i)*Av          ) * (4.0*alphaCj +     betaCj + 2.0*gammaCj)
  +  21.0 * (a(i)*Dv + d(i)*Av          ) * (5.0*alphaCj + 2.0*betaCj + 2.0*gammaCj)
  +  21.0 * (b(i)*Cv + c(i)*Bv          ) * (5.0*alphaCj + 2.0*betaCj + 2.0*gammaCj)
  +  42.0 * (a(i)*Ev + e(i)*Av + b(i)*Bv) * (5.0*alphaCj + 3.0*betaCj +     gammaCj)
  +  42.0 * (a(i)*Fv + f(i)*Av + c(i)*Cv) * (5.0*alphaCj +     betaCj + 3.0*gammaCj)
  +   7.0 * (b(i)*Dv + d(i)*Bv          ) * (6.0*alphaCj + 3.0*betaCj + 2.0*gammaCj)
  +   7.0 * (c(i)*Ev + e(i)*Cv          ) * (6.0*alphaCj + 3.0*betaCj + 2.0*gammaCj)
  +   7.0 * (b(i)*Fv + f(i)*Bv          ) * (6.0*alphaCj + 2.0*betaCj + 3.0*gammaCj)
  +   7.0 * (c(i)*Dv + d(i)*Cv          ) * (6.0*alphaCj + 2.0*betaCj + 3.0*gammaCj)
  +  21.0 * (b(i)*Ev + e(i)*Bv          ) * (6.0*alphaCj + 4.0*betaCj +     gammaCj)
  +  21.0 * (c(i)*Fv + f(i)*Cv          ) * (6.0*alphaCj +     betaCj + 4.0*gammaCj)
  +   3.0 * (d(i)*Ev + e(i)*Dv          ) * (7.0*alphaCj + 4.0*betaCj + 2.0*gammaCj)
  +   3.0 * (d(i)*Fv + f(i)*Dv          ) * (7.0*alphaCj + 2.0*betaCj + 4.0*gammaCj)
  +   2.0 * (e(i)*Fv + f(i)*Ev + d(i)*Dv) * (7.0*alphaCj + 3.0*betaCj + 3.0*gammaCj)
  +  12.0 * (e(i)*Ev                    ) * (7.0*alphaCj + 5.0*betaCj +     gammaCj)
  +  12.0 * (f(i)*Fv                    ) * (7.0*alphaCj +     betaCj + 5.0*gammaCj)
  ;
  return ayij/1260.0;
}
