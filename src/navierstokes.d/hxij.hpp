static double hxij(long i, long j)
{
  return (1.0/12.0) * (12.0*(alphaB(i)*ad(j))
                         +4.0*( betaB(i)*ad(j) + alphaB(i)*bd(j) + gammaB(i)*ad(j) + alphaB(i)*cd(j))
                         +1.0*( betaB(i)*cd(j) + gammaB(i)*bd(j))
                         +2.0*( betaB(i)*bd(j) + gammaB(i)*cd(j))                                     );
}


#if 0
static double hxij(long i, long j)
{
  return (1.0/12.0) * (
  12.0*(alphaB(i)*ad(j))
  +4.0*( betaB(i)*ad(j) + alphaB(i)*bd(j) + gammaB(i)*ad(j) + alphaB(i)*cd(j))
  +1.0*( betaB(i)*cd(j) + gammaB(i)*bd(j))
  +2.0*( betaB(i)*bd(j) + gammaB(i)*cd(j))
  );
}
#endif 
