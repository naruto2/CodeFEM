static double dij(long i, long j)
{
  return (Delta/12.0) * (
  12.0  *(alphaB(i)*alphaB(j) + alphaC(i)*alphaC(j))
  + 4.0 *(alphaB(i)* betaB(j) +  betaB(i)*alphaB(j))
  + 4.0 *(alphaB(i)*gammaB(j) + gammaB(i)*alphaB(j))
  + 4.0 *(alphaC(i)* betaC(j) +  betaC(i)*alphaC(j))
  + 4.0 *(alphaC(i)*gammaC(j) + gammaC(i)*alphaC(j))
  + 1.0 *( betaB(i)*gammaB(j) + gammaB(i)* betaB(j))
  + 1.0 *( betaC(i)*gammaC(j) + gammaC(i)* betaC(j))
  + 2.0 *( betaB(i)* betaB(j) + gammaB(i)*gammaB(j))
  + 2.0 *( betaC(i)* betaC(j) + gammaC(i)*gammaC(j))
  );
}
