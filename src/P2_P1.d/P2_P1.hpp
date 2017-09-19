static double B1=0.0, B2=0.0, C1=0.0, C2=0.0;
static double Delta=1.0;

static double a(long i){
  switch(i) {
  case 1: return  0.0;
  case 2: return  0.0;
  case 3: return  1.0;
  case 4: return  0.0;
  case 5: return  0.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double b(long i){
  switch(i) {
  case 1: return -1.0;
  case 2: return  0.0;
  case 3: return -3.0;
  case 4: return  0.0;
  case 5: return  4.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double c(long i){
  switch(i) {
  case 1: return  0.0;
  case 2: return -1.0;
  case 3: return -3.0;
  case 4: return  4.0;
  case 5: return  0.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double d(long i){
  switch(i) {
  case 1: return  0.0;
  case 2: return  0.0;
  case 3: return  4.0;
  case 4: return -4.0;
  case 5: return -4.0;
  case 6: return  4.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double e(long i){
  switch(i) {
  case 1: return  2.0;
  case 2: return  0.0;
  case 3: return  2.0;
  case 4: return  0.0;
  case 5: return -4.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double f(long i){
  switch(i) {
  case 1: return  0.0;
  case 2: return  2.0;
  case 3: return  2.0;
  case 4: return -4.0;
  case 5: return  0.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double alphaB(long j)
{
  switch(j){
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: return b(j)*B1 + c(j)*B2;
  default: abort();
  }
  abort();
  return NAN;
}


static double betaB(long j)
{
  switch(j){
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: return 2.0*e(j)*B1 + d(j)*B2;
  default: abort();
  }
  abort();
  return NAN;
}


static double gammaB(long j)
{
  switch(j){
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: return d(j)*B1 + 2.0*f(j)*B2;
  default: abort();
  }
  abort();
  return NAN;
}


static double alphaC(long j)
{
  switch(j){
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: return b(j)*C1 + c(j)*C2;
  default: abort();
  }
  abort();
  return NAN;
}


static double betaC(long j)
{
  switch(j){
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: return 2.0*e(j)*C1 + d(j)*C2;
  default: abort();
  }
  abort();
  return NAN;
}


static double gammaC(long j)
{
  switch(j){
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: return d(j)*C1 + 2.0*f(j)*C2;
  default: abort();
  }
  abort();
  return NAN;
}


static double ad(long j){
  switch(j) {
  case 1: return  0.0;
  case 2: return  0.0;
  case 3: return  1.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double bd(long j){
  switch(j) {
  case 1: return  1.0;
  case 2: return  0.0;
  case 3: return -1.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double cd(long j){
  switch(j) {
  case 1: return  0.0;
  case 2: return  1.0;
  case 3: return -1.0;
  default: abort();
  }
  abort();
  return NAN;
}
