//*****************************************************************
// Iterative template routine -- BiCGSTAB
//
// BiCGSTAB solves the unsymmetric linear system Ax = b 
// using the Preconditioned BiConjugate Gradient Stabilized method
//
// BiCGSTAB follows the algorithm described on p. 27 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************
double gp_norm(int n, double *x);
double gp_dot(int n, double *y, double *x);

double cl_norm(int n, double *x);
void   cl_copy(int n, double *y, double *x);
void   cl_init(int argc, char **argv);
double cl_dot(int n, double *y, double *x);
void cl_matrixvector(int n, double *r, double *Aa, int *col_ind,
		     int *row_ptr, double *b, int w);
double cl_phase0(int n, double *r, 
		 double *Aa, int *col_ind, int *row_ptr,
		 double *x, double *rtilde, double *b, int w);
void cl_phase1(int n, double *p, double *r, double *v,
	       double beta, double omega);
double cl_phase2(int n, double *v,
		 double *Aa, int *col_ind, int *row_ptr,
		 double *phat, double *rtilde,int w);
double cl_phase3(int n, double *s, double *r, double *v, double alpha);
void cl_phase4(int n, double *s, double *phat, double alpha);
double  cl_phase5(int n, double *t, 
	       double *Aa, int *col_ind, int *row_ptr,
	       double *shat, double *s, int w);
double cl_phase6(int n, double *x, double *s, double *r, double *t,
	       double *phat, double *shat, double alpha, double omega);
int cl_send_A(int n,int w, double *Aa, int *col_ind, int *row_ptr);

static double norm(int n, double *x)
{
  double tmp=0.0;
  for (int k=1;k<n;k++) tmp += x[k]*x[k];
  return sqrt(tmp);
}


static void copy(int n, double *p, double *q)
{
  for (int k=1;k<n;k++) p[k]=q[k];
}


static void presolve(int n, double *x, double *dinv, double *d)
{
  if ( dinv[1] == 0.0 ) {
    copy(n,x,d);
  }
  else {
    for (int k=1; k<n; k++ ) x[k] = dinv[k]*d[k];
  }
}

static double dot(int n, double *p, double *q)
{
  double tmp=0.0;
  for (int k=1; k<n; k++ ) tmp += p[k]*q[k];
  return tmp;
}

static double phase0(int n, double *r, 
		     double *Aa, int *col_ind, int *row_ptr,
		     double *x, double *rtilde, double *b, int w)
{
  double tmp=0.0;
  for (int k=1;k<n;k++ ) {
    r[k] = 0.0;
    for (int j=row_ptr[k];j<row_ptr[k+1];j++) r[k] += Aa[j] * x[col_ind[j]];
    rtilde[k] = r[k] = b[k] - r[k];
    tmp += r[k]*r[k];
  }
  return sqrt(tmp);
}


static void phase1(int n, double *p, double *r, double *v,
		   double beta, double omega)
{
  for(int k=1;k<n;k++) p[k] = r[k] + beta * (p[k] - omega *v[k]);
}

static double phase2(int n, double *v,
		     double *Aa, int *col_ind, int *row_ptr,
		     double *phat, double *rtilde, int w)
{
  double tmp=0.0;

  for (int k=1;k<n;k++ ) {
    v[k] = 0.0;
    for (int j=row_ptr[k];j<row_ptr[k+1];j++){
      v[k] += Aa[j] * phat[col_ind[j]];
    }
    tmp += rtilde[k]*v[k];
  }
  return tmp;
}


static double phase3(int n, double *s, double * r, double * v, double alpha)
{
  for(int k=1;k<n;k++) s[k] = r[k] - alpha * v[k];
  return norm(n,s);
}


static void phase4(int n, double *x, double *phat, double alpha)
{
  for(int k=1;k<n;k++) x[k] = x[k] + alpha*phat[k];
}


static double phase5(int n, double *t, 
		     double *Aa, int *col_ind, int *row_ptr,
		     double *shat, double *s, int w)
{
  double tmp=0.0, tmq=0.0;
  for (int k=1;k<n;k++ ) {
    t[k] = 0.0;
    for (int j=row_ptr[k];j<row_ptr[k+1];j++)
      t[k] += Aa[j] * shat[col_ind[j]];
    tmp += t[k]*s[k];
    tmq += t[k]*t[k];
  }
  return tmp/tmq;
}

static double phase6(int n, double *x, double *s, double *r, double *t,
		   double *phat, double *shat, double alpha, double omega)
{
  for(int k=1;k<n;k++) {
    x[k] = x[k] + alpha*phat[k] + omega*shat[k];
    r[k] = s[k] - omega * t[k];
  }
  return norm(n,r);
}
    

int sparse__BiCGSTAB(const sparse::matrix<double> &A, double *x, double *b,
		     int &max_iter, double &tol, double *dinv)
	
{
  static int nn, k, ww, w, *row_ptr, *col_ind, n = A.size();
  static double *r,*p,*phat,*s,*shat,*t,*v,*rtilde, *Aa;

  if ( nn < n ) { free(r); free(row_ptr); r=0; }
  if (!r) {
    r = (double*)malloc(8*n*sizeof(double));
    p = &r[n]; s = &r[2*n]; t = &r[3*n]; v = &r[4*n];
    phat = &r[5*n]; shat = &r[6*n]; rtilde = &r[7*n];
    row_ptr = (int*)malloc((n+1)*sizeof(int));
    nn = n;
  }

  w = 1;
  for ( k=1;k<n;k++ ) {
    row_ptr[k] = w;
    for (auto j : A[k] ) if(j.second != 0.0) w++;
  }
  row_ptr[n] = w;

  if( ww < w ) { free(Aa); free(col_ind); Aa=0; }
  if (!Aa){
    Aa      = (double*)malloc((w+1)*sizeof(double));
    col_ind = (int   *)malloc((w+1)*sizeof(int));
    ww = w;
  }
  int ii=1;
  for ( k=1;k<n;k++ ) 
    for (auto j : A[k] ) if(j.second != 0.0) {
	Aa[ii] = j.second;
	col_ind[ii] = j.first;
	ii++;
      }
  cl_send_A(n,w,Aa,col_ind,row_ptr);
  double resid,rho_1,rho_2,alpha,beta,omega, normb = gp_norm(n,b);
  if (normb == 0.0) normb = 1;

  if ((resid = phase0(n,r,Aa,col_ind,row_ptr,x,rtilde,b,w)/normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }
  for (int i = 1; i <= max_iter; i++) {
    rho_1 = gp_dot(n,rtilde,r);
    if (rho_1 == 0) {
      tol = gp_norm(n,r)/normb;
      return 2;
    }
    if (i == 1)
      copy(n,p,r);
    else {
      beta = (rho_1/rho_2) * (alpha/omega);
      phase1(n,p,r,v,beta,omega);
    }
    presolve(n,phat,dinv,p);
    alpha = rho_1/phase2(n,v,Aa,col_ind,row_ptr,phat,rtilde,w);
    
    if ((resid = phase3(n,s,r,v,alpha)/normb) < tol) {
      phase4(n,x,phat,alpha);
      tol = resid;
      return 0;
    }
    presolve(n,shat,dinv,s);
    omega = phase5(n,t,Aa,col_ind,row_ptr,shat,s,w);
    rho_2 = rho_1;
    if ((resid = phase6(n,x,s,r,t,phat,shat,alpha,omega)/normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
    if (omega == 0) {
      tol = gp_norm(n,r)/normb;
      return 3;
    }
  }
  tol = resid;
  return 1;
}
