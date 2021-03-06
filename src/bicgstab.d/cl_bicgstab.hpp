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
void   gp_copy(int n, double *y, double *x);
double gp_dot(int n, double *y, double *x);
void   gp_presolve_pointjacobi(int n, double *x, double *dinv, double *d);
void   gp_presolve(int n, double *x, double *dinv, double *d);

double gp_phase0(int n, double *r, 
		 double *Aa, int *col_ind, int *row_ptr,
		 double *x, double *rtilde, double *b, int w);
void gp_phase1(int n, double *p, double *r, double *v,
	       double beta, double omega);
double gp_phase2(int n, double *v,
		 double *Aa, int *col_ind, int *row_ptr,
		 double *phat, double *rtilde,int w);
double gp_phase3(int n, double *s, double *r, double *v, double alpha);
void   gp_phase4(int n, double *s, double *phat, double alpha);
double gp_phase5(int n, double *t, 
	       double *Aa, int *col_ind, int *row_ptr,
	       double *shat, double *s, int w);
double gp_phase6(int n, double *x, double *s, double *r, double *t,
	       double *phat, double *shat, double alpha, double omega);
double  gp_bicgstab(int n,int w, double*Aa,  int*col_ind,
		    int *row_ptr,  double *x,  double *b,
		    double *r,  double *p,  double *phat,
		    double *s,  double *shat,  double *t,
		    double *v,  double *rtilde,  double *dinv,
		    int max_iter, double tol);

int gp_send_A(int n,int w, double *Aa, int *col_ind, int *row_ptr);


static int sp_send_A(int n,int w,double *Aa, int *col_ind, int *row_ptr)
{
  return 0;
}


static double sp_norm(int n, double *x)
{
  double tmp=0.0;
  for (int k=1;k<n;k++) tmp += x[k]*x[k];
  return sqrt(tmp);
}


static double sp_dot(int n, double *p, double *q)
{
  double tmp=0.0;
  for (int k=1; k<n; k++ ) tmp += p[k]*q[k];
  return tmp;
}


static void sp_copy(int n, double *p, double *q)
{
  for (int k=1;k<n;k++) p[k]=q[k];
}


static void sp_presolve_pointjacobi(int n, double *x, double *dinv, double *d)
{
  for (int k=1; k<n; k++ ) x[k] = dinv[k]*d[k];
}


static void sp_presolve(int n, double *x, double *dinv, double *d)
{
  if ( dinv[1] == 0.0 ) {
    sp_copy(n,x,d);
  }
  else {
    sp_presolve_pointjacobi(n,x,dinv,d);
  }
}


static double sp_phase0(int n, double *r, 
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


static void sp_phase1(int n, double *p, double *r, double *v,
		   double beta, double omega)
{
  for(int k=1;k<n;k++) p[k] = r[k] + beta * (p[k] - omega *v[k]);
}

static double sp_phase2(int n, double *v,
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


static double sp_phase3(int n, double *s, double * r, double * v, double alpha)
{
  for(int k=1;k<n;k++) s[k] = r[k] - alpha * v[k];
  return sp_norm(n,s);
}


static void sp_phase4(int n, double *x, double *phat, double alpha)
{
  for(int k=1;k<n;k++) x[k] = x[k] + alpha*phat[k];
}


static double sp_phase5(int n, double *t, 
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

static double sp_phase6(int n, double *x, double *s, double *r, double *t,
		   double *phat, double *shat, double alpha, double omega)
{
  for(int k=1;k<n;k++) {
    x[k] = x[k] + alpha*phat[k] + omega*shat[k];
    r[k] = s[k] - omega * t[k];
  }
  return sp_norm(n,r);
}

extern size_t np;

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
  double (*norm)(int n, double *x);
  void (*copy)(int n, double *y, double *x);
  double (*dot)(int n, double *y, double *x);
  void   (*presolve_pointjacobi)(int n, double *x, double *dinv, double *d);
  void   (*presolve)(int n, double *x, double *dinv, double *d);

  double (*phase0)(int n, double *r,
                 double *Aa, int *col_ind, int *row_ptr,
                 double *x, double *rtilde, double *b, int w);
  void (*phase1)(int n, double *p, double *r, double *v,
               double beta, double omega);
  double (*phase2)(int n, double *v,
                 double *Aa, int *col_ind, int *row_ptr,
                 double *phat, double *rtilde,int w);
  double (*phase3)(int n, double *s, double *r, double *v, double alpha);
  void   (*phase4)(int n, double *s, double *phat, double alpha);

  double (*phase5)(int n, double *t,
               double *Aa, int *col_ind, int *row_ptr,
               double *shat, double *s, int w);
  double (*phase6)(int n, double *x, double *s, double *r, double *t,
               double *phat, double *shat, double alpha, double omega);
  double  (*bicgstab)(int n,int w, double*Aa,  int*col_ind,
                    int *row_ptr,  double *x,  double *b,
                    double *r,  double *p,  double *phat,
                    double *s,  double *shat,  double *t,
                    double *v,  double *rtilde,  double *dinv,
                    int max_iter, double tol);

  if ( np ) {
    norm = gp_norm;
    copy = gp_copy;
    dot = gp_dot;
    presolve_pointjacobi = gp_presolve_pointjacobi;
    presolve = gp_presolve;
    phase0 = gp_phase0;
    phase1 = gp_phase1;
    phase2 = gp_phase2;
    phase3 = gp_phase3;
    phase4 = gp_phase4;
    phase5 = gp_phase5;
    phase6 = gp_phase6;
    gp_send_A(n,w,Aa,col_ind,row_ptr);
  } else {
    norm = sp_norm;
    copy = sp_copy;
    dot = sp_dot;
    presolve_pointjacobi = sp_presolve_pointjacobi;
    presolve = sp_presolve;
    phase0 = sp_phase0;
    phase1 = sp_phase1;
    phase2 = sp_phase2;
    phase3 = sp_phase3;
    phase4 = sp_phase4;
    phase5 = sp_phase5;
    phase6 = sp_phase6;
  }


  
  double resid,rho_1,rho_2,alpha,beta,omega, normb = norm(n,b);
  if (normb == 0.0) normb = 1;

  if ((resid = phase0(n,r,Aa,col_ind,row_ptr,x,rtilde,b,w)/normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <=max_iter; i++) {
    rho_1 = dot(n,rtilde,r);

    if (rho_1 == 0) {
      tol = norm(n,r)/normb;
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
      tol = norm(n,r)/normb;
      return 3;
    }
  }
  tol = resid;
  return 1;
}
