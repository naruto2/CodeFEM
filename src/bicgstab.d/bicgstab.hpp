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

int 
sparse__BiCGSTAB(const sparse::matrix<double> &A, vector<double> &x,
		 const vector<double> &b,
		 int &max_iter, double &tol)
{
  int k, n = A.size();
  double tmp;
  tmp = 0.0;
  for (k=1;k<n;k++) tmp += b[k]*b[k];
  double tmq,resid,rho_1,rho_2,alpha,beta,omega, normb = sqrt(tmp);
  vector<double> r(n), p(n), phat(n), s(n), shat(n), t(n), v(n), rtilde(n);
  for (k=1;k<n;k++ ) {
    r[k] = 0.0;
    for (auto j : A[k] ) r[k] += j.second * x[j.first];
  }
  for(k=1;k<n;k++) rtilde[k] = r[k] = b[k] - r[k];

  if (normb == 0.0) normb = 1;
  tmp=0.0;
  for (k=1;k<n;k++) tmp += r[k]*r[k];  
  if ((resid = sqrt(tmp) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++) {

    rho_1=0.0;
    for (k=1; k<n; k++ ) rho_1 += rtilde[k]*r[k];

    if (rho_1 == 0) {
      tmp=0.0;
      for (k=1;k<n;k++) tmp += r[k]*r[k];
      tol = sqrt(tmp) / normb;
      return 2;
    }
    if (i == 1)
      for(k=1;k<n;k++) p[k]=r[k];
    else {
      beta = (rho_1/rho_2) * (alpha/omega);
      for(k=1;k<n;k++) p[k] = r[k] + beta * (p[k] - omega *v[k]);
    }
    for(k=1;k<n;k++) phat[k] = p[k];

    for (k=1;k<n;k++ ) {
      v[k] = 0.0;
      for (auto j : A[k] ) v[k] += j.second * phat[j.first];
    }

    tmp=0.0;
    for (k=1;k<n;k++) tmp += rtilde[k]*v[k];
    alpha = rho_1 / tmp;

    for(k=1;k<n;k++) s[k] = r[k] - alpha * v[k];

    tmp=0.0;
    for (k=1;k<n;k++) tmp += s[k]*s[k];
    if ((resid = sqrt(tmp)/normb) < tol) {
      for(k=1;k<n;k++) x[k] = x[k] + alpha*phat[k];
      tol = resid;
      return 0;
    }
    for(k=1;k<n;k++) shat[k] = s[k];

    for (k=1;k<n;k++ ) {
      t[k] = 0.0;
      for (auto j : A[k] ) t[k] += j.second * shat[j.first];
    }
    tmp=0.0;
    for (k=1;k<n;k++) tmp += t[k]*s[k];
    tmq=0.0;
    for (k=1;k<n;k++) tmq += t[k]*t[k];
    omega = tmp / tmq;

    for(k=1;k<n;k++) {
      x[k] = x[k] + alpha*phat[k] + omega*shat[k];
      r[k] = s[k] - omega * t[k];
    }
    rho_2 = rho_1;

    tmp=0.0;
    for (k=1;k<n;k++) tmp += r[k]*r[k];
    if ((resid = sqrt(tmp) / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
    if (omega == 0) {
      tmp=0.0;
      for (k=1;k<n;k++) tmp += r[k]*r[k];
      tol = sqrt(tmp) / normb;
      return 3;
    }
  }

  tol = resid;
  return 1;
}
