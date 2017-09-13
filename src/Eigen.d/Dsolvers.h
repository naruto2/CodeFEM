//*****************************************************************
// Iterative template routine -- BiCG
//
// BiCG solves the unsymmetric linear system Ax = b 
// using the Preconditioned BiConjugate Gradient method
//
// BiCG follows the algorithm described on p. 22 of the 
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

template < class Matrix, class Vector, class Preconditioner, class Real >
int 
BiCG(const Matrix &A, Vector &x, const Vector &b,
     const Preconditioner &M, int &max_iter, Real &tol)
{
  Real resid;
  Vector rho_1(1), rho_2(1), alpha(1), beta(1);
  Vector z, ztilde, p, ptilde, q, qtilde;

  Real normb = norm(b);
  Vector r = b - A * x;
  Vector rtilde = r;

  if (normb == 0.0)
    normb = 1;
  
  if ((resid = norm(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++) {

    z = M.solve(r);
    ztilde = M.trans_solve(rtilde);
    rho_1[0] = dot(z, rtilde);
    if (rho_1[0] == 0) { 
      tol = norm(r) / normb;
      max_iter = i;
      return 2;
    }
    if (i == 1) {
      p = z;
      ptilde = ztilde;
    } else {
      beta[0] = rho_1[0] / rho_2[0];
      p = z + beta[0] * p;
      ptilde = ztilde + beta[0] * ptilde;
    }
    q = A * p;
    qtilde = A.transpose() * ptilde;
    alpha[0] = rho_1[0] / dot(ptilde, q);
    x = x + alpha[0] * p;
    r = r - alpha[0] * q;
    rtilde = rtilde - alpha[0] * qtilde;

    rho_2[0] = rho_1[0];
    if ((resid = norm(r) / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }

  tol = resid;
  return 1;
}

Vector bicg(Smatrix &A, Vector &b){
  Preconditioner M;
  Vector x = b;
  int max_iter = 1000000;
  double tol = 0.0000001;
  BiCG(A, x, b, M, max_iter, tol);
  return x;
}
  
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

template < class Matrix, class Vector, class Preconditioner, class Real >
int 
DBiCGSTAB(const Matrix &A, Vector &x, const Vector &b,
         const Preconditioner &M, int &max_iter, Real &tol)
{
  Real resid;
  Vector rho_1(1), rho_2(1), alpha(1), beta(1), omega(1);
  Vector p, phat, s, shat, t, v;

  Real normb = norm(b);
  Vector r = b - A * x;
  Vector rtilde = r;

  if (normb == 0.0)
    normb = 1;
  
  if ((resid = norm(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++) {
    rho_1[0] = dot(rtilde, r);
    if (rho_1[0] == 0) {
      tol = norm(r) / normb;
      return 2;
    }
    if (i == 1)
      p = r;
    else {
      beta[0] = (rho_1[0]/rho_2[0]) * (alpha[0]/omega[0]);
      p = r + beta[0] * (p - omega[0] * v);
    }
    phat = M.solve(p);
    v = A * phat;
    alpha[0] = rho_1[0] / dot(rtilde, v);
    s = r - alpha[0] * v;
    if ((resid = norm(s)/normb) < tol) {
      x = x + alpha[0] * phat;
      tol = resid;
      return 0;
    }
    shat = M.solve(s);
    t = A * shat;
    omega[0] = dot(t,s) / dot(t,t);


    Vector x1, x2;
    x1 = alpha[0] * phat;
    x2 = omega[0] * shat;
    x = x + x1;
    x = x + x2;
    //x = x + alpha[0] * phat + omega[0] * shat;


    r = s - omega[0] * t;
    rho_2[0] = rho_1[0];
    if ((resid = norm(r) / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
    if (omega[0] == 0) {
      tol = norm(r) / normb;
      return 3;
    }
  }

  tol = resid;
  return 1;
}

Vector bicgstab(Smatrix &A, Vector &b){
  Preconditioner M;
  Vector x = b;
  int max_iter = 1000000;
  double tol = 0.0000001;
  DBiCGSTAB(A, x, b, M, max_iter, tol);
  return x;
}

//*****************************************************************
// Iterative template routine -- CG
//
// CG solves the symmetric positive definite linear
// system Ax=b using the Conjugate Gradient method.
//
// CG follows the algorithm described on p. 15 in the 
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

template < class Matrix, class Vector, class Preconditioner, class Real >
int 
CG(const Matrix &A, Vector &x, const Vector &b,
   const Preconditioner &M, int &max_iter, Real &tol)
{
  Real resid;
  Vector p, z, q;
  Vector alpha(1), beta(1), rho(1), rho_1(1);

  Real normb = norm(b);

  Vector r = b - A*x;

  if (normb == 0.0) 
    normb = 1;
  
  if ((resid = norm(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }
  
  for (int i = 1; i <= max_iter; i++) {
    z = M.solve(r);
    //z = r;
    
    rho[0] = dot(r, z);
    
    if (i == 1)
      p = z;
    else {

      beta[0] = rho[0] / rho_1[0];

      p = z + beta[0] * p;

    }
    q = A*p;

    alpha[0] = rho[0] / dot(p, q);

    x = x+alpha[0] * p;

    r = r-alpha[0] * q;

    if ((resid = norm(r) / normb) <= tol) {
      tol = resid;
      max_iter = i;
      return 0;     
    }
    
    rho_1[0] = rho[0];
  }

  tol = resid;
  return 1;
}

Vector cg(Smatrix &A, Vector &b){
  Preconditioner M;
  Vector x = b;
  int max_iter = 1000000;
  double tol = 0.0000001;
  if ( IsSymmetric(A) )
    CG(A, x, b, M, max_iter, tol);
  return x;
}

//*****************************************************************
// Iterative template routine -- CGS
//
// CGS solves the unsymmetric linear system Ax = b 
// using the Conjugate Gradient Squared method
//
// CGS follows the algorithm described on p. 26 of the 
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

template < class Matrix, class Vector, class Preconditioner, class Real >
int 
CGS(const Matrix &A, Vector &x, const Vector &b,
    const Preconditioner &M, int &max_iter, Real &tol)
{
  Real resid;
  Vector rho_1(1), rho_2(1), alpha(1), beta(1);
  Vector p, phat, q, qhat, vhat, u, uhat;

  Real normb = norm(b);
  Vector r = b - A*x;
  Vector rtilde = r;

  if (normb == 0.0)
    normb = 1;
  
  if ((resid = norm(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++) {

    rho_1[0] = dot(rtilde, r);
    if (rho_1[0] == 0) {
      tol = norm(r) / normb;
      return 2;
    }
    if (i == 1) {
      u = r;
      p = u;
    } else {
      beta[0] = rho_1[0] / rho_2[0];
      u = r + beta[0] * q;
      p = u + beta[0] * (q + beta[0] * p);
    }
    phat = M.solve(p);
    vhat = A*phat;
    alpha[0] = rho_1[0] / dot(rtilde, vhat);
    q = u - alpha[0] * vhat;
    uhat = M.solve(u + q);
    x = x + alpha[0] * uhat;
    qhat = A * uhat;
    r = r - alpha[0] * qhat;
    rho_2[0] = rho_1[0];
    if ((resid = norm(r) / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }
  
  tol = resid;
  return 1;
}

Vector cgs(Smatrix &A, Vector &b){
  Preconditioner M;
  Vector x = b;
  int max_iter = 1000000;
  double tol = 0.0000001;
  CGS(A, x, b, M, max_iter, tol);
  return x;
}

//*****************************************************************
// Iterative template routine -- QMR
//
// QMR.h solves the unsymmetric linear system Ax = b using the
// Quasi-Minimal Residual method following the algorithm as described
// on p. 24 in the SIAM Templates book.
//
//   -------------------------------------------------------------
//   return value     indicates
//   ------------     ---------------------
//        0           convergence within max_iter iterations
//        1           no convergence after max_iter iterations
//                    breakdown in:
//        2             rho
//        3             beta
//        4             gamma
//        5             delta
//        6             ep
//        7             xi
//   -------------------------------------------------------------
//   
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax=b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//*****************************************************************


#include <math.h>

template < class Matrix, class Vector, class Preconditioner1,
           class Preconditioner2, class Real >
int 
QMR(const Matrix &A, Vector &x, const Vector &b, const Preconditioner1 &M1, 
    const Preconditioner2 &M2, int &max_iter, Real &tol)
{
  Real resid;

  Vector rho(1), rho_1(1), xi(1), gamma(1), gamma_1(1), theta(1), theta_1(1);
  Vector eta(1), delta(1), ep(1), beta(1);

  Vector r, v_tld, y, w_tld, z;
  Vector v, w, y_tld, z_tld;
  Vector p, q, p_tld, d, s;

  Real normb = norm(b);

  r = b - A * x;

  if (normb == 0.0)
    normb = 1;

  if ((resid = norm(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  v_tld = r;
  y = M1.solve(v_tld);
  rho[0] = norm(y);

  w_tld = r;
  z = M2.trans_solve(w_tld);
  xi[0] = norm(z);

  gamma[0] = 1.0;
  eta[0] = -1.0;
  theta[0] = 0.0;

  for (int i = 1; i <= max_iter; i++) {

    if (rho[0] == 0.0)
      return 2;                        // return on breakdown

    if (xi[0] == 0.0)
      return 7;                        // return on breakdown

    v = (1. / rho[0]) * v_tld;
    y = (1. / rho[0]) * y;

    w = (1. / xi[0]) * w_tld;
    z = (1. / xi[0]) * z;

    delta[0] = dot(z, y);
    if (delta[0] == 0.0)
      return 5;                        // return on breakdown

    y_tld = M2.solve(y);               // apply preconditioners
    z_tld = M1.trans_solve(z);

    if (i > 1) {
      p = y_tld - (xi[0] * delta[0] / ep[0]) * p;
      q = z_tld - (rho[0] * delta[0] / ep[0]) * q;
    } else {
      p = y_tld;
      q = z_tld;
    }

    p_tld = A * p;
    ep[0] = dot(q, p_tld);
    if (ep[0] == 0.0)
      return 6;                        // return on breakdown

    beta[0] = ep[0] / delta[0];
    if (beta[0] == 0.0)
      return 3;                        // return on breakdown

    v_tld = p_tld - beta[0] * v;
    y = M1.solve(v_tld);

    rho_1[0] = rho[0];
    rho[0] = norm(y);
    w_tld = A.transpose() * q - beta[0] * w;
    z = M2.trans_solve(w_tld);

    xi[0] = norm(z);

    gamma_1[0] = gamma[0];
    theta_1[0] = theta[0];

    theta[0] = rho[0] / (gamma_1[0] * beta[0]);
    gamma[0] = 1.0 / sqrt(1.0 + theta[0] * theta[0]);

    if (gamma[0] == 0.0)
      return 4;                        // return on breakdown

    eta[0] = -eta[0] * rho_1[0] * gamma[0] * gamma[0] / 
      (beta[0] * gamma_1[0] * gamma_1[0]);

    if (i > 1) {
      Vector d1, d2, s1, s2;
      d1 = eta[0] * p;
      d2 = (theta_1[0] * theta_1[0] * gamma[0] * gamma[0]) * d;
      d = d1 + d2;
      //d = eta[0] * p + (theta_1[0] * theta_1[0] * gamma[0] * gamma[0]) * d;
      s1 = eta[0] * p_tld;
      s2 = (theta_1[0] * theta_1[0] * gamma[0] * gamma[0]) * s;
      s = s1 + s2;
      //s = eta[0] * p_tld + (theta_1[0] * theta_1[0] * gamma[0] * gamma[0]) * s;
    } else {
      d = eta[0] * p;
      s = eta[0] * p_tld;
    }

    x = x + d;                            // update approximation vector
    r = r - s;                            // compute residual

    if ((resid = norm(r) / normb) <= tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }

  tol = resid;
  return 1;                            // no convergence
}

Vector qmr(Smatrix &A, Vector &b){
  Preconditioner M;
  Vector x = b;
  int max_iter = 1000000;
  double tol = 0.0000001;
  QMR(A, x, b, M, M, max_iter, tol);
  return x;
}

// ogita 荻田氏の方法(対称行列化)
Vector ogita(Smatrix &A, Vector &b) {

  Smatrix AT = A.transpose();

  Smatrix AA = vcatSmatrix( catSmatrix(A+AT, A-AT),
			    catSmatrix(AT-A,-A-AT));

  Vector  bb = catVector(2*b, -2*b);
  
  return halfVector(cg(AA,bb));
}
