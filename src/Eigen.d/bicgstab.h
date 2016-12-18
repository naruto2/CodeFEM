// Preconditioner M; 前処理クラス
class Preconditioner{
 public:

  Vector &solve(Vector &p) const{
    static Vector q;
    q = p;
    return q;
  }
};

// norm(b); L2ノルム
double norm(const Vector &b) {
  return sqrt(b.dot(b));
}

// dot(a,b); ベクトルの内積
double dot(const Vector &a, const Vector &b) {
  return a.dot(b);
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

/****************************************************************************/
Vector bicgstab(Smatrix &A, Vector &b){
  Preconditioner M;
  Vector x = b;
  int max_iter = 1000000;
  double tol = 0.0000001;
  DBiCGSTAB(A, x, b, M, max_iter, tol);
  return x;
}

