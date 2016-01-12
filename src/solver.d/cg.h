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

  printf("normb = %f\n",norm(b));
  Vector r = b - A*x;
  printf("b = %f %f\n",b[0],b[1]);
  printf("x = %f %f\n",x[0],x[1]);
  printf("r = b - A*x\n");
  printf("r = %f %f\n",r[0],r[1]);
  
  if (normb == 0.0) 
    normb = 1;
  
  if ((resid = norm(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  printf("------------------------------------------------------------\n");
  
  for (int i = 1; i <= max_iter; i++) {

    //z = M.solve(r);
    z = r;
    printf("r = %f %f\n",r[0],r[1]);
    printf("z = %f %f\n",z[0],z[1]);
    
    rho[0] = dot(r, z);

    printf("rho[0]=%f\n",rho[0]);
    
    
    if (i == 1)
      p = z;
    else {

      beta[0] = rho[0] / rho_1[0];

      printf("beta[0]=%f\n",beta[0]);
      
      p = z + beta[0] * p;
      printf("p = %f %f\n",p[0],p[1]);
    }
    
    printf("p = %f %f\n",p[0],p[1]);
    printf("A = %f %f\n %f %f\n",A[0][0],A[0][1],A[1][0],A[1][1]);

    q = A*p;
    printf("q = A*p\n");
    printf("q = %f %f\n",q[0],q[1]);

    alpha[0] = rho[0] / dot(p, q);
    
    printf("alpha[0] = %f\n",alpha[0]);
    

    x = x+alpha[0] * p;
    printf("x = %f %f\n",x[0],x[1]);


    printf("##############################################################\n");
    printf("r = %f %f\n",r[0],r[1]);
    printf("alpha[0] = %f\n",alpha[0]);
    printf("q = %f %f\n",q[0],q[1]);
    printf("r = r-alpha[0] * q\n");
    r = r-alpha[0] * q;
    printf("r = %f %f\n",r[0],r[1]);


    
    
    if ((resid = norm(r) / normb) <= tol) {
      tol = resid;
      max_iter = i;
      return 0;     
    }

    printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    
    rho_1[0] = rho[0];
    printf("rho_1[0] = %f\n",rho_1[0]);

  }

  tol = resid;
  return 1;
}

