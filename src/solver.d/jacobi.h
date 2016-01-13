//*****************************************************************
// Iterative template routine -- Jacobi
//
// Jacobi solves the symmetric positive definite linear
// system Ax=b using the Jacobi method.
//
// Jacobi follows the algorithm described on p. 12 in the 
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
Jacobi(const Matrix &A, Vector &x, const Vector &b,
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



  int n = A.size();
  
  for (int k = 1; k <= max_iter; k++) {


    for ( int i=0; i<n; i++ ) {
      internal_map Ai = A[i];
      internal_map::iterator j = Ai.begin();
      for ( x[i]=0.0; j != Ai.end(); j++ ) if ( i != j->first ) {
	  x[i] += j->second * r[j->first];
	}
      x[i] = (b[i]-x[i])/A[i][i];
    }
    
    r = x;    

    if ((resid = norm(b - A*x) / normb) <= tol) {
      tol = resid;
      max_iter = k;
      return 0;     
    }
    printf("k = %d    resid = %f\n",k,resid);
  }
  tol = resid;
  return 1;
}

