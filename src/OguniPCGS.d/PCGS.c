/* PCGS.f -- translated by f2c (version 12.02.01).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdlib.h> /* For exit() */
#include <f2c.h>

void do_lio(integer* , integer *, char *, integer);

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;


typedef struct {
  doublereal *d__;
  doublereal   *a;
  integer     *ia;
  integer      *m;
  integer      *n;
  integer     *nl;
  integer  a_dim1;
  integer ia_dim1;
  doublereal  *dd;
} Amatrix;


static void cp(doublereal *r, doublereal *x, integer n)
{
  integer i;

  for (i = 1; i <= n; i++)
    x[i] = r[i];
}


static doublereal *incompleteMLUdecomposition(doublereal *dd, doublereal  *s, 
					      doublereal  th, Amatrix A)
{
  integer i__, i__1, i__2, i__3, j, k, nn;
  doublereal ss, sw;

  doublereal *d__ = A.d__;
  doublereal   *a = A.a;
  integer     *ia = A.ia;
  integer      *m = A.m;
  integer      *n = A.n;
  integer     *nl = A.nl;
  integer  a_dim1 = A.a_dim1;
  integer ia_dim1 = A.ia_dim1;


  dd[1] = 1. / (*s * d__[1]);
  i__1 = *n;
  for (i__ = 2; i__ <= i__1; i__++) {
    ss = *s * d__[i__];
    i__2 = m[i__];
    for (k = 1; k <= i__2; ++k) {
      sw = 0.;
      nn = ia[i__ + k * ia_dim1];
      if (nn != 0) {
	i__3 = *nl + m[nn + *n];
	for (j = *nl + 1; j <= i__3; ++j) {
	  if (ia[nn + j * ia_dim1] != i__) {
	    sw += a[nn + j * a_dim1] * th;
	  } else {
	    sw += a[nn + j * a_dim1];
	  }
	}
      }
      ss -= a[i__ + k * a_dim1] * sw * dd[nn];
    }
    dd[i__] = 1. / ss;
  }
  return dd;
}


static doublereal *incompleteLUdecomposition(doublereal *dd, Amatrix A)
{
  integer i__, i__1, i__2, i__3, j, k, nn;
  doublereal ss, sw;

  doublereal *d__ = A.d__;
  doublereal   *a = A.a;
  integer     *ia = A.ia;
  integer      *m = A.m;
  integer      *n = A.n;
  integer     *nl = A.nl;
  integer  a_dim1 = A.a_dim1;
  integer ia_dim1 = A.ia_dim1;

  
  dd[1] = 1. / d__[1];
      i__1 = *n;
      for (i__ = 2; i__ <= i__1; ++i__) {
	ss = d__[i__];
	i__2 = m[i__];
	for (k = 1; k <= i__2; ++k) {
	  nn = ia[i__ + k * ia_dim1];
	  i__3 = *nl + m[nn + *n];
	  for (j = *nl + 1; j <= i__3; ++j) {
	    if (ia[nn + j * ia_dim1] == i__) {
	      ss -= a[i__ + k * a_dim1] * a[nn + j * a_dim1] * dd[nn];
	    }
	  }
	}
	dd[i__] = 1. / ss;
      }
      return dd;
}


static void multiply(doublereal *q,  Amatrix A, doublereal *x)
{    
  integer i, j, k;
  
  doublereal *d__ = A.d__;
  doublereal   *a = A.a;
  integer     *ia = A.ia;
  integer      *m = A.m;
  integer      *n = A.n;
  integer     *nl = A.nl;
  integer  a_dim1 = A.a_dim1;
  integer ia_dim1 = A.ia_dim1;
  
  for (i = 1; i <= *n; i++) {
    q[i] = d__[i] * x[i];
    k = m[i];
    for (j = 1; j <= k; j++) {
      q[i] += a[i + j * a_dim1] * x[ia[i + j * ia_dim1]];
    }
    k = *nl + m[i + *n];
    for (j = *nl + 1; j <= k; j++) {
      q[i] += a[i + j * a_dim1] * x[ia[i + j * ia_dim1]];
    }
  }
}


static doublereal dot(doublereal *p, doublereal *q, integer n)
{
  doublereal c = 0.0;
  integer i;

  for (i = 1; i <= n; i++)
    c += p[i] * q[i];

  return c;
}


static void presolve(doublereal *r__, Amatrix A, doublereal *r__1)
{    

  integer i__, i__1, i__2, j;
  doublereal sw;
  
  doublereal *d__ = A.d__;
  doublereal   *a = A.a;
  integer     *ia = A.ia;
  integer      *m = A.m;
  integer      *n = A.n;
  integer     *nl = A.nl;
  integer  a_dim1 = A.a_dim1;
  integer ia_dim1 = A.ia_dim1;
  doublereal *dd  = A.dd;
  
  i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = m[i__];
	for (j = 1; j <= i__2; ++j) {
	    r__[i__] -= a[i__ + j * a_dim1] * r__1[ia[i__ + j * ia_dim1]];
	}
	r__[i__] *= dd[i__];
    }

    for (i__ = *n; i__ >= 1; --i__) {
	sw = 0.;
	i__1 = *nl + m[i__ + *n];
	for (j = *nl + 1; j <= i__1; ++j) {
	    sw += a[i__ + j * a_dim1] * r__1[ia[i__ + j * ia_dim1]];
	}
	r__[i__] -= dd[i__] * sw;
    }
}


static void alpha_calculation(doublereal *h, doublereal *e,
			      doublereal alpha, doublereal *q, integer n)
{
  integer i;

  for (i = 1; i <= n; i++)
    h[i] = e[i] + alpha * q[i];
}


/* Subroutine */ int pcgs_(doublereal *d__, doublereal *a, integer *ia, 
	integer *n, integer *n1, integer *nl, doublereal *b, doublereal *eps, 
	integer *itr, doublereal *s, doublereal *x, doublereal *dd, 
	doublereal *p, doublereal *q, doublereal *r__, doublereal *r0, 
	doublereal *e, doublereal *h__, doublereal *w, integer *m, integer *
	ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, ia_dim1, ia_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal y, c1, c2, c3, x1, x2, th;
    static integer nn;
    static doublereal ss, sw, res, beta, alpha;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };


/* ********************************************************************** */
/*  CONJUGATE GRADIENT SQURED METHOD WITH INCOMPLETE LU DECOMPOSITION. * */
/*                                                                     * */
/*  PARAMETERS: SAME AS ROUTINE PCG EXCEPT A AND IA.                   * */
/*   ON ENTRY:                                                         * */
/*     A      NON-ZERO ELEMENTS OF THE LOWER TRIANGULAR PART OF THE    * */
/*            MATRIX A INTO 1-ST TO NL-TH POSITION OF THE ARRAY A.     * */
/*            NON-ZERO ELEMENTS OF THE UPPER TRIANGULAR PART OF THE    * */
/*            MATRIX A INTO NL+1-ST TO 2*NL-TH POSITION OF THE ARRAY A.* */
/*     IA     COLUMN INDEX OF CORRESPONDING ELEMENT IN THE ARRAY A.    * */
/*   OTHERS:  WORKING PARAMETERS.                                      * */
/*                                                                     * */
/*  COPYRIGHT:     TSUTOMU OGUNI       FEB. 1 1993      VER. 2         * */
/* ********************************************************************** */

    /* Parameter adjustments */
    --m;
    --h__;
    --e;
    --r0;
    --b;
    --d__;
    ia_dim1 = *n1;
    ia_offset = 1 + ia_dim1;
    ia -= ia_offset;
    a_dim1 = *n1;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *ier = 0;
    if (*n1 < *n || *s < 0.f) {
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "(SUBR. PCGS) INVALID ARGUMENT. ", (ftnlen)31);
	do_lio(&c__3, &c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*n1), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*nl), (ftnlen)sizeof(integer));
	do_lio(&c__5, &c__1, (char *)&(*s), (ftnlen)sizeof(doublereal));
	e_wsle();
	*ier = 2;
	return 0;
    }

    th = 1.;
    if (*s > 0.f && *s < 1.f) {
	th = *s;
	*s = 1.;
    }
    i__1 = *n << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L5: */
	m[i__] = 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dd[i__] = 0.;
	i__2 = *nl;
	for (j = 1; j <= i__2; ++j) {
	    if (ia[i__ + j * ia_dim1] != 0) {
		++m[i__];
	    }
/* L7: */
	}
	i__2 = *nl << 1;
	for (j = *nl + 1; j <= i__2; ++j) {
	    if (ia[i__ + j * ia_dim1] != 0) {
		++m[i__ + *n];
	    }
/* L8: */
	}
/* L6: */
    }
    dd[0] = 0.;
    x[0] = 0.;
    p[0] = 0.;
    q[0] = 0.;
    r__[0] = 0.;
    w[0] = 0.;

    Amatrix A;
    A.d__=d__;
    A.a  =  a;
    A.ia = ia;
    A.m  =  m;
    A.n  =  n;
    A.nl = nl;
    A.a_dim1  = a_dim1;
    A.ia_dim1 = ia_dim1;
    
    if (*s != 0.f) 
      incompleteMLUdecomposition(dd, s, th, A);
    else 
      incompleteLUdecomposition(dd, A);

    A.dd = dd;
    
    multiply(q, A, x);
    alpha_calculation(r__,  b, -1.0, q, *n);
    presolve(r__, A, r__);
    cp(r__, r0, *n);
    cp(r__, p,  *n);
    cp(r__, e,  *n);
    c1 = dot(r__, r__, *n);

/*  ITERATION PHASE */
    i__1 = *itr;
    for (k = 1; k <= i__1; ++k) {

      multiply(q,A,p);
      presolve(q,A,q);

      c2 = dot(q,r0,*n);
      
      if (c2 == 0.f) {
	    *ier = 3;
	    *itr = k;
	    goto L300;
	}
	alpha = c1 / c2;
	c3 = 0.;
	x1 = 0.;
	x2 = 0.;
	alpha_calculation(h__, e, -alpha, q, *n);
	alpha_calculation(w, e, 1.0, h__, *n);
	multiply(q,A,w);
	presolve(q,A,q);
	alpha_calculation(r__, r__, -alpha, q, *n);

	for (i__ = 1; i__ <= *n; ++i__) {
	  y = x[i__];
	  x[i__] += alpha * w[i__];
	  c3 += r__[i__] * r0[i__];
	  x1 += y * y;
	  d__1 = x[i__] - y;
	  x2 += d__1 * d__1;
	}

	if (x1 != 0.f) {
	  res = sqrt(x2 / x1);
	  if (res <= *eps) {
	    *itr = k;
	    *ier = 0;
	    goto L300;
	  }
	}
	if (c1 == 0.f) {
	  *ier = 4;
	  *itr = k;
	  goto L300;
	}
	beta = c3 / c1;
	c1 = c3;

	alpha_calculation(e, r__, beta, h__, *n);
	for (i__ = 1; i__ <= *n; ++i__) {
	  p[i__] = e[i__] + beta * (h__[i__] + beta * p[i__]);
	}
	
    }
    *ier = 1;
    s_wsle(&io___18);
    do_lio(&c__9, &c__1, "(SUBR. MLUCGS) NO CONVERGENCE. ", (ftnlen)31);
    e_wsle();
L300:
    *eps = res;
    if (th != 1.) 
      *s = th;
    return 0;
}
