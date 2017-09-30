/* PCGSB.F -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;


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
    static doublereal ss, sw, res, beta;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal alpha;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);

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
/*  INCOMPLETE CHOLESKY DECOMPOSITION */
    if (*s != 0.f) {
	dd[1] = 1. / (*s * d__[1]);
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
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
/* L17: */
		    }
		}
/* L16: */
		ss -= a[i__ + k * a_dim1] * sw * dd[nn];
	    }
/* L15: */
	    dd[i__] = 1. / ss;
	}
    } else {
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
			ss -= a[i__ + k * a_dim1] * a[nn + j * a_dim1] * dd[
				nn];
		    }
/* L22: */
		}
/* L21: */
	    }
/* L20: */
	    dd[i__] = 1. / ss;
	}
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__] = d__[i__] * x[i__];
	i__2 = m[i__];
	for (j = 1; j <= i__2; ++j) {
/* L32: */
	    q[i__] += a[i__ + j * a_dim1] * x[ia[i__ + j * ia_dim1]];
	}
	i__2 = *nl + m[i__ + *n];
	for (j = *nl + 1; j <= i__2; ++j) {
/* L34: */
	    q[i__] += a[i__ + j * a_dim1] * x[ia[i__ + j * ia_dim1]];
	}
/* L30: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	r__[i__] = b[i__] - q[i__];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = m[i__];
	for (j = 1; j <= i__2; ++j) {
/* L52: */
	    r__[i__] -= a[i__ + j * a_dim1] * r__[ia[i__ + j * ia_dim1]];
	}
/* L50: */
	r__[i__] *= dd[i__];
    }
    for (i__ = *n; i__ >= 1; --i__) {
	sw = 0.;
	i__1 = *nl + m[i__ + *n];
	for (j = *nl + 1; j <= i__1; ++j) {
/* L62: */
	    sw += a[i__ + j * a_dim1] * r__[ia[i__ + j * ia_dim1]];
	}
/* L60: */
	r__[i__] -= dd[i__] * sw;
    }
/*       C1 = 0.0D0 */
/*       DO 70 I=1,N */
/*        R0(I) = R(I) */
/*        P(I) = R(I) */
/*        E(I) = R(I) */
/*  70    C1 = C1 + R(I) * R(I) */
    dcopy_(n, &r__[1], &c__1, &r0[1], &c__1);
    dcopy_(n, &r__[1], &c__1, &p[1], &c__1);
    dcopy_(n, &r__[1], &c__1, &e[1], &c__1);
    c1 = ddot_(n, &r__[1], &c__1, &r__[1], &c__1);
/*  ITERATION PHASE */
    i__1 = *itr;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    q[i__] = d__[i__] * p[i__];
	    i__3 = m[i__];
	    for (j = 1; j <= i__3; ++j) {
/* L85: */
		q[i__] += a[i__ + j * a_dim1] * p[ia[i__ + j * ia_dim1]];
	    }
	    i__3 = *nl + m[*n + i__];
	    for (j = *nl + 1; j <= i__3; ++j) {
/* L87: */
		q[i__] += a[i__ + j * a_dim1] * p[ia[i__ + j * ia_dim1]];
	    }
/* L80: */
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = m[i__];
	    for (j = 1; j <= i__3; ++j) {
/* L95: */
		q[i__] -= a[i__ + j * a_dim1] * q[ia[i__ + j * ia_dim1]];
	    }
/* L90: */
	    q[i__] = dd[i__] * q[i__];
	}
	for (i__ = *n; i__ >= 1; --i__) {
	    sw = 0.;
	    i__2 = *nl + m[i__ + *n];
	    for (j = *nl + 1; j <= i__2; ++j) {
/* L105: */
		sw += a[i__ + j * a_dim1] * q[ia[i__ + j * ia_dim1]];
	    }
/* L100: */
	    q[i__] -= dd[i__] * sw;
	}
/*       C2 = 0.0D0 */
/*       DO 110 I=1,N */
/*  110   C2 = C2 + Q(I) * R0(I) */
	c2 = ddot_(n, &q[1], &c__1, &r0[1], &c__1);
	if (c2 == 0.f) {
	    *ier = 3;
	    *itr = k;
	    goto L300;
	}
	alpha = c1 / c2;
/* 	C3 = 0.0D0 */
/* 	X1 = 0.0D0 */
	x2 = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    h__[i__] = e[i__] - alpha * q[i__];
/* L120: */
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L130: */
	    w[i__] = e[i__] + h__[i__];
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    q[i__] = d__[i__] * w[i__];
	    i__3 = m[i__];
	    for (j = 1; j <= i__3; ++j) {
/* L142: */
		q[i__] += a[i__ + j * a_dim1] * w[ia[i__ + j * ia_dim1]];
	    }
	    i__3 = *nl + m[i__ + *n];
	    for (j = *nl + 1; j <= i__3; ++j) {
/* L144: */
		q[i__] += a[i__ + j * a_dim1] * w[ia[i__ + j * ia_dim1]];
	    }
/* L140: */
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = m[i__];
	    for (j = 1; j <= i__3; ++j) {
/* L155: */
		q[i__] -= a[i__ + j * a_dim1] * q[ia[i__ + j * ia_dim1]];
	    }
/* L150: */
	    q[i__] = dd[i__] * q[i__];
	}
	for (i__ = *n; i__ >= 1; --i__) {
	    sw = 0.;
	    i__2 = *nl + m[i__ + *n];
	    for (j = *nl + 1; j <= i__2; ++j) {
/* L165: */
		sw += a[i__ + j * a_dim1] * q[ia[i__ + j * ia_dim1]];
	    }
/* L160: */
	    q[i__] -= dd[i__] * sw;
	}
	d__1 = -alpha;
	daxpy_(n, &d__1, &q[1], &c__1, &r__[1], &c__1);
	c3 = ddot_(n, &r__[1], &c__1, &r0[1], &c__1);
	x1 = ddot_(n, &x[1], &c__1, &x[1], &c__1);
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y = x[i__];
/*        R(I) = R(I) - ALPHA * Q(I) */
	    x[i__] += alpha * w[i__];
/*        C3 = C3 + R(I) * R0(I) */
/*        X1 = X1 + Y * Y */
/* L170: */
/* Computing 2nd power */
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
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e[i__] = r__[i__] + beta * h__[i__];
	    p[i__] = e[i__] + beta * (h__[i__] + beta * p[i__]);
/* L180: */
	}

/* L200: */
    }
    *ier = 1;
    s_wsle(&io___18);
    do_lio(&c__9, &c__1, "(SUBR. MLUCGS) NO CONVERGENCE. ", (ftnlen)31);
    e_wsle();
L300:
    *eps = res;
    if (th != 1.) {
	*s = th;
    }
    return 0;
/*  END OF PCGS */
} /* pcgs_ */

