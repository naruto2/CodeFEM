#pragma GCC diagnostic ignored "-Wimplicit-function-declaration"

/* GLU1.F -- translated by f2c (version 12.02.01).
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


/* Subroutine */ int glu1_(doublereal *a, integer *n, integer *n1, doublereal 
	*eps, doublereal *wk, integer *ip, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t, w, aik;
    static integer ipk;
    static doublereal amax;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___10 = { 0, 6, 0, 0, 0 };



/*        GLU1 */
/*               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1 */

/*               SOLVES SIMULTANEOUS LINEAR EQUATIONS */
/*               BY GAUSSIAN ELIMINATION METHOD. */

/*        INPUT - - */
/*             A(N1,N)  R *8  : 2-DIM. ARRAY CONTAINING THE COEFFICIENTS. */
/*             N        I *4  : ORDER OF MATRIX. */
/*             N1       I *4  : SIZE OF ARRAY A. */
/*             EPS      R *8  : PARAMETER TO CHECK SINGULARITY OF THE */
/*                              MATRIX. ( STANDARD VALUE 3.52D-15 ) */
/*        OUTPUT - - */
/*             A(N1,N)        : RESULT OF GAUSSIAN ELIMINATION. */
/*             IP(N)    I *4  : PIVOT NUMBER. */
/*             IER      I *4  : = 0,  FOR NORMAL EXECUTION. */
/*                              = 1,  FOR SINGULAR MATRIX. */
/*                              = 3,  FOR INVALID ARGUEMENT. */
/*        WORKING  - */
/*             WK(N)    R *8  : 1-DIM. ARRAY. */

/*             LEFT HAND SIDE */
    /* Parameter adjustments */
    a_dim1 = *n1;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --wk;
    --ip;

    /* Function Body */
    if (*eps < 0.) {
	*eps = 3.52e-15;
    }
    if (*n > *n1 || *n <= 0) {
	*ier = 3;
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "  (SUBR. GLU1)  INVALID ARGUMENT.  N1, N =", (
		ftnlen)42);
	do_lio(&c__3, &c__1, (char *)&(*n1), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	e_wsle();
	return 0;
    }

    *ier = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/*             FIND MAXIMUM ELEMENT IN THE K-TH COLUMN. */
	amax = (d__1 = a[k + k * a_dim1], abs(d__1));
	ipk = k;
	i__2 = *n;
	for (i__ = k + 1; i__ <= i__2; ++i__) {
	    aik = (d__1 = a[i__ + k * a_dim1], abs(d__1));
	    if (aik > amax) {
		ipk = i__;
		amax = aik;
	    }
/* L110: */
	}
	ip[k] = ipk;

	if (amax > *eps) {
	    if (ipk != k) {
		w = a[ipk + k * a_dim1];
		a[ipk + k * a_dim1] = a[k + k * a_dim1];
		a[k + k * a_dim1] = w;
	    }
/*             COMPUTE ALFA */
	    i__2 = *n;
	    for (i__ = k + 1; i__ <= i__2; ++i__) {
		a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] / a[k + k * a_dim1]
			;
/* L120: */
		wk[i__] = a[i__ + k * a_dim1];
	    }

	    i__2 = *n;
	    for (j = k + 1; j <= i__2; ++j) {
		if (ipk != k) {
		    w = a[ipk + j * a_dim1];
		    a[ipk + j * a_dim1] = a[k + j * a_dim1];
		    a[k + j * a_dim1] = w;
		}
/*             GAUSSIAN ELIMINATION */
		t = a[k + j * a_dim1];
		i__3 = *n;
		for (i__ = k + 1; i__ <= i__3; ++i__) {
/* L140: */
		    a[i__ + j * a_dim1] += wk[i__] * t;
		}
/* L130: */
	    }
/*             MATRIX IS SINGULAR. */
	} else {
	    *ier = 1;
	    ip[k] = k;
	    i__2 = *n;
	    for (i__ = k + 1; i__ <= i__2; ++i__) {
/* L150: */
		a[i__ + k * a_dim1] = 0.;
	    }
	    s_wsle(&io___10);
	    do_lio(&c__9, &c__1, "  (SUBR. GLU1)  MATRIX IS SINGULAR AT K =", 
		    (ftnlen)41);
	    do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsle();
	    return 0;
	}
/* L100: */
    }
    return 0;
} /* glu1_ */


/* Subroutine */ int gslv1_(doublereal *a, integer *n, integer *n1, 
	doublereal *b, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k;
    static doublereal t, w;


/*        GSLV1 */
/*               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1 */

/*               SOLVES SIMULTANEOUS LINEAR EQUATIONS */
/*               BY GAUSSIAN ELIMINATION METHOD. */

/*        INPUT - - */
/*             A(N1,N)  R *8  : RESULT OF GAUSSIAN ELIMINATION. */
/*             N        I *4  : ORDER OF MATRIX. */
/*             N1       I *4  : SIZE OF ARRAY A. */
/*             B(N)     R *8  : 1-DIM. ARRAY CONTAINING THE RIGHT HAND */
/*                              SIDE VECTOR. */
/*             IP(N)    I *4  : PIVOT NUMBER. */
/*        OUTPUT - - */
/*             B(N)           : SOLUTION. */

/*             FORWARD ELIMINATION PROCESS */
    /* Parameter adjustments */
    a_dim1 = *n1;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;
    --ip;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (ip[k] != k) {
	    w = b[ip[k]];
	    b[ip[k]] = b[k];
	    b[k] = w;
	}
/*             GAUSSIAN ELIMINATION */
	t = b[k];
	i__2 = *n;
	for (i__ = k + 1; i__ <= i__2; ++i__) {
/* L110: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
/* L100: */
    }
/*             BACKWARD SUBSTITUTION PROCESS */
    b[*n] /= a[*n + *n * a_dim1];
    for (k = *n - 1; k >= 1; --k) {
	t = b[k + 1];
	i__1 = k;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L210: */
	    b[i__] -= a[i__ + (k + 1) * a_dim1] * t;
	}
	b[k] /= a[k + k * a_dim1];
/* L200: */
    }
    return 0;
} /* gslv1_ */

