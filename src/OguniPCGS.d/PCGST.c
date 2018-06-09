/* PCGST.f -- translated by f2c (version 12.02.01).
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


/* Actual main program */
int main(int argc, char **argv)
{
	extern int MAIN__();
	libf2c_init(argc, argv);
	MAIN__();
	libf2c_close();
	exit(0);
	return 0;
}

/* Main program */
int MAIN__(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal a[2000]	/* was [100][20] */, b[100], d__[100], e[100],
	     h__[100];
    static integer i__, j, l, m[200], n;
    static doublereal p[101], q[101], r__[101], s, w[101], x[101], r0[100], 
	    dd[101];
    static integer ia[2000]	/* was [100][20] */, nl, ier;
    static doublereal eps;
    static integer itr;
    extern /* Subroutine */ int pcgs_(doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };


/* ********************************************************************** */
/*   SAMPLE PROGRAM FOR ROUTINE PCGS.                VER. 1            * */
/* ********************************************************************** */

    n = 10;
    nl = 3;
    s = 0.;
    l = 100;
L5:
    itr = 50;
    eps = 1e-7;
    i__1 = nl << 1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ia[i__ + j * 100 - 101] = 0;
/* L10: */
	    a[i__ + j * 100 - 101] = 0.;
	}
    }
    i__2 = n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L15: */
	x[i__] = 0.;
    }
    b[0] = 20.;
    b[1] = 23.;
    b[2] = 23.;
    b[3] = 25.;
    b[4] = 24.;
    b[5] = 24.;
    b[6] = 25.;
    b[7] = 23.;
    b[8] = 23.;
    b[9] = 20.;
    i__2 = n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L20: */
	d__[i__ - 1] = 14.;
    }
    i__2 = n;
    for (i__ = 2; i__ <= i__2; ++i__) {
	a[i__ + 298] = 3.;
	ia[i__ + 298] = i__;
	ia[i__ - 1] = i__ - 1;
/* L30: */
	a[i__ - 1] = 3.;
    }
    i__2 = n;
    for (i__ = 4; i__ <= i__2; ++i__) {
	a[i__ + 396] = 2.;
	ia[i__ + 396] = i__;
	ia[i__ + 99] = i__ - 3;
/* L40: */
	a[i__ + 99] = 2.;
    }
    i__2 = n;
    for (i__ = 7; i__ <= i__2; ++i__) {
	a[i__ + 493] = 1.;
	ia[i__ + 493] = i__;
	ia[i__ + 199] = i__ - 6;
/* L50: */
	a[i__ + 199] = 1.;
    }
    pcgs_(d__, a, ia, &n, &l, &nl, b, &eps, &itr, &s, x, dd, p, q, r__, r0, e,
	     h__, w, m, &ier);
    s_wsle(&io___24);
    do_lio(&c__9, &c__1, "EXAMPLE OF PCGS. ", (ftnlen)17);
    do_lio(&c__3, &c__1, (char *)&itr, (ftnlen)sizeof(integer));
    do_lio(&c__5, &c__1, (char *)&eps, (ftnlen)sizeof(doublereal));
    do_lio(&c__5, &c__1, (char *)&s, (ftnlen)sizeof(doublereal));
    e_wsle();
    s_wsle(&io___25);
    i__2 = n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_lio(&c__5, &c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsle();
    if (s == 0.f) {
	s = .95;
	goto L5;
    }
    if (s < 1.f) {
	s = 1.02;
	goto L5;
    }
    if (s > 1.f) {
	s = 1.;
	goto L5;
    }
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Main program alias */ int pcgsts_ () { MAIN__ (); return 0; }
