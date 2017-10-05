#pragma GCC diagnostic ignored "-Wimplicit-function-declaration"

/* GLU1T.F -- translated by f2c (version 12.02.01).
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

static integer c__51 = 51;
static doublereal c_b7 = 3.52e-15;
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;


/*             GLU1T  : SAMPLE PROGRAM OF GLU1 AND GSLV1 */
/*                           H.HASEGAWA, OCT.  4 1991 */

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
    /* Format strings */
    static char fmt_600[] = "(\002 \002,a,i4)";
    static char fmt_605[] = "(\002 \002,i5,d20.10)";
    static char fmt_610[] = "(\002 \002,a,d10.3)";

    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal a[2601]	/* was [51][51] */, b[51];
    static integer i__, j;
    static doublereal x[51];
    static integer ii, ip[51], ir;
    static doublereal xj, wk1[51], dif, sum;
    extern /* Subroutine */ int glu1_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    static doublereal smax;
    extern /* Subroutine */ int gslv1_(doublereal *, integer *, integer *, 
	    doublereal *, integer *);
    static doublereal snorm;

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_600, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, fmt_605, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, fmt_610, 0 };
    static cilist io___22 = { 0, 6, 0, fmt_610, 0 };


/*             GIVE A SOLUTION X */
    for (j = 1; j <= 51; ++j) {
/* L10: */
	x[j - 1] = 1.;
    }
/*             MAKE A MATRIX AND COMPUTE RIGHT HAND SIDE */
    for (j = 1; j <= 51; ++j) {
/* L100: */
	b[j - 1] = 0.;
    }
    for (j = 1; j <= 51; ++j) {
	xj = x[j - 1];
	for (i__ = 1; i__ <= 51; ++i__) {
	    ii = 52 - i__;
	    a[ii + j * 51 - 52] = (doublereal) (52 - max(i__,j));
	    b[ii - 1] += a[ii + j * 51 - 52] * xj;
/* L120: */
	}
/* L110: */
    }
/*             SOLVES SIMULTANEOUS LINEAR EQUATIONS */
    glu1_(a, &c__51, &c__51, &c_b7, wk1, ip, &ir);
    if (ir == 0) {
	gslv1_(a, &c__51, &c__51, b, ip);
    } else {
	s_wsle(&io___11);
	do_lio(&c__9, &c__1, " IR : ", (ftnlen)6);
	do_lio(&c__3, &c__1, (char *)&ir, (ftnlen)sizeof(integer));
	e_wsle();
    }
    s_wsle(&io___12);
    do_lio(&c__9, &c__1, " ", (ftnlen)1);
    e_wsle();
    s_wsfe(&io___13);
    do_fio(&c__1, " N             : ", (ftnlen)17);
    do_fio(&c__1, (char *)&c__51, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsle(&io___14);
    do_lio(&c__9, &c__1, " ", (ftnlen)1);
    e_wsle();
    for (j = 1; j <= 51; ++j) {
/* L300: */
	s_wsfe(&io___15);
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&b[j - 1], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
/*             CHECK OF COMPUTATION */
    sum = 0.f;
    smax = 0.f;
    for (j = 1; j <= 51; ++j) {
	dif = (d__1 = x[j - 1] - b[j - 1], abs(d__1));
	smax = max(dif,smax);
/* Computing 2nd power */
	d__1 = dif;
	sum += d__1 * d__1;
/* L400: */
    }
    snorm = sqrt(sum);
    s_wsle(&io___20);
    do_lio(&c__9, &c__1, " ", (ftnlen)1);
    e_wsle();
    s_wsfe(&io___21);
    do_fio(&c__1, " MAX. OF ERROR : ", (ftnlen)17);
    do_fio(&c__1, (char *)&smax, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___22);
    do_fio(&c__1, " NORM OF ERROR : ", (ftnlen)17);
    do_fio(&c__1, (char *)&snorm, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Main program alias */ int glu1t_ () { MAIN__ (); return 0; }
