/* FORFUL.F -- translated by f2c (version 12.02.01).
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


/* Subroutine */ int forful_(integer *ia, integer *l, integer *ma, integer *m,
	 integer *ip, integer *jp, integer *ir, integer *ic, integer *kerns, 
	integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, nd, nh, nk, nn, nl, np;

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___9 = { 0, 6, 0, 0, 0 };


/* ********************************************************************** */
/*  FORD-FULKERSON METHOD FOR ORDERING.                                * */
/*                                                                     * */
/*  PARAMETERS:                                                        * */
/*   ON ENTRY:  SAME AS STWART ROUTINE.                                * */
/*   ON RETURN:                                                        * */
/*     IP     COLUMN INDEX IN EACH I-TH ROW. ON THE CASE OF SINGLETON, * */
/*            MINUS SIGN.                                              * */
/*     JP     ROW INDEX IN EACH J-TH COLUMN.                           * */
/*     IER    THE ERROR CODE. IF IER=0, NORMAL RETURN.                 * */
/*   OTHERS:  WORKING PARAMETERS.                                      * */
/*                                                                     * */
/*  COPYRIGHT:    TSUTOMU OGUNI      SEP. 1 1991         VER. 1        * */
/* ********************************************************************** */

    /* Parameter adjustments */
    --ia;
    --ic;
    --ir;
    --jp;
    --ip;

    /* Function Body */
    *ier = 0;
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	jp[k] = 0;
/* L1: */
	ip[k] = 0;
    }
    i__1 = *kerns - 1;
    for (k = 1; k <= i__1; ++k) {
	ip[ir[k]] = -ic[k];
/* L3: */
	jp[ic[k]] = ir[k];
    }
    i__ = 1;
    j = 0;
L2:
    if (jp[i__] == 0) {
	i__1 = ma[i__];
	for (k = ma[i__ - 1] + 1; k <= i__1; ++k) {
	    if (ip[ia[k]] == 0) {
		nk = ia[k];
		ip[nk] = i__;
		jp[i__] = nk;
		goto L20;
	    }
/* L10: */
	}
	goto L50;
    }
L20:
    if (i__ >= *m) {
	return 0;
    }
    ++i__;
    goto L2;
L50:
    nn = 0;
L52:
    ++nn;
    if (nn > ma[i__] - ma[i__ - 1]) {
	s_wsle(&io___6);
	do_lio(&c__9, &c__1, "(SUBR. FORFUL) ERROR STOP. ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(integer));
	e_wsle();
	*ier = 1;
	return 0;
    }
    nh = ia[ma[i__ - 1] + nn];
    j = ip[nh];
    if (j <= 0) {
	goto L52;
    }
    ip[nh] = i__;
    jp[i__] = nh;
L55:
    i__1 = ma[j];
    for (k = ma[j - 1] + 1; k <= i__1; ++k) {
	if (ip[ia[k]] == 0) {
	    nl = ia[k];
	    ip[nl] = j;
	    jp[j] = nl;
	    goto L20;
	}
/* L60: */
    }
    nn = 0;
L62:
    ++nn;
    if (nn > ma[j] - ma[j - 1]) {
	s_wsle(&io___9);
	do_lio(&c__9, &c__1, "(SUBR. FORFUL) ERROR STOP. ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(integer));
	e_wsle();
	*ier = 1;
	return 0;
    }
    np = ia[ma[j - 1] + nn];
    nd = ip[np];
    if (nd <= 0) {
	goto L62;
    }
    ip[np] = j;
    jp[j] = np;
    j = nd;
    goto L55;
/*  END OF FORFUL */
} /* forful_ */

