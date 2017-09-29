/* STWART.F -- translated by f2c (version 12.02.01).
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


/* Subroutine */ int stwart_(integer *ia, integer *l, integer *ma, integer *m,
	 integer *r__, integer *c__, integer *ir, integer *ic, integer *jrow, 
	integer *jcol, integer *ip, integer *jp, integer *kerns, integer *
	mend, integer *iw, integer *lg, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, n, mm, is, min__, kern, irow;
    extern /* Subroutine */ int single_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), forful_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 6, 0, 0, 0 };


/* ********************************************************************** */
/*  STWART METHOD OF BLOCKING FOR NON-SYMMETRIC SPARSE MATRIX.         * */
/*                                                                     * */
/*  PARAMETERS:                                                        * */
/*   ON ENTRY:                                                         * */
/*     IA     THE ARRAY WHICH CONTAINS ROW INDEX OF NON-ZERO ELEMENTS  * */
/*            OF THE MATRIX.                                           * */
/*     L      THE LEADING DIMENSION OF THE ARRAY A.                    * */
/*     MA     ACCUMULATED SUM OF NUMBERS OF NON-ZERO ELEMENTS IN       * */
/*            EACH COLUMN OF THE MATRIX A.                             * */
/*     M      THE ORDER OF THE MATRIX A.                               * */
/*   ON RETURN:                                                        * */
/*     JROW   THE INFORMATION ABOUT CHANGE OF ROWS.                    * */
/*     JCOL   THE INFORMATION ABOUT CHANGE OF COLUMNS.                 * */
/*     IW     BLOCK INDEX OF EACH COLUMN.                              * */
/*     LG     NUMBER OF BLOCKS OF THE KERNEL.                          * */
/*     IER    THE ERROR CODE. IF IER=0, NORMAL RETURN.                 * */
/*   OTHERS:  WORKING PARAMETERS.                                      * */
/*                                                                     * */
/*  COPYRIGHT:     TSUTOMU OGUNI     SEP. 1 1991        VER. 1         * */
/* ********************************************************************** */

    /* Parameter adjustments */
    --ia;
    --iw;
    --jp;
    --ip;
    --jcol;
    --jrow;
    --ic;
    --ir;
    --c__;
    --r__;

    /* Function Body */
    *ier = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__[i__] = 0;
	c__[i__] = 0;
	iw[i__] = 0;
	jrow[i__] = 0;
	jcol[i__] = 0;
	ir[i__] = 0;
/* L1: */
	ic[i__] = 0;
    }

    single_(&ia[1], l, &ir[1], &ic[1], &r__[1], &c__[1], ma, m, mend, kerns, &
	    jcol[1], &jrow[1]);
    forful_(&ia[1], l, ma, m, &ip[1], &jp[1], &ir[1], &ic[1], kerns, ier);

    kern = *kerns;
    mm = *m;
    *lg = 1;
L5:

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__[i__] = 0;
/* L10: */
	c__[i__] = 0;
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (jcol[j] == 0) {
	    i__2 = ma[j];
	    for (n = ma[j - 1] + 1; n <= i__2; ++n) {
		if (jrow[ia[n]] == 0) {
		    ++r__[ia[n]];
		    ++c__[j];
		}
/* L30: */
	    }
	}
/* L20: */
    }
/*  SEARCH OF MINIMUM ROW COUNT */
    min__ = 100000;
    irow = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (jrow[i__] == 0) {
	    if (r__[i__] < min__) {
		min__ = r__[i__];
		irow = i__;
	    }
	}
/* L40: */
    }
    if (irow == 0) {
	s_wsle(&io___8);
	do_lio(&c__9, &c__1, "(SUBR. STWART) STOP AT. ", (ftnlen)24);
	do_lio(&c__3, &c__1, (char *)&kern, (ftnlen)sizeof(integer));
	e_wsle();
	*ier = 1;
	return 0;
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (jcol[j] == 0) {
	    i__2 = ma[j];
	    for (n = ma[j - 1] + 1; n <= i__2; ++n) {
		if (ia[n] == irow) {
		    iw[j] = *lg;
		    ++kern;
		}
/* L55: */
	    }
	}
/* L50: */
    }

L57:
    is = 0;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (iw[j] == *lg) {
	    irow = jp[j];
	    i__2 = *m;
	    for (k = 1; k <= i__2; ++k) {
		if (jcol[k] == 0) {
		    if (iw[k] == 0) {
			i__3 = ma[k];
			for (n = ma[k - 1] + 1; n <= i__3; ++n) {
			    if (ia[n] == irow) {
				iw[k] = *lg;
				++kern;
				++is;
			    }
/* L80: */
			}
		    }
		}
/* L70: */
	    }
	}
/* L60: */
    }

    if (is != 0) {
	goto L57;
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	if (iw[j] == *lg) {
	    jcol[j] = mm;
	    jrow[jp[j]] = mm;
	    --mm;
	}
/* L90: */
    }
    if (*lg < *m) {
	if (kern <= *m) {
	    ++(*lg);
	    goto L5;
	}
    }
    *ier = 0;
    return 0;
/*  END OF STWART */
} /* stwart_ */

