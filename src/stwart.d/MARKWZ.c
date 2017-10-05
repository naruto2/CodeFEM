#pragma GCC diagnostic ignored "-Wimplicit-function-declaration"

/* MARKWZ.F -- translated by f2c (version 12.02.01).
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


/* Subroutine */ int markwz_(doublereal *a, integer *ia, integer *l, integer *
	ma, integer *m, doublereal *e, integer *ie, integer *ir, integer *ic, 
	integer *r__, integer *c__, doublereal *w, doublereal *z__, integer *
	me, integer *nume, integer *mend, integer *kerns, integer *jrow, 
	integer *jcol, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, n, s, t, kk, val, min__, icol;
    extern /* Subroutine */ int single_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };


/* ********************************************************************** */
/*  MARKOWIZ TRIANGULARIZATION MEHOD FOR NON-SYMMETRIC SPARSE MATRIX.  * */
/*                                                                     * */
/*  PARAMETERS:                                                        * */
/*   ON ENTRY:                                                         * */
/*     A      THE ARRAY WHICH CONTAINS NON-ZERO ELEMENTS               * */
/*            OF THE MATRIX IN COLUMN-WISE.                            * */
/*     IA     THE ARRAY WHICH CONTAINS ROW INDEX OF CORRESPONDING      * */
/*            ELEMENTS IN THE ARRAY A.                                 * */
/*     L      THE LEADING DIMENSION OF THE ARRAY A AND IA.             * */
/*     MA     THE ARRAY WHICH CONTAINS ACCUMULATED SUM OF ELEMENTS     * */
/*            IN EACH COLUMN OF THE ARRAY A.                           * */
/*   ON RETURN:                                                        * */
/*     E      THE ARRAY WHICH IS CONSTRUCTED BY THE MATRIX A AND       * */
/*            IT'S FILL-IN ELEMENTS.                                   * */
/*     IE     THE ARRAY WHICH CONTAINS ROW INDEX OF CORRESPOIDING      * */
/*            ELEMENTS IN THE ARRAY E.                                 * */
/*     ME     THE ARRAY WHICH CONTAINS ACCUMULATED SUM OF ELEMENTS     * */
/*            IN EACH COLUMN OF THE ARRAY E.                           * */
/*     NUME   THE NUMBER OF NON-ZERO ELEMENTS IN THE MATRIX E.         * */
/*     JROW   THE INFORMATION ABOUT CHANGE BETWEEN ROWS.               * */
/*     JCOL   THE INFORMATION ABOUT CHANGE BETWEEN COLUMNS.            * */
/*     IER    ERROR CODE. IF IER = 0, NORMAL RETURN.                   * */
/*   OTHER PARAMETERS:  WORKING PARAMETERS.                            * */
/*                                                                     * */
/*  COPYRIGHT:     TSUTOMU OGUNI    SEP. 1 1992      VER. 2            * */
/* ********************************************************************** */

    /* Parameter adjustments */
    --ie;
    --e;
    --ia;
    --a;
    --jcol;
    --jrow;
    --z__;
    --w;
    --c__;
    --r__;
    --ic;
    --ir;

    /* Function Body */
    *ier = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__[i__] = 0;
	c__[i__] = 0;
	jrow[i__] = 0;
	jcol[i__] = 0;
	ir[i__] = 0;
	ic[i__] = 0;
/* L1: */
    }
    me[0] = 0;

    single_(&ia[1], l, &ir[1], &ic[1], &r__[1], &c__[1], ma, m, mend, kerns, &
	    jcol[1], &jrow[1]);

    *nume = 0;
    kk = *kerns;
/*  ROW SINGLETON */
    if (*mend != 0) {
	i__1 = *mend;
	for (k = 1; k <= i__1; ++k) {
	    icol = ic[k];
	    i__2 = ma[icol];
	    for (n = ma[icol - 1] + 1; n <= i__2; ++n) {
		++(*nume);
		e[*nume] = a[n];
/* L14: */
		ie[*nume] = ia[n];
	    }
/* L13: */
	    me[k] = *nume;
	}
    }
/*  COLUMN SINGLETON */
    if (*mend != *kerns - 1) {
	i__1 = *kerns - 1;
	for (k = *mend + 1; k <= i__1; ++k) {
	    j = ic[k];
	    i__2 = ma[j];
	    for (n = ma[j - 1] + 1; n <= i__2; ++n) {
		++(*nume);
		e[*nume] = a[n];
/* L76: */
		ie[*nume] = ia[n];
	    }
/* L75: */
	    me[k] = *nume;
	}
    }

    i__1 = *m;
    for (k = *kerns; k <= i__1; ++k) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (jrow[i__] == 0) {
		r__[i__] = 0;
	    }
/* L11: */
	}
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    if (jcol[j] == 0) {
		c__[j] = 0;
	    }
/* L12: */
	}
	s = 0;
	t = 0;
	min__ = 100000;

	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    if (jcol[j] == 0) {
		i__3 = *m;
		for (i__ = 1; i__ <= i__3; ++i__) {
/* L22: */
		    z__[i__] = 0.;
		}
		i__3 = ma[j];
		for (n = ma[j - 1] + 1; n <= i__3; ++n) {
/* L24: */
		    z__[ia[n]] = 1.;
		}
/*  UPDATING BY NEW COLUMNS */
		if (kk != *kerns) {
		    i__3 = *kerns - 1;
		    for (i__ = kk; i__ <= i__3; ++i__) {
			if (z__[ir[i__]] != 0.f) {
			    i__4 = me[i__];
			    for (n = me[i__ - 1] + 1; n <= i__4; ++n) {
				if (jrow[ie[n]] == 0) {
				    z__[ie[n]] = 1.;
				}
/* L26: */
			    }
			}
/* L25: */
		    }
		}
/*  COMPUTE R AND C */
		i__3 = *m;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    if (z__[i__] != 0.) {
			if (jrow[i__] == 0) {
			    ++r__[i__];
			    ++c__[j];
			}
		    }
/* L27: */
		}
		if (c__[j] == 0) {
		    s_wsle(&io___10);
		    do_lio(&c__9, &c__1, "(SUBR. MARKWZ) ZERO-COLUMN. ", (
			    ftnlen)28);
		    do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(integer));
		    e_wsle();
		    *ier = 1;
		    return 0;
		}
	    }
/* L180: */
	}

	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    if (jcol[j] == 0) {
		i__3 = *m;
		for (i__ = 1; i__ <= i__3; ++i__) {
/* L290: */
		    z__[i__] = 0.;
		}
		i__3 = ma[j];
		for (n = ma[j - 1] + 1; n <= i__3; ++n) {
/* L300: */
		    z__[ia[n]] = 1.;
		}
		if (kk != *kerns) {
		    i__3 = *kerns - 1;
		    for (i__ = kk; i__ <= i__3; ++i__) {
			if (z__[ir[i__]] != 0.f) {
			    i__4 = me[i__];
			    for (n = me[i__ - 1] + 1; n <= i__4; ++n) {
				if (jrow[ie[n]] == 0) {
				    z__[ie[n]] = 1.;
				}
/* L330: */
			    }
			}
/* L320: */
		    }
		}
		i__3 = *m;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    if (jrow[i__] == 0) {
			if (z__[i__] != 0.f) {
			    val = (c__[j] - 1) * (r__[i__] - 1);
			    if (val < min__) {
				min__ = val;
				s = j;
				t = i__;
			    }
			}
		    }
/* L50: */
		}
	    }
/* L280: */
	}
	if (s == 0) {
	    s_wsle(&io___12);
	    do_lio(&c__9, &c__1, "(SUBR. MARKWZ) STOP AT ", (ftnlen)23);
	    do_lio(&c__3, &c__1, (char *)&(*kerns), (ftnlen)sizeof(integer));
	    e_wsle();
	    *ier = 1;
	    return 0;
	}
/*  GENERATE NEW COLUMN */
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w[i__] = 0.;
/* L41: */
	    z__[i__] = 0.;
	}
	i__2 = ma[s];
	for (n = ma[s - 1] + 1; n <= i__2; ++n) {
	    w[ia[n]] = a[n];
/* L42: */
	    z__[ia[n]] = 1.;
	}

	if (kk != *kerns) {
	    i__2 = *kerns - 1;
	    for (i__ = kk; i__ <= i__2; ++i__) {
		if (z__[ir[i__]] != 0.f) {
		    i__3 = me[i__];
		    for (n = me[i__ - 1] + 1; n <= i__3; ++n) {
			if (jrow[ie[n]] == 0) {
			    z__[ie[n]] = 1.;
			}
/* L44: */
		    }
		}
/* L43: */
	    }
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (z__[i__] != 0.f) {
		++(*nume);
		e[*nume] = w[i__];
		ie[*nume] = i__;
	    }
/* L45: */
	}
	me[*kerns] = *nume;
	ir[*kerns] = t;
	ic[*kerns] = s;
	jcol[s] = *kerns;
	jrow[t] = *kerns;
	++(*kerns);

/* L70: */
    }
    *ier = 0;
    return 0;
/*  END OF MARKWZ */
} /* markwz_ */

