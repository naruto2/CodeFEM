/* MARKWZT.F -- translated by f2c (version 12.02.01).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdio.h>
#include <stdlib.h> /* For exit() */
#include <f2c.h>

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__9 = 9;


/* Actual main program */
int acutualmain(int argc, char **argv, integer *ib, integer ll, integer *mb, integer mm)
{
  extern int MAIN__(integer *ib, integer ll, integer *mb, integer mm);
	libf2c_init(argc, argv);
	MAIN__(ib,ll,mb,mm);
	libf2c_close();
	exit(0);
	return 0;
}

/* Main program */
int MAIN__(integer *ib, integer ll, integer *mb, integer mm)
{
    /* Initialized data */

  static integer ia[100] =
    { 
         2,                      10, 11,    13,    15, 16,   //  6
	       4, 5,             10,                   16,   // 10
                       	      9,     11, 12,                 // 13
      1,          5,    7, 8,                                // 17
      1, 2,    4, 5, 6,       9, 10,                         // 24
	       4,             9,         12,                 // 27
         2, 3, 4,                               14,          // 31
      1, 2, 9,                                          16,  // 35
      	       4,                                   15,      // 37
	             6, 7, 8,                   14, 15,      // 42
         2,                                             16,  // 44
	                      9,     11, 12,                 // 47
	             6,                         14, 15,      // 50
	                                            15,      // 51
	       4,    6,                  12,        15,      // 55
         2,    4,                10,        13,              // 59
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0 };
    static integer ma[21] = { 0,6,10,13,17,24,27,31,35,37,42,44,47,50,51,55,
	    59,0,0,0,0 };
    static integer m = 16;
    static integer l = 100;
    static integer ip[100] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	    0,0,0,0,0,0,0,0,0,0,0 };

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal a[100];
    static integer c__[20];
    static doublereal e[200];
    static integer i__, j, r__[20];
    static doublereal w[20], z__[20];
    static integer ic[20], ie[200], me[21], lg, jp[20], ir[20], is, iw[20], 
	    ier, num, mend, jcol[20], nume, jrow[20], kerns;
    extern /* Subroutine */ int markwz_(doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *),
	     stwart_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___27 = { 0, 6, 0, 0, 0 };
    static cilist io___28 = { 0, 6, 0, 0, 0 };
    static cilist io___32 = { 0, 6, 0, 0, 0 };
    static cilist io___33 = { 0, 6, 0, 0, 0 };
    static cilist io___34 = { 0, 6, 0, 0, 0 };
    static cilist io___35 = { 0, 6, 0, 0, 0 };
    static cilist io___36 = { 0, 6, 0, 0, 0 };
    static cilist io___37 = { 0, 6, 0, 0, 0 };
    static cilist io___38 = { 0, 6, 0, 0, 0 };
    static cilist io___39 = { 0, 6, 0, 0, 0 };


/* ********************************************************************** */
/*  SAMPLE PROGRAM FOR MARKWZ AND STWART ROUTINE.          VER. 1      * */
/* ********************************************************************** */

    is = 0;
L5:
    num = 0;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ma[i__];
	for (j = ma[i__ - 1] + 1; j <= i__2; ++j) {
	    ++num;
	    a[num - 1] = (doublereal) num;
/* L20: */
	}
    }
    if (is == 0) {
	markwz_(a, ia, &l, ma, &m, e, ie, ir, ic, r__, c__, w, z__, me, &nume,
		 &mend, &kerns, jrow, jcol, &ier);
	s_wsle(&io___27);
	do_lio(&c__3, &c__1, (char *)&num, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&nume, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&mend, (ftnlen)sizeof(integer));
	i__2 = kerns - 1;
	do_lio(&c__3, &c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___28);
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_lio(&c__3, &c__1, (char *)&me[i__], (ftnlen)sizeof(integer));
	}
	e_wsle();
    } else {
      integer *R    = (integer*)calloc(sizeof(integer),mm);
      integer *C    = (integer*)calloc(sizeof(integer),mm);
      integer *IR   = (integer*)calloc(sizeof(integer),mm);
      integer *IC   = (integer*)calloc(sizeof(integer),mm);
      integer *JROW = (integer*)calloc(sizeof(integer),mm);
      integer *JCOL = (integer*)calloc(sizeof(integer),mm);
      integer *IP   = (integer*)calloc(sizeof(integer),mm);
      integer *JP   = (integer*)calloc(sizeof(integer),mm);
      integer *IW   = (integer*)calloc(sizeof(integer),mm);


      /*
      stwart_(ib, &ll, mb, &mm, r__, c__, ir, ic, jrow, jcol, ip, jp, &kerns, 
      &mend, iw, &lg, &ier); */
      stwart_(ib, &ll, mb, &mm, R, C, IR, IC, JROW, JCOL, IP, JP, &kerns, 
		&mend, IW, &lg, &ier);

      printf("m=%d, l=%d, kerns=%d, mend=%d, lg=%d, ier=%d\n",
	     mm,    ll,    kerns,    mend,    lg,    ier);


      free(R);
      free(C);
      free(IR);
      free(IC);
      free(JROW);
      free(JCOL);
      free(IP);
      free(JP);
      
        s_wsle(&io___32);
	do_lio(&c__9, &c__1, "STEWART METHOD. ", (ftnlen)16);
	do_lio(&c__3, &c__1, (char *)&mend, (ftnlen)sizeof(integer));
	i__2 = kerns - 1;
	do_lio(&c__3, &c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&lg, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___33);
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_lio(&c__3, &c__1, (char *)&ip[i__ - 1], (ftnlen)sizeof(integer)
		    );
	}
	e_wsle();
	s_wsle(&io___34);
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_lio(&c__3, &c__1, (char *)&jp[i__ - 1], (ftnlen)sizeof(integer)
		    );
	}
	e_wsle();
	s_wsle(&io___35);
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    do_lio(&c__3, &c__1, (char *)&iw[i__ - 1], (ftnlen)sizeof(integer)
		    );
	}
	e_wsle();
    }
    s_wsle(&io___36);
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_lio(&c__3, &c__1, (char *)&jcol[i__ - 1], (ftnlen)sizeof(integer));
    }
    e_wsle();
    s_wsle(&io___37);
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_lio(&c__3, &c__1, (char *)&jrow[i__ - 1], (ftnlen)sizeof(integer));
    }
    e_wsle();
    s_wsle(&io___38);
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_lio(&c__3, &c__1, (char *)&ic[i__ - 1], (ftnlen)sizeof(integer));
    }
    e_wsle();
    s_wsle(&io___39);
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	do_lio(&c__3, &c__1, (char *)&ir[i__ - 1], (ftnlen)sizeof(integer));
    }
    e_wsle();
    if (is == 0) {
	is = 1;
	goto L5;
    }
    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

// /* Main program alias */ int ordrng_ () { MAIN__ (); return 0; }
