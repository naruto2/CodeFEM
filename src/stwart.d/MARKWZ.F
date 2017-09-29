*
      SUBROUTINE MARKWZ(A,IA,L,MA,M,E,IE,IR,IC,R,C,W,Z,ME,NUME,
     *                  MEND,KERNS,JROW,JCOL,IER)
***********************************************************************
*  MARKOWIZ TRIANGULARIZATION MEHOD FOR NON-SYMMETRIC SPARSE MATRIX.  *
*                                                                     *
*  PARAMETERS:                                                        *
*   ON ENTRY:                                                         *
*     A      THE ARRAY WHICH CONTAINS NON-ZERO ELEMENTS               *
*            OF THE MATRIX IN COLUMN-WISE.                            *
*     IA     THE ARRAY WHICH CONTAINS ROW INDEX OF CORRESPONDING      *
*            ELEMENTS IN THE ARRAY A.                                 *
*     L      THE LEADING DIMENSION OF THE ARRAY A AND IA.             *
*     MA     THE ARRAY WHICH CONTAINS ACCUMULATED SUM OF ELEMENTS     *
*            IN EACH COLUMN OF THE ARRAY A.                           *
*   ON RETURN:                                                        *
*     E      THE ARRAY WHICH IS CONSTRUCTED BY THE MATRIX A AND       *
*            IT'S FILL-IN ELEMENTS.                                   *
*     IE     THE ARRAY WHICH CONTAINS ROW INDEX OF CORRESPOIDING      *
*            ELEMENTS IN THE ARRAY E.                                 *
*     ME     THE ARRAY WHICH CONTAINS ACCUMULATED SUM OF ELEMENTS     *
*            IN EACH COLUMN OF THE ARRAY E.                           *
*     NUME   THE NUMBER OF NON-ZERO ELEMENTS IN THE MATRIX E.         *
*     JROW   THE INFORMATION ABOUT CHANGE BETWEEN ROWS.               *
*     JCOL   THE INFORMATION ABOUT CHANGE BETWEEN COLUMNS.            *
*     IER    ERROR CODE. IF IER = 0, NORMAL RETURN.                   *
*   OTHER PARAMETERS:  WORKING PARAMETERS.                            *
*                                                                     *
*  COPYRIGHT:     TSUTOMU OGUNI    SEP. 1 1992      VER. 2            *
***********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       INTEGER*4 R, C, S, T, VAL
       DIMENSION A(L), IA(L), E(L), IE(L), IR(M), IC(M), R(M), C(M),
     *           W(M), Z(M), ME(0:M), MA(0:M), JROW(M), JCOL(M)  
C
      IER = 0
      DO 1 I=1,M
       R(I) = 0
       C(I) = 0
       JROW(I) = 0
       JCOL(I) = 0
       IR(I) = 0
       IC(I) = 0
    1 CONTINUE 
      ME(0) = 0
C 
      CALL SINGLE(IA,L,IR,IC,R,C,MA,M,MEND,KERNS,JCOL,JROW)
C
      NUME = 0
      KK = KERNS
C  ROW SINGLETON
      IF (MEND .NE. 0) THEN
       DO 13 K=1,MEND
        ICOL = IC(K)
        DO 14 N=MA(ICOL-1)+1,MA(ICOL)
         NUME = NUME + 1
         E(NUME) = A(N)
   14    IE(NUME) = IA(N)
   13   ME(K) = NUME
      ENDIF
C  COLUMN SINGLETON
      IF (MEND .NE. KERNS-1) THEN
       DO 75 K=MEND+1,KERNS-1
        J = IC(K)
        DO 76 N=MA(J-1)+1,MA(J)
         NUME = NUME + 1
         E(NUME) = A(N)
   76    IE(NUME) = IA(N)       
   75   ME(K) = NUME
      ENDIF
C
      DO 70 K=KERNS,M
       DO 11 I=1,M
        IF (JROW(I) .EQ. 0) R(I) = 0
   11  CONTINUE
       DO 12 J=1,M
        IF (JCOL(J) .EQ. 0) C(J) = 0
   12  CONTINUE
       S = 0
       T = 0
       MIN = 100000
C
       DO 180 J=1,M
        IF (JCOL(J) .EQ. 0) THEN
         DO 22 I=1,M
  22      Z(I) = 0.0D0
         DO 24 N=MA(J-1)+1,MA(J)
  24      Z(IA(N)) = 1.0D0
C  UPDATING BY NEW COLUMNS
         IF (KK .NE. KERNS) THEN
          DO 25 I=KK,KERNS-1
           IF (Z(IR(I)) .NE. 0.0) THEN
            DO 26 N=ME(I-1)+1,ME(I)
             IF (JROW(IE(N)) .EQ. 0) Z(IE(N)) = 1.0D0
  26        CONTINUE
           ENDIF
  25      CONTINUE
         ENDIF
C  COMPUTE R AND C
         DO 27 I=1,M
          IF (Z(I) .NE. 0) THEN
           IF (JROW(I) .EQ. 0) THEN
            R(I) = R(I) + 1
            C(J) = C(J) + 1
           ENDIF
          ENDIF
   27    CONTINUE
         IF (C(J) .EQ. 0) THEN
          WRITE(*,*) '(SUBR. MARKWZ) ZERO-COLUMN. ', J
          IER = 1
          RETURN
         ENDIF
        ENDIF
  180  CONTINUE
C
       DO 280 J=1,M
        IF (JCOL(J) .EQ. 0) THEN
         DO 290 I=1,M
  290     Z(I) = 0.0D0
         DO 300 N=MA(J-1)+1,MA(J)
  300     Z(IA(N)) = 1.0D0
         IF (KK .NE. KERNS) THEN
          DO 320 I=KK,KERNS-1
           IF (Z(IR(I)) .NE. 0.0) THEN
            DO 330 N=ME(I-1)+1,ME(I)
             IF (JROW(IE(N)) .EQ. 0) Z(IE(N)) = 1.0D0
  330       CONTINUE
           ENDIF
  320     CONTINUE
         ENDIF
         DO 50 I=1,M
          IF (JROW(I) .EQ. 0) THEN
           IF (Z(I) .NE. 0.0) THEN
            VAL = (C(J) - 1) * (R(I) - 1)
            IF (VAL .LT. MIN) THEN
             MIN = VAL
             S = J
             T = I
            ENDIF
           ENDIF
          ENDIF
   50    CONTINUE      
        ENDIF
  280  CONTINUE
       IF (S .EQ. 0) THEN
        WRITE(*,*) '(SUBR. MARKWZ) STOP AT ', KERNS
        IER = 1
        RETURN
       ENDIF
C  GENERATE NEW COLUMN
       DO 41 I=1,M
        W(I) = 0.0D0
   41   Z(I) = 0.0D0
       DO 42 N=MA(S-1)+1,MA(S)
        W(IA(N)) = A(N)
   42   Z(IA(N)) = 1.0D0
C
       IF (KK .NE. KERNS) THEN
        DO 43 I=KK,KERNS-1
         IF (Z(IR(I)) .NE. 0.0) THEN
          DO 44 N=ME(I-1)+1,ME(I)
           IF (JROW(IE(N)) .EQ. 0) Z(IE(N)) = 1.0D0
   44     CONTINUE
         ENDIF
   43   CONTINUE
       ENDIF
       DO 45 I=1,M
        IF (Z(I) .NE. 0.0) THEN
         NUME = NUME + 1
         E(NUME) = W(I)
         IE(NUME) = I
        ENDIF
   45  CONTINUE
       ME(KERNS) = NUME
       IR(KERNS) = T
       IC(KERNS) = S
       JCOL(S) = KERNS
       JROW(T) = KERNS
       KERNS = KERNS + 1
C
   70 CONTINUE
      IER = 0
      RETURN
C  END OF MARKWZ
      END
