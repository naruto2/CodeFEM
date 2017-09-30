*
*             GLU1T  : SAMPLE PROGRAM OF GLU1 AND GSLV1
*                           H.HASEGAWA, OCT.  4 1991
*
      PROGRAM GLU1T
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER ( N=51, N1=N, EPS=3.52D-15 )
      DIMENSION A(N1,N1), IP(N1), X(N), B(N1), WK1(N)
*             GIVE A SOLUTION X
      DO 10 J = 1, N
   10   X(J) = 1.0D0
*             MAKE A MATRIX AND COMPUTE RIGHT HAND SIDE
      DO 100 J = 1, N
  100   B(J) = 0.0D0
      DO 110 J = 1, N
        XJ = X(J)
        DO 120 I = 1, N
          II = N+1-I
          A(II,J) = DBLE(N+1-MAX(I,J))
          B(II) = B(II) + A(II,J)*XJ
  120   CONTINUE
  110 CONTINUE
*             SOLVES SIMULTANEOUS LINEAR EQUATIONS
      CALL GLU1( A, N, N, EPS, WK1, IP, IR )
      IF( IR.EQ.0 )  THEN
         CALL GSLV1( A, N, N, B, IP )
      ELSE
         WRITE(*,*) ' IR : ', IR
      END IF
      WRITE(*,*) ' '
      WRITE(*,600) ' N             : ', N
      WRITE(*,*) ' '
      DO 300 J = 1, N
  300   WRITE(*,605)  J, B(J)
  600 FORMAT(' ',A,I4)
  605 FORMAT(' ',I5,D20.10)
*             CHECK OF COMPUTATION
      SUM = 0.0
      SMAX = 0.0
      DO 400 J = 1, N
        DIF = ABS( X(J)-B(J) )
        SMAX = MAX( DIF, SMAX )
        SUM = SUM + DIF**2
  400 CONTINUE
      SNORM = SQRT(SUM)
      WRITE(*,*) ' '
      WRITE(*,610) ' MAX. OF ERROR : ', SMAX
      WRITE(*,610) ' NORM OF ERROR : ', SNORM
  610 FORMAT(' ',A,D10.3)
      STOP
      END
