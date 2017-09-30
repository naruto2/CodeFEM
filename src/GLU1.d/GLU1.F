*
      SUBROUTINE GLU1( A, N, N1, EPS, WK, IP, IER )
*
*        GLU1
*               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1
*
*               SOLVES SIMULTANEOUS LINEAR EQUATIONS
*               BY GAUSSIAN ELIMINATION METHOD.
*
*        INPUT - -
*             A(N1,N)  R *8  : 2-DIM. ARRAY CONTAINING THE COEFFICIENTS.
*             N        I *4  : ORDER OF MATRIX.
*             N1       I *4  : SIZE OF ARRAY A.
*             EPS      R *8  : PARAMETER TO CHECK SINGULARITY OF THE
*                              MATRIX. ( STANDARD VALUE 3.52D-15 )
*        OUTPUT - -
*             A(N1,N)        : RESULT OF GAUSSIAN ELIMINATION.
*             IP(N)    I *4  : PIVOT NUMBER.
*             IER      I *4  : = 0,  FOR NORMAL EXECUTION.
*                              = 1,  FOR SINGULAR MATRIX.
*                              = 3,  FOR INVALID ARGUEMENT.
*        WORKING  -
*             WK(N)    R *8  : 1-DIM. ARRAY.
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N1,*), IP(*), WK(*)
*             LEFT HAND SIDE
      IF( EPS.LT.0.0D0 )  EPS = 3.52D-15
      IF( ( N.GT.N1 ).OR.( N.LE.0 ) )  THEN
         IER = 3
         WRITE(*,*) '  (SUBR. GLU1)  INVALID ARGUMENT.  N1, N =', N1, N
         RETURN
      END IF
*
      IER = 0
      DO 100 K = 1, N
*             FIND MAXIMUM ELEMENT IN THE K-TH COLUMN.
        AMAX = ABS(A(K,K))
        IPK = K
        DO 110 I = K+1, N
          AIK = ABS(A(I,K))
          IF( AIK.GT.AMAX )  THEN
             IPK = I
             AMAX = AIK
          END IF
  110   CONTINUE
        IP(K) = IPK
*
        IF( AMAX.GT.EPS )  THEN
           IF( IPK.NE.K )  THEN
              W = A(IPK,K)
              A(IPK,K) = A(K,K)
              A(K,K) = W
           END IF
*             COMPUTE ALFA
           DO 120 I = K+1, N
             A(I,K) = -A(I,K)/A(K,K)
  120        WK(I) = A(I,K)
*
           DO 130 J = K+1, N
             IF( IPK.NE.K )  THEN
                W = A(IPK,J)
                A(IPK,J) = A(K,J)
                A(K,J) = W
             END IF
*             GAUSSIAN ELIMINATION
             T = A(K,J)
             DO 140 I = K+1, N
  140          A(I,J) = A(I,J)+WK(I)*T
  130      CONTINUE
*             MATRIX IS SINGULAR.
        ELSE
           IER = 1
           IP(K) = K
           DO 150 I = K+1, N
  150        A(I,K) = 0.0D0
           WRITE(*,*) '  (SUBR. GLU1)  MATRIX IS SINGULAR AT K =', K
           RETURN
        END IF
  100 CONTINUE
      RETURN
      END
*
      SUBROUTINE GSLV1( A, N, N1, B, IP )
*
*        GSLV1
*               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1
*
*               SOLVES SIMULTANEOUS LINEAR EQUATIONS
*               BY GAUSSIAN ELIMINATION METHOD.
*
*        INPUT - -
*             A(N1,N)  R *8  : RESULT OF GAUSSIAN ELIMINATION.
*             N        I *4  : ORDER OF MATRIX.
*             N1       I *4  : SIZE OF ARRAY A.
*             B(N)     R *8  : 1-DIM. ARRAY CONTAINING THE RIGHT HAND
*                              SIDE VECTOR.
*             IP(N)    I *4  : PIVOT NUMBER.
*        OUTPUT - -
*             B(N)           : SOLUTION.
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N1,*), B(*), IP(*)
*             FORWARD ELIMINATION PROCESS
      DO 100 K = 1, N
        IF( IP(K).NE.K ) THEN
           W = B(IP(K))
           B(IP(K)) = B(K)
           B(K) = W
        END IF
*             GAUSSIAN ELIMINATION
        T = B(K)
        DO 110 I = K+1, N
  110     B(I) = B(I)+A(I,K)*T
  100 CONTINUE
*             BACKWARD SUBSTITUTION PROCESS
      B(N) = B(N)/A(N,N)
      DO 200 K = N-1, 1, -1
        T = B(K+1)
        DO 210 I = 1, K
  210     B(I) = B(I)-A(I,K+1)*T
        B(K) = B(K)/A(K,K)
  200 CONTINUE
      RETURN
      END
