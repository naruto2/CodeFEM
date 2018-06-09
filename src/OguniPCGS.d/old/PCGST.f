*
      PROGRAM PCGSTS
***********************************************************************
*   SAMPLE PROGRAM FOR ROUTINE PCGS.                VER. 1            *
***********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION A(100,20),B(100),X(0:100),DD(0:100),W(0:100),D(100)
     *  ,P(0:100),Q(0:100),R(0:100),R0(100),E(100),H(100),M(200)
     *  ,IA(100,20) 
C
      N=10
      NL=3
      S=0.0D0
      L=100
    5 ITR=50
      EPS=1.0D-7
      DO 10 J=1,2*NL
       DO 10 I=1,N
        IA(I,J)=0
   10   A(I,J)=0.0D0
      DO 15 I=1,N
   15  X(I)=0.0D0
      B(1)=20.0D0
      B(2)=23.0D0
      B(3)=23.0D0
      B(4)=25.0D0
      B(5)=24.0D0
      B(6)=24.0D0
      B(7)=25.0D0
      B(8)=23.0D0
      B(9)=23.0D0
      B(10)=20.0D0
      DO 20 I=1,N
   20  D(I)=14.0D0
      DO 30 I=2,N
       A(I-1,4)=3.0D0
       IA(I-1,4)=I
       IA(I,1)=I-1
   30  A(I,1)=3.0D0
      DO 40 I=4,N
       A(I-3,5)=2.0D0
       IA(I-3,5)=I
       IA(I,2)=I-3
   40  A(I,2)=2.0D0
      DO 50 I=7,N
       A(I-6,6)=1.0D0
       IA(I-6,6)=I 
       IA(I,3)=I-6
   50  A(I,3)=1.0D0
      CALL PCGS(D,A,IA,N,L,NL,B,EPS,ITR,S,X,DD,P,Q,R,R0,E,H,W,M,IER)
      WRITE(*,*) 'EXAMPLE OF PCGS. ', ITR, EPS, S
      WRITE(*,*) (X(I),I=1,N)
      IF ( S .EQ. 0.0) THEN
       S=0.95D0
       GO TO 5
      ENDIF 
      IF (S .LT. 1.0) THEN
       S=1.02D0
       GO TO 5
      ENDIF
      IF (S .GT. 1.0) THEN
       S = 1.0D0
       GO TO 5
      ENDIF
      STOP
      END
