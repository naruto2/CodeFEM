*
      PROGRAM ORDRNG
***********************************************************************
*  SAMPLE PROGRAM FOR MARKWZ AND STWART ROUTINE.          VER. 1      *
***********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       INTEGER*4 R, C
       DIMENSION A(100),IA(100),IR(20),R(20),C(20),IP(100),JCOL(20),
     *  JROW(20),MA(0:20),IC(20),E(200),IE(200),W(20),Z(20),ME(0:20)
     *  ,JP(20),IW(20)
       DATA IA/2,10,11,13,15,16,4,5,10,16,9,11,12,1,5,7,8,1,2,4,5,6,9,
     *  10,4,9,12,2,3,4,14,1,2,9,16,4,15,6,7,8,14,15,2,16,9,11,12,6,14,
     *  15,15,4,6,12,15,2,4,10,13,41*0/
       DATA MA/0,6,10,13,17,24,27,31,35,37,42,44,47,50,51,55,59,4*0/
       DATA M,L,THETA,IP/16,100,1.0D-12,100*0/
C
      IS=0
    5 CONTINUE
      NUM=0
      DO 20 I=1,M
       DO 20 J=MA(I-1)+1,MA(I)
        NUM=NUM+1
        A(NUM)=DFLOAT(NUM)
   20  CONTINUE
      IF (IS .EQ. 0) THEN
       CALL MARKWZ(A,IA,L,MA,M,E,IE,IR,IC,R,C,W,Z,ME,NUME,MEND,KERNS,
     *             JROW,JCOL,IER)
       WRITE(*,*) NUM,NUME,MEND,KERNS-1
       WRITE(*,*) (ME(I),I=1,M)
      ELSE
       CALL STWART(IA,L,MA,M,R,C,IR,IC,JROW,JCOL,IP,JP,KERNS,MEND,
     *             IW,LG,IER)
       WRITE(*,*) 'STEWART METHOD. ', MEND, KERNS - 1, LG
       WRITE(*,*) (IP(I),I=1,M)
       WRITE(*,*) (JP(I),I=1,M)
       WRITE(*,*) (IW(I),I=1,M)
      ENDIF 
      WRITE(*,*) (JCOL(I),I=1,M)
      WRITE(*,*) (JROW(I),I=1,M)
      WRITE(*,*) (IC(I),I=1,M)
      WRITE(*,*) (IR(I),I=1,M)
      IF (IS .EQ. 0) THEN
       IS=1
       GO TO 5
      ENDIF 
      STOP
      END
