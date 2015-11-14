C
C @(#)blas.f	3.4 (BNP) 12/9/88; from blas.f 2.2 10/13/87
C
C       (standard double precision BLAS)
C
        SUBROUTINE DATX(N,DA,DX,INCX,DY,INCY)
C
C       dy := da*dx
C
        DOUBLE PRECISION DX(*),DY(*),DA
        INTEGER I,INCX,INCY,IX,IY,M,MP1,N
        IF (N.LE.0) RETURN
        IF (DA.EQ.0.0D0) RETURN
        IF (INCX.EQ.1.AND.INCY.EQ.1) GO TO 20
C
C       unequal increments or equal increments .ne. one
C
        IX = 1
        IY = 1
        IF (INCX.LT.0) IX = (-N+1)*INCX+1
        IF (INCY.LT.0) IY = (-N+1)*INCY+1
        DO 10 I = 1,N
          DY(IY) = DA*DX(IX)
          IX = IX+INCX
          IY = IY+INCY
10      CONTINUE
        RETURN
C
C       code for both increments equal to 1
C
20      M = MOD(N,4)
        IF (M.EQ.0) GO TO 40
        DO 30 I = 1,M
          DY(I) = DA*DX(I)
30      CONTINUE
        IF (N.LT.4) RETURN
40      MP1 = M+1
        DO 50 I = MP1,N,4
          DY(I) = DA*DX(I)
          DY(I+1) = DA*DX(I+1)
          DY(I+2) = DA*DX(I+2)
          DY(I+3) = DA*DX(I+3)
50      CONTINUE
        RETURN
        END
C
        SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C
C       dy := dy+da*dx
C
        DOUBLE PRECISION DX(*),DY(*),DA
        INTEGER I,INCX,INCY,IX,IY,M,MP1,N
        IF (N.LE.0) RETURN
        IF (DA.EQ.0.0D0) RETURN
        IF (INCX.EQ.1.AND.INCY.EQ.1) GO TO 20
C
C       unequal increments or equal increments .ne. one
C
        IX = 1
        IY = 1
        IF (INCX.LT.0) IX = (-N+1)*INCX+1
        IF (INCY.LT.0) IY = (-N+1)*INCY+1
        DO 10 I = 1,N
          DY(IY) = DY(IY)+DA*DX(IX)
          IX = IX+INCX
          IY = IY+INCY
10      CONTINUE
        RETURN
C
C       code for both increments equal to 1
C
20      M = MOD(N,4)
        IF (M.EQ.0) GO TO 40
        DO 30 I = 1,M
          DY(I) = DY(I)+DA*DX(I)
30      CONTINUE
        IF (N.LT.4) RETURN
40      MP1 = M+1
        DO 50 I = MP1,N,4
          DY(I) = DY(I)+DA*DX(I)
          DY(I+1) = DY(I+1)+DA*DX(I+1)
          DY(I+2) = DY(I+2)+DA*DX(I+2)
          DY(I+3) = DY(I+3)+DA*DX(I+3)
50      CONTINUE
        RETURN
        END
C
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
C
C     copy double precision dx to double precision dy.
C     for i = 0 to n-1, copy dx(lx+i*incx) to dy(ly+i*incy),
C     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
C     defined in a similar way using incy.
C
      INTEGER N,INCX,INCY,IX,IY,I,M,MP1,NS
      DOUBLE PRECISION DX(*),DY(*)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.INCY) IF (INCX-1) 5,20,60
    5 CONTINUE
C
C        code for unequal or nonpositive increments.
C
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX+1
      IF (INCY.LT.0) IY = (-N+1)*INCY+1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX+INCX
        IY = IY+INCY
   10 CONTINUE
      RETURN
C
C        code for both increments equal to 1
C
C        clean-up loop so remaining vector length is a multiple of 7.
C
   20 M = MOD(N,7)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF (N.LT.7) RETURN
   40 MP1 = M+1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I+1) = DX(I+1)
        DY(I+2) = DX(I+2)
        DY(I+3) = DX(I+3)
        DY(I+4) = DX(I+4)
        DY(I+5) = DX(I+5)
        DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
C
C        code for equal, positive, nonunit increments.
C
   60 CONTINUE
      NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DX(I)
   70 CONTINUE
      RETURN
      END
C
        DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C       forms the dot product of two vectors
C
        DOUBLE PRECISION DX(*),DY(*),STEMP
        INTEGER I,INCX,INCY,IX,IY,M,MP1,N
        STEMP = 0.0D0
        DDOT = 0.0D0
        IF (N.LE.0) RETURN
        IF (INCX.EQ.1.AND.INCY.EQ.1) GO TO 20
C
C       unequal increments or increments .ne. 1
C
        IX = 1
        IY = 1
        IF (INCX.LT.0) IX = (-N+1)*INCX+1
        IF (INCY.LT.0) IY = (-N+1)*INCY+1
        DO 10 I = 1,N
          STEMP = STEMP+DX(IX)*DY(IY)
          IX = IX+INCX
          IY = IY+INCY
10      CONTINUE
        DDOT = STEMP
        RETURN
C
C       both increments equal 1
C
20      M = MOD(N,5)
        IF (M.EQ.0) GO TO 40
        DO 30 I = 1,M
          STEMP = STEMP+DX(I)*DY(I)
30      CONTINUE
        IF (N.LT.5) GO TO 60
40      MP1 = M+1
        DO 50 I = MP1,N,5
          STEMP = STEMP+DX(I)*DY(I)+DX(I+1)*DY(I+1)+
     *      DX(I+2)*DY(I+2)+DX(I+3)*DY(I+3)+DX(I+4)*DY(I+4)
50      CONTINUE
60      DDOT = STEMP
        RETURN
        END
C
        SUBROUTINE DSCAL(N,DA,DX,INCX)
C
C       scales a vector by a scalar
C
        DOUBLE PRECISION DA,DX(*)
        INTEGER I,INCX,M,MP1,N,NINCX
        IF (N.LE.0) RETURN
        IF (INCX.EQ.1) GO TO 20
C
C       increment not equal to 1
C
        NINCX = N*INCX
        DO 10 I = 1,NINCX,INCX
          DX(I) = DA*DX(I)
10      CONTINUE
        RETURN
C
C       increment equal 1
C
20      M = MOD(N,5)
        IF (M.EQ.0) GO TO 40
        DO 30 I = 1,M
          DX(I) = DA*DX(I)
30      CONTINUE
        IF (N.LT.5) RETURN
40      MP1 = M+1
        DO 50 I = MP1,N,5
          DX(I) = DA*DX(I)
          DX(I+1) = DA*DX(I+1)
          DX(I+2) = DA*DX(I+2)
          DX(I+3) = DA*DX(I+3)
          DX(I+4) = DA*DX(I+4)
50      CONTINUE
        RETURN
        END
C
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
C
C     interchanges two vectors.
C     uses unrolled loops for increments equal one.
C     jack dongarra, linpack, 3/11/78.
C
      DOUBLE PRECISION DX(*),DY(*),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1.AND.INCY.EQ.1) GO TO 20
C
C       code for unequal increments or equal increments not equal
C         to 1
C
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX+1
      IF (INCY.LT.0) IY = (-N+1)*INCY+1
      DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX+INCX
        IY = IY+INCY
   10 CONTINUE
      RETURN
C
C       code for both increments equal to 1
C
C       clean-up loop
C
   20 M = MOD(N,3)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
   30 CONTINUE
      IF (N.LT.3) RETURN
   40 MP1 = M+1
      DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I+1)
        DX(I+1) = DY(I+1)
        DY(I+1) = DTEMP
        DTEMP = DX(I+2)
        DX(I+2) = DY(I+2)
        DY(I+2) = DTEMP
   50 CONTINUE
      RETURN
      END
C
        INTEGER FUNCTION IDAMAX(N,DX,INCX)
C
C       Finds the index of the element having maximu. absolute value
C       cf LINPACK
C
        DOUBLE PRECISION DX(*),SMAX
        INTEGER I,INCX,IX,N
        IDAMAX = 0
        IF (N.LT.1) RETURN
        IDAMAX = 1
        IF (N.EQ.1) RETURN
        IF (INCX.EQ.1) GO TO 20
C
C       code for increment not equal to 1
C
        IX = 1
        SMAX = ABS(DX(1))
        IX = IX+INCX
        DO 10 I = 2,N
          IF (ABS(DX(IX)).LE.SMAX) GO TO 5
          IDAMAX = I
          SMAX = ABS(DX(IX))
5         IX = IX+INCX
10      CONTINUE
C
C       code for increment equal to 1
C
20      SMAX = ABS(DX(1))
        DO 30 I = 2,N
          IF (ABS(DX(I)).LE.SMAX) GO TO 30
          IDAMAX = I
          SMAX = ABS(DX(I))
30      CONTINUE
        RETURN
        END
