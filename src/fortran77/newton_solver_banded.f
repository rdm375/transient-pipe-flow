C=======================================================================
C  DAMPED NEWTON SOLVER FOR ONE IMPLICIT TRANSIENT SYSTEM
C=======================================================================

      SUBROUTINE NEWTON_SOLVE_BANDED(N,U,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &                        DX,DT,THETA,D,A,T,Z,RS,MU,EPS,
     &                        RTOL,STOL,MAXIT,VERBOSE,INFO,NITER)
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      IMPLICIT NONE
      INTEGER N,MAXIT,VERBOSE,INFO,NITER
      INTEGER NU,I,LS,LINFO,ITER
      INTEGER KL,KU,LDAB
      PARAMETER (KL=2,KU=1,LDAB=6)
      DOUBLE PRECISION U(*),PO(0:*),MO(0:*)
      DOUBLE PRECISION MINO,MOUTO,MOUTN,PINN,DX,DT,THETA
      DOUBLE PRECISION D,A,T,Z,RS,MU,EPS,RTOL,STOL
      DOUBLE PRECISION R(202),RT(202),AB(LDAB,202)
      DOUBLE PRECISION RHS(202),DU(202),UT(202),WORK(202)
      DOUBLE PRECISION RNORM,RNEW,SNORM,USCALE,LAMBDA

      LOGICAL PIPE_VALID
      EXTERNAL PIPE_VALID
      NITER = 0
      IF (N.LT.2.OR.N.GT.100) THEN
         INFO = 3
         NITER = 0
         RETURN
      ENDIF
      NU = 2*N + 2

      IF (.NOT.PIPE_VALID(N,PO,MO,MINO,DX,DT,THETA,
     & D,A,T,Z,RS,MU,EPS,RTOL,STOL,MAXIT)) THEN
         INFO=3
         RETURN
      ENDIF
      IF (.NOT.IEEE_IS_FINITE(MOUTO).OR.
     &    .NOT.IEEE_IS_FINITE(MOUTN).OR.
     &    .NOT.IEEE_IS_FINITE(PINN)) THEN
         INFO=3
         RETURN
      ENDIF
      IF (PINN.LE.0D0) THEN
         INFO=3
         RETURN
      ENDIF
      DO 4 I=1,NU
         IF (.NOT.IEEE_IS_FINITE(U(I))) THEN
            INFO=3
            RETURN
         ENDIF
    4 CONTINUE
      DO 5 I=1,NU
         WORK(I) = U(I)
    5 CONTINUE

      CALL ASSEMBLE_RESIDUAL(N,WORK,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &     DX,DT,THETA,D,A,T,Z,RS,MU,EPS,R)

      DO 100 ITER=0,MAXIT-1
         RNORM = 0.0D0
         DO 10 I=1,NU
            RNORM = DMAX1(RNORM,DABS(R(I)))
   10    CONTINUE
         IF (RNORM .LE. RTOL) THEN
            NITER = ITER
            INFO = 0
            GOTO 200
         ENDIF

         CALL ASSEMBLE_JACOBIAN_BANDED(N,WORK,DX,DT,THETA,
     &        D,A,T,Z,RS,MU,EPS,AB,LDAB)
         DO 40 I=1,NU
            RHS(I) = -R(I)
   40    CONTINUE
         CALL SOLVE_BANDED(NU,AB,LDAB,KL,KU,RHS,DU,LINFO)
         IF (LINFO .NE. 0) THEN
            NITER = ITER + 1
            INFO = 2
            RETURN
         ENDIF

         SNORM = 0.0D0
         USCALE = 1.0D0
         DO 50 I=1,NU
            SNORM = DMAX1(SNORM,DABS(DU(I)))
            USCALE = DMAX1(USCALE,DABS(WORK(I)))
   50    CONTINUE

         LAMBDA = 1.0D0
         DO 70 LS=1,20
            DO 60 I=1,NU
               UT(I) = WORK(I) + LAMBDA*DU(I)
   60       CONTINUE
            CALL ASSEMBLE_RESIDUAL(N,UT,PO,MO,MINO,MOUTO,MOUTN,
     &           PINN,DX,DT,THETA,D,A,T,Z,RS,MU,EPS,RT)
            RNEW = 0.0D0
            DO 65 I=1,NU
               RNEW = DMAX1(RNEW,DABS(RT(I)))
   65       CONTINUE
            IF (RNEW .LT. RNORM) GOTO 80
            LAMBDA = 0.5D0*LAMBDA
   70    CONTINUE
         NITER = ITER + 1
         INFO = 4
         RETURN

   80    CONTINUE
         IF (VERBOSE .NE. 0) THEN
            WRITE(*,'(I4,3(1X,1PE12.4))') ITER,RNORM,SNORM,
     &           LAMBDA
         ENDIF
         DO 90 I=1,NU
            WORK(I) = UT(I)
            R(I) = RT(I)
   90    CONTINUE
         IF (LAMBDA*SNORM .LE. STOL*USCALE .AND.
     &       RNEW .LE. 10.0D0*RTOL) THEN
            NITER = ITER + 1
            INFO = 0
            GOTO 200
         ENDIF
  100 CONTINUE

      NITER = MAXIT
      INFO = 1
      RETURN

  200 CONTINUE
      DO 210 I=1,NU
         U(I) = WORK(I)
  210 CONTINUE
      RETURN
      END
