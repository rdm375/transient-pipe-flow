C=======================================================================
C  DAMPED NEWTON SOLVER FOR ONE IMPLICIT TRANSIENT SYSTEM
C=======================================================================

      SUBROUTINE NEWTON_SOLVE(N,U,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &                        DX,DT,THETA,D,A,T,Z,RS,MU,EPS,
     &                        RTOL,STOL,MAXIT,VERBOSE,INFO,NITER)
      IMPLICIT NONE
      INTEGER N,MAXIT,VERBOSE,INFO,NITER
      INTEGER NU,LDJ,I,J,LS,LINFO,ITER
      PARAMETER (LDJ=202)
      DOUBLE PRECISION U(*),PO(0:*),MO(0:*)
      DOUBLE PRECISION MINO,MOUTO,MOUTN,PINN,DX,DT,THETA
      DOUBLE PRECISION D,A,T,Z,RS,MU,EPS,RTOL,STOL
      DOUBLE PRECISION R(202),RT(202),JAC(LDJ,202),AC(LDJ,202)
      DOUBLE PRECISION RHS(202),DU(202),UT(202)
      SAVE JAC,AC
      DOUBLE PRECISION RNORM,RNEW,SNORM,USCALE,LAMBDA

      NU = 2*N + 2
      IF (NU .GT. 202) THEN
         INFO = 3
         NITER = 0
         RETURN
      ENDIF

      CALL ASSEMBLE_RESIDUAL(N,U,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &     DX,DT,THETA,D,A,T,Z,RS,MU,EPS,R)

      DO 100 ITER=0,MAXIT-1
         RNORM = 0.0D0
         DO 10 I=1,NU
            RNORM = DMAX1(RNORM,DABS(R(I)))
   10    CONTINUE
         IF (RNORM .LE. RTOL) THEN
            NITER = ITER
            INFO = 0
            RETURN
         ENDIF

         CALL ASSEMBLE_JACOBIAN(N,U,DX,DT,THETA,D,A,T,Z,RS,
     &        MU,EPS,JAC,LDJ)
         DO 30 J=1,NU
            DO 20 I=1,NU
               AC(I,J) = JAC(I,J)
   20       CONTINUE
   30    CONTINUE
         DO 40 I=1,NU
            RHS(I) = -R(I)
   40    CONTINUE
         CALL SOLVE_DENSE(NU,AC,LDJ,RHS,DU,LINFO)
         IF (LINFO .NE. 0) THEN
            INFO = 2
            RETURN
         ENDIF

         SNORM = 0.0D0
         USCALE = 1.0D0
         DO 50 I=1,NU
            SNORM = DMAX1(SNORM,DABS(DU(I)))
            USCALE = DMAX1(USCALE,DABS(U(I)))
   50    CONTINUE

         LAMBDA = 1.0D0
         DO 70 LS=1,20
            DO 60 I=1,NU
               UT(I) = U(I) + LAMBDA*DU(I)
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
         INFO = 4
         RETURN

   80    CONTINUE
         IF (VERBOSE .NE. 0) THEN
            WRITE(*,'(I4,3(1X,1PE12.4))') ITER,RNORM,SNORM,
     &           LAMBDA
         ENDIF
         DO 90 I=1,NU
            U(I) = UT(I)
            R(I) = RT(I)
   90    CONTINUE
         IF (LAMBDA*SNORM .LE. STOL*USCALE .AND.
     &       RNEW .LE. 10.0D0*RTOL) THEN
            NITER = ITER + 1
            INFO = 0
            RETURN
         ENDIF
  100 CONTINUE

      NITER = MAXIT
      INFO = 1
      RETURN
      END
