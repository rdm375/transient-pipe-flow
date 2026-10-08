C=======================================================================
C  ONE IMPLICIT THETA-METHOD TIME STEP
C=======================================================================

      SUBROUTINE TRANSIENT_STEP(N,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &                          DX,DT,THETA,D,A,T,Z,RS,MU,EPS,
     &                          RTOL,STOL,MAXIT,VERBOSE,U,
     &                          INFO,NITER)
      IMPLICIT NONE
      INTEGER N,MAXIT,VERBOSE,INFO,NITER,I
      DOUBLE PRECISION GUESS(202)
      DOUBLE PRECISION PO(0:*),MO(0:*),U(*)
      DOUBLE PRECISION MINO,MOUTO,MOUTN,PINN,DX,DT,THETA
      DOUBLE PRECISION D,A,T,Z,RS,MU,EPS,RTOL,STOL

C     Old state is the initial Newton guess.  The prescribed new inlet
C     pressure is inserted explicitly; all other new-time quantities
C     begin from their old-time values.

      IF (N.LT.2.OR.N.GT.100) THEN
         INFO=3
         NITER=0
         RETURN
      ENDIF

      GUESS(1) = PINN
      GUESS(2) = MINO
      DO 10 I=0,N-1
         GUESS(2*I+3) = MO(I)
         GUESS(2*I+4) = PO(I+1)
   10 CONTINUE

      CALL NEWTON_SOLVE(N,GUESS,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &     DX,DT,THETA,D,A,T,Z,RS,MU,EPS,RTOL,STOL,MAXIT,
     &     VERBOSE,INFO,NITER)

      IF (INFO.EQ.0) THEN
         DO 30 I=1,2*N+2
            U(I)=GUESS(I)
   30    CONTINUE
      ENDIF
      RETURN
      END
