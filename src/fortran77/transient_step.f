C=======================================================================
C  ONE IMPLICIT THETA-METHOD TIME STEP
C=======================================================================

      SUBROUTINE TRANSIENT_STEP(N,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &                          DX,DT,THETA,D,A,T,Z,RS,MU,EPS,
     &                          RTOL,STOL,MAXIT,VERBOSE,U,
     &                          INFO,NITER)
      IMPLICIT NONE
      INTEGER N,MAXIT,VERBOSE,INFO,NITER,I
      DOUBLE PRECISION PO(0:*),MO(0:*),U(*)
      DOUBLE PRECISION MINO,MOUTO,MOUTN,PINN,DX,DT,THETA
      DOUBLE PRECISION D,A,T,Z,RS,MU,EPS,RTOL,STOL

C     Old state is the initial Newton guess.  The prescribed new inlet
C     pressure is inserted explicitly; all other new-time quantities
C     begin from their old-time values.

      U(1) = PINN
      U(2) = MINO
      DO 10 I=0,N-1
         U(2*I+3) = MO(I)
         U(2*I+4) = PO(I+1)
   10 CONTINUE

      CALL NEWTON_SOLVE(N,U,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &     DX,DT,THETA,D,A,T,Z,RS,MU,EPS,RTOL,STOL,MAXIT,
     &     VERBOSE,INFO,NITER)

      RETURN
      END
