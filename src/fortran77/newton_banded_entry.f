C=======================================================================
C  BANDED BACKEND ENTRY POINT
C  Link this file instead of newton_solver.f.
C=======================================================================

      SUBROUTINE NEWTON_SOLVE(N,U,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &                        DX,DT,THETA,D,A,T,Z,RS,MU,EPS,
     &                        RTOL,STOL,MAXIT,VERBOSE,INFO,NITER)
      IMPLICIT NONE
      INTEGER N,MAXIT,VERBOSE,INFO,NITER
      DOUBLE PRECISION U(*),PO(0:*),MO(0:*)
      DOUBLE PRECISION MINO,MOUTO,MOUTN,PINN,DX,DT,THETA
      DOUBLE PRECISION D,A,T,Z,RS,MU,EPS,RTOL,STOL

      CALL NEWTON_SOLVE_BANDED(N,U,PO,MO,MINO,MOUTO,MOUTN,
     &     PINN,DX,DT,THETA,D,A,T,Z,RS,MU,EPS,
     &     RTOL,STOL,MAXIT,VERBOSE,INFO,NITER)

      RETURN
      END
