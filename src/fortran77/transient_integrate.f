C=======================================================================
C  MULTI-STEP TRANSIENT INTEGRATION
C=======================================================================

      SUBROUTINE INTEGRATE_TRANSIENT(N,PO,MO,MINO,PIN,MDOT0,
     &                               MDOT1,TRAMP,DT,NSTEPS,THETA,
     &                               DX,D,A,T,Z,RS,MU,EPS,RTOL,
     &                               STOL,MAXIT,TIME,HMIN,HOUT,
     &                               HPOUT,HLINE,HNIT,HMAXBAL,
     &                               INFO)
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      IMPLICIT NONE
      INTEGER N,NSTEPS,MAXIT,INFO,K,I,NITER
      DOUBLE PRECISION PO(0:*),MO(0:*),TIME(0:*),HMIN(0:*)
      DOUBLE PRECISION HOUT(0:*),HPOUT(0:*),HLINE(0:*)
      DOUBLE PRECISION HMAXBAL
      INTEGER HNIT(0:*)
      DOUBLE PRECISION PIN,MDOT0,MDOT1,TRAMP,DT,THETA,DX
      DOUBLE PRECISION D,A,T,Z,RS,MU,EPS,RTOL,STOL,MINO
      DOUBLE PRECISION U(202),PNEW(0:100),MNEW(0:99)
      DOUBLE PRECISION MOUTO,MOUTN,TNEW,BALERR,DEMAND_RAMP

      LOGICAL PIPE_VALID
      EXTERNAL PIPE_VALID
      IF (N.LT.2.OR.N.GT.100.OR.NSTEPS.LT.0) THEN
         INFO = 3
         RETURN
      ENDIF

      IF (.NOT.PIPE_VALID(N,PO,MO,MINO,DX,DT,THETA,
     & D,A,T,Z,RS,MU,EPS,RTOL,STOL,MAXIT).OR.
     & .NOT.IEEE_IS_FINITE(PIN).OR.
     & .NOT.IEEE_IS_FINITE(MDOT0).OR.
     & .NOT.IEEE_IS_FINITE(MDOT1).OR.
     & .NOT.IEEE_IS_FINITE(TRAMP)) THEN
         INFO=3
         RETURN
      ENDIF
      IF (PIN.LE.0D0.OR.TRAMP.LT.0D0) THEN
         INFO=3
         RETURN
      ENDIF
      MOUTO = MDOT0
      TIME(0) = 0.0D0
      HMIN(0) = MINO
      HOUT(0) = MOUTO
      HPOUT(0) = PO(N)
      CALL LINEPACK_CZ(N,PO,DX,A,T,Z,RS,HLINE(0))
      HNIT(0) = 0
      HMAXBAL = 0.0D0

      DO 100 K=1,NSTEPS
         TNEW = DBLE(K)*DT
         MOUTN = DEMAND_RAMP(TNEW,MDOT0,MDOT1,TRAMP)
         CALL TRANSIENT_STEP(N,PO,MO,MINO,MOUTO,MOUTN,PIN,
     &        DX,DT,THETA,D,A,T,Z,RS,MU,EPS,RTOL,STOL,MAXIT,0,
     &        U,INFO,NITER)
         IF (INFO .NE. 0) RETURN

         PNEW(0) = U(1)
         DO 10 I=0,N-1
            MNEW(I) = U(2*I+3)
            PNEW(I+1) = U(2*I+4)
   10    CONTINUE

         TIME(K) = TNEW
         HMIN(K) = U(2)
         HOUT(K) = MOUTN
         HPOUT(K) = PNEW(N)
         CALL LINEPACK_CZ(N,PNEW,DX,A,T,Z,RS,HLINE(K))
         HNIT(K) = NITER
         BALERR = DABS((HLINE(K)-HLINE(K-1))/DT
     &      -(THETA*(HMIN(K)-HOUT(K))
     &      +(1.0D0-THETA)*(HMIN(K-1)-HOUT(K-1))))
         HMAXBAL = DMAX1(HMAXBAL,BALERR)

         DO 20 I=0,N
            PO(I) = PNEW(I)
   20    CONTINUE
         DO 30 I=0,N-1
            MO(I) = MNEW(I)
   30    CONTINUE
         MINO = HMIN(K)
         MOUTO = MOUTN
  100 CONTINUE

      INFO = 0
      RETURN
      END

C=======================================================================
C  LINEAR DOWNSTREAM DEMAND RAMP FOLLOWED BY A CONSTANT HOLD.
C=======================================================================

      DOUBLE PRECISION FUNCTION DEMAND_RAMP(TIME,MDOT0,MDOT1,TRAMP)
      IMPLICIT NONE
      DOUBLE PRECISION TIME,MDOT0,MDOT1,TRAMP,FRACTION

      IF (TRAMP .LE. 0.0D0 .OR. TIME .GE. TRAMP) THEN
         DEMAND_RAMP = MDOT1
      ELSE IF (TIME .LE. 0.0D0) THEN
         DEMAND_RAMP = MDOT0
      ELSE
         FRACTION = TIME/TRAMP
         DEMAND_RAMP = MDOT0 + FRACTION*(MDOT1-MDOT0)
      ENDIF
      RETURN
      END

C=======================================================================
C  CONSTANT-Z TRAPEZOIDAL LINEPACK DIAGNOSTIC.
C=======================================================================

      SUBROUTINE LINEPACK_CZ(N,P,DX,A,T,Z,RS,MASS)
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION P(0:*),DX,A,T,Z,RS,MASS,W
      DOUBLE PRECISION EOS_CZ_RHO

      MASS = 0.0D0
      DO 10 I=0,N
         W = 1.0D0
         IF (I .EQ. 0 .OR. I .EQ. N) W = 0.5D0
         MASS = MASS + W*EOS_CZ_RHO(P(I),T,Z,RS)
   10 CONTINUE
      MASS = A*DX*MASS
      RETURN
      END
