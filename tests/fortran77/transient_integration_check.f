      PROGRAM TRANSIENTINTEGRATIONCHECK
      IMPLICIT NONE
      INTEGER N,NSTEPS,I,K,INFO,MAXIT
      PARAMETER (N=20,NSTEPS=360)
      DOUBLE PRECISION PO(0:N),MO(0:N-1),P0(0:N),P1(0:N)
      DOUBLE PRECISION TIME(0:NSTEPS),HMIN(0:NSTEPS)
      DOUBLE PRECISION HOUT(0:NSTEPS),HPOUT(0:NSTEPS)
      DOUBLE PRECISION HLINE(0:NSTEPS)
      INTEGER HNIT(0:NSTEPS)
      DOUBLE PRECISION PI,D,A,T,Z,RS,MU,EPS,L,DX,DT,THETA
      DOUBLE PRECISION PIN,MDOT0,MDOT1,TRAMP,MINO,RE,FD0,FD1
      DOUBLE PRECISION COEF0,COEF1,X,RTOL,STOL,HMAXBAL
      DOUBLE PRECISION CUMBOUND,CUMERR,RELERR,FINALERR,MAXPERR
      DOUBLE PRECISION MAXMERR,FINALMASS,EXPECTED,DEMAND_RAMP
      DOUBLE PRECISION REYNOLDS_MASS,FRICTION_SJ

      PI = 4.0D0*DATAN(1.0D0)
      D = 1.0D0
      A = PI*D*D/4.0D0
      T = 288.15D0
      Z = 0.90D0
      RS = 500.0D0
      MU = 1.1D-5
      EPS = 4.5D-5
      L = 100000.0D0
      DX = L/DBLE(N)
      DT = 60.0D0
      THETA = 0.65D0
      PIN = 8.0D6
      MDOT0 = 100.0D0
      MDOT1 = 105.0D0
      TRAMP = 600.0D0
      MINO = MDOT0
      RTOL = 1.0D-11
      STOL = 1.0D-12
      MAXIT = 30

      RE = REYNOLDS_MASS(MDOT0,D,A,MU)
      FD0 = FRICTION_SJ(RE,EPS,D)
      RE = REYNOLDS_MASS(MDOT1,D,A,MU)
      FD1 = FRICTION_SJ(RE,EPS,D)
      COEF0 = FD0*Z*RS*T*MDOT0*MDOT0/(D*A*A)
      COEF1 = FD1*Z*RS*T*MDOT1*MDOT1/(D*A*A)
      DO 10 I=0,N
         X = DBLE(I)*DX
         P0(I) = DSQRT(PIN*PIN-COEF0*X)
         P1(I) = DSQRT(PIN*PIN-COEF1*X)
         PO(I) = P0(I)
   10 CONTINUE
      DO 20 I=0,N-1
         MO(I) = MDOT0
   20 CONTINUE

      WRITE(*,*)
      WRITE(*,*) 'MULTI-STEP TRANSIENT INTEGRATION CHECK'
      WRITE(*,*) '--------------------------------------'
      WRITE(*,'(A,I5,A,F8.1,A,F6.2)') ' steps=',NSTEPS,
     &     ' dt=',DT,' theta=',THETA
      WRITE(*,'(A,F8.1,A,F8.1,A,F8.1)') ' demand ',MDOT0,
     &     ' -> ',MDOT1,' ramp [s]=',TRAMP

      CALL INTEGRATE_TRANSIENT(N,PO,MO,MINO,PIN,MDOT0,MDOT1,
     &     TRAMP,DT,NSTEPS,THETA,DX,D,A,T,Z,RS,MU,EPS,RTOL,
     &     STOL,MAXIT,TIME,HMIN,HOUT,HPOUT,HLINE,HNIT,HMAXBAL,
     &     INFO)
      IF (INFO .NE. 0) THEN
         WRITE(*,*) 'FAIL: integration info=',INFO
         STOP 1
      ENDIF

C     Independent cumulative theta-weighted boundary mass integral.

      CUMBOUND = 0.0D0
      DO 30 K=1,NSTEPS
         CUMBOUND = CUMBOUND + DT*(THETA*(HMIN(K)-HOUT(K))
     &      +(1.0D0-THETA)*(HMIN(K-1)-HOUT(K-1)))
         EXPECTED = DEMAND_RAMP(TIME(K),MDOT0,MDOT1,TRAMP)
         IF (DABS(HOUT(K)-EXPECTED) .GT. 1.0D-12) THEN
            WRITE(*,*) 'FAIL: demand history at step ',K
            STOP 1
         ENDIF
   30 CONTINUE
      CUMERR = (HLINE(NSTEPS)-HLINE(0))-CUMBOUND
      RELERR = DABS(CUMERR)/DMAX1(1.0D0,DABS(CUMBOUND),
     &         DABS(HLINE(NSTEPS)-HLINE(0)))

C     Compare final state with independent analytic 105 kg/s state.

      MAXPERR = 0.0D0
      DO 40 I=0,N
         MAXPERR = DMAX1(MAXPERR,DABS(PO(I)-P1(I))/P1(I))
   40 CONTINUE
      MAXMERR = DABS(MINO-MDOT1)/MDOT1
      DO 50 I=0,N-1
         MAXMERR = DMAX1(MAXMERR,DABS(MO(I)-MDOT1)/MDOT1)
   50 CONTINUE
      FINALERR = DMAX1(MAXPERR,MAXMERR)
      FINALMASS = HLINE(NSTEPS)-HLINE(0)

      WRITE(*,'(A,1PE14.6)') ' initial outlet p [Pa] = ',P0(N)
      WRITE(*,'(A,1PE14.6)') ' final outlet p [Pa]   = ',HPOUT(NSTEPS)
      WRITE(*,'(A,1PE14.6)') ' target outlet p [Pa]  = ',P1(N)
      WRITE(*,'(A,1PE14.6)') ' final inlet flow      = ',HMIN(NSTEPS)
      WRITE(*,'(A,1PE14.6)') ' total linepack change = ',FINALMASS
      WRITE(*,'(A,1PE14.6)') ' boundary mass integral= ',CUMBOUND
      WRITE(*,'(A,1PE14.6)') ' cumulative rel error  = ',RELERR
      WRITE(*,'(A,1PE14.6)') ' max step balance abs  = ',HMAXBAL
      WRITE(*,'(A,1PE14.6)') ' final steady rel error= ',FINALERR
      WRITE(*,'(A,I4)') ' final Newton iterations= ',HNIT(NSTEPS)

      IF (HPOUT(NSTEPS) .GE. HPOUT(0)) THEN
         WRITE(*,*) 'FAIL: outlet pressure did not decrease'
         STOP 1
      ENDIF
      IF (HLINE(NSTEPS) .GE. HLINE(0)) THEN
         WRITE(*,*) 'FAIL: linepack did not decrease'
         STOP 1
      ENDIF
      IF (RELERR .GT. 5.0D-10) THEN
         WRITE(*,*) 'FAIL: cumulative mass conservation'
         STOP 1
      ENDIF
      IF (FINALERR .GT. 5.0D-7) THEN
         WRITE(*,*) 'FAIL: did not relax to final steady state'
         STOP 1
      ENDIF

      WRITE(*,*) 'PASS: multi-step transient integration verified'
      END
