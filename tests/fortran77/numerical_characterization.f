      PROGRAM NUMERICALCHARACTERIZATION
      IMPLICIT NONE
      INTEGER MAXSTEP,NTEMP,NSPACE,NTHETA
      PARAMETER (MAXSTEP=3000,NTEMP=5,NSPACE=5,NTHETA=7)
      INTEGER I,INFO,NSTEP,NREF
      INTEGER NTV(NTEMP),NSV(NSPACE)
      DOUBLE PRECISION DTV(NTEMP),THV(NTHETA)
      DOUBLE PRECISION PREF,MINREF,LREF,PVAL,MINVAL,LVAL
      DOUBLE PRECISION ERRP,ERRM,ERRL,TEND,THETA,DT
      DOUBLE PRECISION ETP(NTEMP),ESP(NSPACE),ETH(NTHETA)
      DOUBLE PRECISION POUT(0:MAXSTEP),HMIN(0:MAXSTEP)
      DOUBLE PRECISION HLINE(0:MAXSTEP)
      DOUBLE PRECISION TIME(0:MAXSTEP)
      CHARACTER*80 FNAME

      DATA DTV /120.0D0,60.0D0,30.0D0,15.0D0,7.5D0/
      DATA NTV /15,30,60,120,240/
      DATA NSV /10,20,40,80,100/
      DATA THV /0.50D0,0.55D0,0.60D0,0.65D0,0.70D0,
     &          0.80D0,1.00D0/

C     Smooth 30-minute demand transition, observed at 30 minutes.
C     The cosine ramp avoids derivative corners inside the interval.

      TEND = 1800.0D0
      NREF = 480
      CALL RUN_CASE(100,3.75D0,NREF,0.50D0,TEND,MAXSTEP,
     &     TIME,POUT,HMIN,HLINE,PREF,MINREF,LREF,INFO)
      IF (INFO .NE. 0) STOP 1

      CALL SYSTEM('mkdir -p benchmarks/results-m6')

C     Temporal refinement at fixed N=40.  Use an independent fine
C     temporal reference on the same grid so spatial error cancels.

      CALL RUN_CASE(40,3.75D0,NREF,0.50D0,TEND,MAXSTEP,
     &     TIME,POUT,HMIN,HLINE,PREF,MINREF,LREF,INFO)
      IF (INFO .NE. 0) STOP 1
      FNAME = 'benchmarks/results-m6/temporal.csv'
      OPEN(10,FILE=FNAME,STATUS='REPLACE')
      WRITE(10,'(A)') 'n,dt,theta,pout,min,linepack,'//
     &                 'relerr_pout,relerr_min,relerr_linepack'
      DO 10 I=1,NTEMP
         CALL RUN_CASE(40,DTV(I),NTV(I),0.50D0,TEND,MAXSTEP,
     &        TIME,POUT,HMIN,HLINE,PVAL,MINVAL,LVAL,INFO)
         IF (INFO .NE. 0) STOP 1
         ERRP = DABS(PVAL-PREF)/DABS(PREF)
         ERRM = DABS(MINVAL-MINREF)/DMAX1(1.0D0,DABS(MINREF))
         ERRL = DABS(LVAL-LREF)/DABS(LREF)
         ETP(I) = ERRP
         WRITE(10,'(I4,A,F10.4,A,F6.2,6(A,1PE16.8))')
     &        40,',',DTV(I),',',0.50D0,',',PVAL,',',MINVAL,',',
     &        LVAL,',',ERRP,',',ERRM,',',ERRL
   10 CONTINUE
      CLOSE(10)

C     Spatial refinement at fixed dt=7.5 s, theta=0.5.
C     Use N=100 at the identical dt as the spatial reference so the
C     comparison does not mix temporal and spatial truncation error.

      CALL RUN_CASE(100,7.5D0,240,0.50D0,TEND,MAXSTEP,
     &     TIME,POUT,HMIN,HLINE,PREF,MINREF,LREF,INFO)
      IF (INFO .NE. 0) STOP 1
      FNAME = 'benchmarks/results-m6/spatial.csv'
      OPEN(11,FILE=FNAME,STATUS='REPLACE')
      WRITE(11,'(A)') 'n,dx,dt,theta,pout,min,linepack,'//
     &                 'relerr_pout,relerr_min,relerr_linepack'
      DO 20 I=1,NSPACE
         CALL RUN_CASE(NSV(I),7.5D0,240,0.50D0,TEND,MAXSTEP,
     &        TIME,POUT,HMIN,HLINE,PVAL,MINVAL,LVAL,INFO)
         IF (INFO .NE. 0) STOP 1
         ERRP = DABS(PVAL-PREF)/DABS(PREF)
         ERRM = DABS(MINVAL-MINREF)/DMAX1(1.0D0,DABS(MINREF))
         ERRL = DABS(LVAL-LREF)/DABS(LREF)
         ESP(I) = ERRP
         WRITE(11,'(I4,A,F10.2,A,F10.4,A,F6.2,6(A,1PE16.8))')
     &        NSV(I),',',100000.0D0/DBLE(NSV(I)),',',7.5D0,',',
     &        0.50D0,',',PVAL,',',MINVAL,',',LVAL,',',ERRP,',',
     &        ERRM,',',ERRL
   20 CONTINUE
      CLOSE(11)

C     Theta study at the established engineering dt=60 s, N=40.
C     Reference is fine-step Crank-Nicolson on the same grid.

      CALL RUN_CASE(40,3.75D0,NREF,0.50D0,TEND,MAXSTEP,
     &     TIME,POUT,HMIN,HLINE,PREF,MINREF,LREF,INFO)
      IF (INFO .NE. 0) STOP 1
      FNAME = 'benchmarks/results-m6/theta.csv'
      OPEN(12,FILE=FNAME,STATUS='REPLACE')
      WRITE(12,'(A)') 'n,dt,theta,pout,min,linepack,'//
     &                 'relerr_pout,relerr_min,relerr_linepack'
      DO 30 I=1,NTHETA
         THETA = THV(I)
         DT = 60.0D0
         NSTEP = 30
         CALL RUN_CASE(40,DT,NSTEP,THETA,TEND,MAXSTEP,
     &        TIME,POUT,HMIN,HLINE,PVAL,MINVAL,LVAL,INFO)
         IF (INFO .NE. 0) STOP 1
         ERRP = DABS(PVAL-PREF)/DABS(PREF)
         ERRM = DABS(MINVAL-MINREF)/DMAX1(1.0D0,DABS(MINREF))
         ERRL = DABS(LVAL-LREF)/DABS(LREF)
         ETH(I) = ERRP
         WRITE(12,'(I4,A,F10.4,A,F6.2,6(A,1PE16.8))')
     &        40,',',DT,',',THETA,',',PVAL,',',MINVAL,',',LVAL,
     &        ',',ERRP,',',ERRM,',',ERRL
   30 CONTINUE
      CLOSE(12)

C     Frozen M6 regression envelopes.  These are intentionally wider
C     than the measured baseline and detect numerical regressions, not
C     compiler-level roundoff changes.

      DO 35 I=1,NTEMP-1
         IF (ETP(I+1) .GE. ETP(I)) THEN
            WRITE(*,*) 'FAIL: temporal refinement is not monotone'
            STOP 1
         ENDIF
   35 CONTINUE
      DO 36 I=1,NSPACE-2
         IF (ESP(I+1) .GE. ESP(I)) THEN
            WRITE(*,*) 'FAIL: spatial refinement is not monotone'
            STOP 1
         ENDIF
   36 CONTINUE
      IF (ETP(2) .GT. 5.0D-7) THEN
         WRITE(*,*) 'FAIL: dt=60 temporal baseline tolerance'
         STOP 1
      ENDIF
      IF (ESP(3) .GT. 2.0D-7) THEN
         WRITE(*,*) 'FAIL: N=40 spatial baseline tolerance'
         STOP 1
      ENDIF
      IF (ETH(4) .GT. 7.0D-7) THEN
         WRITE(*,*) 'FAIL: theta=0.65 baseline tolerance'
         STOP 1
      ENDIF

      WRITE(*,*)
      WRITE(*,*) 'M6 NUMERICAL CHARACTERIZATION'
      WRITE(*,*) '-----------------------------'
      WRITE(*,*) 'reference: N=100 dt=3.75 theta=0.50'
      WRITE(*,'(A,1PE14.6)') ' reference outlet p = ',PREF
      WRITE(*,*) 'wrote benchmarks/results-m6/temporal.csv'
      WRITE(*,*) 'wrote benchmarks/results-m6/spatial.csv'
      WRITE(*,*) 'wrote benchmarks/results-m6/theta.csv'
      WRITE(*,*) 'PASS: M6 characterization completed'
      END

C=======================================================================
C  RUN ONE SMOOTH-RAMP TRANSIENT CASE.
C=======================================================================

      SUBROUTINE RUN_CASE(N,DT,NSTEPS,THETA,TEND,MAXSTEP,
     &                    TIME,HPOUT,HMIN,HLINE,POUTF,MINF,
     &                    LINEF,INFO)
      IMPLICIT NONE
      INTEGER N,NSTEPS,MAXSTEP,INFO,I,K,NITER,MAXIT
      DOUBLE PRECISION DT,THETA,TEND,TIME(0:MAXSTEP)
      DOUBLE PRECISION HPOUT(0:MAXSTEP),HMIN(0:MAXSTEP)
      DOUBLE PRECISION HLINE(0:MAXSTEP)
      DOUBLE PRECISION PO(0:100),MO(0:99),PNEW(0:100)
      DOUBLE PRECISION MNEW(0:99),U(202)
      DOUBLE PRECISION PI,D,A,T,Z,RS,MU,EPS,L,DX,PIN
      DOUBLE PRECISION MDOT0,MDOT1,TRAMP,MINO,MOUTO,MOUTN
      DOUBLE PRECISION RE,FD,COEF,X,RTOL,STOL,POUTF,MINF,LINEF
      DOUBLE PRECISION DEMAND_SMOOTH,REYNOLDS_MASS,FRICTION_SJ

      IF (N .GT. 100 .OR. NSTEPS .GT. MAXSTEP) THEN
         INFO = 3
         RETURN
      ENDIF
      IF (DABS(DBLE(NSTEPS)*DT-TEND) .GT. 1.0D-10*TEND) THEN
         INFO = 5
         RETURN
      ENDIF

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
      PIN = 8.0D6
      MDOT0 = 100.0D0
      MDOT1 = 105.0D0
      TRAMP = TEND
      MINO = MDOT0
      MOUTO = MDOT0
      RTOL = 1.0D-11
      STOL = 1.0D-12
      MAXIT = 30

      RE = REYNOLDS_MASS(MDOT0,D,A,MU)
      FD = FRICTION_SJ(RE,EPS,D)
      COEF = FD*Z*RS*T*MDOT0*MDOT0/(D*A*A)
      DO 10 I=0,N
         X = DBLE(I)*DX
         PO(I) = DSQRT(PIN*PIN-COEF*X)
   10 CONTINUE
      DO 20 I=0,N-1
         MO(I) = MDOT0
   20 CONTINUE

      TIME(0) = 0.0D0
      HMIN(0) = MINO
      HPOUT(0) = PO(N)
      CALL LINEPACK_CZ(N,PO,DX,A,T,Z,RS,HLINE(0))

      DO 100 K=1,NSTEPS
         TIME(K) = DBLE(K)*DT
         MOUTN = DEMAND_SMOOTH(TIME(K),MDOT0,MDOT1,TRAMP)
         CALL TRANSIENT_STEP(N,PO,MO,MINO,MOUTO,MOUTN,PIN,
     &        DX,DT,THETA,D,A,T,Z,RS,MU,EPS,RTOL,STOL,MAXIT,0,
     &        U,INFO,NITER)
         IF (INFO .NE. 0) RETURN
         PNEW(0) = U(1)
         DO 30 I=0,N-1
            MNEW(I) = U(2*I+3)
            PNEW(I+1) = U(2*I+4)
   30    CONTINUE
         HMIN(K) = U(2)
         HPOUT(K) = PNEW(N)
         CALL LINEPACK_CZ(N,PNEW,DX,A,T,Z,RS,HLINE(K))
         DO 40 I=0,N
            PO(I) = PNEW(I)
   40    CONTINUE
         DO 50 I=0,N-1
            MO(I) = MNEW(I)
   50    CONTINUE
         MINO = HMIN(K)
         MOUTO = MOUTN
  100 CONTINUE

      POUTF = HPOUT(NSTEPS)
      MINF = HMIN(NSTEPS)
      LINEF = HLINE(NSTEPS)
      INFO = 0
      RETURN
      END

C=======================================================================
C  C1-SMOOTH HALF-COSINE DEMAND RAMP.
C=======================================================================

      DOUBLE PRECISION FUNCTION DEMAND_SMOOTH(TIME,MDOT0,MDOT1,
     &                                        TRAMP)
      IMPLICIT NONE
      DOUBLE PRECISION TIME,MDOT0,MDOT1,TRAMP,PI,S
      PI = 4.0D0*DATAN(1.0D0)
      IF (TIME .LE. 0.0D0) THEN
         DEMAND_SMOOTH = MDOT0
      ELSE IF (TIME .GE. TRAMP) THEN
         DEMAND_SMOOTH = MDOT1
      ELSE
         S = 0.5D0*(1.0D0-DCOS(PI*TIME/TRAMP))
         DEMAND_SMOOTH = MDOT0 + S*(MDOT1-MDOT0)
      ENDIF
      RETURN
      END
