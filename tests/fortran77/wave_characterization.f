      PROGRAM M7BWAVES
      IMPLICIT NONE
      INTEGER MAXSTEP,NS,INFO,MODE,I,J,K,CASEID,N
      PARAMETER (MAXSTEP=2400)
      DOUBLE PRECISION T(0:MAXSTEP),P(0:MAXSTEP)
      DOUBLE PRECISION M(0:MAXSTEP),L(0:MAXSTEP)
      DOUBLE PRECISION DT,THETA,PF,MF,LF
      DOUBLE PRECISION TH(3),DTS(3)
      INTEGER NN(2)
      DATA TH /0.50D0,0.65D0,1.00D0/
      DATA DTS /15.0D0,60.0D0,120.0D0/
      DATA NN /20,40/
      CALL SYSTEM('mkdir -p benchmarks/results-m7b')
      OPEN(20,FILE='benchmarks/results-m7b/sensors.csv',
     & STATUS='REPLACE')
      WRITE(20,'(A)') 'case_id,forcing,n,dt,theta,time,'//
     & 'p25_pa,p50_pa,p75_pa,pout_pa,min_kg_s,'//
     & 'linepack_kg,step_mass_defect_kg,cum_mass_defect_kg'
      CASEID = 0
      DO 100 MODE=1,3
         DO 90 I=1,2
            N=NN(I)
            DO 80 J=1,3
               DT=DTS(J)
               NS=IDNINT(1800.0D0/DT)
               DO 70 K=1,3
                  THETA=TH(K)
                  CASEID=CASEID+1
                  CALL M7B_CASE(N,DT,NS,THETA,1800.0D0,
     &               MAXSTEP,T,P,M,L,PF,MF,LF,MODE,INFO,20,CASEID)
                  IF (INFO .NE. 0) THEN
                     WRITE(*,*) 'FAIL M7b case',CASEID,INFO
                     STOP 1
                  ENDIF
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      CLOSE(20)
      WRITE(*,*) 'PASS: M7b sensor histories generated'
      END

      SUBROUTINE M7B_CASE(N,DT,NSTEPS,THETA,TEND,MAXSTEP,
     &                    TIME,HPOUT,HMIN,HLINE,POUTF,MINF,
     &                    LINEF,MODE,INFO,UNIT,CASEID)
      IMPLICIT NONE
      INTEGER N,NSTEPS,MAXSTEP,INFO,I,K,NITER,MAXIT,MODE
      INTEGER UNIT,CASEID
      DOUBLE PRECISION DT,THETA,TEND,TIME(0:MAXSTEP)
      DOUBLE PRECISION HPOUT(0:MAXSTEP),HMIN(0:MAXSTEP)
      DOUBLE PRECISION HLINE(0:MAXSTEP)
      DOUBLE PRECISION PO(0:100),MO(0:99),PNEW(0:100)
      DOUBLE PRECISION MNEW(0:99),U(202)
      DOUBLE PRECISION PI,D,A,T,Z,RS,MU,EPS,L,DX,PIN
      DOUBLE PRECISION MDOT0,MDOT1,TRAMP,MINO,MOUTO,MOUTN
      DOUBLE PRECISION RE,FD,COEF,X,RTOL,STOL,POUTF,MINF,LINEF
      DOUBLE PRECISION DEMAND_DYNAMIC,REYNOLDS_MASS,FRICTION_SJ
      DOUBLE PRECISION LOLD,DEFECT,ACCUM,INTEGRAL

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

      ACCUM = 0.0D0
      LOLD = HLINE(0)
      WRITE(UNIT,900) CASEID,MODE,N,DT,THETA,0.0D0,
     & PO(N/4),PO(N/2),PO(3*N/4),HPOUT(0),HMIN(0),
     & HLINE(0),0.0D0,0.0D0
      DO 100 K=1,NSTEPS
         TIME(K) = DBLE(K)*DT
         MOUTN = DEMAND_DYNAMIC(TIME(K),MDOT0,MDOT1,MODE)
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
         INTEGRAL = DT*(THETA*(HMIN(K)-MOUTN)
     &         +(1.0D0-THETA)*(MINO-MOUTO))
         ACCUM = ACCUM + INTEGRAL
         DEFECT = HLINE(K)-LOLD-INTEGRAL
         WRITE(UNIT,900) CASEID,MODE,N,DT,THETA,TIME(K),
     & PNEW(N/4),PNEW(N/2),PNEW(3*N/4),HPOUT(K),
     & HMIN(K),HLINE(K),DEFECT,HLINE(K)-HLINE(0)-ACCUM
         LOLD = HLINE(K)
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
  900 FORMAT(I5,',',I1,',',I3,',',F9.3,',',F5.2,',',
     & F12.3,8(',',1PE20.12))
      END

      DOUBLE PRECISION FUNCTION DEMAND_DYNAMIC(T,M0,M1,MODE)
      IMPLICIT NONE
      INTEGER MODE
      DOUBLE PRECISION T,M0,M1,PI,S
      PI = 4.0D0*DATAN(1.0D0)
      IF (MODE .EQ. 1) THEN
C        Step demand applied at the first positive time.
         DEMAND_DYNAMIC = M1
      ELSE IF (MODE .EQ. 2) THEN
C        Rectangular pulse: 300 <= t < 900 seconds.
         DEMAND_DYNAMIC = M0
         IF (T .GE. 300.0D0 .AND. T .LT. 900.0D0)
     &       DEMAND_DYNAMIC = M1
      ELSE
C        Smooth cosine ramp over 600 seconds.
         S = DMIN1(1.0D0,DMAX1(0.0D0,T/600.0D0))
         DEMAND_DYNAMIC = M0 + 0.5D0*(1.0D0-DCOS(PI*S))*
     &                    (M1-M0)
      ENDIF
      RETURN
      END
