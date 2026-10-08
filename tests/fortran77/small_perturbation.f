C     M7B SMALL-PERTURBATION NONLINEAR REFERENCE.
      PROGRAM SMALLVERIFY
      IMPLICIT NONE
      INTEGER MAXSTEP,NS,INFO,CASEID,I,J,K
      PARAMETER (MAXSTEP=2400)
      DOUBLE PRECISION TIME(0:MAXSTEP),P(0:MAXSTEP)
      DOUBLE PRECISION M(0:MAXSTEP),L(0:MAXSTEP)
      DOUBLE PRECISION AMPS(4),THS(3),DTS(2),PF,MF,LF
      DATA AMPS /1.0D0,0.5D0,0.25D0,0.125D0/
      DATA THS /0.5D0,0.65D0,1.0D0/
      DATA DTS /15.0D0,60.0D0/
      CALL SYSTEM('mkdir -p benchmarks/results-m7b')
      OPEN(20,FILE='benchmarks/results-m7b/small_nonlinear.csv',
     & STATUS='REPLACE')
      WRITE(20,'(A)') 'case_id,forcing,n,dt,theta,time,'//
     & 'p25_pa,p50_pa,p75_pa,pout_pa,min_kg_s,'//
     & 'linepack_kg,step_mass_defect_kg,cum_mass_defect_kg'
      CASEID=0
      DO 30 I=1,4
         DO 20 J=1,3
            DO 10 K=1,2
               CASEID=CASEID+1
               NS=IDNINT(1800.0D0/DTS(K))
               CALL SMALL_CASE(40,DTS(K),NS,THS(J),1800.0D0,
     &              MAXSTEP,TIME,P,M,L,PF,MF,LF,3,INFO,20,
     &              CASEID,AMPS(I))
               IF (INFO .NE. 0) THEN
                  WRITE(*,*) 'FAIL small verification',CASEID,INFO
                  STOP 1
               ENDIF
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      CLOSE(20)
      WRITE(*,*) 'PASS: small nonlinear histories'
      END

      SUBROUTINE SMALL_CASE(N,DT,NSTEPS,THETA,TEND,MAXSTEP,
     &                    TIME,HPOUT,HMIN,HLINE,POUTF,MINF,
     &                    LINEF,MODE,INFO,UNIT,CASEID,AMP)
      IMPLICIT NONE
      INTEGER N,NSTEPS,MAXSTEP,INFO,I,K,NITER,MAXIT,MODE
      INTEGER UNIT,CASEID
      DOUBLE PRECISION DT,THETA,TEND,TIME(0:MAXSTEP),AMP
      DOUBLE PRECISION HPOUT(0:MAXSTEP),HMIN(0:MAXSTEP)
      DOUBLE PRECISION HLINE(0:MAXSTEP)
      DOUBLE PRECISION PO(0:100),MO(0:99),PNEW(0:100)
      DOUBLE PRECISION MNEW(0:99),U(202)
      DOUBLE PRECISION PI,D,A,T,Z,RS,MU,EPS,L,DX,PIN
      DOUBLE PRECISION MDOT0,MDOT1,TRAMP,MINO,MOUTO,MOUTN
      DOUBLE PRECISION RE,FD,COEF,X,RTOL,STOL,POUTF,MINF,LINEF
      DOUBLE PRECISION SMALL_DEMAND,REYNOLDS_MASS,FRICTION_SJ
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
      MDOT1 = MDOT0 + AMP
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
         MOUTN = SMALL_DEMAND(TIME(K),MDOT0,AMP)
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
     & F12.3,8(',',ES24.16))
      END

      DOUBLE PRECISION FUNCTION SMALL_DEMAND(T,M0,AMP)
      IMPLICIT NONE
      DOUBLE PRECISION T,M0,AMP,S,PI
      PI=4.0D0*DATAN(1.0D0)
      S=DMIN1(1.0D0,DMAX1(0.0D0,T/600.0D0))
      SMALL_DEMAND=M0+AMP*0.5D0*(1.0D0-DCOS(PI*S))
      RETURN
      END
