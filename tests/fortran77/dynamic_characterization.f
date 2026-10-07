      PROGRAM DYNAMICCHARACTERIZATION
      IMPLICIT NONE
      INTEGER MAXSTEP,NCASE,NTH,I,J,K,INFO,MODE,NS,NREF
      PARAMETER (MAXSTEP=2400,NCASE=3,NTH=7)
      DOUBLE PRECISION TH(NTH),DT,TEND
      DOUBLE PRECISION T(0:MAXSTEP),P(0:MAXSTEP)
      DOUBLE PRECISION M(0:MAXSTEP),L(0:MAXSTEP)
      DOUBLE PRECISION TR(0:MAXSTEP),PR(0:MAXSTEP)
      DOUBLE PRECISION MR(0:MAXSTEP),LR(0:MAXSTEP)
      DOUBLE PRECISION PF,MF,LF,RP,RM,RL,EP,EM,EL
      DOUBLE PRECISION MAXP,MAXM,MAXL,MINP,MAXPVAL
      DATA TH /0.50D0,0.55D0,0.60D0,0.65D0,0.70D0,
     &         0.80D0,1.00D0/
      CALL SYSTEM('mkdir -p benchmarks/results-m7')
      OPEN(10,FILE='benchmarks/results-m7/dynamics.csv',
     &     STATUS='REPLACE')
      WRITE(10,'(A)') 'forcing,n,dt,theta,max_rel_pout,'//
     & 'max_rel_min,max_rel_linepack,peak_pout,low_pout'
      OPEN(11,FILE='benchmarks/results-m7/history.csv',
     &     STATUS='REPLACE')
      WRITE(11,'(A)') 'forcing,theta,time,pout,min,linepack'
      TEND = 1800.0D0
      NREF = 1200
      DO 100 MODE=1,NCASE
C        Fine time reference on identical spatial grid.
         CALL DYNAMIC_CASE(40,1.5D0,NREF,0.5D0,TEND,
     &      MAXSTEP,TR,PR,MR,LR,RP,RM,RL,MODE,INFO)
         IF (INFO .NE. 0) STOP 1
         DO 90 I=1,NTH
            DT = 60.0D0
            NS = 30
            CALL DYNAMIC_CASE(40,DT,NS,TH(I),TEND,
     &       MAXSTEP,T,P,M,L,PF,MF,LF,MODE,INFO)
            IF (INFO .NE. 0) STOP 1
            MAXP=0.0D0
            MAXM=0.0D0
            MAXL=0.0D0
            MINP=P(0)
            MAXPVAL=P(0)
            DO 20 K=0,NS
               J=40*K
               EP=DABS(P(K)-PR(J))/DABS(PR(J))
               EM=DABS(M(K)-MR(J))/DMAX1(1.0D0,DABS(MR(J)))
               EL=DABS(L(K)-LR(J))/DABS(LR(J))
               MAXP=DMAX1(MAXP,EP)
               MAXM=DMAX1(MAXM,EM)
               MAXL=DMAX1(MAXL,EL)
               MINP=DMIN1(MINP,P(K))
               MAXPVAL=DMAX1(MAXPVAL,P(K))
               WRITE(11,'(I1,A,F5.2,A,F9.2,3(A,1PE16.8))')
     &            MODE,',',TH(I),',',T(K),',',P(K),',',
     &            M(K),',',L(K)
   20       CONTINUE
            WRITE(10,'(I1,A,I4,A,F9.3,A,F5.2,5(A,1PE16.8))')
     &         MODE,',',40,',',DT,',',TH(I),',',MAXP,',',
     &         MAXM,',',MAXL,',',MAXPVAL,',',MINP
C           Loose regression guards, not fitted precision claims.
            IF (MAXP .GT. 2.0D-3 .OR. MAXM .GT. 2.0D-2
     &          .OR. MAXL .GT. 2.0D-3) THEN
               WRITE(*,*) 'FAIL: dynamics regression',MODE,I
               STOP 1
            ENDIF
   90    CONTINUE
  100 CONTINUE
      CLOSE(10)
      CLOSE(11)
      WRITE(*,*) 'PASS: M7 dynamic histories and comparisons'
      END

      SUBROUTINE DYNAMIC_CASE(N,DT,NSTEPS,THETA,TEND,MAXSTEP,
     &                    TIME,HPOUT,HMIN,HLINE,POUTF,MINF,
     &                    LINEF,MODE,INFO)
      IMPLICIT NONE
      INTEGER N,NSTEPS,MAXSTEP,INFO,I,K,NITER,MAXIT,MODE
      DOUBLE PRECISION DT,THETA,TEND,TIME(0:MAXSTEP)
      DOUBLE PRECISION HPOUT(0:MAXSTEP),HMIN(0:MAXSTEP)
      DOUBLE PRECISION HLINE(0:MAXSTEP)
      DOUBLE PRECISION PO(0:100),MO(0:99),PNEW(0:100)
      DOUBLE PRECISION MNEW(0:99),U(202)
      DOUBLE PRECISION PI,D,A,T,Z,RS,MU,EPS,L,DX,PIN
      DOUBLE PRECISION MDOT0,MDOT1,TRAMP,MINO,MOUTO,MOUTN
      DOUBLE PRECISION RE,FD,COEF,X,RTOL,STOL,POUTF,MINF,LINEF
      DOUBLE PRECISION DEMAND_DYNAMIC,REYNOLDS_MASS,FRICTION_SJ

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
