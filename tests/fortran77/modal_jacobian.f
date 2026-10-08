C     M7B LINEARIZED OPERATOR FROM PRODUCTION ANALYTIC JACOBIAN.
C     Inlet pressure fixed; inlet flow algebraic.  Dynamic unknowns
C     are face flows and interior pressures (U(3:2*N+2)).
      PROGRAM MODALJAC
      IMPLICIT NONE
      INTEGER N,I,K,NU,ROW
      DOUBLE PRECISION U(202),J(202,202),DX,DT,THETA
      DOUBLE PRECISION D,A,T,Z,RS,MU,EPS,L,PIN,MDOT
      DOUBLE PRECISION RE,FD,COEF,PI,CRHO,MASS
      DOUBLE PRECISION REYNOLDS_MASS,FRICTION_SJ
      DOUBLE PRECISION EOS_CZ_DRHODP
      CALL SYSTEM('mkdir -p benchmarks/results-m7b')
      OPEN(21,FILE='benchmarks/results-m7b/modal_matrix.csv',
     & STATUS='REPLACE')
      WRITE(21,'(A)') 'n,row,col,mass,operator'
      DO 100 N=20,40,20
         NU=2*N+2
         PI=4.0D0*DATAN(1.0D0)
         D=1.0D0
         A=PI*D*D/4.0D0
         T=288.15D0
         Z=0.90D0
         RS=500.0D0
         MU=1.1D-5
         EPS=4.5D-5
         L=100000.0D0
         DX=L/DBLE(N)
         PIN=8.0D6
         MDOT=100.0D0
         DT=1.0D0
         THETA=1.0D0
         RE=REYNOLDS_MASS(MDOT,D,A,MU)
         FD=FRICTION_SJ(RE,EPS,D)
         COEF=FD*Z*RS*T*MDOT*MDOT/(D*A*A)
         U(1)=PIN
         U(2)=MDOT
         DO 10 I=0,N-1
            U(2*I+3)=MDOT
            U(2*I+4)=DSQRT(PIN*PIN-COEF*DBLE(I+1)*DX)
   10    CONTINUE
         CALL ASSEMBLE_JACOBIAN(N,U,DX,DT,THETA,D,A,T,Z,RS,
     &                          MU,EPS,J,202)
         CRHO=EOS_CZ_DRHODP(T,Z,RS)
         DO 30 I=3,NU
            ROW=I-2
            IF (MOD(I,2) .EQ. 1) THEN
               MASS=1.0D0/A
            ELSE IF (I .EQ. NU) THEN
               MASS=A*DX*CRHO/2.0D0
            ELSE
               MASS=A*DX*CRHO
            ENDIF
            DO 20 K=3,NU
               IF (I .EQ. K) THEN
                  WRITE(21,900) N,ROW,K-2,MASS,
     &                    -(J(I,K)-MASS/DT)/MASS
               ELSE
                  WRITE(21,900) N,ROW,K-2,MASS,-J(I,K)/MASS
               ENDIF
   20       CONTINUE
   30    CONTINUE
  100 CONTINUE
      CLOSE(21)
      WRITE(*,*) 'PASS: production Jacobian modal matrices'
  900 FORMAT(I4,',',I4,',',I4,',',ES24.16,',',ES24.16)
      END
