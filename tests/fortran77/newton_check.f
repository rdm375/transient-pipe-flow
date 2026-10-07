      PROGRAM NEWTONCHECK
      IMPLICIT NONE
      INTEGER N,NU,I,ICASE,INFO,NITER
      PARAMETER (N=4,NU=10)
      DOUBLE PRECISION U(NU),USTAR(NU),PO(0:N),MO(0:N-1)
      DOUBLE PRECISION PI,D,A,T,Z,RS,MU,EPS,DX,DT,THETA
      DOUBLE PRECISION MINO,MOUTO,MOUTN,PINN,MDOT,RE,FD,COEF,X
      DOUBLE PRECISION PERT,MAXERR,PSCALE,MSCALE,ERR
      DOUBLE PRECISION REYNOLDS_MASS,FRICTION_SJ
      DOUBLE PRECISION RTOL,STOL

      PI = 4.0D0*DATAN(1.0D0)
      D = 1.0D0
      A = PI*D*D/4.0D0
      T = 288.15D0
      Z = 0.90D0
      RS = 500.0D0
      MU = 1.1D-5
      EPS = 4.5D-5
      DX = 25000.0D0
      DT = 60.0D0
      THETA = 0.65D0
      MDOT = 100.0D0
      PINN = 8.0D6
      MINO = MDOT
      MOUTO = MDOT
      MOUTN = MDOT
      RTOL = 1.0D-11
      STOL = 1.0D-12

      RE = REYNOLDS_MASS(MDOT,D,A,MU)
      FD = FRICTION_SJ(RE,EPS,D)
      COEF = FD*Z*RS*T*MDOT*MDOT/(D*A*A)

      DO 10 I=0,N
         X = DBLE(I)*DX
         PO(I) = DSQRT(PINN*PINN-COEF*X)
   10 CONTINUE
      DO 20 I=0,N-1
         MO(I) = MDOT
   20 CONTINUE

      USTAR(1) = PO(0)
      USTAR(2) = MDOT
      DO 30 I=0,N-1
         USTAR(2*I+3) = MDOT
         USTAR(2*I+4) = PO(I+1)
   30 CONTINUE

      WRITE(*,*)
      WRITE(*,*) 'DAMPED NEWTON RECOVERY CHECK'
      WRITE(*,*) '----------------------------'
      WRITE(*,'(A,1PE12.4)') ' Swamee-Jain f_D = ',FD

      DO 100 ICASE=1,4
         IF (ICASE .EQ. 1) PERT = 1.0D-4
         IF (ICASE .EQ. 2) PERT = 1.0D-3
         IF (ICASE .EQ. 3) PERT = 1.0D-2
         IF (ICASE .EQ. 4) PERT = 5.0D-2
         DO 40 I=1,NU
            IF (MOD(I,2) .EQ. 1 .OR. I .EQ. NU) THEN
               U(I) = USTAR(I)*(1.0D0+PERT*(-1.0D0)**I)
            ELSE
               U(I) = USTAR(I)*(1.0D0+0.5D0*PERT*(-1.0D0)**I)
            ENDIF
   40    CONTINUE
         CALL NEWTON_SOLVE(N,U,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &        DX,DT,THETA,D,A,T,Z,RS,MU,EPS,RTOL,STOL,30,1,
     &        INFO,NITER)
         IF (INFO .NE. 0) THEN
            WRITE(*,*) 'FAIL: Newton did not converge, info=',INFO
            STOP 1
         ENDIF
         MAXERR = 0.0D0
         DO 50 I=1,NU
            IF (I .EQ. 1 .OR. I .EQ. NU .OR.
     &          (I .GE. 4 .AND. MOD(I,2) .EQ. 0)) THEN
               PSCALE = DMAX1(DABS(USTAR(I)),1.0D0)
               ERR = DABS(U(I)-USTAR(I))/PSCALE
            ELSE
               MSCALE = DMAX1(DABS(USTAR(I)),1.0D0)
               ERR = DABS(U(I)-USTAR(I))/MSCALE
            ENDIF
            MAXERR = DMAX1(MAXERR,ERR)
   50    CONTINUE
         WRITE(*,'(A,1PE9.2,A,I3,A,1PE12.4)') ' perturb=',PERT,
     &        ' iterations=',NITER,' max rel state error=',MAXERR
         IF (MAXERR .GT. 2.0D-10) THEN
            WRITE(*,*) 'FAIL: Newton recovered wrong state'
            STOP 1
         ENDIF
  100 CONTINUE

C     Explicit failure-path test: zero iterations cannot converge
C     from a perturbed state.

      DO 60 I=1,NU
         U(I) = USTAR(I)
   60 CONTINUE
      U(1) = 1.01D0*U(1)
      CALL NEWTON_SOLVE(N,U,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &     DX,DT,THETA,D,A,T,Z,RS,MU,EPS,RTOL,STOL,0,0,
     &     INFO,NITER)
      IF (INFO .NE. 1) THEN
         WRITE(*,*) 'FAIL: Newton failure reporting'
         STOP 1
      ENDIF

      WRITE(*,*) 'PASS: damped Newton recovery and failure paths'
      END
