      PROGRAM TIMESTEPCHECK
      IMPLICIT NONE
      INTEGER N,NU,I,ICASE,INFO,NITER
      PARAMETER (N=4,NU=10)
      DOUBLE PRECISION U(NU),PO(0:N),MO(0:N-1),PSTAR(0:N)
      DOUBLE PRECISION PI,D,A,T,Z,RS,MU,EPS,DX,DT,THETA
      DOUBLE PRECISION MINO,MOUTO,MOUTN,PINN,MDOT,RE,FD,COEF,X
      DOUBLE PRECISION RTOL,STOL,MAXERR,ERR,MOLD,MNEW,BAL,RHS
      DOUBLE PRECISION RERR,MINN,PNEWN
      DOUBLE PRECISION REYNOLDS_MASS,FRICTION_SJ

      PI = 4.0D0*DATAN(1.0D0)
      D = 1.0D0
      A = PI*D*D/4.0D0
      T = 288.15D0
      Z = 0.90D0
      RS = 500.0D0
      MU = 1.1D-5
      EPS = 4.5D-5
      DX = 25000.0D0
      MDOT = 100.0D0
      PINN = 8.0D6
      MINO = MDOT
      MOUTO = MDOT
      RTOL = 1.0D-11
      STOL = 1.0D-12

      RE = REYNOLDS_MASS(MDOT,D,A,MU)
      FD = FRICTION_SJ(RE,EPS,D)
      COEF = FD*Z*RS*T*MDOT*MDOT/(D*A*A)
      DO 10 I=0,N
         X = DBLE(I)*DX
         PSTAR(I) = DSQRT(PINN*PINN-COEF*X)
         PO(I) = PSTAR(I)
   10 CONTINUE
      DO 20 I=0,N-1
         MO(I) = MDOT
   20 CONTINUE

      WRITE(*,*)
      WRITE(*,*) 'SINGLE IMPLICIT TIME-STEP CHECK'
      WRITE(*,*) '-------------------------------'

C     Exact steady-state preservation for several dt/theta pairs.

      DO 100 ICASE=1,4
         IF (ICASE .EQ. 1) THEN
            DT = 10.0D0
            THETA = 0.50D0
         ELSE IF (ICASE .EQ. 2) THEN
            DT = 60.0D0
            THETA = 0.65D0
         ELSE IF (ICASE .EQ. 3) THEN
            DT = 600.0D0
            THETA = 0.80D0
         ELSE
            DT = 3600.0D0
            THETA = 1.00D0
         ENDIF
         MOUTN = MDOT
         CALL TRANSIENT_STEP(N,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &        DX,DT,THETA,D,A,T,Z,RS,MU,EPS,RTOL,STOL,30,0,
     &        U,INFO,NITER)
         IF (INFO .NE. 0) THEN
            WRITE(*,*) 'FAIL: steady time step, info=',INFO
            STOP 1
         ENDIF
         MAXERR = 0.0D0
         MAXERR = DMAX1(MAXERR,DABS(U(1)-PSTAR(0))/PSTAR(0))
         MAXERR = DMAX1(MAXERR,DABS(U(2)-MDOT)/MDOT)
         DO 30 I=0,N-1
            MAXERR = DMAX1(MAXERR,DABS(U(2*I+3)-MDOT)/MDOT)
            ERR = DABS(U(2*I+4)-PSTAR(I+1))/PSTAR(I+1)
            MAXERR = DMAX1(MAXERR,ERR)
   30    CONTINUE
         CALL CHECK_BALANCE(N,PO,U,MINO,MOUTO,MOUTN,DX,DT,
     &        THETA,A,T,Z,RS,MOLD,MNEW,BAL,RHS,RERR)
         WRITE(*,'(A,F7.1,A,F5.2,A,I3,A,1PE11.3,A,1PE11.3)')
     &        ' dt=',DT,' theta=',THETA,' iter=',NITER,
     &        ' state err=',MAXERR,' balance err=',RERR
         IF (MAXERR .GT. 2.0D-10 .OR. RERR .GT. 2.0D-11) THEN
            WRITE(*,*) 'FAIL: steady preservation or balance'
            STOP 1
         ENDIF
  100 CONTINUE

C     A failed step must not overwrite its output buffer.
      DO 101 I=1,NU
         U(I)=-12345D0
  101 CONTINUE
      MOUTN=1.05D0*MDOT
      CALL TRANSIENT_STEP(N,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &     DX,60D0,0.65D0,D,A,T,Z,RS,MU,EPS,RTOL,STOL,0,0,
     &     U,INFO,NITER)
      IF (INFO.NE.1.OR.NITER.NE.0) THEN
         WRITE(*,*) 'FAIL: timestep failure status'
         STOP 1
      ENDIF
      DO 102 I=1,NU
         IF (TRANSFER(U(I),0_8).NE.
     &       TRANSFER(-12345D0,0_8)) THEN
            WRITE(*,*) 'FAIL: timestep changed output on failure'
            STOP 1
         ENDIF
  102 CONTINUE

C     One genuine transient step: increase downstream demand by 5%.
C     The inlet pressure is held fixed.  The outlet pressure and total
C     linepack should fall, while the solved inlet flow is diagnostic.

      DT = 60.0D0
      THETA = 0.65D0
      MOUTN = 1.05D0*MDOT
      CALL TRANSIENT_STEP(N,PO,MO,MINO,MOUTO,MOUTN,PINN,
     &     DX,DT,THETA,D,A,T,Z,RS,MU,EPS,RTOL,STOL,30,0,
     &     U,INFO,NITER)
      IF (INFO .NE. 0) THEN
         WRITE(*,*) 'FAIL: perturbed time step, info=',INFO
         STOP 1
      ENDIF
      CALL CHECK_BALANCE(N,PO,U,MINO,MOUTO,MOUTN,DX,DT,
     &     THETA,A,T,Z,RS,MOLD,MNEW,BAL,RHS,RERR)
      MINN = U(2)
      PNEWN = U(NU)
      WRITE(*,*)
      WRITE(*,*) ' downstream-demand perturbation'
      WRITE(*,'(A,I3)') ' Newton iterations       = ',NITER
      WRITE(*,'(A,1PE14.6)') ' new inlet flow [kg/s]  = ',MINN
      WRITE(*,'(A,1PE14.6)') ' old outlet p [Pa]      = ',PO(N)
      WRITE(*,'(A,1PE14.6)') ' new outlet p [Pa]      = ',PNEWN
      WRITE(*,'(A,1PE14.6)') ' linepack change [kg]   = ',MNEW-MOLD
      WRITE(*,'(A,1PE14.6)') ' balance lhs [kg/s]     = ',BAL
      WRITE(*,'(A,1PE14.6)') ' balance rhs [kg/s]     = ',RHS
      WRITE(*,'(A,1PE14.6)') ' relative balance error = ',RERR

      IF (PNEWN .GE. PO(N)) THEN
         WRITE(*,*) 'FAIL: outlet pressure did not fall'
         STOP 1
      ENDIF
      IF (MNEW .GE. MOLD) THEN
         WRITE(*,*) 'FAIL: linepack did not fall'
         STOP 1
      ENDIF
      IF (RERR .GT. 2.0D-11) THEN
         WRITE(*,*) 'FAIL: discrete global mass balance'
         STOP 1
      ENDIF
      IF (DABS(U(1)-PINN) .GT. 1.0D-8) THEN
         WRITE(*,*) 'FAIL: new inlet pressure history'
         STOP 1
      ENDIF
      IF (DABS(MOUTO-MDOT) .GT. 1.0D-12 .OR.
     &    DABS(MOUTN-1.05D0*MDOT) .GT. 1.0D-12) THEN
         WRITE(*,*) 'FAIL: outlet boundary history'
         STOP 1
      ENDIF

      WRITE(*,*) 'PASS: single implicit time step verified'
      END

C=======================================================================
C  Independent linepack and global balance calculation.
C=======================================================================

      SUBROUTINE CHECK_BALANCE(N,PO,U,MINO,MOUTO,MOUTN,DX,DT,
     &                         THETA,A,T,Z,RS,MOLD,MNEW,BAL,
     &                         RHS,RERR)
      IMPLICIT NONE
      INTEGER N,I,IP
      DOUBLE PRECISION PO(0:*),U(*)
      DOUBLE PRECISION MINO,MOUTO,MOUTN,DX,DT,THETA,A,T,Z,RS
      DOUBLE PRECISION MOLD,MNEW,BAL,RHS,RERR,RHO
      DOUBLE PRECISION EOS_CZ_RHO

      MOLD = 0.0D0
      MNEW = 0.0D0
      DO 10 I=0,N
         IF (I .EQ. 0 .OR. I .EQ. N) THEN
            RHO = 0.5D0
         ELSE
            RHO = 1.0D0
         ENDIF
         MOLD = MOLD + RHO*EOS_CZ_RHO(PO(I),T,Z,RS)
         IF (I .EQ. 0) THEN
            IP = 1
         ELSE
            IP = 2*I + 2
         ENDIF
         MNEW = MNEW + RHO*EOS_CZ_RHO(U(IP),T,Z,RS)
   10 CONTINUE
      MOLD = A*DX*MOLD
      MNEW = A*DX*MNEW
      BAL = (MNEW-MOLD)/DT
      RHS = THETA*(U(2)-MOUTN)
     &    +(1.0D0-THETA)*(MINO-MOUTO)
      RERR = DABS(BAL-RHS)/DMAX1(1.0D0,DABS(BAL),DABS(RHS))

      RETURN
      END
