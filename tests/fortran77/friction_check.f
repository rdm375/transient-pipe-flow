      PROGRAM FRICTIONCHECK
      IMPLICIT NONE

      INTEGER I
      DOUBLE PRECISION PI,D,A,MU,EPS,MDOT,RHOF
      DOUBLE PRECISION RE,FD,SF,DFDRE,DFDM,DSFDM
      DOUBLE PRECISION H,FP,FM,FDP,FDM,FDNUM,SFNUM
      DOUBLE PRECISION ERR,MAXERR
      DOUBLE PRECISION REVAL(5)
      DOUBLE PRECISION REYNOLDS_MASS,DREYNOLDS_DMDOT
      DOUBLE PRECISION FRICTION_SJ,DFRICTION_SJ_DRE
      DOUBLE PRECISION FRICTION_SOURCE
      DOUBLE PRECISION DFRICTION_SOURCE_DMDOT

      DATA REVAL /1.0D4,1.0D5,1.0D6,1.0D7,1.0D8/

      PI = 4.0D0*DATAN(1.0D0)
      D = 1.0D0
      A = PI*D*D/4.0D0
      MU = 1.1D-5
      EPS = 4.5D-5
      MDOT = 100.0D0
      RHOF = 60.0D0

      RE = REYNOLDS_MASS(MDOT,D,A,MU)
      FD = FRICTION_SJ(RE,EPS,D)
      SF = FRICTION_SOURCE(FD,D,A,RHOF,MDOT)

      WRITE(*,*)
      WRITE(*,*) 'SWAMEE-JAIN FRICTION CHECK'
      WRITE(*,*) '---------------------------'
      WRITE(*,'(A,1PE24.16)') ' Reynolds number = ',RE
      WRITE(*,'(A,1PE24.16)') ' Darcy factor    = ',FD
      WRITE(*,'(A,1PE24.16)') ' Friction source = ',SF

      IF (DABS(RE-1.1574904952137843D7) .GT. 1.0D-6) STOP 1
      IF (DABS(FD-1.070235485102057D-2) .GT. 1.0D-14) STOP 1
      IF (DABS(SF-1.4458336816876263D0) .GT. 1.0D-13) STOP 1

C     Check df/dRe over four decades using centered differences.

      MAXERR = 0.0D0
      DO 10 I=1,5
         RE = REVAL(I)
         H = 1.0D-5*RE
         DFDRE = DFRICTION_SJ_DRE(RE,EPS,D)
         FP = FRICTION_SJ(RE+H,EPS,D)
         FM = FRICTION_SJ(RE-H,EPS,D)
         FDNUM = (FP-FM)/(2.0D0*H)
         ERR = DABS(DFDRE-FDNUM)/DMAX1(DABS(FDNUM),1.0D-30)
         MAXERR = DMAX1(MAXERR,ERR)
   10 CONTINUE

      WRITE(*,'(A,1PE12.4)') ' max rel error df/dRe = ',MAXERR
      IF (MAXERR .GT. 1.0D-8) STOP 1

C     Check df/dm and dSf/dm for forward and reverse flow.

      MAXERR = 0.0D0
      DO 20 I=1,2
         IF (I .EQ. 1) THEN
            MDOT = 100.0D0
         ELSE
            MDOT = -100.0D0
         ENDIF

         RE = REYNOLDS_MASS(MDOT,D,A,MU)
         FD = FRICTION_SJ(RE,EPS,D)
         DFDRE = DFRICTION_SJ_DRE(RE,EPS,D)
         DFDM = DFDRE*DREYNOLDS_DMDOT(MDOT,D,A,MU)
         DSFDM = DFRICTION_SOURCE_DMDOT(FD,DFDM,D,A,RHOF,
     &                                  MDOT)

         H = 1.0D-4
         FDP = FRICTION_SJ(REYNOLDS_MASS(MDOT+H,D,A,MU),
     &                     EPS,D)
         FDM = FRICTION_SJ(REYNOLDS_MASS(MDOT-H,D,A,MU),
     &                     EPS,D)
         FDNUM = (FDP-FDM)/(2.0D0*H)
         ERR = DABS(DFDM-FDNUM)/DMAX1(DABS(FDNUM),1.0D-30)
         MAXERR = DMAX1(MAXERR,ERR)

         FP = FRICTION_SOURCE(FDP,D,A,RHOF,MDOT+H)
         FM = FRICTION_SOURCE(FDM,D,A,RHOF,MDOT-H)
         SFNUM = (FP-FM)/(2.0D0*H)
         ERR = DABS(DSFDM-SFNUM)/DMAX1(DABS(SFNUM),1.0D-30)
         MAXERR = DMAX1(MAXERR,ERR)
   20 CONTINUE

      WRITE(*,'(A,1PE12.4)') ' max rel error flow deriv = ',MAXERR
      IF (MAXERR .GT. 1.0D-8) STOP 1

      WRITE(*,*) 'PASS: Swamee-Jain values and derivatives'
      END
