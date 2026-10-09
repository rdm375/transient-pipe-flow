C=======================================================================
C  SWAMEE-JAIN DARCY FRICTION MODEL
C=======================================================================

      DOUBLE PRECISION FUNCTION REYNOLDS_MASS(MDOT,D,A,MU)
      IMPLICIT NONE
      DOUBLE PRECISION MDOT,D,A,MU

      IF (D .LE. 0.0D0) STOP 'ERROR: nonpositive diameter'
      IF (A .LE. 0.0D0) STOP 'ERROR: nonpositive area'
      IF (MU .LE. 0.0D0) STOP 'ERROR: nonpositive viscosity'

      REYNOLDS_MASS = DABS(MDOT)*D/(A*MU)
      RETURN
      END


      DOUBLE PRECISION FUNCTION DREYNOLDS_DMDOT(MDOT,D,A,MU)
      IMPLICIT NONE
      DOUBLE PRECISION MDOT,D,A,MU

      IF (D .LE. 0.0D0) STOP 'ERROR: nonpositive diameter'
      IF (A .LE. 0.0D0) STOP 'ERROR: nonpositive area'
      IF (MU .LE. 0.0D0) STOP 'ERROR: nonpositive viscosity'

      DREYNOLDS_DMDOT = DSIGN(D/(A*MU),MDOT)
      RETURN
      END


      DOUBLE PRECISION FUNCTION FRICTION_SJ(RE,EPS,D)
      IMPLICIT NONE
      DOUBLE PRECISION RE,EPS,D
      DOUBLE PRECISION ARG,LARG

      IF (RE .LE. 0.0D0) STOP 'ERROR: nonpositive Reynolds number'
      IF (EPS .LT. 0.0D0) STOP 'ERROR: negative roughness'
      IF (D .LE. 0.0D0) STOP 'ERROR: nonpositive diameter'

      ARG = EPS/(3.7D0*D) + 5.74D0/(RE**0.9D0)
      LARG = DLOG10(ARG)
      FRICTION_SJ = 0.25D0/(LARG*LARG)
      RETURN
      END


      DOUBLE PRECISION FUNCTION DFRICTION_SJ_DRE(RE,EPS,D)
      IMPLICIT NONE
      DOUBLE PRECISION RE,EPS,D
      DOUBLE PRECISION ARG,LARG,DXDRE,LN10

      IF (RE .LE. 0.0D0) STOP 'ERROR: nonpositive Reynolds number'
      IF (EPS .LT. 0.0D0) STOP 'ERROR: negative roughness'
      IF (D .LE. 0.0D0) STOP 'ERROR: nonpositive diameter'

      ARG = EPS/(3.7D0*D) + 5.74D0/(RE**0.9D0)
      LARG = DLOG10(ARG)
      LN10 = DLOG(10.0D0)
      DXDRE = -0.9D0*5.74D0/(RE**1.9D0)

      DFRICTION_SJ_DRE =
     &   -0.5D0*DXDRE/(ARG*LN10*LARG*LARG*LARG)
      RETURN
      END


      DOUBLE PRECISION FUNCTION FRICTION_SOURCE(FD,D,A,RHOF,
     &                                          MDOT)
      IMPLICIT NONE
      DOUBLE PRECISION FD,D,A,RHOF,MDOT

      IF (D .LE. 0.0D0) STOP 'ERROR: nonpositive diameter'
      IF (A .LE. 0.0D0) STOP 'ERROR: nonpositive area'
      IF (RHOF .LE. 0.0D0) STOP 'ERROR: nonpositive face density'

      FRICTION_SOURCE =
     &   FD*MDOT*DABS(MDOT)/(2.0D0*D*A*A*RHOF)
      RETURN
      END


      DOUBLE PRECISION FUNCTION DFRICTION_SOURCE_DMDOT(FD,DFDM,
     &                                                  D,A,RHOF,
     &                                                  MDOT)
      IMPLICIT NONE
      DOUBLE PRECISION FD,DFDM,D,A,RHOF,MDOT
      DOUBLE PRECISION C

      IF (D .LE. 0.0D0) STOP 'ERROR: nonpositive diameter'
      IF (A .LE. 0.0D0) STOP 'ERROR: nonpositive area'
      IF (RHOF .LE. 0.0D0) STOP 'ERROR: nonpositive face density'

      C = 1.0D0/(2.0D0*D*A*A*RHOF)
      DFRICTION_SOURCE_DMDOT = C*(DFDM*MDOT*DABS(MDOT)
     &                           + 2.0D0*FD*DABS(MDOT))
      RETURN
      END

C=======================================================================
C  COMBINED SWAMEE-JAIN VALUE AND REYNOLDS DERIVATIVE
C=======================================================================

      SUBROUTINE FRICTION_SJ_VALUE_DERIVATIVE(RE,EPS,D,FD,DFDRE)
      IMPLICIT NONE
      DOUBLE PRECISION RE,EPS,D,FD,DFDRE
      DOUBLE PRECISION ARG,LARG,DXDRE,LN10,RE09

      IF (RE .LE. 0.0D0) STOP 'ERROR: nonpositive Reynolds number'
      IF (EPS .LT. 0.0D0) STOP 'ERROR: negative roughness'
      IF (D .LE. 0.0D0) STOP 'ERROR: nonpositive diameter'

      RE09 = RE**0.9D0
      ARG = EPS/(3.7D0*D) + 5.74D0/RE09
      LARG = DLOG10(ARG)

      FD = 0.25D0/(LARG*LARG)

      LN10 = DLOG(10.0D0)
      DXDRE = -0.9D0*5.74D0/(RE*RE09)

      DFDRE = -0.5D0*DXDRE/
     &        (ARG*LN10*LARG*LARG*LARG)

      RETURN
      END
