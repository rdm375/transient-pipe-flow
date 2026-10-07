C=======================================================================
C  CONSTANT-Z ISOTHERMAL EQUATION OF STATE
C
C       RHO     = P / (Z RS T)
C       DRHO/DP = 1 / (Z RS T)
C=======================================================================

      DOUBLE PRECISION FUNCTION EOS_CZ_RHO(P,T,Z,RS)
      IMPLICIT NONE
      DOUBLE PRECISION P,T,Z,RS

      EOS_CZ_RHO = P/(Z*RS*T)

      RETURN
      END


      DOUBLE PRECISION FUNCTION EOS_CZ_DRHODP(T,Z,RS)
      IMPLICIT NONE
      DOUBLE PRECISION T,Z,RS

      EOS_CZ_DRHODP = 1.0D0/(Z*RS*T)

      RETURN
      END
