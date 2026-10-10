C=======================================================================
C  ANALYTIC JACOBIAN OF THE TRANSIENT RESIDUAL
C=======================================================================

      SUBROUTINE ASSEMBLE_JACOBIAN_BANDED(N,U,DX,DT,THETA,
     &                             D,A,T,Z,RS,MU,EPS,AB,LDAB)
      IMPLICIT NONE
      INTEGER N,LDAB,I,K,IP,IL,IR,IM,ROW,NU,KV
      DOUBLE PRECISION U(*),AB(LDAB,*)
      DOUBLE PRECISION DX,DT,THETA,D,A,T,Z,RS,MU,EPS
      DOUBLE PRECISION CRHO,RHOL,RHOR,RHOF,RE,FD,DFDRE,DREDM
      DOUBLE PRECISION DFDM,SF,DSDM,DSDP
      DOUBLE PRECISION EOS_CZ_RHO,EOS_CZ_DRHODP
      DOUBLE PRECISION REYNOLDS_MASS
      DOUBLE PRECISION DREYNOLDS_DMDOT
      DOUBLE PRECISION FRICTION_SOURCE,DFRICTION_SOURCE_DMDOT

      NU = 2*N + 2
      KV = 3
      DO 20 K=1,NU
         DO 10 I=1,LDAB
            AB(I,K) = 0.0D0
   10    CONTINUE
   20 CONTINUE

      CRHO = EOS_CZ_DRHODP(T,Z,RS)

C     Inlet pressure BC and inlet half-cell continuity.

      AB(KV+1+(1)-(1),(1)) = 1.0D0
      AB(KV+1+(2)-(1),(1)) = A*DX*CRHO/(2.0D0*DT)
      AB(KV+1+(2)-(2),(2)) = -THETA
      AB(KV+1+(2)-(3),(3)) = THETA

      DO 100 I=0,N-1
         IM = 2*I + 3
         IF (I .EQ. 0) THEN
            IL = 1
         ELSE
            IL = 2*I + 2
         ENDIF
         IR = 2*I + 4
         ROW = 2*I + 3

         RHOL = EOS_CZ_RHO(U(IL),T,Z,RS)
         RHOR = EOS_CZ_RHO(U(IR),T,Z,RS)
         RHOF = 0.5D0*(RHOL+RHOR)
         RE = REYNOLDS_MASS(U(IM),D,A,MU)
         CALL FRICTION_SJ_VALUE_DERIVATIVE(RE,EPS,D,FD,DFDRE)
         DREDM = DREYNOLDS_DMDOT(U(IM),D,A,MU)
         DFDM = DFDRE*DREDM
         SF = FRICTION_SOURCE(FD,D,A,RHOF,U(IM))
         DSDM = DFRICTION_SOURCE_DMDOT(FD,DFDM,D,A,RHOF,
     &                                  U(IM))
         DSDP = -SF*CRHO/(2.0D0*RHOF)

         AB(KV+1+(ROW)-(IL),(IL)) = THETA*(-1.0D0/DX+DSDP)
         AB(KV+1+(ROW)-(IM),(IM)) = 1.0D0/(A*DT) + THETA*DSDM
         AB(KV+1+(ROW)-(IR),(IR)) = THETA*(1.0D0/DX+DSDP)

         IP = IR
         IF (I .LT. N-1) THEN
            AB(KV+1+(ROW+1)-(IP),(IP)) = A*DX*CRHO/DT
            AB(KV+1+(ROW+1)-(IM),(IM)) = -THETA
            AB(KV+1+(ROW+1)-(IM+2),(IM+2)) = THETA
         ELSE
            AB(KV+1+(ROW+1)-(IP),(IP)) = A*DX*CRHO/(2.0D0*DT)
            AB(KV+1+(ROW+1)-(IM),(IM)) = -THETA
         ENDIF
  100 CONTINUE

      RETURN
      END
