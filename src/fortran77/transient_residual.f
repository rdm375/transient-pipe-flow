C=======================================================================
C  TRANSIENT RESIDUAL ASSEMBLY
C
C  Unknown ordering for N pipe intervals:
C
C    P0, MIN, M1/2, P1, M3/2, P2, ..., M(N-1/2), PN
C
C  Equation ordering:
C
C    inlet pressure BC, inlet continuity,
C    momentum 0, continuity 1, ..., momentum N-1, outlet continuity
C=======================================================================

      SUBROUTINE ASSEMBLE_RESIDUAL(N,U,PO,MO,MINO,MOUTO,MOUTN,
     &                             PINN,DX,DT,THETA,D,A,T,Z,RS,
     &                             MU,EPS,R)
      IMPLICIT NONE
      INTEGER N,I,IP,IL,IR,IM,ROW
      DOUBLE PRECISION U(*),PO(0:*),MO(0:*),R(*)
      DOUBLE PRECISION MINO,MOUTO,MOUTN,PINN,DX,DT,THETA
      DOUBLE PRECISION D,A,T,Z,RS,MU,EPS
      DOUBLE PRECISION COLD,CNEW,RHOL,RHOR,RHOF,RE,FD,SF
      DOUBLE PRECISION RHOL0,RHOR0,RHOF0,RE0,FD0,SF0
      DOUBLE PRECISION EOS_CZ_RHO,REYNOLDS_MASS,FRICTION_SJ
      DOUBLE PRECISION FRICTION_SOURCE

      COLD = 1.0D0 - THETA
      CNEW = THETA

C     Inlet pressure boundary condition.

      R(1) = U(1) - PINN

C     Inlet half control volume.

      R(2) = A*DX/(2.0D0*DT)
     &     *(EOS_CZ_RHO(U(1),T,Z,RS)-EOS_CZ_RHO(PO(0),T,Z,RS))
     &     + CNEW*(U(3)-U(2)) + COLD*(MO(0)-MINO)

      DO 100 I=0,N-1
         IM = 2*I + 3
         IF (I .EQ. 0) THEN
            IL = 1
         ELSE
            IL = 2*I + 2
         ENDIF
         IR = 2*I + 4
         ROW = 2*I + 3

C        New-time friction source.

         RHOL = EOS_CZ_RHO(U(IL),T,Z,RS)
         RHOR = EOS_CZ_RHO(U(IR),T,Z,RS)
         RHOF = 0.5D0*(RHOL+RHOR)
         RE = REYNOLDS_MASS(U(IM),D,A,MU)
         FD = FRICTION_SJ(RE,EPS,D)
         SF = FRICTION_SOURCE(FD,D,A,RHOF,U(IM))

C        Old-time friction source.

         RHOL0 = EOS_CZ_RHO(PO(I),T,Z,RS)
         RHOR0 = EOS_CZ_RHO(PO(I+1),T,Z,RS)
         RHOF0 = 0.5D0*(RHOL0+RHOR0)
         RE0 = REYNOLDS_MASS(MO(I),D,A,MU)
         FD0 = FRICTION_SJ(RE0,EPS,D)
         SF0 = FRICTION_SOURCE(FD0,D,A,RHOF0,MO(I))

C        Reduced momentum equation.

         R(ROW) = (U(IM)-MO(I))/(A*DT)
     &          + CNEW*((U(IR)-U(IL))/DX+SF)
     &          + COLD*((PO(I+1)-PO(I))/DX+SF0)

C        Continuity equation to the right of this face.  The last
C        equation is the outlet half control volume.

         IP = IR
         IF (I .LT. N-1) THEN
            R(ROW+1) = A*DX/DT
     &         *(EOS_CZ_RHO(U(IP),T,Z,RS)
     &          -EOS_CZ_RHO(PO(I+1),T,Z,RS))
     &         + CNEW*(U(IM+2)-U(IM))
     &         + COLD*(MO(I+1)-MO(I))
         ELSE
            R(ROW+1) = A*DX/(2.0D0*DT)
     &         *(EOS_CZ_RHO(U(IP),T,Z,RS)
     &          -EOS_CZ_RHO(PO(N),T,Z,RS))
     &         + CNEW*(MOUTN-U(IM))
     &         + COLD*(MOUTO-MO(I))
         ENDIF
  100 CONTINUE

      RETURN
      END
