      PROGRAM JACOBIANCHECK
      IMPLICIT NONE
      INTEGER N,NU,LDJ
      PARAMETER (N=4,NU=10,LDJ=10)
      INTEGER I,J,K,ICASE,NNZ,IMAX,JMAX
      DOUBLE PRECISION U(NU),UP(NU),UM(NU),R0(NU),RP(NU),RM(NU)
      DOUBLE PRECISION PO(0:N),MO(0:N-1),JA(LDJ,NU),JF(LDJ,NU)
      DOUBLE PRECISION PI,D,A,T,Z,RS,MU,EPS,DX,DT,THETA
      DOUBLE PRECISION MINO,MOUTO,MOUTN,PINN,SGN,H,EP,ERR,DEN
      DOUBLE PRECISION MAXREL,MAXZERO,MAXRES,AV,FV
      DOUBLE PRECISION RELTOL,ZEROTOL

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
      PINN = 8.04D6
      EP = EPSILON(1.0D0)**(1.0D0/3.0D0)
      RELTOL = 2.0D-7
      ZEROTOL = 2.0D-10

      PO(0) = 8.00D6
      PO(1) = 7.96D6
      PO(2) = 7.92D6
      PO(3) = 7.88D6
      PO(4) = 7.84D6

      WRITE(*,*)
      WRITE(*,*) 'TRANSIENT RESIDUAL JACOBIAN CHECK'
      WRITE(*,*) '---------------------------------'
      WRITE(*,*) ' intervals = ',N
      WRITE(*,*) ' unknowns  = ',NU
      WRITE(*,*) ' theta     = ',THETA

      DO 200 ICASE=1,2
         IF (ICASE .EQ. 1) THEN
            SGN = 1.0D0
         ELSE
            SGN = -1.0D0
         ENDIF

         MINO = SGN*101.0D0
         MOUTO = SGN*96.0D0
         MOUTN = SGN*95.0D0
         MO(0) = SGN*100.0D0
         MO(1) = SGN*99.0D0
         MO(2) = SGN*98.0D0
         MO(3) = SGN*97.0D0

         U(1) = 8.05D6
         U(2) = SGN*103.0D0
         U(3) = SGN*98.0D0
         U(4) = 7.99D6
         U(5) = SGN*101.0D0
         U(6) = 7.94D6
         U(7) = SGN*96.0D0
         U(8) = 7.90D6
         U(9) = SGN*99.0D0
         U(10) = 7.85D6

         CALL ASSEMBLE_RESIDUAL(N,U,PO,MO,MINO,MOUTO,MOUTN,
     &        PINN,DX,DT,THETA,D,A,T,Z,RS,MU,EPS,R0)
         CALL ASSEMBLE_JACOBIAN(N,U,DX,DT,THETA,D,A,T,Z,RS,
     &        MU,EPS,JA,LDJ)

         MAXRES = 0.0D0
         DO 20 I=1,NU
            MAXRES = DMAX1(MAXRES,DABS(R0(I)))
   20    CONTINUE

C        Independent centered finite-difference Jacobian.

         DO 60 J=1,NU
            DO 30 K=1,NU
               UP(K) = U(K)
               UM(K) = U(K)
   30       CONTINUE
            H = EP*DMAX1(DABS(U(J)),1.0D0)
            UP(J) = UP(J) + H
            UM(J) = UM(J) - H
            CALL ASSEMBLE_RESIDUAL(N,UP,PO,MO,MINO,MOUTO,MOUTN,
     &           PINN,DX,DT,THETA,D,A,T,Z,RS,MU,EPS,RP)
            CALL ASSEMBLE_RESIDUAL(N,UM,PO,MO,MINO,MOUTO,MOUTN,
     &           PINN,DX,DT,THETA,D,A,T,Z,RS,MU,EPS,RM)
            DO 50 I=1,NU
               JF(I,J) = (RP(I)-RM(I))/(2.0D0*H)
   50       CONTINUE
   60    CONTINUE

         MAXREL = 0.0D0
         MAXZERO = 0.0D0
         IMAX = 0
         JMAX = 0
         NNZ = 0
         DO 100 J=1,NU
            DO 90 I=1,NU
               AV = JA(I,J)
               FV = JF(I,J)
               IF (DABS(AV) .GT. 1.0D-14) THEN
                  NNZ = NNZ + 1
                  DEN = DMAX1(DABS(AV),DABS(FV))
                  ERR = DABS(AV-FV)/DEN
                  IF (ERR .GT. MAXREL) THEN
                     MAXREL = ERR
                     IMAX = I
                     JMAX = J
                  ENDIF
               ELSE
                  MAXZERO = DMAX1(MAXZERO,DABS(FV))
               ENDIF
   90       CONTINUE
  100    CONTINUE

         WRITE(*,*)
         IF (ICASE .EQ. 1) THEN
            WRITE(*,*) ' forward-flow non-equilibrium state'
         ELSE
            WRITE(*,*) ' reverse-flow non-equilibrium state'
         ENDIF
         WRITE(*,'(A,1PE12.4)') ' max residual magnitude = ',MAXRES
         WRITE(*,'(A,I6)') ' analytic nonzeros      = ',NNZ
         WRITE(*,'(A,1PE12.4)') ' max nonzero rel error  = ',MAXREL
         WRITE(*,'(A,1PE12.4)') ' max analytic-zero FD   = ',MAXZERO
         WRITE(*,'(A,2I4)') ' worst nonzero row,col  = ',IMAX,JMAX
         WRITE(*,'(A,1PE12.4)') ' analytic at worst      = ',
     &        JA(IMAX,JMAX)
         WRITE(*,'(A,1PE12.4)') ' finite diff at worst   = ',
     &        JF(IMAX,JMAX)

         IF (NNZ .NE. 6*N+3) THEN
            WRITE(*,*) 'FAIL: unexpected Jacobian sparsity'
            STOP 1
         ENDIF
         IF (MAXREL .GT. RELTOL) THEN
            WRITE(*,*) 'FAIL: analytic/FD Jacobian mismatch'
            STOP 1
         ENDIF
         IF (MAXZERO .GT. ZEROTOL) THEN
            WRITE(*,*) 'FAIL: unexpected off-pattern derivative'
            STOP 1
         ENDIF
  200 CONTINUE

      WRITE(*,*)
      WRITE(*,*) 'PASS: transient residual Jacobian verified'
      END
