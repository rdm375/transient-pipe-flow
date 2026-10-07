C=======================================================================
C  ISOTHERMAL PIPE SIMULATOR
C
C  First numerical validation:
C    Constant-Z EOS
C    Constant Darcy friction factor
C    Reduced momentum equation
C    Staggered spatial grid
C    Analytic steady-state pressure-squared solution
C
C  Purpose:
C    Assemble the discrete residual at the exact analytic steady state.
C    No nonlinear solver is used here.
C=======================================================================

      PROGRAM STEADYRES

      IMPLICIT NONE

      INTEGER NMAX
      PARAMETER (NMAX=1000)

      INTEGER N,I,NARG
      CHARACTER*32 ARG

      DOUBLE PRECISION L,D,A,DX
      DOUBLE PRECISION TEMP,Z,RS,FD
      DOUBLE PRECISION PIN,MDOT0
      DOUBLE PRECISION PI
      DOUBLE PRECISION COEF

      DOUBLE PRECISION P(0:NMAX)
      DOUBLE PRECISION RHO(0:NMAX)
      DOUBLE PRECISION MF(0:NMAX-1)

      DOUBLE PRECISION RC(0:NMAX)
      DOUBLE PRECISION RM(0:NMAX-1)

      DOUBLE PRECISION RHOF,SF,PG
      DOUBLE PRECISION MINLET,MOUTLET
      DOUBLE PRECISION MAXC,MAXM,MAXR
      DOUBLE PRECISION RELM,MAXRELM,DENOM
      DOUBLE PRECISION RELTOL
      DOUBLE PRECISION PSQ

C----- Problem definition.

      PI = 4.0D0*DATAN(1.0D0)

      N    = 100

C----- Optional command-line grid size.

      NARG = COMMAND_ARGUMENT_COUNT()

      IF (NARG .GE. 1) THEN
         CALL GET_COMMAND_ARGUMENT(1,ARG)
         READ(ARG,*,ERR=900) N
      ENDIF

      IF (N .LT. 1 .OR. N .GT. NMAX) THEN
         WRITE(*,*) 'ERROR: N must satisfy 1 <= N <= ',NMAX
         STOP 1
      ENDIF

      L    = 100000.0D0
      D    = 1.0D0
      A    = PI*D*D/4.0D0
      DX   = L/DBLE(N)

      TEMP = 288.15D0
      Z    = 0.90D0
      RS   = 500.0D0
      FD   = 0.010D0

      PIN   = 8.0D6
      MDOT0 = 100.0D0

C----- Analytic steady pressure-squared coefficient:
C
C      p(x)^2 = pin^2 - coef*x
C
C      coef = f_D Z R_s T mdot^2 / (D A^2)

      COEF = FD*Z*RS*TEMP*MDOT0*MDOT0/(D*A*A)

C----- Dimensionless steady momentum validation tolerance.

      RELTOL = 1.0D-10

C----- Construct exact analytic pressure field and EOS density.

      DO 100 I=0,N

         PSQ = PIN*PIN - COEF*DBLE(I)*DX

         IF (PSQ .LE. 0.0D0) THEN
            WRITE(*,*) 'ERROR: non-positive pressure squared at node ',I
            STOP 1
         ENDIF

         P(I)   = DSQRT(PSQ)
         RHO(I) = P(I)/(Z*RS*TEMP)

  100 CONTINUE

C----- Uniform steady mass flow at all interior staggered faces.

      DO 110 I=0,N-1
         MF(I) = MDOT0
  110 CONTINUE

      MINLET  = MDOT0
      MOUTLET = MDOT0

C=======================================================================
C  CONTINUITY RESIDUALS
C
C  These are integral control-volume residuals in kg/s.
C  At steady state all density time derivatives vanish.
C=======================================================================

C----- Inlet half control volume.

      RC(0) = MF(0) - MINLET

C----- Interior control volumes.

      DO 200 I=1,N-1
         RC(I) = MF(I) - MF(I-1)
  200 CONTINUE

C----- Outlet half control volume.

      RC(N) = MOUTLET - MF(N-1)

C=======================================================================
C  REDUCED MOMENTUM RESIDUALS
C
C  R_m = (p_{i+1}-p_i)/dx + S_f
C
C  with arithmetic face density
C
C      rho_f = (rho_i + rho_{i+1})/2
C
C  and
C
C      S_f = f_D mdot |mdot| / (2 D A^2 rho_f)
C=======================================================================

      DO 300 I=0,N-1

         RHOF = 0.5D0*(RHO(I)+RHO(I+1))

         SF = FD*MF(I)*DABS(MF(I))
     &        /(2.0D0*D*A*A*RHOF)

         RM(I) = (P(I+1)-P(I))/DX + SF

  300 CONTINUE

C=======================================================================
C  RESIDUAL NORMS
C=======================================================================

      MAXC = 0.0D0

      DO 400 I=0,N
         MAXC = DMAX1(MAXC,DABS(RC(I)))
  400 CONTINUE

      MAXM    = 0.0D0
      MAXRELM = 0.0D0

      DO 410 I=0,N-1

         MAXM = DMAX1(MAXM,DABS(RM(I)))

         RHOF = 0.5D0*(RHO(I)+RHO(I+1))

         SF = FD*MF(I)*DABS(MF(I))
     &        /(2.0D0*D*A*A*RHOF)

         PG = (P(I+1)-P(I))/DX

         DENOM = DMAX1(DABS(PG),DABS(SF))

         IF (DENOM .GT. 0.0D0) THEN
            RELM = DABS(PG+SF)/DENOM
         ELSE
            RELM = 0.0D0
         ENDIF

         MAXRELM = DMAX1(MAXRELM,RELM)

  410 CONTINUE

      MAXR = DMAX1(MAXC,MAXM)

C=======================================================================
C  REPORT
C=======================================================================

      WRITE(*,*)
      WRITE(*,*) 'ISOTHERMAL PIPE STEADY RESIDUAL TEST'
      WRITE(*,*) '------------------------------------'
      WRITE(*,'(A,I8)')    ' intervals             = ',N
      WRITE(*,'(A,ES24.16)') ' dx [m]                = ',DX
      WRITE(*,'(A,ES24.16)') ' inlet pressure [Pa]   = ',P(0)
      WRITE(*,'(A,ES24.16)') ' outlet pressure [Pa]  = ',P(N)
      WRITE(*,'(A,ES24.16)') ' mass flow [kg/s]      = ',MDOT0
      WRITE(*,*)
      WRITE(*,'(A,ES24.16)') ' max continuity resid  = ',MAXC
      WRITE(*,'(A,ES24.16)') ' max momentum resid    = ',MAXM
      WRITE(*,'(A,ES24.16)') ' max residual          = ',MAXR
      WRITE(*,'(A,ES24.16)') ' max relative mom resid= ',MAXRELM
      WRITE(*,*)

C----- The momentum residual is formed by subtracting two quantities
C     evaluated independently in floating point, so exact bitwise zero
C     is neither required nor expected.  This first executable reports
C     the residual rather than imposing an arbitrary absolute tolerance.
C     Validation instead uses a dimensionless residual relative to the
C     pressure-gradient and friction terms being balanced.

      IF (MAXC .GT. 0.0D0) THEN
         WRITE(*,*) 'FAIL: steady continuity residual is nonzero'
         STOP 2
      ENDIF

      IF (MAXRELM .GT. RELTOL) THEN
         WRITE(*,*) 'FAIL: relative momentum residual exceeds tolerance'
         STOP 3
      ENDIF

      WRITE(*,*) 'PASS: analytic steady state satisfies'
      WRITE(*,*) '      discrete equations'

      STOP

  900 CONTINUE
      WRITE(*,*) 'ERROR: invalid interval count'
      STOP 1

      END
