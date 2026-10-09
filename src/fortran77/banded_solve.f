C=======================================================================
C  BANDED LINEAR SOLVER WITH PARTIAL PIVOTING
C
C  Input matrix A is stored in AB:
C
C    AB(KL+KU+1+I-J,J) = A(I,J)
C
C  Original bandwidth: KL lower, KU upper.
C  Storage accommodates KL additional upper diagonals for pivot fill.
C  LDAB >= 2*KL+KU+1.
C
C  AB and X are overwritten.
C  INFO=0 success; INFO>0 pivot failure; INFO<0 invalid input.
C
C  Uses the same absolute pivot threshold as SOLVE_DENSE.
C=======================================================================

      SUBROUTINE SOLVE_BANDED(N,AB,LDAB,KL,KU,B,X,INFO)
      IMPLICIT NONE
      INTEGER N,LDAB,KL,KU,INFO
      INTEGER I,J,K,IP,IMAX,JMAX,KV
      DOUBLE PRECISION AB(LDAB,*),B(*),X(*)
      DOUBLE PRECISION AMAX,TMP,FACTOR,SUM,PIVTOL
      DOUBLE PRECISION BAND_GET
      EXTERNAL BAND_GET

      INFO = 0
      IF (N.LT.1.OR.KL.LT.0.OR.KU.LT.0) THEN
         INFO = -1
         RETURN
      ENDIF
      IF (LDAB.LT.2*KL+KU+1) THEN
         INFO = -2
         RETURN
      ENDIF

      KV = KL+KU
      PIVTOL = 100.0D0*EPSILON(1.0D0)

      DO 10 I=1,N
         X(I) = B(I)
   10 CONTINUE

      DO 100 K=1,N-1
         IMAX = MIN0(N,K+KL)
         JMAX = MIN0(N,K+KV)

C        Pivot search within the lower bandwidth.

         IP = K
         AMAX = DABS(BAND_GET(AB,LDAB,KV,K,K))
         DO 20 I=K+1,IMAX
            TMP = DABS(BAND_GET(AB,LDAB,KV,I,K))
            IF (TMP.GT.AMAX) THEN
               AMAX = TMP
               IP = I
            ENDIF
   20    CONTINUE

         IF (AMAX.LE.PIVTOL) THEN
            INFO = K
            RETURN
         ENDIF

C        Swap active row segments only.
C        Columns before K have already been eliminated.

         IF (IP.NE.K) THEN
            DO 30 J=K,JMAX
               TMP = BAND_GET(AB,LDAB,KV,K,J)
               CALL BAND_SET(AB,LDAB,KV,K,J,
     &                       BAND_GET(AB,LDAB,KV,IP,J))
               CALL BAND_SET(AB,LDAB,KV,IP,J,TMP)
   30       CONTINUE
            TMP = X(K)
            X(K) = X(IP)
            X(IP) = TMP
         ENDIF

C        Eliminate entries below the pivot.

         DO 50 I=K+1,IMAX
            FACTOR = BAND_GET(AB,LDAB,KV,I,K)/
     &               BAND_GET(AB,LDAB,KV,K,K)
            CALL BAND_SET(AB,LDAB,KV,I,K,0.0D0)

            DO 40 J=K+1,JMAX
               TMP = BAND_GET(AB,LDAB,KV,I,J) -
     &                FACTOR*BAND_GET(AB,LDAB,KV,K,J)
               CALL BAND_SET(AB,LDAB,KV,I,J,TMP)
   40       CONTINUE

            X(I) = X(I)-FACTOR*X(K)
   50    CONTINUE
  100 CONTINUE

      IF (DABS(BAND_GET(AB,LDAB,KV,N,N)).LE.PIVTOL) THEN
         INFO = N
         RETURN
      ENDIF

C     Back substitution over the filled upper band.

      DO 120 I=N,1,-1
         SUM = X(I)
         DO 110 J=I+1,MIN0(N,I+KV)
            SUM = SUM-BAND_GET(AB,LDAB,KV,I,J)*X(J)
  110    CONTINUE
         X(I) = SUM/BAND_GET(AB,LDAB,KV,I,I)
  120 CONTINUE

      RETURN
      END

C=======================================================================
C  BAND ACCESS HELPERS
C=======================================================================

      DOUBLE PRECISION FUNCTION BAND_GET(AB,LDAB,KV,I,J)
      IMPLICIT NONE
      INTEGER LDAB,KV,I,J,ROW
      DOUBLE PRECISION AB(LDAB,*)

      ROW = KV+1+I-J
      IF (ROW.LT.1.OR.ROW.GT.LDAB) THEN
         BAND_GET = 0.0D0
      ELSE
         BAND_GET = AB(ROW,J)
      ENDIF
      RETURN
      END

      SUBROUTINE BAND_SET(AB,LDAB,KV,I,J,VALUE)
      IMPLICIT NONE
      INTEGER LDAB,KV,I,J,ROW
      DOUBLE PRECISION AB(LDAB,*),VALUE

      ROW = KV+1+I-J
      IF (ROW.GE.1.AND.ROW.LE.LDAB) THEN
         AB(ROW,J) = VALUE
      ENDIF
      RETURN
      END
