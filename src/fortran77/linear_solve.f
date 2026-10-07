C=======================================================================
C  DENSE REFERENCE LINEAR SOLVER WITH PARTIAL PIVOTING
C=======================================================================

      SUBROUTINE SOLVE_DENSE(N,A,LDA,B,X,INFO)
      IMPLICIT NONE
      INTEGER N,LDA,INFO,I,J,K,IP
      DOUBLE PRECISION A(LDA,*),B(*),X(*)
      DOUBLE PRECISION AMAX,TMP,FACTOR,SUM,PIVTOL

      INFO = 0
      PIVTOL = 100.0D0*EPSILON(1.0D0)

      DO 10 I=1,N
         X(I) = B(I)
   10 CONTINUE

      DO 50 K=1,N-1
         IP = K
         AMAX = DABS(A(K,K))
         DO 20 I=K+1,N
            IF (DABS(A(I,K)) .GT. AMAX) THEN
               AMAX = DABS(A(I,K))
               IP = I
            ENDIF
   20    CONTINUE
         IF (AMAX .LE. PIVTOL) THEN
            INFO = K
            RETURN
         ENDIF
         IF (IP .NE. K) THEN
            DO 30 J=K,N
               TMP = A(K,J)
               A(K,J) = A(IP,J)
               A(IP,J) = TMP
   30       CONTINUE
            TMP = X(K)
            X(K) = X(IP)
            X(IP) = TMP
         ENDIF
         DO 40 I=K+1,N
            FACTOR = A(I,K)/A(K,K)
            A(I,K) = 0.0D0
            DO 35 J=K+1,N
               A(I,J) = A(I,J) - FACTOR*A(K,J)
   35       CONTINUE
            X(I) = X(I) - FACTOR*X(K)
   40    CONTINUE
   50 CONTINUE

      IF (DABS(A(N,N)) .LE. PIVTOL) THEN
         INFO = N
         RETURN
      ENDIF

      DO 70 I=N,1,-1
         SUM = X(I)
         DO 60 J=I+1,N
            SUM = SUM - A(I,J)*X(J)
   60    CONTINUE
         X(I) = SUM/A(I,I)
   70 CONTINUE

      RETURN
      END
