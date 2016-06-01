      SUBROUTINE TRED1(NM,N,A,D,E,E2)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      DOUBLE PRECISION A(NM,N),D(N),E(N),E2(N)
      DOUBLE PRECISION F,G,H,SCALE
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      DO 100 I = 1, N
         D(I) = A(N,I)
         A(N,I) = A(I,I)
  100 CONTINUE
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(D(K))
C
         IF (SCALE .NE. 0.0D0) GO TO 140
C
         DO 125 J = 1, L
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = 0.0D0
  125    CONTINUE
C
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 300
C
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         IF (L .EQ. 1) GO TO 285
C     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0D0
C
         DO 240 J = 1, L
            F = D(J)
            G = E(J) + A(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + A(K,J) * D(K)
               E(K) = E(K) + A(K,J) * F
  200       CONTINUE
C
  220       E(J) = G
  240    CONTINUE
C     .......... FORM P ..........
         F = 0.0D0
C
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
C
         H = F / (H + H)
C     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - H * D(J)
C     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
C
            DO 260 K = J, L
  260       A(K,J) = A(K,J) - F * E(K) - G * D(K)
C
  280    CONTINUE
C
  285    DO 290 J = 1, L
            F = D(J)
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = F * SCALE
  290    CONTINUE
C
  300 CONTINUE
C
      RETURN
      END
