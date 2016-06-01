      SUBROUTINE COMQR2(NM,N,LOW,IGH,ORTR,ORTI,HR,HI,WR,WI,ZR,ZI,IERR)
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,NM,NN,IGH,IP1,
     X        ITN,ITS,LOW,LP1,ENM1,IEND,IERR
      DOUBLE PRECISION HR(NM,N),HI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N),
     X       ORTR(IGH),ORTI(IGH)
      DOUBLE PRECISION SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,TST1,TST2,
     X       PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE
C     ALGOL PROCEDURE  COMLR2, NUM. MATH. 16, 181-204(1970) BY PETERS
C     AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C     THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
C     (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A COMPLEX UPPER HESSENBERG MATRIX BY THE QR
C     METHOD.  THE EIGENVECTORS OF A COMPLEX GENERAL MATRIX
C     CAN ALSO BE FOUND IF  CORTH  HAS BEEN USED TO REDUCE
C     THIS GENERAL MATRIX TO HESSENBERG FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        ORTR AND ORTI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  CORTH, IF PERFORMED.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, SET ORTR(J) AND
C          ORTI(J) TO 0.0D0 FOR THESE ELEMENTS.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN FURTHER
C          INFORMATION ABOUT THE TRANSFORMATIONS WHICH WERE USED IN THE
C          REDUCTION BY  CORTH, IF PERFORMED.  IF THE EIGENVECTORS OF
C          THE HESSENBERG MATRIX ARE DESIRED, THESE ELEMENTS MAY BE
C          ARBITRARY.
C
C     ON OUTPUT
C
C        ORTR, ORTI, AND THE UPPER HESSENBERG PORTIONS OF HR AND HI
C          HAVE BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS.  THE EIGENVECTORS
C          ARE UNNORMALIZED.  IF AN ERROR EXIT IS MADE, NONE OF
C          THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS CSROOT FOR COMPLEX SQUARE ROOT.
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
C     .......... INITIALIZE EIGENVECTOR MATRIX ..........
      DO 101 J = 1, N
C
         DO 100 I = 1, N
            ZR(I,J) = 0.0D0
            ZI(I,J) = 0.0D0
  100    CONTINUE
         ZR(J,J) = 1.0D0
  101 CONTINUE
C     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
C                FROM THE INFORMATION LEFT BY CORTH ..........
      IEND = IGH - LOW - 1
      IF (IEND) 180, 150, 105
C     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
  105 DO 140 II = 1, IEND
         I = IGH - II
         IF (ORTR(I) .EQ. 0.0D0 .AND. ORTI(I) .EQ. 0.0D0) GO TO 140
         IF (HR(I,I-1) .EQ. 0.0D0 .AND. HI(I,I-1) .EQ. 0.0D0) GO TO 140
C     .......... NORM BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
         NORM = HR(I,I-1) * ORTR(I) + HI(I,I-1) * ORTI(I)
         IP1 = I + 1
C
         DO 110 K = IP1, IGH
            ORTR(K) = HR(K,I-1)
            ORTI(K) = HI(K,I-1)
  110    CONTINUE
C
         DO 130 J = I, IGH
            SR = 0.0D0
            SI = 0.0D0
C
            DO 115 K = I, IGH
               SR = SR + ORTR(K) * ZR(K,J) + ORTI(K) * ZI(K,J)
               SI = SI + ORTR(K) * ZI(K,J) - ORTI(K) * ZR(K,J)
  115       CONTINUE
C
            SR = SR / NORM
            SI = SI / NORM
C
            DO 120 K = I, IGH
               ZR(K,J) = ZR(K,J) + SR * ORTR(K) - SI * ORTI(K)
               ZI(K,J) = ZI(K,J) + SR * ORTI(K) + SI * ORTR(K)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
  150 L = LOW + 1
C
      DO 170 I = L, IGH
         LL = MIN0(I+1,IGH)
         IF (HI(I,I-1) .EQ. 0.0D0) GO TO 170
         NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
         YR = HR(I,I-1) / NORM
         YI = HI(I,I-1) / NORM
         HR(I,I-1) = NORM
         HI(I,I-1) = 0.0D0
C
         DO 155 J = I, N
            SI = YR * HI(I,J) - YI * HR(I,J)
            HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
            HI(I,J) = SI
  155    CONTINUE
C
         DO 160 J = 1, LL
            SI = YR * HI(J,I) + YI * HR(J,I)
            HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
            HI(J,I) = SI
  160    CONTINUE
C
         DO 165 J = LOW, IGH
            SI = YR * ZI(J,I) + YI * ZR(J,I)
            ZR(J,I) = YR * ZR(J,I) - YI * ZI(J,I)
            ZI(J,I) = SI
  165    CONTINUE
C
  170 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0D0
      TI = 0.0D0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 IF (EN .LT. LOW) GO TO 680
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         TST1 = DABS(HR(L-1,L-1)) + DABS(HI(L-1,L-1))
     X            + DABS(HR(L,L)) + DABS(HI(L,L))
         TST2 = TST1 + DABS(HR(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1)
      XI = HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0D0 .AND. XI .EQ. 0.0D0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0D0
      YI = (HI(ENM1,ENM1) - SI) / 2.0D0
      CALL CSROOT(YR**2-YI**2+XR,2.0D0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0D0) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
C     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = DABS(HR(EN,ENM1)) + DABS(HR(ENM1,EN-2))
      SI = 0.0D0
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
      LP1 = L + 1
C
      DO 500 I = LP1, EN
         SR = HR(I,I-1)
         HR(I,I-1) = 0.0D0
         NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR)
         XR = HR(I-1,I-1) / NORM
         WR(I-1) = XR
         XI = HI(I-1,I-1) / NORM
         WI(I-1) = XI
         HR(I-1,I-1) = NORM
         HI(I-1,I-1) = 0.0D0
         HI(I,I-1) = SR / NORM
C
         DO 490 J = I, N
            YR = HR(I-1,J)
            YI = HI(I-1,J)
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
            HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
            HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
            HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  490    CONTINUE
C
  500 CONTINUE
C
      SI = HI(EN,EN)
      IF (SI .EQ. 0.0D0) GO TO 540
      NORM = PYTHAG(HR(EN,EN),SI)
      SR = HR(EN,EN) / NORM
      SI = SI / NORM
      HR(EN,EN) = NORM
      HI(EN,EN) = 0.0D0
      IF (EN .EQ. N) GO TO 540
      IP1 = EN + 1
C
      DO 520 J = IP1, N
         YR = HR(EN,J)
         YI = HI(EN,J)
         HR(EN,J) = SR * YR + SI * YI
         HI(EN,J) = SR * YI - SI * YR
  520 CONTINUE
C     .......... INVERSE OPERATION (COLUMNS) ..........
  540 DO 600 J = LP1, EN
         XR = WR(J-1)
         XI = WI(J-1)
C
         DO 580 I = 1, J
            YR = HR(I,J-1)
            YI = 0.0D0
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            IF (I .EQ. J) GO TO 560
            YI = HI(I,J-1)
            HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  560       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  580    CONTINUE
C
         DO 590 I = LOW, IGH
            YR = ZR(I,J-1)
            YI = ZI(I,J-1)
            ZZR = ZR(I,J)
            ZZI = ZI(I,J)
            ZR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            ZI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
            ZR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            ZI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  590    CONTINUE
C
  600 CONTINUE
C
      IF (SI .EQ. 0.0D0) GO TO 240
C
      DO 630 I = 1, EN
         YR = HR(I,EN)
         YI = HI(I,EN)
         HR(I,EN) = SR * YR - SI * YI
         HI(I,EN) = SR * YI + SI * YR
  630 CONTINUE
C
      DO 640 I = LOW, IGH
         YR = ZR(I,EN)
         YI = ZI(I,EN)
         ZR(I,EN) = SR * YR - SI * YI
         ZI(I,EN) = SR * YI + SI * YR
  640 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 HR(EN,EN) = HR(EN,EN) + TR
      WR(EN) = HR(EN,EN)
      HI(EN,EN) = HI(EN,EN) + TI
      WI(EN) = HI(EN,EN)
      EN = ENM1
      GO TO 220
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  680 NORM = 0.0D0
C
      DO 720 I = 1, N
C
         DO 720 J = I, N
            TR = DABS(HR(I,J)) + DABS(HI(I,J))
            IF (TR .GT. NORM) NORM = TR
  720 CONTINUE
C
      IF (N .EQ. 1 .OR. NORM .EQ. 0.0D0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 2 DO -- ..........
      DO 800 NN = 2, N
         EN = N + 2 - NN
         XR = WR(EN)
         XI = WI(EN)
         HR(EN,EN) = 1.0D0
         HI(EN,EN) = 0.0D0
         ENM1 = EN - 1
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 780 II = 1, ENM1
            I = EN - II
            ZZR = 0.0D0
            ZZI = 0.0D0
            IP1 = I + 1
C
            DO 740 J = IP1, EN
               ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)
               ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)
  740       CONTINUE
C
            YR = XR - WR(I)
            YI = XI - WI(I)
            IF (YR .NE. 0.0D0 .OR. YI .NE. 0.0D0) GO TO 765
               TST1 = NORM
               YR = TST1
  760          YR = 0.01D0 * YR
               TST2 = NORM + YR
               IF (TST2 .GT. TST1) GO TO 760
  765       CONTINUE
            CALL CDIV(ZZR,ZZI,YR,YI,HR(I,EN),HI(I,EN))
C     .......... OVERFLOW CONTROL ..........
            TR = DABS(HR(I,EN)) + DABS(HI(I,EN))
            IF (TR .EQ. 0.0D0) GO TO 780
            TST1 = TR
            TST2 = TST1 + 1.0D0/TST1
            IF (TST2 .GT. TST1) GO TO 780
            DO 770 J = I, EN
               HR(J,EN) = HR(J,EN)/TR
               HI(J,EN) = HI(J,EN)/TR
  770       CONTINUE
C
  780    CONTINUE
C
  800 CONTINUE
C     .......... END BACKSUBSTITUTION ..........
      ENM1 = N - 1
C     .......... VECTORS OF ISOLATED ROOTS ..........
      DO  840 I = 1, ENM1
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
         IP1 = I + 1
C
         DO 820 J = IP1, N
            ZR(I,J) = HR(I,J)
            ZI(I,J) = HI(I,J)
  820    CONTINUE
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW+1 DO -- ..........
      DO 880 JJ = LOW, ENM1
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
         DO 880 I = LOW, IGH
            ZZR = 0.0D0
            ZZI = 0.0D0
C
            DO 860 K = LOW, M
               ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)
               ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)
  860       CONTINUE
C
            ZR(I,J) = ZZR
            ZI(I,J) = ZZI
  880 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
