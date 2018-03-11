! compile with
!gfortran -std=f2003  -I ../../../Prog_8/ -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ test1.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a
!
Program OPEXPMULTCTTEST

      Use Operator_mod
      Use MyMats

      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: U, Uold, mytmp, matnew, matold
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: V
      Complex (Kind=Kind(0.D0)), Dimension (:), Allocatable :: Z
      Complex (Kind=Kind(0.D0)) :: tmp, lexp
      Real(Kind = Kind(0.D0)), Dimension(:), allocatable :: E
      Integer, Dimension (:), Allocatable :: P
      Integer :: i, j, n, m, opn, Ndim

      Do opn = 1, 8
      Do Ndim = opn +1, 15
         lwork = 2* opn
         Allocate (U(opn, opn), V(opn, Ndim), P(opn), Z(opn), Uold(opn, opn), E(opn), mytmp(opn, opn))
         Allocate(matnew(Ndim, Ndim), matold(Ndim, Ndim))
         Do i = 1, Ndim
            Do j = 1, Ndim
               matnew (i, j) = CMPLX (i, j, kind(0.D0))
               matold (i, j) = CMPLX (i, j, kind(0.D0))
            End Do
         End Do
!
         Do i = 1, opn
            Do j = 1, Ndim
               V (i, j) = CMPLX (i, j, kind(0.D0))
            End Do
            P (i) = i
            Z (i) = Exp (CMPLX(i, j, kind(0.D0)))
         End Do
        if(opn > 1) then
         Do i = 1, opn
            Do j = 1, opn
               U (i, j) = CMPLX (i + j, j - i, kind(0.D0))
            End Do
         End Do
         CALL DIAG(U, mytmp, E)
         U = mytmp
         else
         U(1, 1) = 1.D0
         Uold(1, 1) = 1.D0
         endif
         if(opn > 1) then
         lexp = DET_C(mytmp, opn)
         DO I = 1, opn
            U(I, 1) = U(I, 1) / lexp
         ENDDO
         Uold = U
         endif

         Call opexpmultct (V, U, P, matnew, Z, opn, Ndim)

! check against old version
         Do n = 1, opn
            lexp = Z (n)
            Do i = 1, Ndim
!             Mat(I,Op%P(n))  =  ExpMOp(n) * zdotu(Op%N,Op%U(1,n),1,VH(1,I),1)
               tmp = CMPLX (0.d0, 0.d0, kind(0.D0))
               Do m = 1, opn
                  tmp = tmp + V (m, i) * conjg (Uold(m, n))
               End Do
               matold (P(n), i) = lexp * tmp
            End Do
         End Do

!          write (*, *) "opn = ", opn
! write (*, *) (matold)
! write (*,*) "================================"
! write (*, *) (matnew)
         Do i = 1, Ndim
            Do j = 1, Ndim
               tmp = matold (i, j) - matnew (i, j)
               IF(ABS(AIMAG(tmp)) > 1D-11) THEN
               If (Abs(Aimag(tmp)) > Max(Abs(Aimag(matnew(i, j))), Abs(Aimag(matold(i, j))))*1D-13) Then
                  Write (*,*) "ERROR ", opn, matold (i, j), matnew (i, j), tmp
                  Stop 2
               End If
               ENDIF
               IF(ABS(DBLE(tmp)) > 1D-11) THEN
               If (Abs(DBLE(tmp)) > Abs(DBLE(matnew(i, j)))*1D-13) Then
                  Write (*,*) "ERROR ", matold (i, j), matnew (i, j)
                  Stop 3
               End If
               ENDIF
            End Do
         End Do
         Deallocate (U, V, P, Z, uold, E, mytmp, matnew, matold)
         ENDDO
      End Do
End Program OPEXPMULTCTTEST
