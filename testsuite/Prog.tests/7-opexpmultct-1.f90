! compile with
!gfortran -std=f2003  -I ../../../Prog_8/ -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ test1.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a
!
Program OPEXPMULTCTTEST

      Use Operator_mod

      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: U
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: V
      Complex (Kind=Kind(0.D0)), Dimension (:), Allocatable :: Z
      Complex (Kind=Kind(0.D0)), Dimension (5, 5) :: matnew, matold
      Complex (Kind=Kind(0.D0)) :: tmp, lexp
      Integer, Dimension (:), Allocatable :: P
      Integer :: i, j, n, m, opn, Ndim

      Ndim = 5
      Do opn = 1, 4
         Allocate (U(opn, opn), V(opn, Ndim), Z(opn), P(opn))
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
!
         Do i = 1, opn
            Do j = 1, opn
               U (i, j) = CMPLX (i, j, kind(0.D0))
            End Do
         End Do

         Call opexpmultct (V, U, P, matnew, Z, opn, Ndim)

! check against old version
         Do n = 1, opn
            lexp = Z (n)
            Do i = 1, Ndim
!             Mat(I,Op%P(n))  =  ExpMOp(n) * zdotu(Op%N,Op%U(1,n),1,VH(1,I),1)
               tmp = CMPLX (0.d0, 0.d0, kind(0.D0))
               Do m = 1, opn
                  tmp = tmp + V (m, i) * conjg (U(m, n))
               End Do
               matold (P(n), i) = lexp * tmp
            End Do
         End Do

         Do i = 1, Ndim
            Do j = 1, Ndim
               tmp = matold (i, j) - matnew (i, j)
               If (Aimag(tmp) > Abs(Aimag(matnew(i, j)))*1.D-15) Then
                  Write (*,*) "ERROR", matold (i, j), matnew (i, j)
                  Stop 2
               End If
               If (Real(tmp) > Abs(Real(matnew(i, j)))*1.D-15) Then
                  Write (*,*) "ERROR", matold (i, j), matnew (i, j)
                  Stop 2
               End If

            End Do
         End Do
         Deallocate (U, V, Z, P)
      End Do
write (*,*) "success"
End Program OPEXPMULTCTTEST
