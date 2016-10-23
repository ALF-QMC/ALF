! compile with
! gfortran -std=f2003  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a
!
Program OPMULTTEST
!
      Use Operator_mod
!
!
      Complex (Kind=Kind(0.D0)), Dimension (:, :), allocatable :: U
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: V
      Complex (Kind=Kind(0.D0)), Dimension (5, 5) :: matnew, matold
      Complex (Kind=Kind(0.D0)) :: tmp
      Integer, Dimension (:), Allocatable :: P
      Integer :: i, j, n, opn, Ndim
!
      Do opn = 1, 4
         Allocate (V(opn, 5), P(opn), U(opn, opn))
         Ndim = 5
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
         End Do
!
         Do i = 1, opn
            Do j = 1, opn
               U (i, j) = CMPLX (i, j, kind(0.D0))
            End Do
         End Do
!
         Call opmult (V, U, P, matnew, opn, Ndim)
!
! check against old version
         Do n = 1, opn
            Do i = 1, Ndim
               tmp = CMPLX (0.d0, 0.d0, kind(0.D0))
               Do m = 1, opn
                  tmp = tmp + U (n, m) * V (m, i)
               End Do
               matold (P(n), i) = tmp
            End Do
         End Do
!
         Do i = 1, Ndim
            Do j = 1, Ndim
               If (matold(i, j) .Ne. matnew(i, j)) Then
                  Stop 2
               End If
            End Do
         End Do
         Deallocate (V, P, U)
      End Do
End Program OPMULTTEST
