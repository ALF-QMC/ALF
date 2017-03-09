! compile with
! gfortran -std=f2003  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a
!
Program CSR
!
      Use Operator_mod
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: Vnew, &
     & Vold
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: mat
      Integer, Dimension (:), Allocatable :: P
      Integer :: i, j, n, opn, NDim
!
      opn = 3
      NDim = 5
      Allocate (Vold(opn, NDim), Vnew(opn, NDim), P(opn), mat(NDim, &
     & NDim))
      Do i = 1, NDim
         Do j = 1, NDim
            mat (i, j) = CMPLX (i, j, kind(0.D0))
         End Do
      End Do
!
      Do i = 1, opn
         P (i) = i
      End Do
!
      Call copy_select_rows (Vnew, mat, P, opn, NDim)
!
! check old version
      Do n = 1, opn
         Do i = 1, NDim
            Vold (n, i) = mat (i, P(n))
         End Do
      End Do
!
      Do i = 1, 3
         Do j = 1, 3
            If (Vold(i, j) .Ne. Vnew(i, j)) Then
               Write (*,*) "ERROR", Vold (i, j), Vnew (i, j)
               Stop 2
            End If
         End Do
      End Do
      Deallocate (Vold, Vnew, P, mat)
End Program CSR
