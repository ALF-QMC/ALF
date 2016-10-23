! compile with
! gfortran -std=f2003  -I ../../../Prog_8/ -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ test2.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a

!
!
Program OPEXPMULTTEST
!
      Use Operator_mod
      Implicit None
!
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: U
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: V
      Complex (Kind=Kind(0.D0)), Dimension (:), Allocatable :: Z
      Complex (Kind=Kind(0.D0)), Dimension (5, 5) :: matnew, matold
      Complex (Kind=Kind(0.D0)) :: tmp, Z1, g
      Real (Kind=Kind(0.D0)), Dimension (:), Allocatable :: E
      Real :: spin
      Integer, Dimension (:), Allocatable :: P
      Integer :: i, j, n, m, opn, Ndim
!
      Ndim = 5
      Do opn = 1, 4
         Allocate (U(opn, opn), V(opn, Ndim), Z(opn), P(opn), E(opn))
         spin = - 1.D0
         g = 2.D0
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
            E (i) = 3.5D1 * i
         End Do
!
         Do i = 1, opn
            Do j = 1, opn
               U (i, j) = CMPLX (i, j, kind(0.D0))
            End Do
         End Do
!
         Z = Exp (g*spin*E)
         Call opexpmult (V, U, P, matnew, Z, opn, Ndim)
!
! check against old version
         Do n = 1, opn
            Z1 = Exp (g*CMPLX(E(n)*spin, 0.d0, kind(0.D0)))
            Do i = 1, Ndim
               tmp = CMPLX (0.d0, 0.d0, kind(0.D0))
               Do m = 1, opn
                  tmp = tmp + V (m, i) * U (m, n)
               End Do
               matold (i, P(n)) = tmp * Z1
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
         Deallocate (U, V, Z, E, P)
      End Do
End Program OPEXPMULTTEST
