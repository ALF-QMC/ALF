! compile with
!  gfortran -std=f2003  -I ../../../Prog/ -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ 7-opexpmultct-2.f90
!
!
Program OPEXPMULTCTTEST
!
      Use Operator_mod
      Use MyMats
      Implicit None
!
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: U, V, Uold, mytmp
      Complex (Kind=Kind(0.D0)), Dimension (:), Allocatable :: Z
      Complex (Kind=Kind(0.D0)), Dimension (5, 5) :: matnew, matold
      Complex (Kind=Kind(0.D0)) :: tmp, Z1, g, lexp
      Real(Kind = Kind(0.D0)) :: spin
      Real(Kind = Kind(0.D0)), Dimension(:), allocatable :: E
      Integer, Dimension (:), Allocatable :: P
      Integer :: i, j, n, m, opn, Ndim
!
      Ndim = 5
      Do opn = 1, 4
         Allocate (U(opn, opn), V(opn, Ndim), P(opn), Z(opn), Uold(opn, opn), E(opn), mytmp(opn, opn))
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
            E(1) = 1.D0
         endif
         if(opn > 1) then
         lexp = DET_C(mytmp, opn)
         DO I = 1, opn
            U(I, 1) = U(I, 1) / lexp
         ENDDO
         Uold = U
         endif
!
         Z = Exp (-g*spin*E)
         Call opexpmultct (V, U, P, matnew, Z, opn, Ndim)
!
! check against old version
         Do n = 1, opn
            Z1 = Exp (-g*CMPLX(E(n)*spin, 0.d0, kind(0.D0)))
            Do i = 1, Ndim
               tmp = CMPLX (0.d0, 0.d0, kind(0.D0))
               Do m = 1, opn
                  tmp = tmp + V (m, i) * conjg (Uold(m, n))
               End Do
               matold (P(n), i) = tmp * Z1
            End Do
         End Do
!
!          write (*, *) "opn = ", opn
! write (*, *) (matold)
! write (*,*) "================================"
! write (*, *) (matnew)
         Do i = 1, Ndim
            Do j = 1, Ndim
               tmp = matold (i, j) - matnew (i, j)
               IF(ABS(AIMAG(tmp)) > 1.D-14) THEN
               If (Abs(Aimag(tmp)) > Max(Abs(Aimag(matnew(i, j))), Abs(Aimag(matold(i, j))))*1D-13) Then
                  Write (*,*) "ERROR ", opn, matold (i, j), matnew (i, j), tmp
                  Stop 2
               End If
               ENDIF
            End Do
         End Do
         Deallocate (U, V, P, Z, uold, E, mytmp)
      End Do
      write (*,*) "success"
End Program OPEXPMULTCTTEST
