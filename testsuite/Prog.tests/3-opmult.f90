! compile with
! gfortran -std=f2003  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas
!
Program OPMULTTEST
!
      Use Operator_mod
      Use MyMats
!
!
      Complex (Kind=Kind(0.D0)), Dimension (:, :), allocatable :: U, Uold, mytmp
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: V
      Complex (Kind=Kind(0.D0)), Dimension (5, 5) :: matnew, matold
      Complex (Kind=Kind(0.D0)) :: tmp, Z
      REAL(Kind= Kind(0.D0)), Dimension(:), allocatable :: E
      Integer, Dimension (:), Allocatable :: P
      Integer :: i, j, n, opn, Ndim
!
      Do opn = 1, 4
         Allocate (V(opn, 5), P(opn), U(opn, opn), Uold(opn, opn), E(opn), myTMP(opn, opn))
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
        if(opn > 1) then
         Do i = 1, opn
            Do j = 1, opn
               U (i, j) = CMPLX (i+j, j-i, kind(0.D0))
            End Do
         End Do
         CALL DIAG(U, mytmp, E)
         U = mytmp

         else
            U(1, 1) = 1.D0
            Uold(1, 1) = 1.D0
         endif
         if (opn > 1) then
         Z = Det_C(mytmp, opn)
         
         DO I = 1, opn
            U(I, 1) = U(I, 1) / Z
         ENDDO
         Uold = U
        endif
         Call opmult (V, U, P, matnew, opn, Ndim)
!
! check against old version
!write (*,*) Uold
         Do n = 1, opn
            Do i = 1, Ndim
               tmp = CMPLX (0.d0, 0.d0, kind(0.D0))
               Do m = 1, opn
                  tmp = tmp + Uold (n, m) * V (m, i)
               End Do
               matold (P(n), i) = tmp
            End Do
         End Do
! write (*, *) "opn = ", opn
! write (*, *) (matold)
! write (*,*) "================================"
! write (*, *) (matnew)
         Do i = 1, Ndim
            Do j = 1, Ndim
               If (Abs(matold(i, j) - matnew(i, j)) > MAX(ABS(matold(i, j)), ABS(matnew(i, j)))*1D-14) Then
               write (*,*) "ERROR"
                  Stop 2
               End If
            End Do
         End Do
         Deallocate (V, P, U, uold, E, mytmp)
      End Do
      write (*,*) "success"
End Program OPMULTTEST

