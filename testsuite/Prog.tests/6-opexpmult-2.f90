! compile with
! gfortran -std=f2003  -I ../../../Prog_8/ -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ test2.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a

!
!
Program OPEXPMULTTEST
!
      Use Operator_mod
      Use MyMats
      Implicit None
!
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: U, Uold, mytmp
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: V
      Complex (Kind=Kind(0.D0)), Dimension (:), Allocatable :: Z
      Complex (Kind=Kind(0.D0)), Dimension (5, 5) :: matnew, matold
      Complex (Kind=Kind(0.D0)) :: tmp, lexp, g, Z1
      Real(Kind = Kind(0.D0)), Dimension(:), allocatable :: E
      Real(Kind = Kind(0.D0)) :: spin
      Complex(Kind=Kind(0.D0)), allocatable, dimension(:) :: work, TAU
      Integer, Dimension (:), Allocatable :: P
      Integer :: i, j, n, m, opn, Ndim, info, lwork
!
      Ndim = 5
      Do opn = 1, 4
         lwork = 2* opn
         Allocate (U(opn, opn), V(opn, 5), P(opn), Z(opn), work(lwork), Tau(opn), Uold(opn, opn), E(opn), mytmp(opn, opn))
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
         TAU = 0.D0
         if(opn > 2) then
            CALL ZGEQRF(opn, opn, U, opn, TAU, WORK, LWORK, INFO)
            DO I = 1, opn-1
                U(I, opn) = TAU(I)
            ENDDO
            U(opn, opn) = 1.D0 - U(opn,opn)*(1.D0 - TAU(opn))! absorb last sign (stored in U(opn, opn)) into redefinition of tau
         endif
         endif
         Z = Exp(g*spin*E)
         Call opexpmult (V, U, P, matnew, Z, opn, Ndim)
!
! check against old version
         Do n = 1, opn
            Z1 = Exp (g*CMPLX(E(n)*spin, 0.d0, kind(0.D0)))
            Do i = 1, Ndim
               tmp = CMPLX (0.d0, 0.d0, kind(0.D0))
               Do m = 1, opn
                  tmp = tmp + V (m, i) * Uold (m, n)
               End Do
               matold (i, P(n)) = tmp * Z1
            End Do
         End Do
!
!     write (*, *) "opn = ", opn
!     DO I = 1, Ndim
!         write (*, *) (matold(I, :))
!     ENDDO
! write (*,*) "================================"
!     DO I = 1, Ndim
!         write (*, *) (matnew(I, :))
!     ENDDO
         Do i = 1, Ndim
            Do j = 1, Ndim
            if (ABS(matold(i, j) - matnew(i, j)) > 1.D-14) then
               If (ABS(matold(i, j) - matnew(i, j)) > ABS(matold(i, j))*1D-15) Then
               write (*,*) "ERROR", matold(i, j), matnew(i, j)
                  Stop 2
               End If
               End if
            End Do
         End Do
         Deallocate (U, V, P, Z, work, uold, E, mytmp, tau)
      End Do
      write (*,*) "success"
End Program OPEXPMULTTEST
