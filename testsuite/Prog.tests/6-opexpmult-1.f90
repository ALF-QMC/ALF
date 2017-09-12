! compile with
! gfortran -std=f2003  -I ../../../Prog/ -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ 6-opexpmult-1.f90 ../../../Prog/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas

!
!
! This test seems to be a bit sensitive with respect to
! the precise Implementation
Program OPEXPMULTTEST
      Use Operator_mod
      Implicit None

!
!
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: U, Uold, mytmp, matnew, matold
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: V
      Complex (Kind=Kind(0.D0)), Dimension (:), Allocatable :: Z
      Complex (Kind=Kind(0.D0)) :: tmp, lexp
      Real(Kind = Kind(0.D0)), Dimension(:), allocatable :: E
      Integer, Dimension (:), Allocatable :: P
      Integer :: i, j, n, m, opn, Ndim
!
!
      Do opn = 1, 8
      Do Ndim = opn+1, 15
         Allocate (U(opn, opn), V(opn, Ndim), P(opn), Z(opn), Uold(opn, opn), E(opn), mytmp(opn, opn))
         Allocate(matnew(Ndim, Ndim), matold(ndim, ndim))
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
        if(opn > 1) then
         Do i = 1, opn
            Do j = 1, opn
               U (i, j) = CMPLX (i + j, j - i, kind(0.D0))
            End Do
         End Do
         CALL DIAG(U, mytmp, E)
         Do i = 1, opn
              Z (i) = Exp (CMPLX(E(i), i, kind(0.D0)))
         enddo
         U = mytmp
         else
            U(1, 1) = 1.D0
            Uold(1, 1) = 1.D0
            Z(1) = Exp (2.D0)
         endif
         if(opn > 1) then
            lexp = DET_C(mytmp, opn)
            DO I = 1, opn
                U(I, 1) = U(I, 1) / lexp
            ENDDO
            Uold = U
!          if(opn > 2) then
!             CALL ZGEQRF(opn, opn, U, opn, TAU, WORK, LWORK, INFO)
!             DO I = 1, opn-1
!                 U(I, opn) = TAU(I)
!             ENDDO
!             U(opn, opn) = 1.D0 - U(opn,opn)*(1.D0 - TAU(opn))! absorb last sign (stored in U(opn, opn)) into redefinition of tau
!          endif
         endif
         Call opexpmult (V, U, P, matnew, Z, opn, Ndim)
!
! check against old version
         Do n = 1, opn
            lexp = Z (n)
            Do i = 1, Ndim
!             Mat(I,Op%P(n))  =  ExpMOp(n) * zdotu(Op%N,Op%U(1,n),1,VH(1,I),1)
               tmp = CMPLX (0.d0, 0.d0, kind(0.D0))
               Do m = 1, opn
                  tmp = tmp + V (m, i) * Uold (m, n)
               End Do
               matold (i, P(n)) = lexp * tmp
            End Do
         End Do

!    write (*, *) "opn = ", opn
!     DO I = 1, Ndim
!         write (*, *) (matold(I, :))
!     ENDDO
! write (*,*) "================================"
!     DO I = 1, Ndim
!         write (*, *) (matnew(I, :))
!     ENDDO
         Do i = 1, Ndim
            Do j = 1, Ndim
               tmp = matold (i, j) - matnew (i, j)
               if(Abs(AIMAG(TMP)) > 1.D-14) THEN
               If (Abs(Aimag(tmp)) > Abs(Aimag(matnew(i, j)))*1.D-14) Then
                  Write (*,*) "ERROR in imag", matold (i, j), matnew (i, j)
                  Stop 2
               End If
               ENDIF
               if(Abs(AIMAG(TMP)) > 1.D-14) THEN
               If (Abs(Real(tmp)) > Abs(Real(matnew(i, j)))*1.D-14) Then
                  Write (*,*) "ERROR in real", matold (i, j), matnew (i, j)
                  Stop 3
               End If
               ENDIF
            End Do
         End Do
         Deallocate (U, V, P, Z, uold, E, mytmp, matnew, matold)
      End Do
      Enddo
      write (*,*) "success"
End Program OPEXPMULTTEST

