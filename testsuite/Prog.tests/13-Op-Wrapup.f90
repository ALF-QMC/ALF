! compile with
! gfortran -Wall -std=f2003 -I ../../../Prog_8/  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a
!
Program Wrapup
!
      Use Operator_mod
      Implicit None
!
      Interface
         Subroutine Op_WrapupFFA (Mat, Op, spin, Ndim, N_Type)
            Use Operator_mod
            Type (Operator), Intent (In) :: Op
            Complex (Kind=8), Intent (Inout) :: Mat (Ndim, Ndim)
            Real (Kind=8), Intent (In) :: spin
            Integer, Intent (In) :: N_Type, Ndim
         End Subroutine
      End Interface
!
      Complex (Kind=Kind(0.D0)) :: Zre, Zim
      Real (Kind=Kind(0.D0)) :: spin
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: VH, &
     & matnew, matold
      Complex (Kind=Kind(0.D0)), Dimension (:), Allocatable :: Expop, &
     & ExpMop
      Integer :: i, n, j, Ndim, N_Type, opn
      Type (Operator) :: Op
!
! setup some test data
      Ndim = 5
!
      Do opn = 1, 4
         Do N_Type = 1, 2
            Allocate (VH(opn, Ndim), matold(Ndim, Ndim), matnew(Ndim, &
           & Ndim), Expop(opn), ExpMop(opn))
            Call Op_seths ()
            Call Op_make (Op, opn)
!
            Do i = 1, Op%n
               Op%P (i) = i
               Do n = 1, Op%n
                  Op%O (i, n) = CMPLX (n+i, n - i, kind(0.D0))
               End Do
            End Do
!
            Op%g = 2.D0
            Op%alpha = 0.D0
            Call Op_set (Op)
!
            spin = - 1.0
!
            Do i = 1, Ndim
               Do n = 1, Ndim
                  matnew (i, n) = CMPLX (i, n, kind(0.D0))
                  matold (i, n) = CMPLX (i, n, kind(0.D0))
               End Do
            End Do
!
!
! check against old version from Operator_FFA.f90
!
            Call Op_WrapupFFA (matold, Op, spin, Ndim, N_Type)
!
            Call Op_Wrapup (matnew, Op, spin, Ndim, N_Type)
!
!
            Do i = 1, Ndim
               Do j = 1, Ndim
                  Zre = real (matnew(i, j)-matold(i, j))
                  Zim = aimag (matnew(i, j)-matold(i, j))
                  If (Abs(Zre) > Max(Abs(real(matnew(i, j))), &
                 & Abs(real(matold(i, j))))*1D-14) Then
                     Write (*,*) "opn: ", opn, "N_type", N_Type, "i = ", i, "j = ", j
                     Write (*,*) "ERROR in real part", real (matnew(i, &
                    & j)), real (matold(i, j))
                     Stop 2
                  End If
                  If (Abs(Zim) > Max(Abs(aimag(matnew(i, j))), &
                 & Abs(aimag(matold(i, j))))*1D-14) Then
                     Write (*,*) "ERROR in imag part", aimag (matnew(i, &
                    & j)), aimag (matold(i, j))
                     Stop 3
                  End If
               End Do
            End Do
!
            Deallocate (VH, matnew, matold, Expop, ExpMop)
            call Op_clear(Op, opn)
         End Do
      End Do

End Program Wrapup
!
Subroutine Op_WrapupFFA (Mat, Op, spin, Ndim, N_Type)
!
      Use Operator_mod
      Implicit None
!
      Integer, Intent (In) :: Ndim
      Type (Operator), Intent (In) :: Op
      Complex (Kind=8), Intent (Inout) :: Mat (Ndim, Ndim)
      Real (Kind=8), Intent (In) :: spin
      Integer, Intent (In) :: N_Type
!
    ! Local
      Complex (Kind=8) :: VH (Ndim, Op%n), Z, Z1
      Integer :: n, i, m, lwork, info, j
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: Uold
      Complex (Kind=Kind(0.D0)), Dimension (:), Allocatable :: R, WORK, TAU
!
!
!
!
    !!!!! N_Type ==1
    !    exp(Op%g*spin*Op%E)*(Op%U^{dagger})*Mat*Op%U*exp(-Op%g*spin*Op%E)
    !
    !!!!!
    !!!!! N_Type == 2
    !    Op%U * Mat * (Op%U^{dagger})
    !!!!!
    
    lwork = 2* Op%N
    Allocate (Uold(Op%n, Op%n), R(Op%n-1), TAU(Op%n), WORK(LWORK))
    Uold = Op%U
    if (Op%N > 2) then
        TAU = Op%U(:, Op%N)
        Do i = 1, Op%N - 1
            R(i) = Op%U(i,i)
        ENDDO
        CALL ZUNGQR(Op%N, Op%N, Op%N, Uold, Op%N, TAU, WORK, LWORK, INFO)
        DO i = 1, Op%N
            Do j = 1, Op%N-1
                Uold(i, j) = Uold(i,j) * R(j)
            ENDDO
        ENDDO
    endif
    
      If (N_Type == 1) Then
         VH = 0.D0
         Do n = 1, Op%n
            Z = CMPLX (1.d0, 0.d0, kind(0.D0))
            If (n <= Op%N_non_Zero) Z = Exp (-Op%g*CMPLX(Op%E(n)*spin, &
           & 0.d0, kind(0.D0)))
            Do m = 1, Op%n
               Z1 = Uold (m, n) * Z
               Do i = 1, Ndim
                  VH (i, n) = VH (i, n) + Mat (i, Op%P(m)) * Z1
               End Do
            End Do
         End Do
         Do n = 1, Op%n
            Do i = 1, Ndim
               Mat (i, Op%P(n)) = VH (i, n)
            End Do
         End Do
!
         VH = 0.D0
         Do n = 1, Op%n
            Z = CMPLX (1.d0, 0.d0, kind(0.D0))
            If (n <= Op%N_non_Zero) Z = Exp (Op%g*CMPLX(Op%E(n)*spin, &
           & 0.d0, kind(0.D0)))
            Do m = 1, Op%n
               Z1 = Z * conjg (Uold(m, n))
               Do i = 1, Ndim
                  VH (i, n) = VH (i, n) + Z1 * Mat (Op%P(m), i)
               End Do
            End Do
         End Do
         Do n = 1, Op%n
            Do i = 1, Ndim
               Mat (Op%P(n), i) = VH (i, n)
            End Do
         End Do
      Else If (N_Type == 2) Then
         VH = 0.D0
         Do n = 1, Op%n
            Do m = 1, Op%n
               Z1 = conjg (Uold(n, m))
               Do i = 1, Ndim
                  VH (i, n) = VH (i, n) + Mat (i, Op%P(m)) * Z1
               End Do
            End Do
         End Do
         Do n = 1, Op%n
            Do i = 1, Ndim
               Mat (i, Op%P(n)) = VH (i, n)
            End Do
         End Do
!
         VH = 0.D0
         Do n = 1, Op%n
            Do m = 1, Op%n
               Z1 = Uold (n, m)
               Do i = 1, Ndim
                  VH (i, n) = VH (i, n) + Z1 * Mat (Op%P(m), i)
               End Do
            End Do
         End Do
         Do n = 1, Op%n
            Do i = 1, Ndim
               Mat (Op%P(n), i) = VH (i, n)
            End Do
         End Do
      End If
    Deallocate (Uold, R, TAU, WORK)
End Subroutine Op_WrapupFFA
