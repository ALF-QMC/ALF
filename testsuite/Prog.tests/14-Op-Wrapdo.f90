! compile with
!gfortran  -Wall -std=f2003 -I ../../../Prog_8/  -I ../../../Libraries/Modules/ -L ../../../Libraries/Modules/ main.f90 ../../../Prog_8/Operator.o ../../../Prog_8/UDV_WRAP.o ../../../Libraries/Modules/modules_90.a -llapack -lblas ../../../Libraries/MyNag/libnag.a

Program OPMULTTEST
!
      Use Operator_mod
      Implicit None
!
      Interface
         Subroutine Op_WrapdoFFA (Mat, Op, spin, Ndim, N_Type)
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
      Ndim = 30
!
      Do opn = 1, 4
         Do N_Type = 1, 2
         write (*,*) "opn = ", opn, "N_Type = ", N_type
            Allocate (VH(opn, Ndim), matold(Ndim, Ndim), matnew(Ndim, &
           & Ndim), Expop(opn), ExpMop(opn))
            Call Op_seths ()
            Call Op_make (Op, opn)
!
            Do i = 1, Op%n
               Op%P (i) = i
               Do n = 1, Op%n
                  Op%O (i, n) = CMPLX (i+n, n-i, kind(0.D0))
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
            Call Op_Wrapdo (matnew, Op, spin, Ndim, N_Type)
!
! check against old version from Operator_FFA.f90
            Call Op_WrapdoFFA (matold, Op, spin, Ndim, N_Type)
!
    write (*, *) "opn = ", op%N
    DO I = 1, Ndim
        write (*, *) (matold(I, :))
    ENDDO
write (*,*) "================================"
    DO I = 1, Ndim
        write (*, *) (matnew(I, :))
    ENDDO
            Do i = 1, Ndim
               Do j = 1, Ndim
                  Zre = DBLE (matnew(i, j)-matold(i, j))
                  Zim = aimag (matnew(i, j)-matold(i, j))
                  If(Abs(Zre) > 1D-13) then
                  If (Abs(Zre) > Max(Abs(DBLE(matnew(i, j))), &
                 & Abs(DBLE(matold(i, j))))*1D-13) Then
                     Write (*,*) "opn: ", opn, "N_type", N_Type, "i = ", i, "j = ", j
                     Write (*,*) "ERROR in real part", DBLE (matnew(i, &
                    & j)), DBLE (matold(i, j)), DBLE(Zre)
                     Stop 2
                  End If
                  endif
                  If(Abs(Zim) > 1D-13) then
                  If (Abs(Zim) > Max(Abs(aimag(matnew(i, j))), &
                 & Abs(aimag(matold(i, j))))*1D-13) Then
                     Write (*,*) "ERROR in imag part", aimag (matnew(i, &
                    & j)), aimag (matold(i, j))
                     Stop 3
                  End If
                  Endif
               End Do
            End Do
!
            Deallocate (VH, matnew, matold, Expop, ExpMop)
         End Do
      End Do
      write (*,*) "success"
End Program OPMULTTEST
!
Subroutine Op_WrapdoFFA (Mat, Op, spin, Ndim, N_Type)
!
      Use Operator_mod
      Implicit None
!
      Integer, Intent (In) :: Ndim
      Type (Operator), Intent (In) :: Op
      Complex (Kind=8), Intent (Inout) :: Mat (Ndim, Ndim)
      Real (Kind=8), Intent (In) :: spin
      Integer, Intent (In) :: N_Type
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: Uold
      Complex (Kind=Kind(0.D0)), Dimension (:), Allocatable :: R, WORK, TAU
!
    ! Local
      Complex (Kind=8) :: VH (Ndim, Op%n), Z, Z1
      Integer :: n, i, j, m, lwork, info
!
    !!!!! N_Type == 1
    !    Op%U*exp(-Op%g*spin*Op%E)*Mat*exp(Op%g*spin*Op%E)*(Op%U^{dagger})
    !
    !!!!!
    !!!!! N_Type == 2
    !    (Op%U^{dagger}) * Mat * Op%U
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
         VH = CMPLX (0.d0, 0.d0, kind(0.D0))
         Do m = 1, Op%n
            Z = CMPLX (1.d0, 0.d0, kind(0.D0))
            If (m <= Op%N_non_Zero) Z = Exp (Op%g*CMPLX(Op%E(m)*spin, &
           & 0.d0, kind(0.D0)))
            Do n = 1, Op%n
               Z1 = Z * conjg (Uold(n, m))
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
         VH = CMPLX (0.d0, 0.d0, kind(0.D0))
         Do m = 1, Op%n
            Z = CMPLX (1.d0, 0.d0, kind(0.D0))
            If (m <= Op%N_non_Zero) Z = Exp (-Op%g*CMPLX(Op%E(m)*spin, &
           & 0.d0, kind(0.D0)))
            Do n = 1, Op%n
               Z1 = Z * Uold(n, m)
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
         VH = CMPLX (0.d0, 0.d0, kind(0.D0))
         Do n = 1, Op%n
            Do m = 1, Op%n
               Z1 = Uold (m, n)
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
         VH = CMPLX (0.d0, 0.d0, kind(0.D0))
         Do n = 1, Op%n
            Do m = 1, Op%n
               Z1 = conjg (Uold(m, n))
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
!
End Subroutine Op_WrapdoFFA
!
