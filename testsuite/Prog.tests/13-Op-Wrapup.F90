! compile with
!gfortran -std=f2003 -I ../../Prog/ -I ../../Libraries/Modules/ -L ../../Libraries/Modules/  13-Op-Wrapup.F90 ../../Prog/Operator_mod.o ../../Prog/Fields_mod.o ../../Libraries/Modules/modules_90.a -llapack -lblas
!
Program Wrapup
!
      Use Operator_mod
      Use Fields_mod
        
      Implicit None
!
      Interface
         Subroutine Op_WrapupFFA (Mat, Op, spin, Ndim, N_Type)
            Use Operator_mod
            Type (Operator), Intent (In) :: Op
            Complex (Kind=kind(0.D0)), allocatable, Intent (Inout) :: Mat (:, :)
            Complex (Kind=kind(0.D0)), Intent (In) :: spin
            Integer, Intent (In) :: N_Type, Ndim
         End Subroutine
      End Interface
!
      Complex (Kind=Kind(0.D0)) :: Zre, Zim
      Complex (Kind=Kind(0.D0)) :: spin
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: VH, &
     & matnew, matold
      Integer :: i, n, j, Ndim, N_Type, opn, nt
      Type (Operator) :: Op
      Type  (Fields)  :: nsigma_single
          
!
! setup some test data
      Ndim = 5
      !

      Call nsigma_single%make(1,1)
      Do nt = 1,4
      Do opn = 1, 4
         Do N_Type = 1, 2
            Allocate (VH(opn, Ndim), matold(Ndim, Ndim), matnew(Ndim,  Ndim))
            Call Op_make (Op, opn)
!
            Do i = 1, Op%n
               Op%P (i) = i
               Do n = 1, Op%n
                 Op%O (i, n) = CMPLX (0.d1*dble(n+i), 0.d1*dble(n-i), kind(0.D0))
               End Do
            End Do
!
            Op%type = nt
            Op%g = 2.D0
            Op%alpha = 0.D0
            Call Op_set (Op)
            !
            Select case(nt)
            case (1) 
               nsigma_single%f(1,1) = cmplx(real(1,kind=kind(0.d0)), 0.d0, kind(0.d0))
            case (2) 
               nsigma_single%f(1,1) = cmplx(real(2,kind=kind(0.d0)), 0.d0, kind(0.d0))
            case (3) 
               nsigma_single%f(1,1) = cmplx(3.14159267d0, 0.d0, kind(0.d0))
            case (4) 
               nsigma_single%f(1,1) = cmplx(-1.d0, 0.5d0, kind(0.d0))
            end Select
            nsigma_single%t(1)   = Op%type 
!
            Do i = 1, Ndim
               Do n = 1, Ndim
                  matnew (i, n) = CMPLX (i, n, kind(0.D0))
                  matold (i, n) = CMPLX (i, n, kind(0.D0))
               End Do
            End Do
!
!
! check against old version from Operator_FFA.F90
!
            Call Op_WrapupFFA (matold, Op, nsigma_single%phi(1,1), Ndim, N_Type)
!
            Call Op_Wrapup (matnew, Op, nsigma_single%f(1,1), Ndim, N_Type, 1)
!
!
            Do i = 1, Ndim
               Do j = 1, Ndim
                  Zre = real (matnew(i, j)-matold(i, j))
                  Zim = aimag (matnew(i, j)-matold(i, j))
                  If (Abs(Zre) > Max(Abs(real(matnew(i, j))), &
                 & Abs(real(matold(i, j))))*5D-14 .and. abs(Zre) > 1D-15) Then
                     Write (*,*) "opn: ", opn, "N_type", N_Type
                     Write (*,*) "ERROR in real part", real (matnew(i, &
                    & j)), real (matold(i, j))
                     Stop 2
                  End If
                  If (Abs(Zim) > Max(Abs(aimag(matnew(i, j))), &
                 & Abs(aimag(matold(i, j))))*5D-14 .and. abs(Zim) > 1D-15) Then
                     Write (*,*) "ERROR in imag part", aimag (matnew(i, &
                    & j)), aimag (matold(i, j))
                     Stop 3
                  End If
               End Do
            End Do
!
            Deallocate (VH, matnew, matold)
            call Op_clear(Op, opn)
         End Do
      End Do
   Enddo
   Call nsigma_single%clear() 
      
   write (*,*) "SUCCESS"
End Program Wrapup
!
Subroutine Op_WrapupFFA (Mat, Op, spin, Ndim, N_Type)
!
      Use Operator_mod
      Implicit None
!
      Integer, Intent (In) :: Ndim
      Type (Operator), Intent (In) :: Op
      Complex (Kind=kind(0.D0)), allocatable, Intent (Inout) :: Mat (:, :)
      Complex (Kind=kind(0.D0)), Intent (In) :: spin
      Integer, Intent (In) :: N_Type
!
    ! Local
      Complex (Kind=kind(0.D0)) :: VH (Ndim, Op%n), Z, Z1
      Integer :: n, i, m
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
      If (N_Type == 1) Then
         VH = 0.D0
         Do n = 1, Op%n
            Z = CMPLX (1.d0, 0.d0, kind(0.D0))
            If (n <= Op%N_non_Zero) Z = Exp (-Op%g*Op%E(n)*spin)
            Do m = 1, Op%n
               Z1 = Op%U (m, n) * Z
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
            If (n <= Op%N_non_Zero) Z = Exp (Op%g*Op%E(n)*spin)
            Do m = 1, Op%n
               Z1 = Z * conjg (Op%U(m, n))
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
               Z1 = conjg (Op%U(n, m))
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
               Z1 = Op%U (n, m)
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
End Subroutine Op_WrapupFFA
