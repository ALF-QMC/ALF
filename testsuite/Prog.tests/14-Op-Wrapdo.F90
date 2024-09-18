! compile with
!gfortran -std=f2003 -I ../../Prog/ -I ../../Libraries/Modules/ -L ../../Libraries/Modules/  14-Op-Wrapdo.F90 ../../Prog/Operator_mod.o ../../Prog/Fields_mod.o ../../Libraries/Modules/modules_90.a -llapack -lblas


!
!
Program OPWRAPDO
!
      Use Operator_mod
      Use Fields_mod
      Implicit None
!
      Interface
         Subroutine Op_WrapdoFFA (Mat, Op, spin, Ndim, N_Type)
            Use Operator_mod
            Type (Operator), Intent (In) :: Op
            Complex (Kind=kind(0.D0)), allocatable, Intent (Inout) :: Mat (:,:)
            Complex (Kind=kind(0.D0)), Intent (In) :: spin
            Integer, Intent (In) :: N_Type, Ndim
         End Subroutine
      End Interface
!
      Complex (Kind=Kind(0.D0)) :: Zre, Zim
      Complex (Kind=Kind(0.D0)) :: spin
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: VH, &
     & matnew, matold
      Complex (Kind=Kind(0.D0)), Dimension (:), Allocatable :: Expop, &
     & ExpMop
      Integer :: i, n, j, Ndim, N_Type, opn
      Type (Operator) :: Op
      Type (Fields) :: nsigma_single
    
      Call nsigma_single%make(1,1)
      
!
! setup some test data
      Ndim = 30
!
      Do opn = 1, 4
         Do N_Type = 1, 2
            Allocate (VH(opn, Ndim), matold(Ndim, Ndim), matnew(Ndim, &
           & Ndim), Expop(opn), ExpMop(opn))
            Call Op_make (Op, opn)
!
            Do i = 1, Op%n
               Op%P (i) = i
               Do n = 1, Op%n
                  Op%O (i, n) = CMPLX (0.25d0*(n+i), 0.25d0*(n-i), kind(0.D0))
               End Do
            End Do
!
            Op%type=1
            Op%g = 2.D0
            Op%alpha = 0.D0
            Call Op_set (Op)
!
            spin = -1.d0
            nsigma_single%f(1,1) = spin
            nsigma_single%t(1) = Op%type 
!
            Do i = 1, Ndim
               Do n = 1, Ndim
                  matnew (i, n) = CMPLX (i, n, kind(0.D0))
                  matold (i, n) = CMPLX (i, n, kind(0.D0))
               End Do
            End Do
!
            Call Op_Wrapdo (matnew, Op, spin, Ndim, N_Type,1)
!
! check against old version from Operator_FFA.F90
            spin = nsigma_single%Phi(1,1) 
            Call Op_WrapdoFFA (matold, Op, spin, Ndim, N_Type)
!
            Do i = 1, 3
               Do j = 1, 3
                  Zre = real (matnew(i, j)-matold(i, j))
                  Zim = aimag (matnew(i, j)-matold(i, j))
                  if(Abs(Zre)+Abs(Zim) > 1.D-14) THEN
                  If (Abs(Zre) > Max(Abs(DBLE(matnew(i, j))), &
                 & Abs(real(matold(i, j))))*1D-12 .and. abs(Zre) > 1D-15) Then
                     Write (*,*) "opn: ", opn, "N_type", N_Type
                     Write (*,*) "ERROR in real part", DBLE(matnew(i, &
                    & j)), DBLE (matold(i, j))
                     Stop 2
                  End If
                  If (Abs(Zim) > Max(Abs(aimag(matnew(i, j))), &
                 & Abs(aimag(matold(i, j))))*1D-12 .and. abs(Zim) > 1D-15) Then
                     Write (*,*) "ERROR in imag part", aimag (matnew(i, &
                    & j)), aimag (matold(i, j))
                     Stop 3
                  End If
                  endif
               End Do
            End Do
!
            call Op_clear(Op,opn)
            Deallocate (VH, matnew, matold, Expop, ExpMop)
         End Do
      End Do
      
      Call nsigma_single%clear() 

      write (*,*) "success"
End Program OPWRAPDO
!
Subroutine Op_WrapdoFFA (Mat, Op, spin, Ndim, N_Type)
!
      Use Operator_mod
      Implicit None
!
      Integer, Intent (In) :: Ndim
      Type (Operator), Intent (In) :: Op
      Complex (Kind=kind(0.D0)), allocatable, Intent (Inout) :: Mat (:,:)
      Complex (Kind=kind(0.D0)), Intent (In) :: spin
      Integer, Intent (In) :: N_Type
!
    ! Local
      Complex (Kind=kind(0.D0)) :: VH (Ndim, Op%n), Z, Z1
      Integer :: n, i, m
!
    !!!!! N_Type == 1
    !    Op%U*exp(-Op%g*spin*Op%E)*Mat*exp(Op%g*spin*Op%E)*(Op%U^{dagger})
    !
    !!!!!
    !!!!! N_Type == 2
    !    (Op%U^{dagger}) * Mat * Op%U
    !!!!!
      If (N_Type == 1) Then
         VH = CMPLX (0.d0, 0.d0, kind(0.D0))
         Do m = 1, Op%n
            Z = CMPLX (1.d0, 0.d0, kind(0.D0))
            If (m <= Op%N_non_Zero) Z = Exp (Op%g*CMPLX(Op%E(m)*Real(spin,kind(0.d0)), &
           & 0.d0, kind(0.D0)))
            Do n = 1, Op%n
               Z1 = Z * conjg (Op%U(n, m))
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
            If (m <= Op%N_non_Zero) Z = Exp (-Op%g*CMPLX(Op%E(m)*Real(spin,kind(0.d0)), &
           & 0.d0, kind(0.D0)))
            Do n = 1, Op%n
               Z1 = Z * Op%U(n, m)
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
               Z1 = Op%U (m, n)
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
               Z1 = conjg (Op%U(m, n))
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
!
End Subroutine Op_WrapdoFFA
!
