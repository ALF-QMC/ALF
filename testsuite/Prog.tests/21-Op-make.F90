! compile with
! gfortran -Wall -std=f2003 -I ../../Prog/  -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ 21-Op-make.F90 ../../Prog/Operator.o ../../Libraries/Modules/modules_90.a -llapack -lblas
!
Program TESTOPMAKE
!
      Use Operator_mod
      Implicit None
!
      Interface
        subroutine Op_exp_FFA(g,Op,Mat)
            Use Operator_mod
            Type (Operator), Intent(IN)  :: Op
            Complex (Kind=8), allocatable, INTENT(INOUT) :: Mat(:,:)
            Complex (Kind=8), INTENT(IN) :: g
         End Subroutine
      End Interface
!
      Real (Kind=Kind(0.D0)) :: Zre, Zim
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: matnew, matold
      Integer :: i, n, j,  opn
      Type (Operator) :: Op
!
      Do opn = 1, 16 ! overflow occurs for our particular O
            Allocate (matold(opn, opn), matnew(opn, opn))
            ! setup some test data
            Call Op_make (Op, opn)
            Op%type = 2 ! Op_make resets the type to 0
!

            Do i = 1, Op%n
               Op%P (i) = i
               Do n = 1, Op%n
                  Op%O (i, n) = CMPLX (n+i, i-n, kind(0.D0))
               End Do
            End Do
!
            Op%g = 2.D0
            Op%alpha = 0.D0
            Call Op_set (Op)
! check against old version from Operator_FFA.F90
!
            Call Op_exp_FFA (Op%g, Op, matold)
!
            Call Op_exp (Op%g, Op, matnew)
!
!
            Do i = 1, opn
               Do j = 1, opn
                  Zre = real (matnew(i, j)-matold(i, j))
                  Zim = aimag (matnew(i, j)-matold(i, j))
                  If (Abs(Zre) > Max(Abs(real(matnew(i, j))), &
                 & Abs(real(matold(i, j))))*1D-15) Then
                     Write (*,*) "opn: ", opn
                     Write (*,*) "ERROR in real part", real (matnew(i, &
                    & j)), real (matold(i, j)), Zre
!                     Stop 2
                  End If
                  If (Abs(Zim) > Max(Abs(aimag(matnew(i, j))), &
                 & Abs(aimag(matold(i, j))))*1D-15) Then
                     Write (*,*) "opn: ", opn
                     Write (*,*) "ERROR in imag part", aimag (matnew(i, &
                    & j)), aimag (matold(i, j)), Zim
!                     Stop 3
                  End If
            End Do
         End Do
        Deallocate (matnew, matold)
        call Op_clear(Op, opn)
     End Do
     
      write(*,*) "success"

End Program TESTOPMAKE
!
  subroutine Op_exp_FFA(g,Op,Mat)
    Use Operator_mod
    Implicit none 
    Type (Operator), Intent(IN)  :: Op
    Complex (Kind=8), allocatable, INTENT(INOUT) :: Mat(:,:)
    Complex (Kind=8), INTENT(IN) :: g
    Complex (Kind=8) :: Z, Z1

    Integer :: n, i,j
      
    Mat = cmplx(0.d0,0.d0, kind(0.D0))
    Do n = 1,Op%N
       Z = exp(g*cmplx(Op%E(n),0.d0, kind(0.D0)))
       do J = 1,Op%N
          Z1 = Z*conjg(Op%U(J,n))
          Do I = 1,Op%N
             Mat(I,J) = Mat(I,J) + Op%U(I,n)*Z1
          enddo
       enddo
    enddo
  end subroutine Op_exp_FFA
