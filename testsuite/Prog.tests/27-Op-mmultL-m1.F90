! compile with
!gfortran -std=f2003 -I ../../Prog/ -I ../../Libraries/Modules/ -L ../../Libraries/Modules/  27-Op-mmultL-m1.F90 ../../Prog/Operator_mod.o ../../Prog/Fields_mod.o ../../Libraries/Modules/modules_90.a -llapack -lblas
!
Program Mmult_m1
!
      Use Operator_mod
      Use Fields_mod
        
      Implicit None
!
   
!
      Complex (Kind=Kind(0.D0)) :: Zre, Zim
      Complex (Kind=Kind(0.D0)) :: spin
      Complex (Kind=Kind(0.D0)), Dimension (:, :), Allocatable :: VH, &
     & matnew, matold
      Integer :: i, n, j, Ndim, opn, nt, n_cop
      Type (Operator) :: Op
      Type  (Fields)  :: nsigma_single
      Character       :: Cop
          
!
! setup some test data
      Ndim = 5
      !
      do  n_cop  =  1,3
         Select case(n_cop)
         case(1)
            cop  ='N'
         case(2)
            cop  ='C'
         case(3)
            cop = 'T'
         end Select
         do opn = 1,5
            Call nsigma_single%make(1,1)
            Do nt = 1,4
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
               Call Op_mmultL(Matnew,Op,nsigma_single%f(1,1),cop, 1, 1)
               Call Op_mmultL(Matnew,Op,nsigma_single%f(1,1),cop, 1,-1)
               
               Do i = 1, Ndim
                  Do j = 1, Ndim
                     Zre = real (matnew(i, j)-matold(i, j))
                     Zim = aimag (matnew(i, j)-matold(i, j))
                     If (Abs(Zre) > Max(Abs(real(matnew(i, j))), &
                          & Abs(real(matold(i, j))))*5D-14 .and. abs(Zre) > 1D-15) Then
                        Write (*,*) "opn: ", opn, " OP\%type", Nt, " Cop ", Cop
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
            Call nsigma_single%clear() 
         enddo
      enddo
   write (*,*) "SUCCESS"
 end Program MMult_m1

