! compile with
!gfortran -std=f2003 -I ../../Prog/ -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ 11-Op-mmultR.F90  ../../Prog/Operator_mod.o ../../Prog/Fields_mod.o ../../Libraries/Modules/modules_90.a -llapack -lblas


Program Opmulttest

        Use Operator_mod
        Use Fields_mod
        implicit none

        Complex (Kind=Kind(0.D0)) :: Matnew(3,3), matold(3,3), VH(3,3), Z, Z1, Zre, Zim, nspin
        Complex    (KIND = KIND(0.D0)) :: spin
        Integer :: i, n, m, j, ndim , nt, n_type
        Type(Operator) :: Op
        Type (Fields) :: nsigma_single
    
        
        
! setup some test data
        Ndim = 3

        Do n_type = 1,4

           Call Op_make(Op, 3)
           do i = 1, Op%N
              !             Op%E(i) = 2*i-3
              Op%P(i) = i
              do n = 1,Op%N
                 Op%O(i,n) = CMPLX(n+i, n-i, kind(0.D0))
              enddo
           enddo
           
           Op%type=n_type
           Op%N_non_zero = 2
           Op%g = 0.02D0
           call Op_set(Op)
           Call nsigma_single%make(1,1)
           nsigma_single%t(1)   = Op%type
           Select case(n_type)
           case (1) 
              nsigma_single%f(1,1) = cmplx(real(1,kind=kind(0.d0)), 0.d0, kind(0.d0))
           case (2) 
              nsigma_single%f(1,1) = cmplx(real(2,kind=kind(0.d0)), 0.d0, kind(0.d0))
           case (3) 
              nsigma_single%f(1,1) = cmplx(3.14159267d0, 0.d0, kind(0.d0))
           case (4) 
              nsigma_single%f(1,1) = cmplx(-1.d0, 0.5d0, kind(0.d0))
           end Select
           do i = 1,Ndim
              do n = 1,Ndim
                 matnew(i,n) = CMPLX(i,n, kind(0.D0))
                 matold(i,n) = CMPLX(i,n, kind(0.D0))
              enddo
           enddo
           Call Op_mmultR(matnew, Op, nsigma_single%f(1,1), 'n',1)
           
           ! check against old version from Operator_FFA.F90
           
           VH = 0.d0
           do n = 1,Op%N
              Z1 = exp(Op%g*Op%E(n)*nsigma_single%phi(1,1))
              Do m = 1,Op%N
                 Z =  conjg(Op%U(m,n))* Z1 
                 DO I = 1,Ndim
                    VH(I,n)  = VH(I,n) + Z* Matold(Op%P(m),I) 
                 Enddo
              enddo
           Enddo
           Do n = 1,Op%N
              Do I = 1,Ndim
                 Matold(Op%P(n),I) =   VH(I,n) 
              Enddo
           Enddo
           VH = 0.d0
           do n = 1,Op%N
              Do m = 1,Op%N
                 Z =  Op%U(n,m)
                 DO I = 1,Ndim
                    VH(I,n)  = VH(I,n) + Z* Matold(Op%P(m),I) 
                 Enddo
              enddo
           Enddo
           Do n = 1,Op%N
              Do I = 1,Ndim
                 Matold(Op%P(n),I) =   VH(I,n)
              Enddo
           Enddo
           
           do i=1,3
              do j=1,3
                 Zre = real(matnew(i,j)-matold(i,j))
                 Zim = aimag(matnew(i,j)-matold(i,j))
                 if (Abs(Zre) > MAX(ABS(real(matnew(i,j))), ABS(real(matold(i,j))) )*1D-14) then
                    write (*,*) "ERROR in real part", real(matnew(i,j)), real(matold(i,j))
                    STOP 2
                 endif
                 if (Abs(Zim) > MAX(ABS(aimag(matnew(i,j))), ABS(aimag(matold(i,j))) )*1D-14) then
                    write (*,*) "ERROR in imag part", aimag(matnew(i,j)), aimag(matold(i,j))
                    STOP 3
                 endif
              enddo
           enddo

           Call Op_clear(Op, 3)
           Call Op_make(Op, 3)
           !Repeat test for diagonal Operator
           Op%O = CMPLX(0.d0, 0.d0, kind(0.D0))
           do i = 1, Op%N
              Op%O(i,i) = 2*i-3
              Op%P(i) = i
              !             Op%U(i,i) = CMPLX(1.d0, 0.d0, kind(0.D0))
           enddo
           ! the following line is neccessary as we circumvent the Op_set routine
           call Op_set(Op)
           
           do i = 1,Ndim
              do n = 1,Ndim
                 matnew(i,n) = CMPLX(i,n, kind(0.D0))
                 matold(i,n) = CMPLX(i,n, kind(0.D0))
              enddo
           enddo
           
           Call Op_mmultR(matnew, Op, nsigma_single%f(1,1), 'n',1)
           
           ! check against old version from Operator_FFA.F90
           
           VH = 0.d0
           do n = 1,Op%N
              Z1 = exp(Op%g*Op%E(n)*nsigma_single%phi(1,1))
              Do m = 1,Op%N
                 Z =  conjg(Op%U(m,n))* Z1 
                 DO I = 1,Ndim
                    VH(I,n)  = VH(I,n) + Z* Matold(Op%P(m),I) 
                 Enddo
              enddo
           Enddo
           Do n = 1,Op%N
              Do I = 1,Ndim
                 Matold(Op%P(n),I) =   VH(I,n) 
              Enddo
           Enddo
           VH = 0.d0
           do n = 1,Op%N
              Do m = 1,Op%N
                 Z =  Op%U(n,m)
                 DO I = 1,Ndim
                    VH(I,n)  = VH(I,n) + Z* Matold(Op%P(m),I) 
                 Enddo
              enddo
           Enddo
           Do n = 1,Op%N
              Do I = 1,Ndim
                 Matold(Op%P(n),I) =   VH(I,n)
              Enddo
           Enddo
           
           do i=1,3
              do j=1,3
                 Zre = real(matnew(i,j)-matold(i,j))
                 Zim = aimag(matnew(i,j)-matold(i,j))
                 if (Abs(Zre) > MAX(ABS(real(matnew(i,j))), ABS(real(matold(i,j))) )*1D-14) then
                    write (*,*) "ERROR in real part", real(matnew(i,j)), real(matold(i,j))
                    STOP 2
                 endif
                 if (Abs(Zim) > MAX(ABS(aimag(matnew(i,j))), ABS(aimag(matold(i,j))) )*1D-14) then
                    write (*,*) "ERROR in imag part", aimag(matnew(i,j)), aimag(matold(i,j))
                    STOP 3
                 endif
              enddo
           enddo
           call Op_clear(Op, 3)
           Call nsigma_single%clear()
           
        enddo
        write (*,*) "success"
        
      end Program OPMULTTEST
