Module Operator_mod

  !!!!!! This version of  Operator.f90 contains optimization carried out by Johannes Hofmann
  !!!!!! The original version of this module can be found in Operator_FFA.f90 
  !!!!!! Both versions  must give the same results

  Use MyMats

  Implicit none

  Real (Kind=8) :: Phi(-2:2,2),  Gaml(-2:2,2)
  Integer ::  NFLIPL(-2:2,3)
  

  ! What information should the operator contain
  Type Operator
     Integer          :: N, N_non_zero
     complex (kind=8), pointer :: O(:,:), U (:,:)
     Real    (kind=8), pointer :: E(:)
     Integer, pointer :: P(:)
     complex (kind=8) :: g
     complex (kind=8) :: alpha
     Integer          :: Type 
     ! P is an N X Ndim matrix such that  P.T*O*P*  =  A  
     ! P has only one non-zero entry per column which is specified by P
     ! All in all.   g * Phi(s,type) * ( c^{dagger} A c  + alpha )
     ! The variable Type allows you to define the type of HS. 
     ! The first N_non_zero elemets of  diagonal matrix E are non-zero. The rest vanish.
  end type Operator

  
Contains

  Subroutine Op_SetHS
    Implicit none
    Integer ::  n 
    Phi = 0.d0
    do n = -2,2
       Phi(n,1) = real(n,kind=8)
    enddo
    Phi(-2,2) = - SQRT(2.D0 * ( 3.D0 + SQRT(6.D0) ) )
    Phi(-1,2) = - SQRT(2.D0 * ( 3.D0 - SQRT(6.D0) ) )
    Phi( 1,2) =   SQRT(2.D0 * ( 3.D0 - SQRT(6.D0) ) )
    Phi( 2,2) =   SQRT(2.D0 * ( 3.D0 + SQRT(6.D0) ) )
    
    Do n = -2,2
       gaml(n,1) = 1.d0
    Enddo
    GAML(-2,2) = 1.D0 - SQRT(6.D0)/3.D0
    GAML( 2,2) = 1.D0 - SQRT(6.D0)/3.D0
    GAML(-1,2) = 1.D0 + SQRT(6.D0)/3.D0
    GAML( 1,2) = 1.D0 + SQRT(6.D0)/3.D0
    
    NFLIPL(-2,1) = -1
    NFLIPL(-2,2) =  1
    NFLIPL(-2,3) =  2
    
    NFLIPL(-1,1) =  1
    NFLIPL(-1,2) =  2
    NFLIPL(-1,3) = -2
    
    NFLIPL( 1,1) =  2
    NFLIPL( 1,2) = -2
    NFLIPL( 1,3) = -1
    
    NFLIPL( 2,1) = -2
    NFLIPL( 2,2) = -1
    NFLIPL( 2,3) =  1 
    
  end Subroutine Op_SetHS
  
  Subroutine  Op_phase(Phase,OP_V,Nsigma,N_SUN) ! This also goes in Operator  (Input is nsigma, Op_V).
    Implicit none

    Complex  (Kind=8), Intent(Inout) :: Phase
    Integer,           Intent(IN)    :: N_SUN
    Integer,           dimension(:,:), Intent(In) :: Nsigma
    Type (Operator),   dimension(:,:), Intent(In) :: Op_V
    
    Integer :: n, nf, nt
    
    do nf = 1,Size(Op_V,2)
       do n = 1,size(Op_V,1)
          do nt = 1,size(nsigma,2)
             Phase = Phase*exp(cmplx(0.d0, Aimag( Op_V(n,nf)%g * Op_V(n,nf)%alpha ) * Phi(nsigma(n,nt),Op_V(n,nf)%type) ) )
          enddo
       enddo
    enddo
    Phase = Phase**dble(N_SUN)
    
  end Subroutine Op_phase
  

  subroutine Op_make(Op,N)
    Implicit none
    Type (Operator), intent(INOUT) :: Op
    Integer, Intent(IN) :: N
    Allocate (Op%O(N,N), Op%U(N,N), Op%E(N), Op%P(N))
    Op%O = cmplx(0.d0,0.d0)
    Op%U = cmplx(0.d0,0.d0)
    Op%E = 0.d0
    Op%P = 0
    Op%N = N
    Op%N_non_zero = N
    Op%g     = cmplx(0.d0,0.d0)
    Op%alpha = cmplx(0.d0,0.d0)
  end subroutine Op_make

  subroutine Op_clear(Op,N)
    Implicit none
    Type (Operator), intent(INOUT) :: Op
    Integer, Intent(IN) :: N
    Deallocate (Op%O, Op%U, Op%E, Op%P)
  end subroutine Op_clear 

!==========================================================================
  subroutine Op_set(Op)
    Implicit none
    Type (Operator), intent(INOUT) :: Op

    Complex (Kind=8), allocatable :: U(:,:)
    Real    (Kind=8), allocatable :: E(:)
    Real    (Kind=8) :: Zero = 1.E-9
    Integer :: N, I,J,np,nz

    If (Op%N > 1) then
       !Write(6,*) 'Calling diag', Op%O(1,2), Size(Op%O,1), Size(Op%U,1), Size(Op%E,1)
       N = Op%N
       Allocate (U(N,N), E(N))
       Call Diag(Op%O,U, E)  
       Np = 0
       Nz = 0
       do I = 1,N
          if ( abs(E(I)) > Zero ) then
             np = np + 1
             do j = 1, N
                Op%U(j,np) = U(j,i)
             enddo
             Op%E(np)   = E(I)
          else
             do j = 1, N
                Op%U(j,N-nz) = U(j,i)
             enddo
             Op%E(N-nz)   = E(I)
             nz = nz + 1
          endif
       enddo
       Op%N_non_zero = np
       !Write(6,*) "Op_set", np,N
       deallocate (U, E)
       ! Op%U,Op%E)
       !Write(6,*) 'Calling diag 1'
    else
       Op%E(1)   = Op%O(1,1)
       Op%U(1,1) = cmplx(1.d0,0.d0)
       Op%N_non_zero = 1
    endif
!==========================================================================
  end subroutine Op_set


  subroutine Op_exp(g,Op,Mat)
    Implicit none 
    Type (Operator), Intent(IN)  :: Op
    Complex (Kind=8), Dimension(:,:), INTENT(OUT) :: Mat
    Complex (Kind=8), INTENT(IN) :: g
    Complex (Kind=8) :: Z, Z1

    Integer :: n, i,j
      
    Mat = cmplx(0.d0,0.d0)
    Do n = 1,Op%N
       Z = exp(g*cmplx(Op%E(n),0.d0))
       do J = 1,Op%N
          Z1 = Z*conjg(Op%U(J,n))
          Do I = 1,Op%N
             Mat(I,J) = Mat(I,J) + Op%U(I,n)*Z1
          enddo
       enddo
    enddo
  end subroutine Op_exp

  subroutine Op_mmultL(Mat,Op,spin,Ndim)
    Implicit none 
    Integer :: Ndim
    Type (Operator) , INTENT(IN )   :: Op
    Complex (Kind=8), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind=8), INTENT(IN )   :: spin


    ! Local 
    Complex (Kind=8) :: VH(Op%N,Ndim), Z, Z1, tmp
    Integer :: n, i, m, m1
    

    ! In  Mat
    ! Out Mat = Mat*exp(spin*Op)

    Do n = 1,Op%N
       Do I = 1,Ndim
          VH(n,I) = Mat(I,Op%P(n)) 
       Enddo
    Enddo
      !$OMP PARALLEL DO PRIVATE(tmp)
    do n = 1,Op%N
       Z = exp(Op%g*cmplx(Op%E(n)*spin,0.d0))
       Do I = 1,Ndim
          tmp = cmplx(0.d0,0.d0)
          DO m = 1,Op%N
             tmp  = tmp + VH(m,I) * Op%U(m,n)
          Enddo
	  Mat(I,Op%P(n))  = tmp * Z
       enddo
    Enddo


    Do n = 1,Op%N
       Do I = 1,Ndim
          VH(n,I) = Mat(I,Op%P(n))
       Enddo
    Enddo
      !$OMP PARALLEL DO PRIVATE(tmp)
    do n = 1,Op%N
       Do I = 1,Ndim
          TMP = cmplx(0.d0,0.d0)
          DO m = 1,Op%N
             tmp  = tmp + VH(m,I) * conjg(Op%U(n,m))
          Enddo
	  Mat(I,Op%P(n))  = tmp
       enddo
    Enddo

    
    
  end subroutine Op_mmultL

  subroutine Op_mmultR(Mat,Op,spin,Ndim)
    Implicit none 
    Integer :: Ndim
    Type (Operator) , INTENT(IN )   :: Op
    Complex (Kind=8), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind=8), INTENT(IN )   :: spin


    ! Local 
    Complex (Kind=8) :: VH(Op%N,Ndim), Z, Z1, tmp
    Integer :: n, i, m, m1
    
    ! In  Mat
    ! Out Mat = exp(spin*Op)*Mat
    Do n = 1,Op%N
       Do I = 1,Ndim
          VH(n,I) = Mat(Op%P(n),I)
       Enddo
    Enddo
    !$OMP PARALLEL DO PRIVATE(tmp, Z1)
    do n = 1,Op%N
       Z1 = exp(Op%g*cmplx(Op%E(n)*spin,0.d0))
       Do I = 1,Ndim
          tmp = cmplx(0.d0,0.d0)
          DO m = 1,Op%N
             tmp  = tmp + conjg(Op%U(m,n)) * VH(m,I) 
          Enddo
	  Mat(Op%P(n),I)  = tmp * Z1
       enddo
    Enddo
    !$OMP END PARALLEL DO
    Do n = 1,Op%N
       Do I = 1,Ndim
          VH(n,I) = Mat(Op%P(n),I) 
       Enddo
    Enddo
    !$OMP PARALLEL DO PRIVATE(tmp)
    do n = 1,Op%N
       Do I = 1,Ndim
          tmp = cmplx(0.d0,0.d0)
          DO m = 1,Op%N
             tmp = tmp + Op%U(n,m)* VH(m,I)
          Enddo
	  Mat(Op%P(n),I)  = tmp
       enddo
    Enddo
    !$OMP END PARALLEL DO


  end subroutine Op_mmultR

  Subroutine Op_Wrapup(Mat,Op,spin,Ndim,N_Type)

    Implicit none 

    Integer :: Ndim
    Type (Operator) , INTENT(IN )   :: Op
    Complex (Kind=8), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind=8), INTENT(IN )   :: spin
    Integer, INTENT(IN) :: N_Type

    ! Local 
    Complex (Kind=8) :: Z, Z1, tmp, ExpOp(Op%N), ExpMop(Op%N), VH(Op%N,Ndim), ExpHere !, zdotu, zdotc
    Integer :: n, i, m, m1 !, nop
    
    
!     nop=size(Op%U,1)
    
    !!!!! N_Type ==1
    !    exp(Op%g*spin*Op%E)*(Op%U^{dagger})*Mat*Op%U*exp(-Op%g*spin*Op%E)
    !    
    !!!!!
    !!!!! N_Type == 2
    !    Op%U * Mat * (Op%U^{dagger})
    !!!!!
    If (N_type == 1) then
       do n= 1,Op%N
          ExpOp(n) = cmplx(1.d0,0.d0)
          If ( n <= OP%N_non_Zero) ExpOp(n) = exp(Op%g*cmplx(Op%E(n)*spin,0.d0))
          ExpMOp(n) = cmplx(1.d0,0.d0)
          If ( n <= OP%N_non_Zero) ExpMOp(n) = exp(-Op%g*cmplx(Op%E(n)*spin,0.d0))
       enddo
       
       Do n = 1,Op%N
          Do I = 1,Ndim
             VH(n,I) = Mat(I,Op%P(n))
          Enddo
       Enddo
       
      !$OMP PARALLEL DO PRIVATE(tmp)
       do n = 1,Op%N
	  ExpHere=ExpMOp(n)
	  DO I = 1,Ndim
!              Mat(I,Op%P(n))  =  ExpMOp(n) * zdotu(Op%N,Op%U(1,n),1,VH(1,I),1) 
	     tmp=cmplx(0.d0,0.d0)
	     Do m = 1,Op%N
		tmp = tmp + VH(m,I) * Op%U(m,n)
             Enddo
             Mat(I,Op%P(n))  =  ExpHere * tmp
          enddo
       Enddo
      !$OMP END PARALLEL DO

       Do n = 1,Op%N
          Do I = 1,Ndim
             VH(n,I) =    Mat(Op%P(n),I)
          Enddo
       Enddo
      !$OMP PARALLEL DO PRIVATE(tmp)
       do n = 1,Op%N
          ExpHere=ExpOp(n)
          DO I = 1,Ndim
!               Mat(Op%P(n),I)  = zdotc(Op%N,Op%U(1,n),1,VH(1,I),1) * ExpOp(n)
             tmp=cmplx(0.d0,0.d0)
	     Do m = 1,Op%N
                tmp = tmp + conjg(Op%U(m,n))* VH(m,I)
             Enddo
              Mat(Op%P(n),I)  = ExpHere * tmp
          enddo
       enddo
      !$OMP END PARALLEL DO
    elseif (N_Type == 2) then
       Do n = 1,Op%N
          Do I = 1,Ndim
             VH(n,I) = Mat(I,Op%P(n))
          Enddo
       Enddo
      !$OMP PARALLEL DO PRIVATE(tmp)
       do n = 1,Op%N
          DO I = 1,Ndim
!              Mat(I,Op%P(n))  = zdotc(Op%N,Op%U(n,1),nop,VH(1,I),1)
             tmp=cmplx(0.d0,0.d0)
	     Do m = 1,Op%N
                tmp = tmp +  VH(m,I) * conjg(Op%U(n,m)) 
             Enddo
             Mat(I,Op%P(n))  =  tmp
          enddo
       Enddo
      !$OMP END PARALLEL DO
       
       Do n = 1,Op%N
          Do I = 1,Ndim
             VH(n,I) = Mat(Op%P(n),I)
          Enddo
       Enddo
      !$OMP PARALLEL DO PRIVATE(tmp)
       do n = 1,Op%N
          DO I = 1,Ndim
! 	     Mat(Op%P(n),I)  = zdotu(Op%N,Op%U(n,1),nop,VH(1,I),1)
             tmp=cmplx(0.d0,0.d0)
	     Do m = 1,Op%N
	        tmp = tmp + Op%U(n,m) * VH(m,I)
             Enddo
	     Mat(Op%P(n),I)  = tmp
          enddo
       enddo
      !$OMP END PARALLEL DO
    endif
  end Subroutine Op_Wrapup

  Subroutine Op_Wrapdo(Mat,Op,spin,Ndim,N_Type)

    Implicit none 

    Integer :: Ndim
    Type (Operator) , INTENT(IN )   :: Op
    Complex (Kind=8), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind=8), INTENT(IN )   :: spin
    Integer, INTENT(IN) :: N_Type

    ! Local 
    Complex (Kind=8) :: VH(Op%N,Ndim), Z, Z1, tmp, ExpOp(Op%N), ExpMop(Op%N), ExpHere !, zdotu, zdotc
!     Complex (Kind=8) :: alpha, beta, tmp2(Op%N,Ndim)
    Integer :: n, i, m, m1 !, nop
    
!     nop=size(Op%U,1)
    
    !!!!! N_Type == 1
    !    Op%U*exp(-Op%g*spin*Op%E)*Mat*exp(Op%g*spin*Op%E)*(Op%U^{dagger})
    !    
    !!!!!
    !!!!! N_Type == 2
    !    (Op%U^{dagger}) * Mat * Op%U
    !!!!!
    If (N_type == 1) then
       do n= 1,Op%N
          ExpOp(n) = cmplx(1.d0,0.d0)
          If ( n <= OP%N_non_Zero) ExpOp(n) = exp(Op%g*cmplx(Op%E(n)*spin,0.d0))
          ExpMOp(n) = cmplx(1.d0,0.d0)
          If ( n <= OP%N_non_Zero) ExpMOp(n) = exp(-Op%g*cmplx(Op%E(n)*spin,0.d0))
       enddo
       
       Do n = 1,Op%N
	  expHere=ExpOp(n)
          Do I = 1,Ndim
             VH(n,I) = Mat(I,Op%P(n)) * ExpHere
          Enddo
       Enddo
       
       ! ZGEMM might be the better multiplication, but the additional copy precess seem to be to expensive
!        alpha = cmplx (1.0d0,0.0d0)
!        beta = cmplx (0.0d0,0.0d0)
!        CALL ZGEMM('T','C',Ndim,Op%N,Op%N,alpha,VH,Op%N,Op%U,nop,beta,tmp,Ndim)
       !$OMP PARALLEL DO PRIVATE(tmp)
       do n = 1,Op%N
          DO I = 1,Ndim
!              Mat(I,Op%P(n))  = zdotc(Op%N,Op%U(n,1),nop,VH(1,I),1)
             tmp=cmplx(0.d0,0.d0)
	     Do m = 1,Op%N
                tmp = tmp + VH(m,I) * conjg(Op%U(n,m))
             Enddo
             Mat(I,Op%P(n))  = tmp
          enddo
       Enddo
       !$OMP END PARALLEL DO

       Do n = 1,Op%N
	  ExpHere=ExpMOp(n)
          Do I = 1,Ndim
             VH(n,I) = ExpHere * Mat(Op%P(n),I)
          Enddo
       Enddo
       
       ! ZGEMM might be the better multiplication, but the additional copy precess seem to be to expensive
!        alpha = cmplx (1.0d0,0.0d0)
!        beta = cmplx (0.0d0,0.0d0)
!        CALL ZGEMM('T','T',Ndim,Op%N,Op%N,alpha,VH,Op%N,Op%U,nop,beta,tmp,Ndim)
      !$OMP PARALLEL DO PRIVATE(tmp)
       do n = 1,Op%N
	 DO I = 1,Ndim
! 	   Mat(Op%P(n),I)  = zdotu(Op%N,Op%U(n,1),nop,VH(1,I),1)
           tmp=cmplx(0.d0,0.d0)
	   Do m = 1,Op%N
	      tmp = tmp + Op%U(n,m)*VH(m,I)
	   Enddo
	   Mat(Op%P(n),I)  = tmp 
         enddo
       enddo
      !$OMP END PARALLEL DO
    elseif (N_Type == 2) then
       Do n = 1,Op%N
          Do I = 1,Ndim
             VH(n,I) = Mat(I,Op%P(n)) 
          Enddo
       Enddo
       
       ! ZGEMM might be the better multiplication, but the additional copy precess seem to be to expensive
!        alpha = cmplx (1.0d0,0.0d0)
!        beta = cmplx (0.0d0,0.0d0)
!        CALL ZGEMM('T','N',Ndim,Op%N,Op%N,alpha,VH,Op%N,Op%U,nop,beta,tmp,Ndim)
       !$OMP PARALLEL DO PRIVATE(tmp)
       do n = 1,Op%N
          DO I = 1,Ndim
!              Mat(I,Op%P(n))  =  zdotu(Op%N,Op%U(1,n),1,VH(1,I),1)
             tmp=cmplx(0.d0,0.d0)
             Do m = 1,Op%N
	        tmp = tmp + VH(m,I) * Op%U(m,n)
             Enddo
             Mat(I,Op%P(n))  =  tmp
          enddo
       Enddo
      !$OMP END PARALLEL DO
       
       Do n = 1,Op%N
          Do I = 1,Ndim
             VH(n,I) = Mat(Op%P(n),I) 
          Enddo
       Enddo
       
       ! ZGEMM might be the better multiplication, but the additional copy precess seem to be to expensive
!        alpha = cmplx (1.0d0,0.0d0)
!        beta = cmplx (0.0d0,0.0d0)
!        CALL ZGEMM('C','N',Op%N,Ndim,Op%N,alpha,Op%U,nop,VH,Op%N,beta,tmp2,Op%N)
      !$OMP PARALLEL DO PRIVATE(tmp)
       do n = 1,Op%N
          DO I = 1,Ndim
!              Mat(Op%P(n),I)  = zdotc(Op%N,Op%U(1,n),1,VH(1,I),1)
             tmp=cmplx(0.d0,0.d0)
	     Do m = 1,Op%N
                tmp = tmp + conjg(Op%U(m,n))* VH(m,I) 
             Enddo
             Mat(Op%P(n),I)  = tmp
          enddo
       enddo
      !$OMP END PARALLEL DO
    endif
    
  end Subroutine Op_Wrapdo


end Module Operator_mod
