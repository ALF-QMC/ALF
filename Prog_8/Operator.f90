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

!--------------------------------------------------------------------
!> @author
!> 
!
!> @brief 
!> calculate the phase of a given set of operators and HS fields.
!
!> @param[inout] Phase
!> @param[in] Op_V An array of Operators
!> @param[in] Nsigma
!> @param[in] N_SUN 
!--------------------------------------------------------------------
  Subroutine  Op_phase(Phase,OP_V,Nsigma,N_SUN) ! This also goes in Operator  (Input is nsigma, Op_V).
    Implicit none

    Complex  (Kind=8), Intent(Inout) :: Phase
    Integer,           Intent(IN)    :: N_SUN
    Integer,           dimension(:,:), Intent(In) :: Nsigma
    Type (Operator),   dimension(:,:), Intent(In) :: Op_V
    Real  (Kind=8)                   :: angle
    
    Integer :: n, nf, nt
    
    do nf = 1,Size(Op_V,2)
       do n = 1,size(Op_V,1)
          do nt = 1,size(nsigma,2)
             angle = Aimag( Op_V(n,nf)%g * Op_V(n,nf)%alpha ) * Phi(nsigma(n,nt),Op_V(n,nf)%type)
             Phase = Phase*CMPLX(cos(angle),sin(angle), Kind(0.D0))
          enddo
       enddo
    enddo
    Phase = Phase**N_SUN
    
  end Subroutine Op_phase
  

  subroutine Op_make(Op,N)
    Implicit none
    Type (Operator), intent(INOUT) :: Op
    Integer, Intent(IN) :: N
    Allocate (Op%O(N,N), Op%U(N,N), Op%E(N), Op%P(N))
    Op%O = cmplx(0.d0, 0.d0, kind(0.D0))
    Op%U = cmplx(0.d0, 0.d0, kind(0.D0))
    Op%E = 0.d0
    Op%P = 0
    Op%N = N
    Op%N_non_zero = N
    Op%g     = cmplx(0.d0,0.d0, kind(0.D0))
    Op%alpha = cmplx(0.d0,0.d0, kind(0.D0))
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
    Integer :: N, I, np,nz

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
             Op%U(:, np) = U(:, i)
             Op%E(np)   = E(I)
          else
             Op%U(:, N-nz) = U(:, i)
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
       Op%E(1)   = REAL(Op%O(1,1), kind(0.D0))
       Op%U(1,1) = cmplx(1.d0, 0.d0 , kind(0.D0))
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

    Integer :: n, j, iters
    
    iters = Op%N
    Mat = cmplx(0.d0, 0.d0, kind(0.D0))
    Do n = 1, iters
       Z = exp(g*Op%E(n))
       do J = 1, iters
          Z1 = Z*conjg(Op%U(J,n))
          Mat(1:iters, J) = Mat(1:iters, J) + Z1 * Op%U(1:iters, n)
       enddo
    enddo
  end subroutine Op_exp

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function copies select rows to the destination matrix V from
!> the source matrix Mat. The decision on which rows to copy is determined
!> by the vector P.
!
!> @param[inout] V storage for the result matrix
!> @param[in] Mat Where to read those entries
!> @param[in] P A vector with which rows to copy
!> @param[in] opn The length of the vector P
!> @param[in] Ndim Mat is an Ndim x Ndim matrix
!--------------------------------------------------------------------

  subroutine copy_select_rows(V, Mat, P, opn, Ndim)
    Implicit none
    Integer, INTENT(IN) :: opn, Ndim
    Complex (Kind = Kind(0.D0)), INTENT(INOUT) :: V(opn, Ndim)
    Complex (Kind = Kind(0.D0)), Dimension(:,:), INTENT(IN) :: Mat
    Integer , Dimension(:), INTENT(IN) :: P
    Integer :: n
    
    Do n = 1, opn
       V(n,:) = Mat(:, P(n))
    Enddo
    
  end subroutine

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function copies select columns to the destination matrix V from
!> the source matrix Mat. The decision on which columns to copy is determined
!> by the vector P.
!>
!> @note 
!> The code can be replaced by an array expression, but it is not
!> immediately clear that it is faster. The generated code looks longer.
!>
!> @param[inout] V storage for the result matrix
!> @param[in] Mat Where to read those columns
!> @param[in] P A vector with which columns to copy
!> @param[in] opn The length of the vector P
!> @param[in] Ndim Mat is an Ndim x Ndim matrix
!--------------------------------------------------------------------

  subroutine copy_select_columns(V, Mat, P, opn, Ndim)
    Implicit none
    Integer, INTENT(IN) :: opn, Ndim
    Complex (Kind = Kind(0.D0)), INTENT(INOUT) :: V(opn, Ndim)
    Complex (Kind = Kind(0.D0)), Dimension(:,:), INTENT(IN) :: Mat
    Integer , Dimension(:), INTENT(IN) :: P
    Integer :: n
    
    Do n = 1, opn
       V(n,:) = Mat(P(n), :)
    Enddo
    
  end subroutine
  
!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function performs a matrix multiplication of U and V and
!> writes the result to columns in mat specified by P.
!
!> @param[in] V
!> @param[in] U
!> @param[in] P a vector with the columns that we write to
!> @param[inout] Mat The Matrix that we update
!> @param[in] opn The length of the vector P
!> @param[in] Ndim Mat is an Ndim x Ndim matrix
!--------------------------------------------------------------------
  
  subroutine opmult(V, U, P, Mat, opn, Ndim)
    Implicit none
    Integer, INTENT(IN) :: opn, Ndim
    Complex (Kind = Kind(0.D0)), INTENT(IN) :: V(opn, Ndim)
    Complex (Kind = Kind(0.D0)), Dimension(:, :), INTENT(IN) :: U
    Complex (Kind = Kind(0.D0)), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Integer, INTENT(IN) :: P(opn)
    Integer :: n,i

    select case (opn)
    case (1)
        DO I = 1, Ndim
            Mat(P(1), I) = U(1,1) * V(1, I)
        enddo
    case (2)
        DO I = 1, Ndim
            Mat(P(1), I) = U(1,1) * V(1, I) + U(1,2) * V(2, I)
            Mat(P(2), I) = U(2,1) * V(1, I) + U(2,2) * V(2, I)
        enddo
    case default
    do n = 1, opn
        DO I = 1, Ndim
            Mat(P(n), I) = Sum(U(n,:) * V(:, I))
        enddo
    enddo
    end select
    
  end subroutine

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief
!> This function performs a matrix multiplication of conjg(U) and V and
!> writes the result to rows in mat specified by P.
!
!> @param[in] V
!> @param[in] U
!> @param[in] P a vector wiht the rows that we write to
!> @param[inout] Mat The Matrix that we update
!> @param[in] opn The length of the vector P
!> @param[in] Ndim Mat is an Ndim x Ndim matrix
!--------------------------------------------------------------------

  subroutine opmultct(V, U, P, Mat, opn, Ndim)
    Implicit none
    Integer, INTENT(IN) :: opn, Ndim
    Complex (Kind = Kind(0.D0)), INTENT(IN) :: V(opn, Ndim)
    Complex (Kind = Kind(0.D0)), Dimension(:, :), INTENT(IN) :: U
    Complex (Kind = Kind(0.D0)), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Integer, INTENT(IN) :: P(opn)
    Integer :: n, i

    select case (opn)
    case (1)
        DO I = 1, Ndim
            Mat(I, P(1)) = conjg(U(1,1)) * V(1, I)
        enddo
    case (2)
        DO I = 1, Ndim
            Mat(I, P(1)) = conjg(U(1,1)) * V(1, I) + conjg(U(1,2)) * V(2, I)
            Mat(I, P(2)) = conjg(U(2,1)) * V(1, I) + conjg(U(2,2)) * V(2, I)
        enddo
    case default
    do n = 1, opn
        DO I = 1, Ndim 
            Mat(I, P(n)) = Dot_Product(U(n, :), V(:, I))!this involves a conjg
        enddo
    enddo
    end select
        
  end subroutine

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function performs a matrix multiplication of U and V and
!> writes the result scaled by an entry given by Z to rows in mat
!> specified by P.
!
!> @param[in] V
!> @param[in] U
!> @param[in] P a vector wiht the rows that we write to
!> @param[inout] Mat The Matrix that we update
!> @param[in] Z A vector that usually contains exponentials
!> @param[in] opn The length of the vector P
!> @param[in] Ndim Mat is an Ndim x Ndim matrix
!--------------------------------------------------------------------

  subroutine opexpmult(V, U, P, Mat, Z, opn, Ndim)
    Implicit none
    Integer, INTENT(IN) :: opn, Ndim
    Complex (Kind = Kind(0.D0)), INTENT(IN) :: Z(opn)
    Complex (Kind = Kind(0.D0)), INTENT(IN) :: V(opn, Ndim)
    Complex (Kind = Kind(0.D0)), Dimension(:, :), INTENT(IN) :: U
    Complex (Kind = Kind(0.D0)), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Integer, INTENT(IN) :: P(opn)
    Integer :: n, i
    Complex (Kind = Kind(0.D0)) lexp

    select case (opn)
    case (1)
        DO I = 1, Ndim
            Mat(I, P(1)) = Z(1) * U(1, 1) * V(1, I)
        enddo
    case (2)
        DO I = 1, Ndim
            Mat(I, P(1)) = Z(1) * (U(1, 1) * V(1, I) + U(2, 1) * V(2, I))
            Mat(I, P(2)) = Z(2) * (U(1, 2) * V(1, I) + U(2, 2) * V(2, I))
        enddo
    case default
    do n = 1, opn
        lexp = Z(n)
        DO I = 1, Ndim
            Mat(I, P(n)) = lexp * Sum(U(:, n) * V(:, I))
        enddo
    Enddo
    end select

  end subroutine

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief
!> This function performs a matrix multiplication of conjg(U) and V and
!> writes the result scaled by an entry given by Z to columns in mat
!> specified by P.
!
!> @param[in] V
!> @param[in] U
!> @param[in] P a vector wiht the rows that we write to
!> @param[inout] Mat The Matrix that we update
!> @param[in] Z A vector that usually contains exponentials
!> @param[in] opn The length of the vector P
!> @param[in] Ndim Mat is an Ndim x Ndim matrix
!--------------------------------------------------------------------

  subroutine opexpmultct(V, U, P, Mat, Z, opn, Ndim)
    Implicit none
    Integer, INTENT(IN) :: opn, Ndim
    Complex (Kind = Kind(0.D0)), INTENT(IN) :: Z(opn)
    Complex (Kind = Kind(0.D0)), INTENT(IN) :: V(opn, Ndim)
    Complex (Kind = Kind(0.D0)), Dimension(:, :), INTENT(IN) :: U
    Complex (Kind = Kind(0.D0)), INTENT(INOUT) :: Mat (Ndim, Ndim)
    Integer, INTENT(IN) :: P(opn)
    Integer :: n, i
    Complex (Kind = Kind(0.D0)) lexp

    select case (opn)
    case (1)
        DO I = 1, Ndim
            Mat(P(1), I) = Z(1) * conjg(U(1, 1)) * V(1, I)
        enddo
    case (2)
        DO I = 1, Ndim
            Mat(P(1), I) = Z(1) * (conjg(U(1, 1)) * V(1, I) + conjg(U(2, 1)) * V(2, I))
            Mat(P(2), I) = Z(2) * (conjg(U(1, 2)) * V(1, I) + conjg(U(2, 2)) * V(2, I))
        enddo
    case default
    do n = 1, opn
        lexp = Z(n)
        DO I = 1, Ndim
            Mat(P(n), I) = lexp * Dot_Product(U(:, n), V(:, I))
        enddo
    Enddo
    end select

  end subroutine

!--------------------------------------------------------------------
!> @author
!> 
!
!> @brief 
!> Out Mat = Mat*exp(spin*Op)
!
!> @param[inout] Mat
!> @param[in] Op The Operator that we exponentiate
!> @param[in] spin The spin direction that we consider
!> @param[in] Ndim The dimension of the matrix Mat
!--------------------------------------------------------------------
  subroutine Op_mmultL(Mat,Op,spin,Ndim)
    Implicit none 
    Integer :: Ndim
    Type (Operator) , INTENT(IN)   :: Op
    Complex (Kind=8), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind=8), INTENT(IN)   :: spin

    ! Local 
    Complex (Kind=8) :: VH(Op%N,Ndim), Z(Op%N)

    ! In  Mat
    ! Out Mat = Mat*exp(spin*Op)

    call copy_select_rows(VH, Mat, Op%P, Op%N, Ndim)
    Z = exp(Op%g * spin * Op%E)
    call opexpmult(VH, Op%U, Op%P, Mat, Z, Op%N, Ndim)
    call copy_select_rows(VH, Mat, Op%P, Op%N, Ndim)
    call opmultct(VH, Op%U, Op%P, Mat, Op%N, Ndim)

  end subroutine Op_mmultL

  !--------------------------------------------------------------------
!> @author
!>
!
!> @brief 
!> Out Mat = exp(spin*Op)*Mat
!
!> @param[inout] Mat
!> @param[in] Op The Operator that we exponentiate
!> @param[in] spin The spin direction that we consider
!> @param[in] Ndim The dimension of the matrix Mat
!--------------------------------------------------------------------
  subroutine Op_mmultR(Mat,Op,spin,Ndim)
    Implicit none
    Integer :: Ndim
    Type (Operator) , INTENT(IN )   :: Op
    Complex (Kind=8), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind=8), INTENT(IN )   :: spin

    ! Local 
    Complex (Kind=8) :: VH(Op%N,Ndim), Z(Op%N)
    
    ! In  Mat
    ! Out Mat = exp(spin*Op)*Mat
    call copy_select_columns(VH, Mat, Op%P, Op%N, Ndim)
    Z = exp(Op%g * spin * Op%E)
    call opexpmultct(VH, Op%U, Op%P, Mat, Z, Op%N, Ndim)    
    call copy_select_columns(VH, Mat, Op%P, Op%N, Ndim)
    call opmult(VH, Op%U, Op%P, Mat, Op%N, Ndim)

  end subroutine Op_mmultR

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function fills the arrays ExpOp nd ExpMop according to the data in Op
!
!> @param[inout] ExpOp
!> @param[inout] ExpMop
!> @param[in] Op The Operator whose eigenvalues we exponentiate
!> @param[in] spin The spin direction that we consider
!--------------------------------------------------------------------
subroutine FillExpOps(ExpOp, ExpMop, Op, spin)
    Implicit none
    Type (Operator) , INTENT(IN) :: Op
    Complex(kind = kind(0.D0)), INTENT(INOUT) :: ExpOp(Op%N), ExpMop(Op%N)
    Real(kind = kind(0.D0)), Intent(in) :: spin
    Integer :: n

       do n = 1, Op%N
          ExpOp(n) = cmplx(1.d0, 0.d0, kind(0.D0))
          ExpMOp(n) = cmplx(1.d0, 0.d0, kind(0.D0))
          if ( n <= OP%N_non_Zero) then
            ExpOp(n) = exp(Op%g*(Op%E(n)*spin))
            ExpMop(n) = 1.D0/ExpOp(n)
          endif
       enddo

end subroutine

  Subroutine Op_Wrapup(Mat,Op,spin,Ndim,N_Type)

    Implicit none 

    Integer :: Ndim
    Type (Operator) , INTENT(IN )   :: Op
    Complex (Kind=8), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind=8), INTENT(IN )   :: spin
    Integer, INTENT(IN) :: N_Type

    ! Local 
    Complex (Kind=8) :: ExpOp(Op%N), ExpMop(Op%N), VH(Op%N,Ndim) !, zdotu, zdotc
    
!     nop=size(Op%U,1)
    
    !!!!! N_Type ==1
    !    exp(Op%g*spin*Op%E)*(Op%U^{dagger})*Mat*Op%U*exp(-Op%g*spin*Op%E)
    !    
    !!!!!
    !!!!! N_Type == 2
    !    Op%U * Mat * (Op%U^{dagger})
    !!!!!
    If (N_type == 1) then
        call FillExpOps(ExpOp, ExpMop, Op, spin)
        call copy_select_rows(VH, Mat, Op%P, Op%N, Ndim)
        call opexpmult(VH, Op%U, Op%P, Mat, ExpMOp, Op%N, Ndim)
        call copy_select_columns(VH, Mat, Op%P, Op%N, Ndim)
        call opexpmultct(VH, Op%U, Op%P, Mat, ExpOp, Op%N, Ndim)
    elseif (N_Type == 2) then
        call copy_select_rows(VH, Mat, Op%P, Op%N, Ndim)
        call opmultct(VH, Op%U, Op%P, Mat, Op%N, Ndim)
        call copy_select_columns(VH, Mat, Op%P, Op%N, Ndim)
        call opmult(VH, Op%U, Op%P, Mat, Op%N, Ndim)
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
    Complex (Kind=8) :: VH(Op%N,Ndim), ExpOp(Op%N), ExpMop(Op%N), ExpHere !, zdotu, zdotc
!     Complex (Kind=8) :: alpha, beta, tmp2(Op%N,Ndim)
    Integer :: n, i
    
!     nop=size(Op%U,1)
    
    !!!!! N_Type == 1
    !    Op%U*exp(-Op%g*spin*Op%E)*Mat*exp(Op%g*spin*Op%E)*(Op%U^{dagger})
    !    
    !!!!!
    !!!!! N_Type == 2
    !    (Op%U^{dagger}) * Mat * Op%U
    !!!!!
    If (N_type == 1) then
       call FillExpOps(ExpOp, ExpMop, Op, spin)
       
       Do n = 1,Op%N
          expHere=ExpOp(n)
          VH(n, :) = ExpHere * Mat(:, Op%P(n))
       Enddo
       
       ! ZGEMM might be the better multiplication, but the additional copy precess seem to be to expensive
!        alpha = cmplx (1.0d0,0.0d0)
!        beta = cmplx (0.0d0,0.0d0)
!        CALL ZGEMM('T','C',Ndim,Op%N,Op%N,alpha,VH,Op%N,Op%U,nop,beta,tmp,Ndim)
       call opmultct(VH, Op%U, Op%P, Mat, Op%N, Ndim)

       Do n = 1,Op%N
          ExpHere=ExpMOp(n)
          VH(n, :) = ExpHere * Mat(Op%P(n), :)
       Enddo
       
       ! ZGEMM might be the better multiplication, but the additional copy precess seem to be to expensive
!        alpha = cmplx (1.0d0,0.0d0)
!        beta = cmplx (0.0d0,0.0d0)
!        CALL ZGEMM('T','T',Ndim,Op%N,Op%N,alpha,VH,Op%N,Op%U,nop,beta,tmp,Ndim)
       call opmult(VH, Op%U, Op%P, Mat, Op%N, Ndim)
    elseif (N_Type == 2) then
    call copy_select_rows(VH, Mat, Op%P, Op%N, Ndim)
       
       ! ZGEMM might be the better multiplication, but the additional copy precess seem to be to expensive
!        alpha = cmplx (1.0d0,0.0d0)
!        beta = cmplx (0.0d0,0.0d0)
!        CALL ZGEMM('T','N',Ndim,Op%N,Op%N,alpha,VH,Op%N,Op%U,nop,beta,tmp,Ndim)
    select case (Op%N)
    case (1)
        DO I = 1, Ndim
            Mat(I, Op%P(1)) = Op%U(1,1) * VH(1, I)
        enddo
    case (2)
        DO I = 1, Ndim
            Mat(I, Op%P(1)) = Op%U(1, 1) * VH(1, I) + Op%U(2, 1) * VH(2, I)
            Mat(I, Op%P(2)) = Op%U(1, 2) * VH(1, I) + Op%U(2, 2) * VH(2, I)
        enddo
    case default
       do n = 1,Op%N
          DO I = 1,Ndim
             Mat(I, Op%P(n))  =  Sum(VH(:, I) * Op%U(:, n))
          enddo
       Enddo
    end select
      call copy_select_columns(VH, Mat, Op%P, Op%N, Ndim)
       
       ! ZGEMM might be the better multiplication, but the additional copy precess seem to be to expensive
!        alpha = cmplx (1.0d0,0.0d0)
!        beta = cmplx (0.0d0,0.0d0)
!        CALL ZGEMM('C','N',Op%N,Ndim,Op%N,alpha,Op%U,nop,VH,Op%N,beta,tmp2,Op%N)
        select case (Op%N)
        case (1)
            DO I = 1, Ndim
                Mat(Op%P(1), I) = conjg(Op%U(1,1)) * VH(1, I)
            enddo
        case (2)
            DO I = 1, Ndim
                Mat(Op%P(1), I) = conjg(Op%U(1,1)) * VH(1, I) + conjg(Op%U(2,1)) * VH(2, I)
                Mat(Op%P(2), I) = conjg(Op%U(1,2)) * VH(1, I) + conjg(Op%U(2,2)) * VH(2, I)
            enddo
        case default
        do n = 1, Op%N
            DO I = 1, Ndim 
                Mat(Op%P(n), I) = Dot_Product(Op%U(n, :), VH(:, I))!this involves a conjg
            enddo
        enddo
        end select
    endif
    
  end Subroutine Op_Wrapdo

end Module Operator_mod
