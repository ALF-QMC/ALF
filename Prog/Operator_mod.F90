!  Copyright (C) 2016 - 2020 The ALF project
! 
!  This file is part of the ALF project.
! 
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!     
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!     
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Defines the Operator type, and provides a number of operations on this type. 
!
!--------------------------------------------------------------------
Module Operator_mod


  Use MyMats
  Use Fields_mod

  Implicit none

  

  Type Operator
     Integer          :: N, N_non_zero
     logical          :: diag
     complex (Kind=Kind(0.d0)), pointer :: O(:,:), U (:,:), M_exp(:,:,:), E_exp(:,:)
     Real    (Kind=Kind(0.d0)), pointer :: E(:)
     Integer, pointer :: P(:)
     complex (Kind=Kind(0.d0)) :: g
     complex (Kind=Kind(0.d0)) :: alpha
     Integer          :: Type 
     ! P is an N X Ndim matrix such that  P.T*O*P*  =  A  
     ! P has only one non-zero entry per column which is specified by P
     ! All in all.   g * Phi(s,type) * ( c^{dagger} A c  + alpha )
     ! The variable Type allows you to define the type of HS. 
     ! The first N_non_zero elements of diagonal matrix E are non-zero. The rest vanish.

     ! !!!!! M_exp and E_exp  are for storage   !!!!!
     ! If Type =1   then the Ising field  takes the values  s = +/- 1 
     !              and M_exp   has dimensions  M_exp(N,N,-1:1)
     ! If Type =2   then the Ising field  takes the values  s = +/- 1,  +/- 2
     !              and M_exp   has dimensions  M_exp(N,N,-2:2)
     ! M_exp(:,:,s) =  e^{g * Phi(s,type) *  O(:,:) }
     !
     ! 
     ! E_exp(:,s) = e^{g * Phi(s,type) *E(:) } and has dimensions E_exp(N,-1:1)  for Type = 1
     !              and dimensions E_exp(N,-2:2)   for Type = 2
     !
     ! !!! If Type .neq. 1,2  then  E_exp  and  M_exp  are  not allocated !!!
     
  end type Operator

  
Contains

!!$  Real (Kind=Kind(0.d0)) :: Phi(-2:2,2),  Gaml(-2:2,2)
!!$  Integer ::  NFLIPL(-2:2,3)
!!$      Subroutine Op_setHS
!!$
!!$        Implicit none
!!$        !Local
!!$        Integer :: n
!!$        
!!$
!!$        Phi = 0.d0
!!$        do n = -2,2
!!$           Phi(n,1) = real(n,Kind=Kind(0.d0))
!!$        enddo
!!$        Phi(-2,2) = - SQRT(2.D0 * ( 3.D0 + SQRT(6.D0) ) )
!!$        Phi(-1,2) = - SQRT(2.D0 * ( 3.D0 - SQRT(6.D0) ) )
!!$        Phi( 1,2) =   SQRT(2.D0 * ( 3.D0 - SQRT(6.D0) ) )
!!$        Phi( 2,2) =   SQRT(2.D0 * ( 3.D0 + SQRT(6.D0) ) )
!!$        
!!$        Do n = -2,2
!!$           gaml(n,1) = 1.d0
!!$        Enddo
!!$        GAML(-2,2) = 1.D0 - SQRT(6.D0)/3.D0
!!$        GAML( 2,2) = 1.D0 - SQRT(6.D0)/3.D0
!!$        GAML(-1,2) = 1.D0 + SQRT(6.D0)/3.D0
!!$        GAML( 1,2) = 1.D0 + SQRT(6.D0)/3.D0
!!$        
!!$        NFLIPL(-2,1) = -1
!!$        NFLIPL(-2,2) =  1
!!$        NFLIPL(-2,3) =  2
!!$        
!!$        NFLIPL(-1,1) =  1
!!$        NFLIPL(-1,2) =  2
!!$        NFLIPL(-1,3) = -2
!!$        
!!$        NFLIPL( 1,1) =  2
!!$        NFLIPL( 1,2) = -2
!!$        NFLIPL( 1,3) = -1
!!$        
!!$        NFLIPL( 2,1) = -2
!!$        NFLIPL( 2,2) = -1
!!$        NFLIPL( 2,3) =  1 
!!$
!!$      end Subroutine Op_setHS

  

!--------------------------------------------------------------------
!> @author
!> 
!
!> @brief 
!> calculate the phase of a given set of operators and HS fields.
!
!> @param[inout] Phase  Complex
!> * On entry: phase of \f$ \prod_{f} \det M_f(C) \f$, f is the flavor index.
!> * On exit:  phase of  \f$ W(C)  =   \left[ \left( \prod_{n,\tau,f}  \exp \left[ g_f(n) \alpha_f(n) \phi(\sigma(n,\tau)) \right] \right) \det(M_f(C))\right]^{N_{SUN}} \prod_{n,\tau }\gamma(\sigma(n,\tau)) \f$ 
!> @param[in] Op_V  Dimension(:,:)  Type(Operator)
!> * List of interaction operators. OP_V(n,f) has no tau index.  
!> @param[in] Nsigma 
!> Type(Fields)
!> * Fields
!> @param[in] N_SUN
!> Integer
!> * Number of colors
!--------------------------------------------------------------------
  Subroutine  Op_phase(Phase,OP_V,Nsigma,N_SUN) 
    Implicit none

    Complex  (Kind=Kind(0.d0)), Intent(Inout) :: Phase
    Integer,           Intent(IN)    :: N_SUN
    Type  (Fields),    Intent(IN)    :: Nsigma
    Type (Operator),   dimension(:,:), Intent(In) :: Op_V
    Real  (Kind=Kind(0.d0))                       :: angle
    
    Integer :: n, nf, nt
    
    do nf = 1,Size(Op_V,2)
       do n = 1,size(Op_V,1)
          do nt = 1,size(nsigma%f,2)
             angle = Aimag( Op_V(n,nf)%g * Op_V(n,nf)%alpha ) * nsigma%Phi(n,nt) 
             Phase = Phase*CMPLX(cos(angle),sin(angle), Kind(0.D0))
          enddo
       enddo
    enddo
    Phase = Phase**N_SUN
    
  end Subroutine Op_phase
  
!--------------------------------------------------------------------
!> @author
!> 
!> @brief 
!> Set up _core_ data of the operator required for type=0 operators
!
!> @param[inout] Op
!> @param[in] N 
!--------------------------------------------------------------------

  Pure subroutine Op_make(Op,N)
    Implicit none
    Type (Operator), intent(INOUT) :: Op
    Integer, Intent(IN) :: N
    Allocate (Op%O(N,N), Op%P(N) )
    ! F.F.A  Op%M_exp and Op%E_exp are allocated  in Op_set once the type is available.
    
    Op%O = cmplx(0.d0, 0.d0, kind(0.D0))
    Op%P = 0
    Op%N = N
    Op%N_non_zero = N
    Op%g     = cmplx(0.d0,0.d0, kind(0.D0))
    Op%alpha = cmplx(0.d0,0.d0, kind(0.D0))
    Op%diag  = .false.
    Op%type = 0
  end subroutine Op_make

!--------------------------------------------------------------------
!> @author
!> 
!> @brief 
!> Deallocate the operator
!
!> @param[inout] Op
!> @param[in] N 
!--------------------------------------------------------------------

  Pure subroutine Op_clear(Op,N)
    Implicit none
    Type (Operator), intent(INOUT) :: Op
    Integer, Intent(IN) :: N
    Deallocate (Op%O, Op%P )

    If ( Op%Type >= 1)   deallocate(OP%M_exp,OP%E_exp,  OP%U, OP%E)

  end subroutine Op_clear

!--------------------------------------------------------------------
!> @author
!> 
!> @brief 
!> Setup storage for type=1,2 and 3 vertices.  Setup the exponential of the operator for
!> type 1 and 2
!
!> @param[inout] Op  Type (Operator)
!--------------------------------------------------------------------
  subroutine Op_set(Op)
    Implicit none
    Type (Operator), intent(INOUT) :: Op

    Complex (Kind=Kind(0.d0)), allocatable :: U(:,:), TMP(:, :)
    Real    (Kind=Kind(0.d0)), allocatable :: E(:)
    Real    (Kind=Kind(0.d0)) :: Zero = 1.D-9 !, Phi(-2:2)
    Integer :: N, I, J, np,nz
    Complex (Kind=Kind(0.d0)) :: Z
    Type  (Fields)   :: nsigma_single
    
    
    Call nsigma_single%make(1,1)

    N = OP%N
    Allocate ( Op%U(N,N), Op%E(N) )
    Op%U = cmplx(0.d0, 0.d0, kind(0.D0))
    Op%E = 0.d0
    
    If (Op%N > 1) then
       N = Op%N
       Op%diag = .true.
       do I=1,N
          do J=i+1,N
             ! Binary comparison is OK here as Op%O was initialized to zero during Op_make.
             if (Op%O(i,j) .ne. cmplx(0.d0,0.d0, kind(0.D0)) .or. Op%O(j,i) .ne. cmplx(0.d0,0.d0, kind(0.D0))) Op%diag=.false.
          enddo
       enddo
       if (Op%diag) then
          do I=1,N
             Op%E(I)=DBLE(Op%O(I,I))
             Op%U(I,I)=cmplx(1.d0,0.d0, kind(0.D0))
          enddo
          Op%N_non_zero = N
          ! FFA Why do we assume that Op%N_non_zero = N for a diagonal operator? 
       else
          Allocate (U(N,N), E(N), TMP(N, N))
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
          ! Write(6,*) "Op_set", np,N
          TMP = Op%U ! that way we have the changes to the determinant due to the permutation
          Z = Det_C(TMP, N)
          ! Scale Op%U to be in SU(N) 
          DO I = 1, N
             Op%U(I,1) = Op%U(I, 1)/Z 
          ENDDO
          deallocate (U, E, TMP)
          ! Op%U,Op%E)
          ! Write(6,*) 'Calling diag 1'
       endif
    else
       Op%E(1)   = REAL(Op%O(1,1), kind(0.D0))
       Op%U(1,1) = cmplx(1.d0, 0.d0 , kind(0.D0))
       Op%N_non_zero = 1
       Op%diag = .true.
    endif
    select case(OP%type)
    case(1)
       Allocate(Op%E_exp(Op%N, -Op%type : Op%type), Op%M_exp(Op%N, Op%N, -Op%type : Op%type))
       nsigma_single%t(1) = 1
       Do I=1,Op%type
          nsigma_single%f(1,1) = real(I,kind=kind(0.d0))
          do n = 1, Op%N
             Op%E_exp(n,I) = cmplx(1.d0, 0.d0, kind(0.D0))
             Op%E_exp(n,-I) = cmplx(1.d0, 0.d0, kind(0.D0))
             if ( n <= Op%N_non_Zero) then
                !Op%E_exp(n,I) = exp(Op%g*Op%E(n)*Phi_st(I,1))
                Op%E_exp(n,I) = exp(Op%g*Op%E(n)*nsigma_single%Phi(1,1))
                Op%E_exp(n,-I) = 1.D0/Op%E_exp(n,I)
             endif
          enddo
          !call Op_exp(Op%g*Phi_st( I,1),Op,Op%M_exp(:,:,I))
          !call Op_exp(Op%g*Phi_st(-I,1),Op,Op%M_exp(:,:,-I))
          call Op_exp( Op%g*nsigma_single%Phi(1,1),Op,Op%M_exp(:,:,I))
          call Op_exp(-Op%g*nsigma_single%Phi(1,1),Op,Op%M_exp(:,:,-I))
       enddo
    case(2)
       Allocate(Op%E_exp(Op%N, -Op%type : Op%type), Op%M_exp(Op%N, Op%N, -Op%type : Op%type))
       nsigma_single%t(1) = 2
       Do I=1,Op%type
          nsigma_single%f(1,1) = real(I,kind=kind(0.d0))
          do n = 1, Op%N
             Op%E_exp(n,I) = cmplx(1.d0, 0.d0, kind(0.D0))
             Op%E_exp(n,-I) = cmplx(1.d0, 0.d0, kind(0.D0))
             if ( n <= Op%N_non_Zero) then
                !Op%E_exp(n,I)  = exp(Op%g*Op%E(n)*Phi_st(I,2))
                Op%E_exp(n,I) = exp(Op%g*Op%E(n)*nsigma_single%Phi(1,1))
                Op%E_exp(n,-I) = 1.D0/Op%E_exp(n,I)
             endif
          enddo
          call Op_exp( Op%g*nsigma_single%Phi(1,1),Op,Op%M_exp(:,:,I))
          call Op_exp(-Op%g*nsigma_single%Phi(1,1),Op,Op%M_exp(:,:,-I))
          !call Op_exp(Op%g*Phi_st( I,2),Op,Op%M_exp(:,:,I))
          !call Op_exp(Op%g*Phi_st(-I,2),Op,Op%M_exp(:,:,-I))
       enddo
    case default
    end select

    Call nsigma_single%clear() 

  end subroutine Op_set

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!
!> @brief
!> Calculate the exponentiated operator and returns the full matrix.
!>
!> @details
!> It is assumed that  the U and E arrays of the
!> operator are present.
!> @param[in]  g  Complex
!> @param[in]  Op Type(Operator)
!> @param[out]
!> Mat   Complex, Dimension(:,:)
!> * On output  \f$ M = e^{g O} \f$
!
!--------------------------------------------------------------------
  Subroutine Op_exp(g,Op,Mat)

    Implicit none 

    Type (Operator), Intent(IN)  :: Op
    Complex (Kind=Kind(0.d0)), Dimension(:,:), INTENT(OUT) :: Mat
    Complex (Kind=Kind(0.d0)), INTENT(IN) :: g
    Complex (Kind=Kind(0.d0)) :: Z, Z1, y, t
    Complex (Kind=Kind(0.d0)), allocatable, dimension(:,:) :: c

    Integer :: n, j, I, iters
    
    iters = Op%N
    Mat = cmplx(0.d0, 0.d0, kind(0.D0))
    if (Op%diag) then
      Do n = 1, iters
        Mat(n,n)=exp(g*Op%E(n))
      enddo
    else
      Allocate (c(iters, iters))
      c = 0.D0
      Do n = 1, iters
        Z = exp(g*Op%E(n))
        do J = 1, iters
            Z1 = Z*conjg(Op%U(J,n))
            do I = 1, iters
               ! This performs Kahan summation so as to improve precision.
               y = Z1 * Op%U(I, n) - c(I, J)
               t = Mat(I, J) + y
               c(I, J) = (t - Mat(I,J)) - y
               Mat(I, J) = t
               !  Mat(I, J) = Mat(I, J) + Z1 * Op%U(I, n)
            enddo
            ! Mat(1:iters, J) = Mat(1:iters, J) + Z1 * Op%U(1:iters, n)
         enddo
      enddo
      Deallocate(c)
    endif
  end subroutine Op_exp

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!
!> @brief 
!> Out Mat = Mat* Op ( exp(spin*g* P^T O T) )
!
!> @param[inout] Mat Complex Dimension(:,:)
!> * On exit Mat = Mat*Op ( exp(spin*g* P^T O T) )
!> @param[in] Op Type(Operator)
!> * The Operator containing g and the sparse matrix P^T O P 
!> @param[in] spin Real
!> * The field
!> @param[in] cop  Character
!> * cop = N,  Op = None
!> * cop = T,  Op = Transposed
!> * cop = C,  Op = Transposed + Complex conjugation
!> 
  
!--------------------------------------------------------------------
  subroutine Op_mmultL(Mat,Op,spin,cop)
    Implicit none 
    Integer :: Ndim
    Type (Operator)          , INTENT(IN)    :: Op
    Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Mat (:,:)
    Real    (Kind=Kind(0.d0)), INTENT(IN)    :: spin
    Character                , Intent(IN)    :: cop

    ! Local 
    Integer :: I, N1, N2, sp
    Complex (Kind=Kind(0.d0)) :: ExpMat (Op%n,Op%n)
    
    N1=size(Mat,1)
    N2=size(Mat,2)

    ! In  Mat
    ! Out Mat = Mat*exp(spin*Op)
    
    ! quick return if possible
    if ( abs(OP%g) < 1.D-12 ) return

    if ( op%type < 3 ) then
       sp = nint(spin)
       if ( Op%diag ) then
          do I=1,Op%N
             if ( cop == 'c' .or. cop =='C' ) then
                call ZSCAL(N1,conjg(Op%E_exp(I,sp)),Mat(1,Op%P(I)),1)
             else
                call ZSCAL(N1,Op%E_exp(I,sp),Mat(1,Op%P(I)),1)
             endif
          enddo
       else
          call ZSLGEMM('r',cop,Op%N,N1,N2,Op%M_exp(:,:,sp),Op%P,Mat)
       endif
    else
       if ( Op%diag ) then
          do I=1,Op%N
             if ( cop == 'c' .or. cop =='C' ) then
                call ZSCAL(N1,conjg(exp(spin*Op%g*Op%E(I))),Mat(1,Op%P(I)),1)
             else
                call ZSCAL(N1,exp(spin*Op%g*Op%E(I)),Mat(1,Op%P(I)),1)
             endif
          enddo
       else
          call Op_exp(Op%g*spin,Op,expmat)
          call ZSLGEMM('r',cop,Op%N,N1,N2,expmat,Op%P,Mat)
       endif
    endif
  end subroutine Op_mmultL

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief 
!> Out Mat =  Op ( exp(spin*g* P^T O T) ) * Mat
!
!> @param[inout] Mat Complex Dimension(:,:)
!> * On exit Mat = Op ( exp(spin*g* P^T O T) )* Mat
!> @param[in] Op Type(Operator)
!> * The Operator containing g and the sparse matrix P^T O P 
!> @param[in] spin Real
!> * The field
!> @param[in] cop  Character
!> * cop = N,  Op = None
!> * cop = T,  Op = Transposed
!> * cop = C,  Op = Transposed + Complex conjugation
!> 
!--------------------------------------------------------------------
  subroutine Op_mmultR(Mat,Op,spin,cop)
    Implicit none
    Type (Operator)          , INTENT(IN )   :: Op
    Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Mat (:,:)
    Real    (Kind=Kind(0.d0)), INTENT(IN )   :: spin
    Character                , Intent(IN)    :: cop

    ! Local 
    Integer :: I, N1, N2, sp
    Complex (Kind=Kind(0.d0)) :: ExpMat (Op%n,Op%n)

    N1=size(Mat,1)
    N2=size(Mat,2)
    
    ! In  Mat
    ! Out Mat = exp(spin*Op)*Mat
    
    ! quick return if possible
    if ( abs(OP%g) < 1.D-12 ) return

    if ( op%type < 3 ) then
       sp = nint(spin)
       if ( Op%diag ) then
          do I=1,Op%N
             if ( cop == 'c' .or. cop =='C' ) then
                call ZSCAL(N2,conjg(Op%E_exp(I,sp)),Mat(Op%P(I),1),N1)
             else
                call ZSCAL(N2,Op%E_exp(I,sp),Mat(Op%P(I),1),N1)
             endif
          enddo
       else
          call ZSLGEMM('L',cop,Op%N,N1,N2,Op%M_exp(:,:,sp),Op%P,Mat)
       endif
    else
       if ( Op%diag ) then
          do I=1,Op%N
             if ( cop == 'c' .or. cop =='C' ) then
                call ZSCAL(N2,conjg(exp(spin*Op%g*Op%E(I))),Mat(Op%P(I),1),N1)
             else
                call ZSCAL(N2,exp(spin*Op%g*Op%E(I)),Mat(Op%P(I),1),N1)
             endif
          enddo
       else
          call Op_exp(Op%g*spin,Op,expmat)
          call ZSLGEMM('L',cop,Op%N,N1,N2,expmat,Op%P,Mat)
       endif
    endif
  end subroutine Op_mmultR

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief 
!> Wrapup the Green function
!>
!> @param[inout] Mat(Ndim,Ndim)  Complex
!> \verbatim
!>  N_type = 1, Mat = exp(Op%g*spin*Op%E)*(Op%U^{dagger}) * Mat * Op%U*exp(-Op%g*spin*Op%E)
!>  N_type = 2, Mat = Op%U * Mat * (Op%U^{dagger})
!> \endverbatim
!>  **Do not** mix up  \p N_type  and    \p OP\%type 
!> @param[in] Op Type(Operator)
!> \verbatim
!>  The operator containing g, U, P
!> \endverbatim
!> @param[in] spin Real
!> \verbatim
!>  The field
!> \endverbatim
!> @param[in] Ndim Integer
!> 
!--------------------------------------------------------------------
  Subroutine Op_Wrapup(Mat,Op,spin,Ndim,N_Type)

    Implicit none 

    Integer :: Ndim
    Type (Operator)           , INTENT(IN )   :: Op
    Complex (Kind=Kind(0.d0)) , INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind=Kind(0.d0)) , INTENT(IN )   :: spin
    Integer                   , INTENT(IN)    :: N_Type

    ! Local 
    Complex (Kind=Kind(0.d0)) :: VH1(Op%N,Op%N)
    Integer :: I,sp

    if ( op%type < 3 ) then
    !if ( op%type == 2 ) then
       sp = nint(spin)
       If (N_type == 1) then
          if(Op%diag) then
             do I=1,Op%N
                call ZSCAL(Ndim,Op%E_Exp(I,sp),Mat(Op%P(I),1),Ndim)
             enddo
             do I=1,Op%N
                call ZSCAL(Ndim,Op%E_Exp(I,-sp),Mat(1,Op%P(I)),1)
             enddo
          else
             Do i = 1,Op%N
                VH1(:,i)=Op%U(:,i)*Op%E_Exp(I,-sp)
             Enddo
             call ZSLGEMM('r','n',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
             Do i = 1,Op%N
                VH1(:,i)=Op%E_Exp(I, sp)*conjg(Op%U(:,i))
             Enddo
             call ZSLGEMM('l','T',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
          endif
       elseif (N_Type == 2 .and. .not. Op%diag) then
          call ZSLGEMM('l','n',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
          call ZSLGEMM('r','c',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
       endif
    else
       If (N_type == 1) then
          if(Op%diag) then
             do I=1,Op%N
                call ZSCAL(Ndim,exp( spin*Op%g*Op%E(I)),Mat(Op%P(I),1),Ndim)
             enddo
             do I=1,Op%N
                call ZSCAL(Ndim,exp(-spin*Op%g*Op%E(I)),Mat(1,Op%P(I)),1)
             enddo
          else
             Do i = 1,Op%N
                VH1(:,i)=Op%U(:,i)*exp(-spin*Op%g*Op%E(I)) 
             Enddo
             call ZSLGEMM('r','n',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
             Do i = 1,Op%N
                VH1(:,i)=exp(spin*Op%g*Op%E(I))*conjg(Op%U(:,i))
             Enddo
             call ZSLGEMM('l','T',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
          endif
       elseif (N_Type == 2 .and. .not. Op%diag) then
          call ZSLGEMM('l','n',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
          call ZSLGEMM('r','c',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
       endif
    endif
  end Subroutine Op_Wrapup

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief 
!> Wrapup the Green function
!>
!> @param[inout] Mat(Ndim,Ndim)  Complex
!> \verbatim
!>  N_type = 1, Mat = Op%U*exp(-Op%g*spin*Op%E)*Mat*exp(Op%g*spin*Op%E)*(Op%U^{dagger})
!>  N_type = 2, Mat = (Op%U^{dagger}) * Mat * Op%U
!> \endverbatim
!>  **Do not** mix up  \p N_type  and    \p OP\%type 
!> @param[in] Op Type(Operator)
!> \verbatim
!>  The operator containing g, U, P
!> \endverbatim
!> @param[in] spin Real
!> \verbatim
!>  The field
!> \endverbatim
!> @param[in] Ndim Integer
!> 
!--------------------------------------------------------------------

  Subroutine Op_Wrapdo(Mat,Op,spin,Ndim,N_Type)
    Implicit none 
    
    Integer :: Ndim
    Type (Operator) , INTENT(IN )   :: Op
    Complex (Kind = Kind(0.D0)), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real (Kind=Kind(0.d0)), INTENT(IN )   :: spin
    Integer, INTENT(IN) :: N_Type

    ! Local 
    Integer :: n, i, sp
    Complex (Kind = Kind(0.D0)) :: VH1(Op%N,OP%N)

    if ( op%type < 3 ) then
       sp = nint(spin)
       If (N_type == 1) then
          if(Op%diag) then
             do I=1,Op%N
                call ZSCAL(Ndim,Op%E_Exp(I,-sp),Mat(Op%P(I),1),Ndim)
             enddo
             do I=1,Op%N
                call ZSCAL(Ndim,Op%E_Exp(I, sp),Mat(1,Op%P(I)),1)
             enddo
          else
             Do n = 1,Op%N
                VH1(:,n)=Op%U(:,n)*Op%E_Exp(n,-sp)
             Enddo
             call ZSLGEMM('l','n',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
             Do n = 1,Op%N
                VH1(:,n)=Op%E_Exp(n, sp)*conjg(Op%U(:,n))
             Enddo
             call ZSLGEMM('r','T',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
          endif
       elseif (N_Type == 2 .and. .not. Op%diag) then
          call ZSLGEMM('r','n',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
          call ZSLGEMM('l','c',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
       endif
    else
       If (N_type == 1) then
          if(Op%diag) then
             do I=1,Op%N
                call ZSCAL(Ndim,exp(-spin*Op%g*Op%E(I)),Mat(Op%P(I),1),Ndim)
             enddo
             do I=1,Op%N
                call ZSCAL(Ndim,exp( spin*Op%g*Op%E(I)),Mat(1,Op%P(I)),1)
             enddo
          else
             Do n = 1,Op%N
                VH1(:,n)=Op%U(:,n)*exp(-spin*Op%g*Op%E(n)) 
             Enddo
             call ZSLGEMM('l','n',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
             Do n = 1,Op%N
                VH1(:,n)=exp(spin*Op%g*Op%E(n))*conjg(Op%U(:,n))
             Enddo
             call ZSLGEMM('r','T',Op%n,Ndim,Ndim,VH1,Op%P,Mat)
          endif
       elseif (N_Type == 2 .and. .not. Op%diag) then
          call ZSLGEMM('r','n',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
          call ZSLGEMM('l','c',Op%n,Ndim,Ndim,Op%U,Op%P,Mat)
       endif
    endif
  end Subroutine Op_Wrapdo

end Module Operator_mod
