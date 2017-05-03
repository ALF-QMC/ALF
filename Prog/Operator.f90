!  Copyright (C) 2016, 2017 The ALF project
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
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
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

Module Operator_mod

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Defines the Operator type, and provides a number of operations on this
!> type. 
!
!--------------------------------------------------------------------

  Use MyMats

  Implicit none

  Real (Kind=Kind(0.d0)) :: Phi(-2:2,2),  Gaml(-2:2,2)
  Integer ::  NFLIPL(-2:2,3)
  

  Type Operator
     Integer          :: N, N_non_zero
     complex (Kind=Kind(0.d0)), pointer :: O(:,:), U (:,:)
     Real    (Kind=Kind(0.d0)), pointer :: E(:)
     Integer, pointer :: P(:)
     complex (Kind=Kind(0.d0)) :: g
     complex (Kind=Kind(0.d0)) :: alpha
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
       Phi(n,1) = real(n,Kind=Kind(0.d0))
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
  Pure Subroutine  Op_phase(Phase,OP_V,Nsigma,N_SUN) 
    Implicit none

    Complex  (Kind=Kind(0.d0)), Intent(Inout) :: Phase
    Integer,           Intent(IN)    :: N_SUN
    Integer,           dimension(:,:), Intent(In) :: Nsigma
    Type (Operator),   dimension(:,:), Intent(In) :: Op_V
    Real  (Kind=Kind(0.d0))                   :: angle
    
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
  
!--------------------------------------------------------------------

  Pure subroutine Op_make(Op,N)
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

!--------------------------------------------------------------------

  Pure subroutine Op_clear(Op,N)
    Implicit none
    Type (Operator), intent(INOUT) :: Op
    Integer, Intent(IN) :: N
    Deallocate (Op%O, Op%U, Op%E, Op%P)
  end subroutine Op_clear 

!--------------------------------------------------------------------

  subroutine Op_set(Op)
    Implicit none
    Type (Operator), intent(INOUT) :: Op

    Complex (Kind=Kind(0.d0)), allocatable :: U(:,:), TMP(:,:), TMP2(:, :), TAU(:), WORK(:)
    Real    (Kind=Kind(0.d0)), allocatable :: E(:)
    Real    (Kind=Kind(0.d0)) :: Zero = 1.E-9, nm
    Integer :: N, I, J, np,nz, LWORK, INFO
    Complex (Kind=Kind(0.d0)) :: Z
    Complex (Kind=Kind(0.d0)) :: ZLADIV

    If (Op%N > 1) then
       N = Op%N
       LWORK= 2*N
       Allocate (U(N,N), E(N), TMP(N,N), TMP2(N,N), WORK(LWORK), TAU(N))
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
           write (*,*) N
       Op%N_non_zero = np
       !Write(6,*) "Op_set", np,N
       TMP = Op%U ! that way we have the changes to the determinant due to the permutation
!        DO I = 1, N
!        write (*,*) DBLE(Op%U(I,:))
!        ENDDO
!        STOP -2
       Z = Det_C(TMP, N)
!       write (*,*) Z, N, Z**N
!       write (*,*) nm,  Abs(nm)**(1.D0/N), Abs(nm)**(1.D0/N)*CMPLX(Cos(4*ATan(1.D0)/Op%N), SIN(4*ATan(1.D0)/Op%N), Kind(0.D0))
       ! Scale Op%U to be in SU(N)
       DO I = 1, N
            Op%U(I,1) = zladiv(Op%U(I,1),Z)
       ENDDO
       
       if(op%N > 2) then
           TMP = Op%U
           CALL ZGEQRF(N, N, TMP, N, TAU, WORK, LWORK, INFO)
    !        write (*, *) TAU
    !        Do I = 1, N
    !        write (*, *) DBLE(TMP(I, :))
    !        ENDDO
    !       CALL ZUNGQR(N, N, N, TMP, N, TAU, WORK, LWORK, INFO)
            DO I = 1, N-1
                TMP(I, N) = TAU(I)
            ENDDO
            TMP(N, N) = 1.D0 - TMP(N,N)*(1.D0 - TAU(N))! absorb last sign (stored in U(opn, opn)) into redefinition of tau
            Op%U = TMP
       endif
       deallocate (U, E, TMP, TAU, WORK)
       ! Op%U,Op%E)
       !Write(6,*) 'Calling diag 1'
    else
       Op%E(1)   = REAL(Op%O(1,1), kind(0.D0))
       Op%U(1,1) = cmplx(1.d0, 0.d0 , kind(0.D0))
       Op%N_non_zero = 1
    endif
  end subroutine Op_set

!--------------------------------------------------------------------


  Pure subroutine Op_exp(g,Op,Mat)
    Implicit none 
    Type (Operator), Intent(IN)  :: Op
    Complex (Kind=Kind(0.d0)), Dimension(:,:), INTENT(OUT) :: Mat
    Complex (Kind=Kind(0.d0)), INTENT(IN) :: g
    Complex (Kind=Kind(0.d0)) :: Z, Z1, y, t
    Complex (Kind=Kind(0.d0)), allocatable, dimension(:,:) :: c

    Integer :: n, j, I, iters
    
    iters = Op%N
    Mat = cmplx(0.d0, 0.d0, kind(0.D0))
    Allocate (c(iters, iters))
    c = 0.D0
    Do n = 1, iters
       Z = exp(g*Op%E(n))
       do J = 1, iters
          Z1 = Z*conjg(Op%U(J,n))
          do I = 1, iters
            y = Z1 * Op%U(I, n) - c(I, J)
            t = Mat(I, J) + y
            c(I, J) = (t - Mat(I,J)) - y
            Mat(I, J) = t
!            Mat(I, J) = Mat(I, J) + Z1 * Op%U(I, n)
          enddo
          
!          Mat(1:iters, J) = Mat(1:iters, J) + Z1 * Op%U(1:iters, n)
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
       call zcopy(Ndim, Mat(1, P(n)), 1, V(n, 1), opn)
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
!> @param[inout] V storage for the result matrix
!> @param[in] Mat Where to read those columns
!> @param[in] P A vector with which columns to copy
!> @param[in] opn The length of the vector P
!> @param[in] Ndim Mat is an Ndim x Ndim matrix
!--------------------------------------------------------------------

  subroutine copy_select_columns(V, Mat, P, opn, Ndim)
    Implicit none
    Integer, INTENT(IN) :: opn, Ndim
    Complex (Kind = Kind(0.D0)), INTENT(INOUT), Dimension(:, :) :: V
    Complex (Kind = Kind(0.D0)), Dimension(:, :), INTENT(IN) :: Mat
    Integer , Dimension(:), INTENT(IN) :: P
    Integer :: n
    
    Do n = 1, opn
        call zcopy(Ndim, Mat(P(n), 1), Ndim, V(n, 1), opn)
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
    Integer :: n,i, j, lwork, info
    Complex (Kind = Kind(0.D0)) :: alpha, beta
    Complex (Kind = Kind(0.D0)), Dimension(:,:), allocatable :: tmp
    Complex (Kind = Kind(0.D0)), Dimension(:), allocatable :: work

    alpha = 1.D0
    beta = 0.D0
    select case (opn)
    case (1)
        DO I = 1, Ndim
            Mat(P(1), I) = V(1, I)
        enddo
    case (2)
        DO I = 1, Ndim
            Mat(P(1), I) = U(1,1) * V(1, I) + U(1,2) * V(2, I)
            Mat(P(2), I) = -conjg(U(1,2)) * V(1, I) + conjg(U(1,1)) * V(2, I)
        enddo
    case default
        lwork = 2*Ndim
        Allocate(tmp(opn, Ndim), work(lwork))
        tmp = V
!        CALL ZGEMM('N','N', opn, Ndim, opn, alpha, U, opn, V, opn, beta, tmp, opn)
        DO I = 1, opn-1
        DO J = 1, Ndim
        tmp(I, J) = tmp(I, J) * U(I, I)
        ENDDO
        ENDDO
        CALL ZUNMQR('L', 'N', opn, Ndim, opn, U, opn, U(1, opn), tmp, opn, work, lwork, info)
        Mat((P), :) = tmp
        Deallocate(tmp, work)
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
    Integer :: n, i, j, lwork, info
    Complex (Kind = Kind(0.D0)) :: alpha, beta
    Complex (Kind = Kind(0.D0)), Dimension(:,:), allocatable :: tmp
    Complex (Kind = Kind(0.D0)), Dimension(:), allocatable :: work

    alpha = 1.D0
    beta = 0.D0
    select case (opn)
    case (1)
        DO I = 1, Ndim
            Mat(I, P(1)) = V(1, I)
        enddo
    case (2)
        DO I = 1, Ndim
            Mat(I, P(1)) = conjg(U(1,1)) * V(1, I) + conjg(U(1,2)) * V(2, I)
            Mat(I, P(2)) = -U(1,2) * V(1, I) + U(1,1) * V(2, I)
        enddo
    case default
        lwork = 2*Ndim
        Allocate(tmp(Ndim, opn), work(lwork))
        DO I = 1, Ndim
        DO J = 1, opn - 1
        tmp(I, J) = V(J, I) * Conjg(U(J, J))
        ENDDO
        tmp(I, opn) = V(opn, I) 
        ENDDO
        CALL ZUNMQR('R', 'C', Ndim, opn, opn, U, opn, U(1, opn), tmp, Ndim, work, lwork, info)
!        CALL ZGEMM('T','C', Ndim, opn, opn, alpha, V, opn, U, opn, beta, tmp, Ndim)
        Mat(:, (P)) = tmp
        Deallocate(tmp, work)
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
!> @param[in] P a vector with the rows that we write to
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
    Integer :: n, i, j, lwork, info
    Complex (Kind = Kind(0.D0)) :: beta
    Complex (Kind = Kind(0.D0)), Dimension(:, :), allocatable :: tmp
    Complex (Kind = Kind(0.D0)), Dimension(:), allocatable :: work

    beta = 0.D0
    select case (opn)
    case (1)
        DO I = 1, Ndim
            Mat(I, P(1)) = Z(1) * V(1, I)
        enddo
    case (2)
        DO I = 1, Ndim
            Mat(I, P(1)) = Z(1) * (U(1, 1) * V(1, I) - conjg(U(1, 2)) * V(2, I))
            Mat(I, P(2)) = Z(2) * (U(1, 2) * V(1, I) + conjg(U(1, 1)) * V(2, I))
        enddo
    case default
        lwork = 2 * Ndim
        allocate(tmp(Ndim, opn), work(lwork))
        tmp = Transpose(V)
        ! multiply with Q
        CALL ZUNMQR('R', 'N', Ndim ,opn, opn, U, opn, U(1, opn), tmp, Ndim, work, lwork, info)
        ! multiply with R
        DO I = 1, Ndim
            DO J = 1, opn - 1
                tmp(I, J) = tmp(I, J) * U(J, J)
            ENDDO
        ENDDO
        DO I = 1, Ndim
        DO J = 1, opn
        tmp(I, J) = Z(J) * tmp(I, J)
        ENDDO
        ENDDO
        Mat(:, P) = tmp
!         do n = 1, opn
!             call zgemv('T', opn, Ndim, Z(n), V, opn, U(:, n), 1, beta, Mat(:, P(n)), 1)
!         Enddo
        deallocate(tmp, work)
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
!> @param[in] P a vector with the rows that we write to
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
    Integer :: n, i, j, lwork, info
    Complex (Kind = Kind(0.D0)) :: beta
    Complex (Kind = Kind(0.D0)), Dimension(:, :), allocatable :: tmp
    Complex (Kind = Kind(0.D0)), Dimension(:), allocatable :: work

    beta = 0.D0
    select case (opn)
    case (1)
        DO I = 1, Ndim
            Mat(P(1), I) = Z(1) * V(1, I)
        enddo
    case (2)
        DO I = 1, Ndim
            Mat(P(1), I) = Z(1) * (conjg(U(1, 1)) * V(1, I) - U(1, 2) * V(2, I))
            Mat(P(2), I) = Z(2) * (conjg(U(1, 2)) * V(1, I) + U(1, 1) * V(2, I))
        enddo
    case default
        lwork = 2 * Ndim
        allocate(tmp(opn, Ndim), work(lwork))
        tmp = V
        ! multiply with Q
        CALL ZUNMQR('L', 'C', opn, Ndim, opn, U, opn, U(1, opn), tmp, opn, work, lwork, info)
        ! multiply with R
        DO I = 1, opn -1
            DO J = 1, Ndim
                tmp(I, J) = tmp(I, J) * U(I, I)
            ENDDO
        ENDDO
        ! exponentials
        DO I = 1, opn
        Do J = 1, Ndim
        tmp(I, J) = tmp(I, J) * Z(I)
        ENDDO
        ENDDO
        Mat(P, :) = tmp
!         do n = 1, opn
!             call zgemv('T', opn, Ndim, Z(n), V, opn, conjg(U(:, n)), 1, beta, Mat(P(n), 1), size(Mat, 1))
!         Enddo
        deallocate(tmp, work)
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
    Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind=Kind(0.d0)), INTENT(IN)   :: spin

    ! Local 
    Complex (Kind=Kind(0.d0)), Dimension(:, :), allocatable :: VH
    Complex (Kind=Kind(0.d0)), Dimension(:), allocatable :: Z

    ! In  Mat
    ! Out Mat = Mat*exp(spin*Op)
    allocate(VH(Op%N,Ndim), Z(Op%N))
    call copy_select_rows(VH, Mat, Op%P, Op%N, Ndim)
    Z = exp(Op%g * spin * Op%E)
    call opexpmult(VH, Op%U, Op%P, Mat, Z, Op%N, Ndim)
    call copy_select_rows(VH, Mat, Op%P, Op%N, Ndim)
    call opmultct(VH, Op%U, Op%P, Mat, Op%N, Ndim)
    deallocate(VH, Z)
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
    Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind=Kind(0.d0)), INTENT(IN )   :: spin

    ! Local 
    Complex (Kind=Kind(0.d0)), Dimension(:, :), allocatable :: VH
    Complex (Kind=Kind(0.d0)), Dimension(:), allocatable :: Z
    
    ! In  Mat
    ! Out Mat = exp(spin*Op)*Mat
    allocate(VH(Op%N,Ndim), Z(Op%N))
    call copy_select_columns(VH, Mat, Op%P, Op%N, Ndim)
    Z = exp(Op%g * spin * Op%E)
    call opexpmultct(VH, Op%U, Op%P, Mat, Z, Op%N, Ndim)    
    call copy_select_columns(VH, Mat, Op%P, Op%N, Ndim)
    call opmult(VH, Op%U, Op%P, Mat, Op%N, Ndim)
    deallocate(VH, Z)
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
  Pure subroutine FillExpOps(ExpOp, ExpMop, Op, spin)
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
    
  end subroutine FillExpOps

!--------------------------------------------------------------------

  Subroutine Op_Wrapup(Mat,Op,spin,Ndim,N_Type)

    Implicit none 

    Integer :: Ndim
    Type (Operator) , INTENT(IN )   :: Op
    Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind=Kind(0.d0)), INTENT(IN )   :: spin
    Integer, INTENT(IN) :: N_Type

    ! Local 
    Complex (Kind=Kind(0.d0)) :: ExpOp(Op%N), ExpMop(Op%N), VH(Op%N,Ndim)
    
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

!--------------------------------------------------------------------

  Subroutine Op_Wrapdo(Mat,Op,spin,Ndim,N_Type)
    Implicit none 
    
    Integer :: Ndim
    Type (Operator) , INTENT(IN )   :: Op
    Complex (Kind = Kind(0.D0)), INTENT(INOUT) :: Mat (Ndim,Ndim)
    Real    (Kind = Kind(0.D0)), INTENT(IN )   :: spin
    Integer, INTENT(IN) :: N_Type

    ! Local 
    Complex (Kind = Kind(0.D0)), Dimension(:), allocatable :: ExpOp, ExpMop, work
    Integer :: n, i, j, lwork, info
    Complex (Kind = Kind(0.D0)) :: alpha, beta, ExpHere
    Complex (Kind = Kind(0.D0)), Dimension(:, :), allocatable :: VH, tmp, tmp2

    alpha = 1.D0
    beta  = 0.D0
    !!!!! N_Type == 1
    !    Op%U*exp(-Op%g*spin*Op%E)*Mat*exp(Op%g*spin*Op%E)*(Op%U^{dagger})
    !    
    !!!!!
    !!!!! N_Type == 2
    !    (Op%U^{dagger}) * Mat * Op%U
    !!!!!
    
    Allocate(VH(Op%N,Ndim), ExpOp(Op%N), ExpMop(Op%N))
    If (N_type == 1) then
       call FillExpOps(ExpOp, ExpMop, Op, spin)
       
       if (Op%N == 1) then
            DO I = 1, Ndim
                Mat(I, Op%P(1)) = ExpOp(1) * Mat(I, Op%P(1))
            enddo
            DO I = 1, Ndim
                Mat(Op%P(1), I) = ExpMOp(1) * Mat(Op%P(1), I)
            enddo
       else
            Do n = 1,Op%N
                expHere=ExpOp(n)
                VH(n, :) = ExpHere * Mat(:, Op%P(n))
            Enddo
            write(*, *) "opmultct"
            call opmultct(VH, Op%U, Op%P, Mat, Op%N, Ndim)
            Do n = 1,Op%N
                ExpHere=ExpMOp(n)
                VH(n, :) = ExpHere * Mat(Op%P(n), :)
            Enddo
            write (*,*) "opmult"
            call opmult(VH, Op%U, Op%P, Mat, Op%N, Ndim)
       endif
    elseif (N_Type == 2) then
        if(Op%N > 1) then
            select case (Op%N)
            case (2)
                call copy_select_rows(VH, Mat, Op%P, 2, Ndim)
                DO I = 1, Ndim
                    Mat(I, Op%P(1)) = Op%U(1, 1) * VH(1, I) + Op%U(2, 1) * VH(2, I)
                    Mat(I, Op%P(2)) = -conjg(Op%U(2, 1)) * VH(1, I) + conjg(Op%U(1, 1)) * VH(2, I)
                enddo
                call copy_select_columns(VH, Mat, Op%P, 2, Ndim)
                DO I = 1, Ndim
                    Mat(Op%P(1), I) = conjg(Op%U(1,1)) * VH(1, I) + conjg(Op%U(2,1)) * VH(2, I)
                    Mat(Op%P(2), I) = -Op%U(2,1) * VH(1, I) + Op%U(1,1) * VH(2, I)
                enddo
            case default
                call copy_select_rows(VH, Mat, Op%P, Op%N, Ndim)
                lwork = 2 * Ndim
                Allocate(tmp(Ndim, Op%N), work(lwork))
                tmp = TRANSPOSE(VH)
                CALL ZUNMQR('R', 'N', Ndim, Op%N, Op%N, Op%U, op%n, Op%U(1, op%N), tmp, Ndim, work, lwork, info)
!                CALL ZGEMM('T','N', Ndim, op%N, op%N, alpha, VH, op%n, Op%U, op%n, beta, tmp, Ndim)
                DO I = 1, Ndim
                    DO J = 1, Op%N -1
                        tmp(I, J) = tmp(I, J) * Op%U(J, J)
                    ENDDO
                ENDDO
                Mat(:, (Op%P)) = tmp
                Deallocate(tmp)
                call copy_select_columns(VH, Mat, Op%P, Op%N, Ndim)
                Allocate(tmp2(Op%N, Ndim))
                tmp2 = VH
                CALL ZUNMQR('L', 'C', Op%N, Ndim, op%N, Op%U, Op%n, Op%U(1, op%N), tmp2, Op%N, work, lwork, info)
                DO I = 1, Op%N -1
                    DO J = 1, Ndim
                        tmp2(I, J) = tmp2(I, J) * conjg(Op%U(I, I))
                    ENDDO
                ENDDO
!                CALL ZGEMM('C','N', op%N, Ndim, op%N, alpha, Op%U, op%n, VH, op%n, beta, tmp2, op%n)
                Mat(Op%P, :) = tmp2
                Deallocate(tmp2)
            end select
        endif
    endif
  end Subroutine Op_Wrapdo

end Module Operator_mod
