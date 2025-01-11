!  Copyright (C) 2016 - 2023 The ALF project
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
module Operator_mod

   use mpi_shared_memory
   use mat_subroutines
   use MyMats
   use Fields_mod

   implicit none

   type operator
      integer          :: N, N_non_zero                      !> dimension of Operator (P and O), number of non-zero eigenvalues
      integer, private :: win_M_exp, win_U                   !> MPI_windows which can be used for fences (memory synch.) and dealloc.
      logical          :: diag                               !> encodes if Operator is diagonal
      logical, private :: U_alloc, M_exp_alloc, g_t_alloc    !> logical to track if memory is allocated
      complex(Kind=kind(0.d0)), pointer :: O(:, :), U(:, :)  !> Storage for operator matrix O and it's eigenvectors U
      complex(Kind=kind(0.d0)), pointer, private :: M_exp(:, :, :), E_exp(:, :)  !>internal storage for exp(O) and exp(E)
      real(Kind=kind(0.d0)), pointer :: E(:)             !> Eigenvalues of O
      integer, pointer :: P(:)                               !> Projector P encoding DoFs that contribute in Operator
      complex(Kind=kind(0.d0)) :: g                         !> coupling constant
      complex(Kind=kind(0.d0)), allocatable :: g_t(:)       !> time dependent  coupling constant
      complex(Kind=kind(0.d0)) :: alpha                     !> operator shift
      integer          :: type                               !> Type of the operator: 1=Ising; 2=discrete HS; 3=continuous scalar HS, 4=  Complex for  three  bondy term.
      integer          :: Flip_protocol = 1                    !> Flip protocol  for  local  updates.  Only  relevant  for  type =3  fields.
      ! P is an N X Ndim matrix such that  P.T*O*P*  =  A
      ! P has only one non-zero entry per column which is specified by P
      ! All in all.   g * Phi(s,type) * ( c^{dagger} A c  + alpha )
      ! The variable Type allows you to define the type of HS.
      ! The first N_non_zero elements of diagonal matrix E are non-zero. The rest vanish.

      ! !!!!! M_exp and E_exp  are for storage   !!!!!
      ! If Type =1   then the Ising field  takes the values  s = +/- 1
      !              and M_exp   has dimensions  M_exp(N,N,3)   last index=1 rep field -1 and index=3 rep field 1
      ! If Type =2   then the Ising field  takes the values  s = +/- 1,  +/- 2
      !              and M_exp   has dimensions  M_exp(N,N,5)
      ! M_exp(:,:,s) =  e^{g * Phi(s,type) *  O(:,:) }
      !
      !
      ! E_exp(:,s) = e^{g * Phi(s,type) *E(:) } and has dimensions E_exp(N,-1:1)  for Type = 1
      !              and dimensions E_exp(N,-2:2)   for Type = 2
      !
      ! !!! If Type .neq. 1,2  then  E_exp  and  M_exp  are  not allocated !!!
   contains
      procedure  :: get_g_t_alloc => operator_get_g_t_alloc     !>    External access to   private  logical  variable g_t_alloc
   end type operator

contains

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
!> @param[in] nf
!> Integer
!> * flavor index
!--------------------------------------------------------------------
   subroutine Op_phase(Phase, OP_V, Nsigma, nf)
      implicit none

      complex(Kind=kind(0.d0)), intent(Inout) :: Phase
      integer, intent(IN)    :: nf
      type(Fields), intent(IN)    :: Nsigma
      type(operator), dimension(:, :), intent(In) :: Op_V
      real(Kind=kind(0.d0))                       :: angle

      integer :: n, nt
      complex(kind=kind(0.d0)) :: g_loc

      do n = 1, size(Op_V, 1)
         do nt = 1, size(nsigma%f, 2)
            g_loc = Op_V(n, nf)%g
            if (op_v(n, nf)%g_t_alloc) g_loc = Op_V(n, nf)%g_t(nt)
            angle = aimag(g_loc*Op_V(n, nf)%alpha*nsigma%Phi(n, nt))
            Phase = Phase*cmplx(cos(angle), sin(angle), kind(0.d0))
         end do
      end do

   end subroutine Op_phase

!--------------------------------------------------------------------
!> @author
!>
!> @brief
!> Set up _core_ data of the operator required for type=0 operators
!
!> @param[inout] Op
!> @param[in] N
!--------------------------------------------------------------------

   pure subroutine Op_make(Op, N)
      implicit none
      type(operator), intent(INOUT) :: Op
      integer, intent(IN) :: N
      allocate (Op%O(N, N), Op%P(N))
      ! F.F.A  Op%M_exp and Op%E_exp are allocated  in Op_set once the type is available.

      Op%O = cmplx(0.d0, 0.d0, kind(0.d0))
      Op%P = 0
      Op%N = N
      Op%N_non_zero = N
      Op%g = cmplx(0.d0, 0.d0, kind(0.d0))
      Op%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
      Op%diag = .false.
      Op%type = 0
      OP%flip_protocol = 1
      Op%U_alloc = .false.
      Op%M_exp_alloc = .false.
      Op%g_t_alloc = .false.
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

   subroutine Op_clear(Op, N)
      implicit none
      type(operator), intent(INOUT) :: Op
      integer, intent(IN) :: N

!     If ( associated(OP%O) )   deallocate(OP%O)
!     If ( associated(OP%P) )   deallocate(OP%P)
!     If ( associated(OP%E_exp) )   deallocate(OP%E_exp)
!     If ( associated(OP%E) )   deallocate(OP%E)
      if (use_mpi_shm) then
         if (Op%U_alloc) then
            deallocate (OP%O, OP%P, OP%E)
            !call deallocate_shared_memory(OP%win_U)
         end if
         if (Op%M_exp_alloc) then
            deallocate (OP%E_exp)
            !call deallocate_shared_memory(OP%win_M_exp)
         end if
      else
         if (Op%M_exp_alloc) deallocate (OP%M_exp, OP%E_exp)
         if (Op%U_alloc) deallocate (OP%U, OP%O, OP%P, OP%E)
         if (Op%g_t_alloc) deallocate (Op%g_t)
      end if

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
      implicit none
      type(operator), intent(INOUT) :: Op

      complex(Kind=kind(0.d0)), allocatable :: U(:, :), TMP(:, :)
      real(Kind=kind(0.d0)), allocatable :: E(:)
      real(Kind=kind(0.d0)) :: Zero = 1.d-9 !, Phi(-2:2)
      integer :: N, I, J, np, nz, noderank, arrayshape2d(2), arrayshape(3), ierr
      complex(Kind=kind(0.d0)) :: Z
      type(Fields)   :: nsigma_single

      if (allocated(OP%g_t)) Op%g_t_alloc = .true.

      call nsigma_single%make(1, 1)
      noderank = 0

      N = OP%N
      allocate (Op%E(N))
      if (use_mpi_shm) then
         arrayshape2d = (/Op%N, Op%N/)
         call allocate_shared_memory(Op%U, Op%win_U, noderank, arrayshape2d)
         if (noderank == 0) Op%U = cmplx(0.d0, 0.d0, kind(0.d0))
      else
         allocate (Op%U(N, N))
         Op%U = cmplx(0.d0, 0.d0, kind(0.d0))
      end if
      Op%U_alloc = .true.
      Op%E = 0.d0

      if (Op%N > 1) then
         N = Op%N
         Op%diag = .true.
         do I = 1, N
            do J = i + 1, N
               ! Binary comparison is OK here as Op%O was initialized to zero during Op_make.
             if (Op%O(i, j) .ne. cmplx(0.d0, 0.d0, kind(0.d0)) .or. Op%O(j, i) .ne. cmplx(0.d0, 0.d0, kind(0.d0))) Op%diag = .false.
            end do
         end do
         if (Op%diag) then
            do I = 1, N
               Op%E(I) = dble(Op%O(I, I))
               if (noderank == 0) Op%U(I, I) = cmplx(1.d0, 0.d0, kind(0.d0))
            end do
            Op%N_non_zero = N
            ! FFA Why do we assume that Op%N_non_zero = N for a diagonal operator?
         else
            allocate (U(N, N), E(N), TMP(N, N))
            call Diag(Op%O, U, E)
            Np = 0
            Nz = 0
            do I = 1, N
               if (abs(E(I)) > Zero) then
                  np = np + 1
                  if (noderank == 0) Op%U(:, np) = U(:, i)
                  Op%E(np) = E(I)
               else
                  if (noderank == 0) Op%U(:, N - nz) = U(:, i)
                  Op%E(N - nz) = E(I)
                  nz = nz + 1
               end if
            end do
            Op%N_non_zero = np
            ! Write(6,*) "Op_set", np,N
            if (noderank == 0) then
               TMP = Op%U ! that way we have the changes to the determinant due to the permutation
               Z = Det_C(TMP, N)
               ! Scale Op%U to be in SU(N)
               do I = 1, N
                  Op%U(I, 1) = Op%U(I, 1)/Z
               end do
            end if
            deallocate (U, E, TMP)
            ! Op%U,Op%E)
            ! Write(6,*) 'Calling diag 1'
         end if
      else
         Op%E(1) = real(Op%O(1, 1), kind(0.d0))
         if (noderank == 0) Op%U(1, 1) = cmplx(1.d0, 0.d0, kind(0.d0))
         Op%N_non_zero = 1
         Op%diag = .true.
      end if
#ifdef MPI
      if (use_mpi_shm) call MPI_WIN_FENCE(0, Op%win_U, ierr)
#endif
      select case (OP%type)
      case (1)
         if (use_mpi_shm) then
            allocate (Op%E_exp(Op%N, -Op%type:Op%type))
            arrayshape = (/Op%N, Op%N, 3/)
            call allocate_shared_memory(Op%M_exp, Op%win_M_exp, noderank, arrayshape)
         else
            allocate (Op%E_exp(Op%N, -Op%type:Op%type), Op%M_exp(Op%N, Op%N, 3))
         end if
         Op%M_exp_alloc = .true.
         nsigma_single%t(1) = 1
         do I = 1, Op%type
            nsigma_single%f(1, 1) = real(I, kind=kind(0.d0))
            do n = 1, Op%N
               Op%E_exp(n, I) = cmplx(1.d0, 0.d0, kind(0.d0))
               Op%E_exp(n, -I) = cmplx(1.d0, 0.d0, kind(0.d0))
               if (n <= Op%N_non_Zero) then
                  !Op%E_exp(n,I) = exp(Op%g*Op%E(n)*Phi_st(I,1))
                  Op%E_exp(n, I) = exp(Op%g*Op%E(n)*nsigma_single%Phi(1, 1))
                  Op%E_exp(n, -I) = 1.d0/Op%E_exp(n, I)
               end if
            end do
            !call Op_exp(Op%g*Phi_st( I,1),Op,Op%M_exp(:,:,I))
            !call Op_exp(Op%g*Phi_st(-I,1),Op,Op%M_exp(:,:,-I))
            if (noderank == 0) then
               call Op_exp(Op%g*nsigma_single%Phi(1, 1), Op, Op%M_exp(:, :, I + 2))
               call Op_exp(-Op%g*nsigma_single%Phi(1, 1), Op, Op%M_exp(:, :, -I + 2))
            end if
#ifdef MPI
            if (use_mpi_shm) call MPI_WIN_FENCE(0, Op%win_M_exp, ierr)
#endif
         end do
      case (2)
         if (use_mpi_shm) then
            allocate (Op%E_exp(Op%N, -Op%type:Op%type))
            arrayshape = (/Op%N, Op%N, 5/)
            call allocate_shared_memory(Op%M_exp, Op%win_M_exp, noderank, arrayshape)
         else
            allocate (Op%E_exp(Op%N, -Op%type:Op%type), Op%M_exp(Op%N, Op%N, 5))
         end if
         Op%M_exp_alloc = .true.
         nsigma_single%t(1) = 2
         do I = 1, Op%type  ! Here Op%type = 2
            nsigma_single%f(1, 1) = real(I, kind=kind(0.d0))
            do n = 1, Op%N
               Op%E_exp(n, I) = cmplx(1.d0, 0.d0, kind(0.d0))
               Op%E_exp(n, -I) = cmplx(1.d0, 0.d0, kind(0.d0))
               if (n <= Op%N_non_Zero) then
                  !Op%E_exp(n,I)  = exp(Op%g*Op%E(n)*Phi_st(I,2))
                  Op%E_exp(n, I) = exp(Op%g*Op%E(n)*nsigma_single%Phi(1, 1))
                  Op%E_exp(n, -I) = 1.d0/Op%E_exp(n, I)
               end if
            end do
            if (noderank == 0) then
               call Op_exp(Op%g*nsigma_single%Phi(1, 1), Op, Op%M_exp(:, :, I + 3))
               call Op_exp(-Op%g*nsigma_single%Phi(1, 1), Op, Op%M_exp(:, :, -I + 3))
            end if
#ifdef MPI
            if (use_mpi_shm) call MPI_WIN_FENCE(0, Op%win_M_exp, ierr)
#endif
            !call Op_exp(Op%g*Phi_st( I,2),Op,Op%M_exp(:,:,I))
            !call Op_exp(Op%g*Phi_st(-I,2),Op,Op%M_exp(:,:,-I))
         end do
      case default
      end select

      call nsigma_single%clear()

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
   subroutine Op_exp(g, Op, Mat)

      implicit none

      type(operator), intent(IN)  :: Op
      complex(Kind=kind(0.d0)), dimension(:, :), intent(OUT) :: Mat
      complex(Kind=kind(0.d0)), intent(IN) :: g
      complex(Kind=kind(0.d0)) :: Z, Z1, y, t
      complex(Kind=kind(0.d0)), allocatable, dimension(:, :) :: c

      integer :: n, j, I, iters

      iters = Op%N
      Mat = cmplx(0.d0, 0.d0, kind(0.d0))
      if (Op%diag) then
         do n = 1, iters
            Mat(n, n) = exp(g*Op%E(n))
         end do
      else
         allocate (c(iters, iters))
         c = 0.d0
         do n = 1, iters
            Z = exp(g*Op%E(n))
            do J = 1, iters
               Z1 = Z*conjg(Op%U(J, n))
               do I = 1, iters
                  ! This performs Kahan summation so as to improve precision.
                  y = Z1*Op%U(I, n) - c(I, J)
                  t = Mat(I, J) + y
                  c(I, J) = (t - Mat(I, J)) - y
                  Mat(I, J) = t
                  !  Mat(I, J) = Mat(I, J) + Z1 * Op%U(I, n)
               end do
               ! Mat(1:iters, J) = Mat(1:iters, J) + Z1 * Op%U(1:iters, n)
            end do
         end do
         deallocate (c)
      end if
   end subroutine Op_exp

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!
!> @brief
!> Out   For  nsigma_single%f(1,1) =  HS_Field   the  routine  computes
!>       Op%Mat = Mat* Op ( exp( sign * nsigma_single%phi(1,1)*g* P^T O T) )
!>       For  Op%type = 1,2  the  exponential is  stored. Otherwise  it is
!>       computed  on the fly
!>
!> @param[inout] Mat Complex Dimension(:,:)
!> * On exit Mat = Mat*Op ( exp(nsigma_single%phi(1,1)*g* P^T O T) )
!> @param[in] Op Type(Operator)
!> * The Operator containing g and the sparse matrix P^T O P
!> @param[in]  HS_Field Complex
!> * The field
!> @param[in] cop  Character
!> * cop = N,  Op = None
!> * cop = T,  Op = Transposed
!> * cop = C,  Op = Transposed + Complex conjugation
!>

!--------------------------------------------------------------------
   subroutine Op_mmultL(Mat, Op, HS_Field, cop, nt, sign)
      implicit none
      type(operator), intent(IN)    :: Op
      complex(Kind=kind(0.d0)), intent(INOUT) :: Mat(:, :)
      integer, intent(IN)    :: sign
      complex(Kind=kind(0.d0)), intent(IN)    :: Hs_Field
      character, intent(IN)    :: cop
      integer, intent(in)    :: nt

      ! Local
      integer :: I, N1, N2, sp
      complex(Kind=kind(0.d0)) :: ExpMat(Op%n, Op%n), g_loc
      type(Fields)            :: nsigma_single

      call nsigma_single%make(1, 1)
      nsigma_single%f(1, 1) = HS_field
      nsigma_single%t(1) = op%type

      N1 = size(Mat, 1)
      N2 = size(Mat, 2)

      g_loc = OP%g
      if (op%g_t_alloc) g_loc = Op%g_t(nt)
      if (abs(g_loc) < 1.d-12) return

      if (op%type < 3) then
         sp = sign*nint(real(HS_Field))
         if (Op%diag) then
            if (Op%g_t_alloc) then
               do I = 1, Op%N
                  if (cop == 'c' .or. cop == 'C') then
                     call ZSCAL(N1, conjg(exp(sign*nsigma_single%phi(1, 1)*Op%g_t(nt)*Op%E(I))), Mat(1, Op%P(I)), 1)
                  else
                     call ZSCAL(N1, exp(sign*nsigma_single%phi(1, 1)*Op%g_t(nt)*Op%E(I)), Mat(1, Op%P(I)), 1)
                  end if
               end do
            else
               do I = 1, Op%N
                  if (cop == 'c' .or. cop == 'C') then
                     call ZSCAL(N1, conjg(Op%E_exp(I, sp)), Mat(1, Op%P(I)), 1)
                  else
                     call ZSCAL(N1, Op%E_exp(I, sp), Mat(1, Op%P(I)), 1)
                  end if
               end do
            end if
         else
            if (Op%g_t_alloc) then
               call Op_exp(sign*nsigma_single%phi(1, 1)*Op%g_t(nt), Op, expmat)
               call ZSLGEMM('r', cop, Op%N, N1, N2, expmat, Op%P, Mat)
            else
               call ZSLGEMM('r', cop, Op%N, N1, N2, Op%M_exp(:, :, sp + op%type + 1), Op%P, Mat)
            end if
         end if
      else
         if (Op%diag) then
            do I = 1, Op%N
               if (cop == 'c' .or. cop == 'C') then
                  call ZSCAL(N1, conjg(exp(sign*nsigma_single%phi(1, 1)*g_loc*Op%E(I))), Mat(1, Op%P(I)), 1)
               else
                  call ZSCAL(N1, exp(sign*nsigma_single%phi(1, 1)*g_loc*Op%E(I)), Mat(1, Op%P(I)), 1)
               end if
            end do
         else
            call Op_exp(sign*g_loc*nsigma_single%phi(1, 1), Op, expmat)
            call ZSLGEMM('r', cop, Op%N, N1, N2, expmat, Op%P, Mat)
         end if
      end if

   end subroutine Op_mmultL

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief
!> Out   For  nsigma_single%f(1,1) =  HS_Field   the  routine  computes
!>       Op%Mat =  Op ( exp( nsigma_single%phi(1,1)*g* P^T O T) ) * Mat
!>       For  Op%type = 1,2  the  exponential is  stored. Otherwise  it is
!>       computed  on the fly
!
!> @param[inout] Mat Complex Dimension(:,:)
!> * On exit Mat = Op ( exp(nsigma_single%phi(1,1)*g* P^T O T) )* Mat
!> @param[in] Op Type(Operator)
!> * The Operator containing g and the sparse matrix P^T O P
!> @param[in] Hs_Field Complex
!> * The field
!> @param[in] cop  Character
!> * cop = N,  Op = None
!> * cop = T,  Op = Transposed
!> * cop = C,  Op = Transposed + Complex conjugation
!>
!--------------------------------------------------------------------
   subroutine Op_mmultR(Mat, Op, HS_Field, cop, nt)
      implicit none
      type(operator), intent(IN)   :: Op
      complex(Kind=kind(0.d0)), intent(INOUT) :: Mat(:, :)
      complex(Kind=kind(0.d0)), intent(IN)   :: Hs_Field
      character, intent(IN)    :: cop
      integer, intent(in)    :: nt

      ! Local
      integer :: I, N1, N2, sp
      complex(Kind=kind(0.d0)) :: ExpMat(Op%n, Op%n), g_loc
      type(Fields)   :: nsigma_single

      call nsigma_single%make(1, 1)
      nsigma_single%f(1, 1) = Hs_Field
      nsigma_single%t(1) = op%type

      N1 = size(Mat, 1)
      N2 = size(Mat, 2)

      ! quick return if possible
      g_loc = Op%g
      if (op%g_t_alloc) g_loc = Op%g_t(nt)
      if (abs(g_loc) < 1.d-12) return

      if (op%type < 3) then
         sp = nint(real(Hs_Field))
         if (Op%diag) then
            if (op%g_t_alloc) then
               do I = 1, Op%N
                  if (cop == 'c' .or. cop == 'C') then
                     call ZSCAL(N2, conjg(exp(nsigma_single%phi(1, 1)*Op%g_t(nt)*Op%E(I))), Mat(Op%P(I), 1), N1)
                  else
                     call ZSCAL(N2, exp(nsigma_single%phi(1, 1)*Op%g_t(nt)*Op%E(I)), Mat(Op%P(I), 1), N1)
                  end if
               end do
            else
               do I = 1, Op%N
                  if (cop == 'c' .or. cop == 'C') then
                     call ZSCAL(N2, conjg(Op%E_exp(I, sp)), Mat(Op%P(I), 1), N1)
                  else
                     call ZSCAL(N2, Op%E_exp(I, sp), Mat(Op%P(I), 1), N1)
                  end if

               end do
            end if
         else
            if (op%g_t_alloc) then
               call Op_exp(Op%g_t(nt)*nsigma_single%phi(1, 1), Op, expmat)
               call ZSLGEMM('L', cop, Op%N, N1, N2, expmat, Op%P, Mat)
            else
               call ZSLGEMM('L', cop, Op%N, N1, N2, Op%M_exp(:, :, sp + op%type + 1), Op%P, Mat)
            end if
         end if
      else
         if (Op%diag) then
            do I = 1, Op%N
               if (cop == 'c' .or. cop == 'C') then
                  call ZSCAL(N2, conjg(exp(nsigma_single%phi(1, 1)*g_loc*Op%E(I))), Mat(Op%P(I), 1), N1)
               else
                  call ZSCAL(N2, exp(nsigma_single%phi(1, 1)*g_loc*Op%E(I)), Mat(Op%P(I), 1), N1)
               end if
            end do
         else
            call Op_exp(g_loc*nsigma_single%phi(1, 1), Op, expmat)
            call ZSLGEMM('L', cop, Op%N, N1, N2, expmat, Op%P, Mat)
         end if
      end if
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
!>  N_type = 1, Mat = exp(Op%g*nsigma%phi(1,1)*Op%E)*(Op%U^{dagger}) * Mat * Op%U*exp(-Op%g*nsigma%phi(1,1)*Op%E)
!>  N_type = 2, Mat = Op%U * Mat * (Op%U^{dagger})
!> \endverbatim
!>  **Do not** mix up  \p N_type  and    \p OP\%type
!> @param[in] Op Type(Operator)
!> \verbatim
!>  The operator containing g, U, P
!> \endverbatim
!> @param[in] HS_Field Complex
!> \verbatim
!>  The field
!> \endverbatim
!> @param[in] Ndim Integer
!>
!--------------------------------------------------------------------
   subroutine Op_Wrapup(Mat, Op, HS_Field, Ndim, N_Type, nt)

      implicit none

      integer :: Ndim
      type(operator), intent(IN)   :: Op
      complex(Kind=kind(0.d0)), intent(INOUT) :: Mat(Ndim, Ndim)
      complex(Kind=kind(0.d0)), intent(IN)   :: HS_Field
      integer, intent(IN)    :: N_Type
      integer, intent(IN)    :: nt

      ! Local
      complex(Kind=kind(0.d0)) :: VH1(Op%N, Op%N)
      integer :: I, sp
      complex(kind=kind(0.d0)) :: g_loc
      type(Fields)   :: nsigma_single

      call nsigma_single%make(1, 1)
      nsigma_single%f(1, 1) = HS_Field
      nsigma_single%t(1) = op%type

      if (op%type < 3) then
         sp = nint(real(HS_Field))
         if (N_type == 1) then
            if (Op%diag) then
               if (op%g_t_alloc) then
                  do I = 1, Op%N
                     call ZSCAL(Ndim, exp(nsigma_single%phi(1, 1)*Op%g_t(nt)*Op%E(I)), Mat(Op%P(I), 1), Ndim)
                  end do
                  do I = 1, Op%N
                     call ZSCAL(Ndim, exp(-nsigma_single%phi(1, 1)*Op%g_t(nt)*Op%E(I)), Mat(1, Op%P(I)), 1)
                  end do
               else
                  do I = 1, Op%N
                     call ZSCAL(Ndim, Op%E_Exp(I, sp), Mat(Op%P(I), 1), Ndim)
                  end do
                  do I = 1, Op%N
                     call ZSCAL(Ndim, Op%E_Exp(I, -sp), Mat(1, Op%P(I)), 1)
                  end do
               end if
            else
               if (op%g_t_alloc) then
                  do i = 1, Op%N
                     VH1(:, i) = Op%U(:, i)*exp(-nsigma_single%phi(1, 1)*Op%g_t(nt)*Op%E(I))
                  end do
                  call ZSLGEMM('r', 'n', Op%n, Ndim, Ndim, VH1, Op%P, Mat)
                  do i = 1, Op%N
                     VH1(:, i) = exp(nsigma_single%phi(1, 1)*Op%g_t(nt)*Op%E(I))*conjg(Op%U(:, i))
                  end do
                  call ZSLGEMM('l', 'T', Op%n, Ndim, Ndim, VH1, Op%P, Mat)

               else
                  do i = 1, Op%N
                     VH1(:, i) = Op%U(:, i)*Op%E_Exp(I, -sp)
                  end do
                  call ZSLGEMM('r', 'n', Op%n, Ndim, Ndim, VH1, Op%P, Mat)
                  do i = 1, Op%N
                     VH1(:, i) = Op%E_Exp(I, sp)*conjg(Op%U(:, i))
                  end do
                  call ZSLGEMM('l', 'T', Op%n, Ndim, Ndim, VH1, Op%P, Mat)
               end if
            end if
         elseif (N_Type == 2 .and. .not. Op%diag) then
            call ZSLGEMM('l', 'n', Op%n, Ndim, Ndim, Op%U, Op%P, Mat)
            call ZSLGEMM('r', 'c', Op%n, Ndim, Ndim, Op%U, Op%P, Mat)
         end if
      else
         if (N_type == 1) then
            g_loc = Op%g
            if (op%g_t_alloc) g_loc = Op%g_t(nt)
            if (Op%diag) then
               do I = 1, Op%N
                  call ZSCAL(Ndim, exp(nsigma_single%phi(1, 1)*g_loc*Op%E(I)), Mat(Op%P(I), 1), Ndim)
               end do
               do I = 1, Op%N
                  call ZSCAL(Ndim, exp(-nsigma_single%phi(1, 1)*g_loc*Op%E(I)), Mat(1, Op%P(I)), 1)
               end do
            else
               do i = 1, Op%N
                  VH1(:, i) = Op%U(:, i)*exp(-nsigma_single%phi(1, 1)*g_loc*Op%E(I))
               end do
               call ZSLGEMM('r', 'n', Op%n, Ndim, Ndim, VH1, Op%P, Mat)
               do i = 1, Op%N
                  VH1(:, i) = exp(nsigma_single%phi(1, 1)*g_loc*Op%E(I))*conjg(Op%U(:, i))
               end do
               call ZSLGEMM('l', 'T', Op%n, Ndim, Ndim, VH1, Op%P, Mat)
            end if
         elseif (N_Type == 2 .and. .not. Op%diag) then
            call ZSLGEMM('l', 'n', Op%n, Ndim, Ndim, Op%U, Op%P, Mat)
            call ZSLGEMM('r', 'c', Op%n, Ndim, Ndim, Op%U, Op%P, Mat)
         end if
      end if
   end subroutine Op_Wrapup

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!>
!> @brief
!> Wrapup the Green function
!>
!> @param[inout] Mat(Ndim,Ndim)  Complex
!> \verbatim
!>  N_type = 1, Mat = Op%U*exp(-Op%g*nsigma%phi(1,1)*Op%E)*Mat*exp(Op%g*nsigma%phi(1,1)*Op%E)*(Op%U^{dagger})
!>  N_type = 2, Mat = (Op%U^{dagger}) * Mat * Op%U
!> \endverbatim
!>  **Do not** mix up  \p N_type  and    \p OP\%type
!> @param[in] Op Type(Operator)
!> \verbatim
!>  The operator containing g, U, P
!> \endverbatim
!> @param[in] HS_Field Complex
!> \verbatim
!>  The field
!> \endverbatim
!> @param[in] Ndim Integer
!>
!--------------------------------------------------------------------

   subroutine Op_Wrapdo(Mat, Op, HS_Field, Ndim, N_Type, nt)
      implicit none

      integer :: Ndim
      type(operator), intent(IN)   :: Op
      complex(Kind=kind(0.d0)), intent(INOUT) :: Mat(Ndim, Ndim)
      complex(Kind=kind(0.d0)), intent(IN)   :: HS_Field
      integer, intent(IN) :: N_Type, nt

      ! Local
      integer :: n, i, sp
      complex(Kind=kind(0.d0)) :: VH1(Op%N, OP%N)
      complex(Kind=kind(0.d0)) :: g_loc
      type(Fields)   :: nsigma_single

      call nsigma_single%make(1, 1)
      nsigma_single%f(1, 1) = HS_Field
      nsigma_single%t(1) = op%type

      if (op%type < 3) then
         sp = nint(real(HS_Field))
         if (N_type == 1) then
            if (Op%diag) then
               if (op%g_t_alloc) then
                  do I = 1, Op%N
                     call ZSCAL(Ndim, exp(-nsigma_single%phi(1, 1)*Op%g_t(nt)*Op%E(I)), Mat(Op%P(I), 1), Ndim)
                  end do
                  do I = 1, Op%N
                     call ZSCAL(Ndim, exp(nsigma_single%phi(1, 1)*Op%g_t(nt)*Op%E(I)), Mat(1, Op%P(I)), 1)
                  end do
               else
                  do I = 1, Op%N
                     call ZSCAL(Ndim, Op%E_Exp(I, -sp), Mat(Op%P(I), 1), Ndim)
                  end do
                  do I = 1, Op%N
                     call ZSCAL(Ndim, Op%E_Exp(I, sp), Mat(1, Op%P(I)), 1)
                  end do
               end if
            else
               if (op%g_t_alloc) then
                  do n = 1, Op%N
                     VH1(:, n) = Op%U(:, n)*exp(-nsigma_single%phi(1, 1)*Op%g_t(nt)*Op%E(n))
                  end do
                  call ZSLGEMM('l', 'n', Op%n, Ndim, Ndim, VH1, Op%P, Mat)
                  do n = 1, Op%N
                     VH1(:, n) = exp(nsigma_single%phi(1, 1)*Op%g_t(nt)*Op%E(n))*conjg(Op%U(:, n))
                  end do
                  call ZSLGEMM('r', 'T', Op%n, Ndim, Ndim, VH1, Op%P, Mat)
               else
                  do n = 1, Op%N
                     VH1(:, n) = Op%U(:, n)*Op%E_Exp(n, -sp)
                  end do
                  call ZSLGEMM('l', 'n', Op%n, Ndim, Ndim, VH1, Op%P, Mat)
                  do n = 1, Op%N
                     VH1(:, n) = Op%E_Exp(n, sp)*conjg(Op%U(:, n))
                  end do
                  call ZSLGEMM('r', 'T', Op%n, Ndim, Ndim, VH1, Op%P, Mat)
               end if
            end if
         elseif (N_Type == 2 .and. .not. Op%diag) then
            call ZSLGEMM('r', 'n', Op%n, Ndim, Ndim, Op%U, Op%P, Mat)
            call ZSLGEMM('l', 'c', Op%n, Ndim, Ndim, Op%U, Op%P, Mat)
         end if
      else
         if (N_type == 1) then
            g_loc = Op%g
            if (op%g_t_alloc) g_loc = Op%g_t(nt)
            if (Op%diag) then
               do I = 1, Op%N
                  call ZSCAL(Ndim, exp(-nsigma_single%phi(1, 1)*g_loc*Op%E(I)), Mat(Op%P(I), 1), Ndim)
               end do
               do I = 1, Op%N
                  call ZSCAL(Ndim, exp(nsigma_single%phi(1, 1)*g_loc*Op%E(I)), Mat(1, Op%P(I)), 1)
               end do
            else
               do n = 1, Op%N
                  VH1(:, n) = Op%U(:, n)*exp(-nsigma_single%phi(1, 1)*g_loc*Op%E(n))
               end do
               call ZSLGEMM('l', 'n', Op%n, Ndim, Ndim, VH1, Op%P, Mat)
               do n = 1, Op%N
                  VH1(:, n) = exp(nsigma_single%phi(1, 1)*g_loc*Op%E(n))*conjg(Op%U(:, n))
               end do
               call ZSLGEMM('r', 'T', Op%n, Ndim, Ndim, VH1, Op%P, Mat)
            end if
         elseif (N_Type == 2 .and. .not. Op%diag) then
            call ZSLGEMM('r', 'n', Op%n, Ndim, Ndim, Op%U, Op%P, Mat)
            call ZSLGEMM('l', 'c', Op%n, Ndim, Ndim, Op%U, Op%P, Mat)
         end if
      end if
   end subroutine Op_Wrapdo

   function Op_is_real(Op) result(retval)
      implicit none

      type(operator), intent(IN)   :: Op
      logical :: retval
      real(Kind=kind(0.d0)) :: myzero
      integer :: i, j

      retval = (abs(aimag(Op%g)) < abs(Op%g)*epsilon(1.d0))
      ! calculate a matrix scale
      myzero = maxval(abs(Op%E))*epsilon(Op%E)

      do i = 1, Op%N
         do j = 1, Op%N
            retval = retval .and. (abs(aimag(Op%O(i, j))) < myzero)
         end do
      end do
   end function Op_is_real

   logical function operator_get_g_t_alloc(this)
      implicit none

      class(operator), intent(IN)   :: this

      operator_get_g_t_alloc = this%g_t_alloc

   end function operator_get_g_t_alloc
end module Operator_mod
