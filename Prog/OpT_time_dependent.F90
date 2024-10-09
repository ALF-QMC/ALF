!  Copyright (C) 2022 The ALF project
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
!       https://alf.physik.uni-wuerzburg.de .
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage https://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

module OpT_time_dependent_mod
   use ContainerElementBase_mod
   use Operator_mod
   use mat_subroutines
   implicit none

   private
   public :: OpT_time_dependent

   !--------------------------------------------------------------------
   !> @author
   !> ALF-project
   !> @brief
   !> Encapsulates operations for complex exponentiated OpTs.
   !>
   !--------------------------------------------------------------------
   type, extends(ContainerElementBase) :: OpT_time_dependent
      complex(kind=kind(0.d0)), allocatable, dimension(:, :)  :: U(:, :)  !> We  store    the unitary  transformation
      real(kind=kind(0.d0)), allocatable, dimension(:)  :: E(:)  !> We  store  the  real  eigenvaules
      complex(kind=kind(0.d0)), allocatable  :: g_t(:)   !>
      real(kind=kind(0.d0)) :: Zero
      integer, allocatable :: P(:)
      integer :: Ndim_hop
      logical :: symm

      ! Assumption: T = U * E * U^\dagger
   contains
      procedure :: init => OpT_time_dependent_init ! initialize and allocate matrices
      procedure :: dealloc => OpT_time_dependent_dealloc ! dealloc matrices
      procedure :: rmult => OpT_time_dependent_rmult ! right multiplication with Op_T
      procedure :: lmult => OpT_time_dependent_lmult
      procedure :: rmultinv => OpT_time_dependent_rmultinv ! right multiplication with Op_T inverse
      procedure :: lmultinv => OpT_time_dependent_lmultinv
      procedure :: adjointaction => OpT_time_dependent_adjointaction
      procedure :: dump => OpT_time_dependent_dump ! dump matrices for debugging to screen
   end type OpT_time_dependent

contains

   subroutine OpT_time_dependent_init(this, Op_T, symm)
      class(OpT_time_dependent) :: this
      type(operator), intent(in) :: Op_T
      logical, intent(in) :: symm

      this%Zero = 1.e-12
      this%Ndim_hop = Op_T%N
      this%P = Op_T%P ! copy all data locally to be consistent and less error prone
      this%U = Op_T%U
      this%E = Op_T%E
      this%g_t = Op_T%g_t
      this%symm = symm

   end subroutine

   subroutine OpT_time_dependent_adjointaction(this, arg, t1, t2)
      class(OpT_time_dependent), intent(in) :: this
      complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
      integer, intent(in) :: t1
      integer, optional, intent(in) :: t2
      integer :: n1, n2, i, j, nt1, nt2
      complex(kind=kind(0.d0)) :: te

      n1 = size(arg, 1)
      n2 = size(arg, 2)
      nt1 = t1
      if (t1 == 0) nt1 = size(this%g_t, 1)

      if (present(t2)) then
         nt2 = t2
         if (t2 == 0) nt2 = size(this%g_t, 1)
      else
         nt2 = nt1
      end if

      if (dble(this%g_t(nt1)*conjg(this%g_t(nt1))) > this%Zero) then

         ! innermost unitary transform
         call ZSLGEMM('L', 'C', this%Ndim_hop, n1, n2, this%U, this%P, arg)

         ! apply both diagonals
         do i = 1, size(this%e, 1)
            te = exp(0.5d0*this%g_t(nt1)*this%E(i))
            do j = 1, n2
               arg(this%P(i), j) = arg(this%P(i), j)*te
            end do
         end do

         ! outermost unitary transform
         call ZSLGEMM('L', 'N', this%Ndim_hop, n1, n2, this%U, this%P, arg)
      end if

      if (dble(this%g_t(nt2)*conjg(this%g_t(nt2))) > this%Zero) then

         ! innermost unitary transform
         call ZSLGEMM('R', 'N', this%Ndim_hop, n1, n2, this%U, this%P, arg)

         ! apply both diagonals
         do i = 1, n1
            do j = 1, size(this%E, 1)
               arg(i, this%P(j)) = arg(i, this%P(j))*exp(-0.5d0*this%g_t(nt2)*this%E(j))
            end do
         end do

         ! outermost unitary transform
         call ZSLGEMM('R', 'C', this%Ndim_hop, n1, n2, this%U, this%P, arg)
      end if

   end subroutine

   subroutine OpT_time_dependent_rmult(this, arg, t)
      class(OpT_time_dependent), intent(in) :: this
      complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
      integer, intent(in) :: t
      integer :: n1, n2, i, j, t1
      complex(kind=kind(0.d0)) :: g

      ! taken from mmthl
      n1 = size(arg, 1)
      n2 = size(arg, 2)
      g = this%g_t(t)
      if (this%symm) then
         t1 = t - 1
         if (t == 1) t1 = size(this%g_t, 1)
         g = 0.5d0*(this%g_t(t) + this%g_t(t1))
      end if

      if (dble(g*conjg(g)) > this%Zero) then
         call ZSLGEMM('R', 'N', this%Ndim_hop, n1, n2, this%U, this%P, arg)
         do i = 1, n1
            do j = 1, size(this%E, 1)
               arg(i, this%P(j)) = arg(i, this%P(j))*exp(g*this%E(j))
            end do
         end do
         call ZSLGEMM('R', 'C', this%Ndim_hop, n1, n2, this%U, this%P, arg)
      end if

   end subroutine

   subroutine OpT_time_dependent_rmultinv(this, arg, t)
      class(OpT_time_dependent), intent(in) :: this
      complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
      integer, intent(in) :: t
      integer :: n1, n2, i, j, t1
      complex(kind=kind(0.d0)) :: g

      ! taken from mmthl_m1
      n1 = size(arg, 1)
      n2 = size(arg, 2)
      g = this%g_t(t)
      if (this%symm) then
         t1 = t - 1
         if (t == 1) t1 = size(this%g_t, 1)
         g = 0.5d0*(this%g_t(t) + this%g_t(t1))
      end if

      if (dble(g*conjg(g)) > this%Zero) then
         call ZSLGEMM('R', 'N', this%Ndim_hop, n1, n2, this%U, this%P, arg)
         do i = 1, n1
            do j = 1, size(this%E, 1)
               arg(i, this%P(j)) = arg(i, this%P(j))*exp(-g*this%E(j))
            end do
         end do
         call ZSLGEMM('R', 'C', this%Ndim_hop, n1, n2, this%U, this%P, arg)
      end if

   end subroutine

   subroutine OpT_time_dependent_lmult(this, arg, t)
      class(OpT_time_dependent), intent(in) :: this
      complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
      integer, intent(in) :: t
      integer :: n1, n2, i, j, t1
      complex(kind=kind(0.d0)) :: te, g

      ! taken from mmthr
      n1 = size(arg, 1)
      n2 = size(arg, 2)
      g = this%g_t(t)
      if (this%symm) then
         t1 = t - 1
         if (t == 1) t1 = size(this%g_t, 1)
         g = 0.5d0*(this%g_t(t) + this%g_t(t1))
      end if

      if (dble(g*conjg(g)) > this%Zero) then
         call ZSLGEMM('L', 'C', this%Ndim_hop, n1, n2, this%U, this%P, arg)
         do i = 1, size(this%e, 1)
            te = exp(g*this%E(i))
            do j = 1, n2
               arg(this%P(i), j) = arg(this%P(i), j)*te
            end do
         end do
         call ZSLGEMM('L', 'N', this%Ndim_hop, n1, n2, this%U, this%P, arg)
      end if

   end subroutine

   subroutine OpT_time_dependent_lmultinv(this, arg, t)
      class(OpT_time_dependent), intent(in) :: this
      complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
      integer, intent(in) :: t
      integer :: n1, n2, i, j, t1
      complex(kind=kind(0.d0)) :: te, g

      n1 = size(arg, 1)
      n2 = size(arg, 2)
      g = this%g_t(t)
      if (this%symm) then
         t1 = t - 1
         if (t == 1) t1 = size(this%g_t, 1)
         g = 0.5d0*(this%g_t(t) + this%g_t(t1))
      end if

      if (dble(g*conjg(g)) > this%Zero) then
         call ZSLGEMM('L', 'C', this%Ndim_hop, n1, n2, this%U, this%P, arg)
         do i = 1, size(this%e, 1)
            te = exp(-g*this%E(i))
            do j = 1, n2
               arg(this%P(i), j) = arg(this%P(i), j)*te
            end do
         end do
         call ZSLGEMM('L', 'N', this%Ndim_hop, n1, n2, this%U, this%P, arg)
      end if

   end subroutine

   subroutine OpT_time_dependent_dump(this)
      class(OpT_time_dependent), intent(in) :: this
      integer :: i, j

      do i = 1, size(this%U, 1)
         write (*, *) (dble(this%U(i, j)), j=1, size(this%U, 2))
      end do
      write (*, *) "------E--------"
      do i = 1, size(this%E, 1)
         write (*, *) this%E(i)
      end do

      write (*, *) "------g_t--------"
      do i = 1, size(this%g_t, 1)
         write (*, *) this%g_t(i)
      end do

   end subroutine

   subroutine OpT_time_dependent_dealloc(this)
      class(OpT_time_dependent), intent(inout) :: this

      deallocate (this%U, this%E, this%g_t, this%P)
   end subroutine

end module OpT_time_dependent_mod
