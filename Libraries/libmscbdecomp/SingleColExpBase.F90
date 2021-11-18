! MIT License
! 
! Copyright (c) 2021 Florian Goth
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights 
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in 
! all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.


! We require a common base class to distinguish for optimization purposes
!between matrices that have no chemical potential, that have a uniform chemical potential and those that have arbitrary.
module SingleColExpBase_mod
    implicit none

    ! Base for defining the interface
    type, abstract :: SingleColExpBase
    integer :: nrofentries
    integer, allocatable :: x(:), y(:) ! the y array is still around but often(?) unused
    complex (kind=kind(0.d0)), allocatable :: s(:), s2(:)
    real (kind=kind(0.d0)), allocatable :: c(:), c2(:)
    contains
    procedure(rmultinterface), deferred :: rmult
    procedure(lmultinterface), deferred :: lmult
    procedure(rmultinvinterface), deferred :: rmultinv
    procedure(lmultinvinterface), deferred :: lmultinv
    procedure(adjointactioninterface), deferred :: adjointaction
    procedure(adjointactionovertwointerface), deferred :: adjoint_over_two
    procedure(initinterface), deferred :: init
    procedure(deallocinterface), deferred :: dealloc
    end type SingleColExpBase

    abstract interface
    
    !--------------------------------------------------------------------
    !> @brief 
    !> multiplies this with arg from the right.
    !
    !> @param[in] this
    !> @param[inout] arg
    !--------------------------------------------------------------------
      subroutine rmultinterface(this, mat)
         import SingleColExpBase
         class(SingleColExpBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout),  dimension(:,:) :: mat
      end subroutine

    !--------------------------------------------------------------------
    !> @brief 
    !> multiplies this with arg from the left.
    !
    !> @param[in] this
    !--------------------------------------------------------------------
      subroutine lmultinterface(this, mat)
         import SingleColExpBase
         class(SingleColExpBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout),  dimension(:,:), contiguous :: mat
      end subroutine

    !--------------------------------------------------------------------
    !> @brief 
    !> multiplies this^-1 with arg from the right.
    !
    !> @param[in] this
    !--------------------------------------------------------------------
      subroutine rmultinvinterface(this, mat)
         import SingleColExpBase
         class(SingleColExpBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout),  dimension(:,:) :: mat
      end subroutine

    !--------------------------------------------------------------------
    !> @brief 
    !> multiplies this^-1 with arg from the left.
    !
    !> @param[in] this
    !--------------------------------------------------------------------
      subroutine lmultinvinterface(this, mat)
         import SingleColExpBase
         class(SingleColExpBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout),  dimension(:,:) :: mat
      end subroutine

    !--------------------------------------------------------------------
    !> @brief 
    !> This inititializes an object.
    !
    !> @param[in] this
    !> @param nodes the array of Nodes
    !> @param nredges how many nodes are there
    !> @param mys an array with the diagonal entries
    !> @param a possible prefactor
    !--------------------------------------------------------------------
      subroutine initinterface(this, nodes, nredges, mys, weight)
        Use Node_mod
        import SingleColExpBase
        class(SingleColExpBase), intent(inout) :: this
        type(node), dimension(:), intent(in) :: nodes
        real(kind=kind(0.D0)), intent(in), allocatable, dimension(:) :: mys
        integer, intent(in) :: nredges
        real(kind=kind(0.D0)), intent(in) :: weight
      end subroutine
      
    !--------------------------------------------------------------------
    !> @brief 
    !> Free the used memory
    !
    !> @param[in] this
    !--------------------------------------------------------------------
      subroutine deallocinterface(this)
         import SingleColExpBase
         class(SingleColExpBase), intent(inout) :: this
      end subroutine
      
    !--------------------------------------------------------------------
    !> @brief 
    !> Perform the similarity transform e^{-T} arg e^{T}
    !
    !> @param[in] this
    !> @param[inout] mat the matrix that we intend to transform.
    !--------------------------------------------------------------------
      subroutine adjointactioninterface(this, mat)
         import SingleColExpBase
         class(SingleColExpBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout), dimension(:,:) :: mat
      end subroutine
      
    !--------------------------------------------------------------------
    !> @brief 
    !> Perform the similarity transform e^{-T/2} arg e^{T/2}
    !
    !> @param[in] this
    !> @param[inout] mat the matrix that we intend to transform.
    !--------------------------------------------------------------------
      subroutine adjointactionovertwointerface(this, mat)
         import SingleColExpBase

         class(SingleColExpBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout), dimension(:,:) :: mat
      end subroutine
    end interface
    
    contains
    
!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Perform the multiplication of an checkerboard exponential with a matrix.
!> This is an internal helper function that finds reuse in multiple places.
!
!> @param[in] c the diagonal data
!> @param[in] s the off-diagonal data
!> @param[in] x the used matrix positions
!> @param[in] nrofentries how many vertices are in this family.
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
pure subroutine rmultbase(c, s, x, nrofentries, mat)
    real (kind=kind(0.d0)), allocatable, intent(in) :: c(:)
    complex (kind=kind(0.d0)), allocatable, intent(in) :: s(:)
    integer, allocatable, intent(in) :: x(:)
    integer, intent(in) ::nrofentries
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, ndim
    complex(kind=kind(0.D0)) :: t1, t2

    ndim = size(mat,1)
    do i = 1, nrofentries! for every matrix
        do j = 1, ndim
        t1 = mat(j, x(2*i-1))
        t2 = mat(j, x(2*i))
        mat(j, x(2*i-1)) = c(i) * t1 + s(i)* t2
        mat(j, x(2*i)) = c(i) * t2 + conjg(s(i))* t1
        enddo
    enddo
end subroutine

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Perform the multiplication of this exponential with a matrix: out = this*mat
!
!> Notes: unifying x and y into one array gave some speedup.
!> Unifying c and s did not...
!> FIXME: ndim divisible by two...
!> This is an internal helper function that finds reuse in multiple places.
!
!> @param[in] c the diagonal data
!> @param[in] s the off-diagonal data
!> @param[in] x the used matrix positions
!> @param[in] nrofentries how many vertices are in this family.
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------

pure subroutine lmultbase(c, s, x, nrofentries, mat)
    real (kind=kind(0.d0)), allocatable, intent(in) :: c(:)
    complex (kind=kind(0.d0)), allocatable, intent(in) :: s(:)
    integer, allocatable, intent(in) :: x(:)
    integer, intent(in) ::nrofentries
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2 ! determined to be fastest on 6x6 hubbard
    complex(kind=kind(0.D0)) :: t1(step), t2(step)
    integer, allocatable, dimension(:) :: xyarray
    complex(kind=kind(0.D0)), allocatable, dimension(:) :: snh
    real(kind=kind(0.D0)), allocatable, dimension(:) :: csh

! The intel compiler is really helped by using these temporary arrays
    allocate(xyarray(nrofentries), csh(nrofentries), snh(nrofentries) )
    xyarray = x
    csh = c
    snh = s

    ndim = size(mat,1)
    loopend = (ndim/step)*step

! ifort 2017
!DIR$ UNROLL_AND_JAM(4)
    do j = 1, loopend, step
        do i = 1, nrofentries! for every matrix
            do k = 1,step
                t1(k) = mat(xyarray(2*i-1), j+k-1)
                t2(k) = mat(xyarray(2*i), j+k-1)
            enddo
            do k = 1, step
                mat(xyarray(2*i-1), j+k-1) = csh(i) * t1(k) + snh(i) * t2(k)
                mat(xyarray(2*i), j+k-1) = csh(i) * t2(k) + conjg(snh(i)) * t1(k)
            enddo
        enddo
    enddo
    deallocate(xyarray, csh, snh)
end subroutine
end module SingleColExpBase_mod
