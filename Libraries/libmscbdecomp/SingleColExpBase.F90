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


!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @file SingleColExpBase.F90
!> @brief This file contains the interface class for all checkerboard objects.
!--------------------------------------------------------------------
module SingleColExpBase_mod
    implicit none
!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @class SingleColExpBase
!> @brief This defines the interface for all four checkerboard types.
!>
!> We require a common base class to distinguish for optimization purposes
!> between matrices where the diagonal is zero, 
!> those that have the same value on the entire diagonal(i.e. homogeneous)
!> and those that have an arbitrary diagonal.
!--------------------------------------------------------------------
    type, abstract :: SingleColExpBase
        integer :: nrofentries
        integer, allocatable :: xy(:)
        complex (kind=kind(0.d0)), allocatable :: s(:), s2(:)
        real (kind=kind(0.d0)), allocatable :: c(:), c2(:)
    contains
        procedure(vecmultinterface), deferred :: vecmult
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
    !> Multiplies this with mat from the right.
    !
    !> @param[in] this
    !> @param[inout] mat a complex matrix.
    !--------------------------------------------------------------------
      subroutine rmultinterface(this, mat)
         import SingleColExpBase
         class(SingleColExpBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout),  dimension(:,:) :: mat
      end subroutine

    !--------------------------------------------------------------------
    !> @brief 
    !> Multiplies this with mat from the left.
    !
    !> @param[in] this
    !> @param[inout] mat
    !--------------------------------------------------------------------
      subroutine lmultinterface(this, mat)
         import SingleColExpBase
         class(SingleColExpBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout),  dimension(:,:), contiguous :: mat
      end subroutine

    !--------------------------------------------------------------------
    !> @brief 
    !> Multiplies this with a vector from the left.
    !
    !> @param[in] this
    !> @param[inout] vec
    !--------------------------------------------------------------------
      subroutine vecmultinterface(this, vec)
         import SingleColExpBase
         class(SingleColExpBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout),  dimension(:) :: vec
      end subroutine

    !--------------------------------------------------------------------
    !> @brief 
    !> multiplies this^-1 with mat from the right.
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
    !> multiplies this^-1 with mat from the left.
    !
    !> @param[in] this
    !--------------------------------------------------------------------
      subroutine lmultinvinterface(this, mat)
         import SingleColExpBase
         class(SingleColExpBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout),  dimension(:,:), contiguous :: mat
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
        real(kind=kind(0.D0)), intent(in), dimension(:) :: mys
        integer, intent(in) :: nredges
        real(kind=kind(0.D0)), intent(in) :: weight
      end subroutine

    !--------------------------------------------------------------------
    !> @brief 
    !> Free the used memory.
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
!> @param[in] xy the used matrix positions
!> @param[in] nrofentries how many vertices are in this family.
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
pure subroutine rmultbase(c, s, xy, nrofentries, mat)
    real (kind=kind(0.d0)), allocatable, intent(in) :: c(:)
    complex (kind=kind(0.d0)), allocatable, intent(in) :: s(:)
    integer, allocatable, intent(in) :: xy(:)
    integer, intent(in) ::nrofentries
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, ndim
    complex(kind=kind(0.D0)) :: t1, t2

    ndim = size(mat,1)
    do i = 1, nrofentries! for every matrix
        do j = 1, ndim
        t1 = mat(j, xy(2*i-1))
        t2 = mat(j, xy(2*i))
        mat(j, xy(2*i-1)) = c(i) * t1 + s(i)* t2
        mat(j, xy(2*i)) = c(i) * t2 + conjg(s(i))* t1
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

pure subroutine lmultbase(c, s, xy, nrofentries, mat)
    real (kind=kind(0.d0)), allocatable, intent(in) :: c(:)
    complex (kind=kind(0.d0)), allocatable, intent(in) :: s(:)
    integer, allocatable, intent(in) :: xy(:)
    integer, intent(in) ::nrofentries
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout), contiguous :: mat
    
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2 ! determined to be fastest on 6x6 hubbard
    complex(kind=kind(0.D0)) :: t1(step), t2(step)
    integer, allocatable, dimension(:) :: xyarray
    complex(kind=kind(0.D0)), allocatable, dimension(:) :: snh
    real(kind=kind(0.D0)), allocatable, dimension(:) :: csh

! The intel compiler is really helped by using these temporary arrays
    allocate(xyarray(size(xy)), csh(nrofentries), snh(nrofentries) )
    xyarray = xy
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
    
    ! remainder loop
    if ((ndim - loopend) .ne. 0) then
        do i = 1, nrofentries! for every matrix
            t1(1) = mat(xyarray(2*i-1), ndim)
            t2(1) = mat(xyarray(2*i), ndim)
            mat(xyarray(2*i-1), ndim) = csh(i) * t1(1) + snh(i) * t2(1)
            mat(xyarray(2*i), ndim) = csh(i) * t2(1) + conjg(snh(i)) * t1(1)
        enddo
    endif
    deallocate(xyarray, csh, snh)
end subroutine

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> A function to calculate the Frobenius norm of hermitian 2x2 matrices.
!
!> @param[in] d1 first diagonal entry
!> @param[in] d2 second diagonal entry
!> @param[in] o off-diagonal entry
!> @return The value of the frobenius norm
!--------------------------------------------------------------------
function frobnorm(d1, d2, o) result(fn)
    real (kind=kind(0.d0)) :: fn, d1, d2
    complex(kind=kind(0.D0)), intent(in) :: o
    fn = sqrt(d1*d1+d2*d2 + 2*dble(o * conjg(o)))
end function
end module SingleColExpBase_mod
