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


! We require a common base class to distinguish between matrices that have no chemical potential and those that have.
module SingleColExpBase_mod
    implicit none

    ! Base for defining the interface
    type, abstract :: SingleColExpBase
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
end module SingleColExpBase_mod
