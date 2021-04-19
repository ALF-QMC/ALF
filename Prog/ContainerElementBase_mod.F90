!  Copyright (C) 2020 The ALF project
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


! Declare a common base class for the interface: it multiplies with complex matrices
module ContainerElementBase_mod
    implicit none

    private
    public :: ContainerElementBase

    ! Base for defining the interface
    type, abstract :: ContainerElementBase
    contains
    procedure(rmultinterface), deferred :: rmult
    procedure(lmultinterface), deferred :: lmult
    procedure(rmultinvinterface), deferred :: rmultinv
    procedure(lmultinvinterface), deferred :: lmultinv
    procedure(adjointactioninterface), deferred :: adjointaction
    procedure(dump), deferred :: dump
    procedure(dealloc), deferred :: dealloc
    end type ContainerElementBase

    abstract interface
    
    
    !--------------------------------------------------------------------
    !> @brief 
    !> multiplies this with arg from the right.
    !
    !> @param[in] this
    !--------------------------------------------------------------------
      subroutine rmultinterface(this, arg)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout),  dimension(:,:) :: arg
      end subroutine

    !--------------------------------------------------------------------
    !> @brief 
    !> multiplies this with arg from the left.
    !
    !> @param[in] this
    !--------------------------------------------------------------------
      subroutine lmultinterface(this, arg)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout),  dimension(:,:) :: arg
      end subroutine

    !--------------------------------------------------------------------
    !> @brief 
    !> multiplies this^-1 with arg from the right.
    !
    !> @param[in] this
    !--------------------------------------------------------------------
      subroutine rmultinvinterface(this, arg)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout),  dimension(:,:) :: arg
      end subroutine

    !--------------------------------------------------------------------
    !> @brief 
    !> multiplies this^-1 with arg from the left.
    !
    !> @param[in] this
    !--------------------------------------------------------------------
      subroutine lmultinvinterface(this, arg)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout),  dimension(:,:) :: arg
      end subroutine

    !--------------------------------------------------------------------
    !> @brief 
    !> This dumps the content to the screen.
    !
    !> @param[in] this
    !--------------------------------------------------------------------
      subroutine dump(this)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
      end subroutine
      
    !--------------------------------------------------------------------
    !> @brief 
    !> Free the used memory
    !
    !> @param[in] this
    !--------------------------------------------------------------------
      subroutine dealloc(this)
         import ContainerElementBase
         class(ContainerElementBase), intent(inout) :: this
      end subroutine
      
    !--------------------------------------------------------------------
    !> @brief 
    !> Perform the similarity transform e^{-T/2} arg e^{T/2}
    !
    !> @param[in] this
    !> @param[inout] the matrix that we intend to transform.
    !--------------------------------------------------------------------
      subroutine adjointactioninterface(this, arg)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         Complex(kind=kind(0.d0)), intent(inout), dimension(:,:) :: arg
      end subroutine
    end interface
end module ContainerElementBase_mod
