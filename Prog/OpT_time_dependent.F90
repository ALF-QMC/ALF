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
        Complex(kind=kind(0.d0)) :: g
        Real(kind=kind(0.d0)) :: Zero
        Integer :: m, n, Ndim_hop
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
    
    subroutine OpT_time_dependent_init(this, Op_T)
        class(OpT_time_dependent) :: this
        Type(Operator), intent(in) :: Op_T
        Integer :: i, j
        
        this%Zero = 1.E-12
        this%Ndim_hop = Op_T%N
        
    end subroutine

    subroutine OpT_time_dependent_adjointaction(this, arg, t)
        class(OpT_time_dependent), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer, intent(in) :: t
        Integer :: n1, n2

        n1 = size(arg,1)
        n2 = size(arg,2)
        
    end subroutine
    
    subroutine OpT_time_dependent_rmult(this, arg, t)
        class(OpT_time_dependent), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer, intent(in) :: t
        Integer :: n1, n2
        
        ! taken from mmthl
        n1 = size(arg,1)
        n2 = size(arg,2)

    end subroutine
    
        subroutine OpT_time_dependent_rmultinv(this, arg, t)
        class(OpT_time_dependent), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer, intent(in) :: t
        Integer :: n1, n2
        
        ! taken from mmthl_m1
        n1 = size(arg,1)
        n2 = size(arg,2)

    end subroutine
    
    subroutine OpT_time_dependent_lmult(this, arg, t)
        class(OpT_time_dependent), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer, intent(in) :: t
        integer :: n1, n2
        
        ! taken from mmthr
        n1 = size(arg,1)
        n2 = size(arg,2)

    end subroutine
    
    subroutine OpT_time_dependent_lmultinv(this, arg, t)
        class(OpT_time_dependent), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer, intent(in) :: t
        integer :: n1, n2

    end subroutine

    subroutine OpT_time_dependent_dump(this)
        class(OpT_time_dependent), intent(in) :: this
        integer :: i,j

    end subroutine
    
    subroutine OpT_time_dependent_dealloc(this)
        class(OpT_time_dependent), intent(inout) :: this
        
    end subroutine

end module OpT_time_dependent_mod
