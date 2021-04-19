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

! helpful Fortran docs:
! https://materials.prace-ri.eu/400/1/advFortranIntro.pdf
! http://www.chem.helsinki.fi/~manninen/fortran2014/7_Object_oriented_features.pdf

module DynamicMatrixArray_mod
    Use ContainerElementBase_mod
    implicit none

    private
    public :: DynamicMatrixArray

    type :: OpTBasePtrWrapper
        class(ContainerElementBase), pointer :: dat => null()
    end type

    type :: DynamicMatrixArray
        integer :: avamem ! amount of available space
        integer :: tail ! last valid Fortran index
        Type(OpTbasePtrWrapper), allocatable, dimension(:) :: data ! actual effective array of the pointers
    contains
        procedure :: init => DynamicMatrixArray_init
        procedure :: dealloc => DynamicMatrixArray_dealloc
        procedure :: pushback => DynamicMatrixArray_pushback
        procedure :: at => DynamicMatrixArray_at
        procedure :: back => DynamicMatrixArray_back
        procedure :: length => DynamicMatrixArray_length
        ! FIXME: do we need insert?
    end type DynamicMatrixArray

contains

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!
!> @brief
!> set up initial state of the vector.
!>
!> @param[inout] this the vector
!--------------------------------------------------------------------
subroutine DynamicMatrixArray_init(this)
    class(DynamicMatrixArray) :: this
    type(OpTbasePtrWrapper) :: temp
    this%tail = 0 !when the vector has no content this is invalid memory
    this%avamem = 4096/(STORAGE_SIZE(temp)/8) ! allocate a page of memory ! Note STORAGE_SIZE: F2008, SIZEOF: GCC Extension
    allocate(this%data(this%avamem))
end subroutine DynamicMatrixArray_init

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!
!> @brief
!> Deallocates internal storage. Pointees are not deleted!
!>
!> @param[inout] this the vector
!--------------------------------------------------------------------
subroutine DynamicMatrixArray_dealloc(this)
    class(DynamicMatrixArray) :: this
    deallocate(this%data)
end subroutine

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!
!> @brief
!> Attach a pointer to the object given by itm at the end of the vector.
!> If out of space the vector grows.
!>
!> @param[inout] this the vector
!> @param[in] itm the object that we like to store a pointer to.
!
!--------------------------------------------------------------------
subroutine DynamicMatrixArray_pushback(this, itm)
    class(DynamicMatrixArray) :: this
    class(ContainerElementBase), intent(in), target :: itm !Type(...) has to match exactly, class(...) allows for polymorphism
    type(OpTbasePtrWrapper), allocatable, dimension(:) :: temp
    integer :: i

    if (this%tail == this%avamem) then
        ! reallocate the memory
        ! write (*,*) "not enough space -> growing."
        call MOVE_ALLOC(this%data, temp)
        allocate(this%data(2*this%avamem))
        do i = 1, this%avamem
            this%data(i) = temp(i)
        enddo
        deallocate(temp)
        this%avamem = 2*this%avamem
    endif
    this%tail = this%tail + 1
    this%data(this%tail)%dat => itm ! let the pointer point to the object
end subroutine

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!
!> @brief
!> return a pointer to the object stored at position i.
!>
!> @param[in]  i the index.
!> @param[out] the content stored at the position i.
!
!--------------------------------------------------------------------
function DynamicMatrixArray_at(this, pos) result(itm)
    class(DynamicMatrixArray), intent(in) :: this
    integer, intent(in) :: pos
    class(ContainerElementBase), pointer :: itm
    itm => this%data(pos)%dat
end function

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!
!> @brief
!> returns the pointer to the last element
!>
!> @param[inout] this the vector.
!> @return itm the element at the end of the vector.
!
!--------------------------------------------------------------------
function DynamicMatrixArray_back(this) result(itm)
    class(DynamicMatrixArray), intent(in) :: this
    class(ContainerElementBase), pointer :: itm
    itm => this%data(this%tail)%dat
end function

!--------------------------------------------------------------------
!> @author
!> The ALF Project contributors
!
!> @brief
!> Inquire the length of the vector.
!>
!> @param[inout] this the vector
!> @return the current length of the vector
!
!--------------------------------------------------------------------
function DynamicMatrixArray_length(this) result(l)
    class(DynamicMatrixArray) :: this
    integer :: l
    l = this%tail
end function

end module DynamicMatrixArray_mod
