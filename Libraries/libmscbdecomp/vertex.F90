! MIT License
! 
! Copyright (c) 2018 Florian Goth
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

module vertex_mod
    implicit none
    
    type :: Vertex
        integer :: degree
        integer, allocatable, dimension(:) :: nbrs !< the index of the neighbour in an associated array of Vertex classes
    contains
        procedure, pass(this) :: vertex_init
        generic, public :: init => vertex_init
        procedure :: destruct => vertex_destruct
    end type Vertex

contains

subroutine vertex_init(this, deg)
    class(vertex), intent(inout) :: this
    integer, intent(in) :: deg
    this%degree = deg
    allocate(this%nbrs(deg))
end subroutine vertex_init

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Tidy up allocated space.
!
!> @param[in] this The vertex that we consider
!--------------------------------------------------------------------
subroutine vertex_destruct(this)
    class(vertex) :: this

    deallocate(this%nbrs)
end subroutine vertex_destruct

end module vertex_mod
