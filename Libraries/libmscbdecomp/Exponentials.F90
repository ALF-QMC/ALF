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

module Exponentials_mod
    implicit none
    
    type :: Simplenode
        integer :: x,y
        complex (kind=kind(0.d0)) :: s,c
    end type Simplenode

    type :: node
        integer :: x,y, col
        complex (kind=kind(0.d0)) :: axy
    end type node

    type :: SingleColExp
        integer :: nrofentries
        type(Simplenode), dimension(:), allocatable :: nodes
    contains
        procedure :: init => SingleColExp_init
        procedure :: dealloc => SingleColExp_dealloc
        procedure :: vecmult => SingleColExp_vecmult
        procedure :: matmult => SingleColExp_matmult
    end type


!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This holds together a set of exponentials that, if applied in the
!> correct order approximate e^A.
!> It provides functions for matrix-matrix and matrix-vector
!>  multiplications.
!--------------------------------------------------------------------
    type :: FullExp
        integer :: nrofcols
        type(SingleColExp), dimension(:), allocatable :: singleexps
    contains
        procedure :: init => FullExp_init
        procedure :: dealloc => FullExp_dealloc
        procedure :: vecmult => FullExp_vecmult
        procedure :: matmult => FullExp_matmult
    end type FullExp
    
contains

subroutine SingleColExp_vecmult(this, vec)
    class(SingleColExp) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i
    complex(kind=kind(0.D0)) :: t1,t2
    do i = 1, this%nrofentries! for every matrix
        t1 = vec(this%nodes(i)%x)
        t2 = vec(this%nodes(i)%y)
        vec(this%nodes(i)%x) = this%nodes(i)%c * t1 + this%nodes(i)%s* t2
        vec(this%nodes(i)%y) = this%nodes(i)%c * t2 + this%nodes(i)%s* t1
    enddo
end subroutine SingleColExp_vecmult

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Perform the multiplication of this exponential with a matrix.
!
!> @param[in] this The vertex that we consider
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
subroutine SingleColExp_matmult(this, mat)
    class(SingleColExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i, j, k, ndim
    integer, parameter :: step = 2
    complex(kind=kind(0.D0)) :: t1(step), t2(step)
    
    ndim = size(mat,1)
    do j = 1, ndim, step
        do i = 1, this%nrofentries! for every matrix
            do k = 1,step
                t1(k) = mat(this%nodes(i)%x, j+k)
                t2(k) = mat(this%nodes(i)%y, j+k)
            enddo
            do k = 1, step
                mat(this%nodes(i)%x, j+k) = this%nodes(i)%c * t1(k) + this%nodes(i)%s* t2(k)
                mat(this%nodes(i)%y, j+k) = this%nodes(i)%c * t2(k) + this%nodes(i)%s* t1(k)
            enddo
        enddo
    enddo
end subroutine SingleColExp_matmult

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This sets up the data to perform the exponentiation of a 
!> strictly sparse matrix.
!
!> @param[inout] this the SingleColExp object.
!> @param[in] nodes The nodes that belng to this color.
!> @param[in] nredges how many nodes of this color.
!--------------------------------------------------------------------
subroutine SingleColExp_init(this, nodes, nredges)
    class(SingleColExp) :: this
    type(node), dimension(:), intent(in) :: nodes
    integer, intent(in) :: nredges
    integer :: i
    allocate(this%nodes(nredges))
    this%nrofentries = nredges
#ifndef NDEBUG
    write(*,*) "Setting up strict. sparse matrix with ", nredges, "edges"
#endif
    do i = 1, nredges
        this%nodes(i)%x = nodes(i)%x
        this%nodes(i)%y = nodes(i)%y
        this%nodes(i)%c = cosh(nodes(i)%axy)
        this%nodes(i)%s = sinh(nodes(i)%axy)
    enddo
! All nodes that we have been passed are now from a single color.
! They constitute now a strictly sparse matrix.
! Further processing of the entries could be done here.
end subroutine SingleColExp_init

subroutine SingleColExp_dealloc(this)
    class(SingleColExp) :: this
    deallocate(this%nodes)
end subroutine SingleColExp_dealloc

subroutine FullExp_dealloc(this)
    class(FullExp) :: this
    deallocate(this%singleexps)
end subroutine FullExp_dealloc

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function multiplies this full exponential with a vector
!
!> @param[in] this The exponential opbject
!> @param[in] vec The vector that we multiply
!--------------------------------------------------------------------
subroutine FullExp_vecmult(this, vec)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%vecmult(vec)
    enddo
end subroutine FullExp_vecmult

subroutine FullExp_matmult(this, mat)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%matmult(mat)
    enddo
end subroutine FullExp_matmult

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function creates an exponential object from an array of nodes.
!
!> @param this The exponential opbject
!> @param[in] nodes The array of nodes
!> @param[in] usedcolors the number of used colors/terms in 
!>                       the decomposition.
!--------------------------------------------------------------------
subroutine FullExp_init(this, nodes, usedcolors)
    class(FullExp) :: this
    type(node), dimension(:), intent(in) :: nodes
    integer, intent(in) :: usedcolors
    integer, dimension(:), allocatable :: nredges, edgectr
    integer :: i, maxedges,
    type(node), dimension(:, :), allocatable :: colsepnodes! An array of nodes separated by color
#ifndef NDEBUG
    write(*,*) "Setting up Full Checkerboard exponential."
#endif
    ! Determine the number of matrix entries in each family
    allocate (nredges(usedcolors), edgectr(usedcolors))
    nredges = 0
    this%nrofcols = usedcolors
    do i = 1, size(nodes)
        nredges(nodes(i)%col) = nredges(nodes(i)%col) + 1
    enddo
    maxedges = maxval(nredges)
    edgectr = 1
    allocate(colsepnodes(usedcolors, maxedges))
    do i = 1, size(nodes)
        colsepnodes(nodes(i)%col, edgectr(nodes(i)%col)) = nodes(i)
        edgectr(nodes(i)%col) = edgectr(nodes(i)%col) + 1
    enddo

!     do i = 1, usedcolors
!     write (*,*) edgectr(i), nredges(i)
!     enddo
    ! Now that we have properly separated which entry of a matrix belongs to
    ! which color we can create an exponential for each color that exploits
    ! the structure that the color decomposition creates strictly sparse matrices.
    allocate(this%singleexps(usedcolors))
    do i = 1, usedcolors
        call this%singleexps(i)%init(colsepnodes(i, :), nredges(i))
    enddo
    deallocate(colsepnodes)
    deallocate(nredges, edgectr)
end subroutine FullExp_init

end module Exponentials_mod
