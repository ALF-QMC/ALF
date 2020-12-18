! MIT License
! 
! Copyright (c) 2018-2020 Florian Goth
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
        complex (kind=kind(0.d0)) :: s,c ! evaluated sinh and cosh
    end type Simplenode

    type :: node
        integer :: x,y, col
        complex (kind=kind(0.d0)) :: axy ! the value of the matrix A_{x,y}
    end type node

    type :: SingleColExp
        integer :: nrofentries
        type(Simplenode), dimension(:), allocatable :: nodes
    contains
        procedure :: init => SingleColExp_init
        procedure :: dealloc => SingleColExp_dealloc
        procedure :: vecmult => SingleColExp_vecmult
        procedure :: lmult => SingleColExp_lmult
        procedure :: lmultinv => SingleColExp_lmultinv
        procedure :: rmult => SingleColExp_rmult
        procedure :: rmultinv => SingleColExp_rmultinv
    end type

    
!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This holds together a set of exponentials that, if applied in the
!> correct order approximate e^A to first order.
!> It provides functions for matrix-matrix and matrix-vector
!>  multiplications in transposed and non-transposed manner.
!--------------------------------------------------------------------
    type :: EulerExp
        integer :: nrofcols
        type(SingleColExp), allocatable :: singleexps(:)
    contains
         procedure :: init => EulerExp_init
        procedure :: dealloc => EulerExp_dealloc
        procedure :: vecmult => EulerExp_vecmult
        procedure :: vecmult_T => EulerExp_vecmult_T
        procedure :: lmult => EulerExp_lmult
        procedure :: lmultinv => EulerExp_lmultinv
        procedure :: rmult => EulerExp_rmult
        procedure :: rmultinv => EulerExp_rmultinv
        procedure :: rmult_T => EulerExp_rmult_T
        procedure :: lmult_T => EulerExp_lmult_T
        procedure :: rmultinv_T => EulerExp_rmultinv_T
        procedure :: lmultinv_T => EulerExp_lmultinv_T
    end type EulerExp

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This holds together a set of Euler Exponentials
!> and applies them in the correct order to obtain higher order
!> approximations.
!--------------------------------------------------------------------
    type :: FullExp
        integer :: method
        integer :: evals
        type(EulerExp), allocatable :: stages(:)
    contains
        procedure :: init => FullExp_init
        procedure :: dealloc => FullExp_dealloc
        procedure :: vecmult => FullExp_vecmult
        procedure :: vecmult_T => FullExp_vecmult_T
        procedure :: lmult => FullExp_lmult
        procedure :: lmultinv => FullExp_lmultinv
        procedure :: rmult => FullExp_rmult
        procedure :: rmultinv => FullExp_rmultinv
        procedure :: lmult_T => FullExp_lmult_T
    end type FullExp
    
contains

subroutine FullExp_init(this, nodes, usedcolors, method, weight)
    class(FullExp) :: this
    type(node), dimension(:), intent(in) :: nodes
    integer, intent(in) :: usedcolors, method
    complex (kind=kind(0.d0)), intent(in) :: weight
    complex (kind=kind(0.d0)) :: tmp
    integer, dimension(:), allocatable :: nredges, edgectr
    integer :: i, maxedges, k
    type(node), dimension(:, :), allocatable :: colsepnodes! An array of nodes separated by color
    character(len=64) :: filename
    type(EulerExp) :: dummy
    
#ifndef NDEBUG
    write(*,*) "Setting up Full Checkerboard exponential."
#endif
    select case (method)
        case (1)! Euler
            this%evals = 2
            allocate(this%stages(this%evals))
            call this%stages(1)%init(nodes, usedcolors, weight)
            tmp = 0.D0 ! cheat to get Euler method in there
            call this%stages(2)%init(nodes, usedcolors, tmp)
        case (2)! Strang
            this%evals = 2
            allocate(this%stages(this%evals))
            tmp = 1.D0/2.D0*weight
            call this%stages(1)%init(nodes, usedcolors, tmp)
            call this%stages(2)%init(nodes, usedcolors, tmp)
        case (3)! SE_2 2
            this%evals = 4
            allocate(this%stages(this%evals))
            tmp = 0.21178*weight
            call this%stages(1)%init(nodes, usedcolors, tmp)
            tmp = 0.28822*weight
            call this%stages(2)%init(nodes, usedcolors, tmp)
            tmp = 0.28822*weight
            call this%stages(3)%init(nodes, usedcolors, tmp)
            tmp = 0.21178*weight
            call this%stages(4)%init(nodes, usedcolors, tmp)
        case (4)! SE_3 4, Yoshida, Neri
            this%evals = 6
            allocate(this%stages(this%evals))
            tmp = 0.675604*weight
            call this%stages(1)%init(nodes, usedcolors, tmp)
            tmp = 0.675604*weight
            call this%stages(2)%init(nodes, usedcolors, tmp)
            tmp = -0.851207*weight
            call this%stages(3)%init(nodes, usedcolors, tmp)
            tmp = -0.851207*weight
            call this%stages(4)%init(nodes, usedcolors, tmp)
            tmp = 0.675604*weight
            call this%stages(5)%init(nodes, usedcolors, tmp)
            tmp = 0.675604*weight
            call this%stages(6)%init(nodes, usedcolors, tmp)
        case (5)! SE_6 4, Blanes
            this%evals = 12
            allocate(this%stages(this%evals))
            tmp = 0.0792037*weight
            call this%stages(1)%init(nodes, usedcolors, tmp)
            tmp = 0.130311*weight
            call this%stages(2)%init(nodes, usedcolors, tmp)
            tmp = 0.222861*weight
            call this%stages(3)%init(nodes, usedcolors, tmp)
            tmp = -0.366713*weight
            call this%stages(4)%init(nodes, usedcolors, tmp)
            tmp = 0.324648*weight
            call this%stages(5)%init(nodes, usedcolors, tmp)
            tmp = 0.109688*weight
            call this%stages(6)%init(nodes, usedcolors, tmp)
            tmp = 0.109688*weight
            call this%stages(7)%init(nodes, usedcolors, tmp)
            tmp = 0.324648*weight
            call this%stages(8)%init(nodes, usedcolors, tmp)
            tmp = -0.366713*weight
            call this%stages(9)%init(nodes, usedcolors, tmp)
            tmp = 0.222861*weight
            call this%stages(10)%init(nodes, usedcolors, tmp)
            tmp = 0.130311*weight
            call this%stages(11)%init(nodes, usedcolors, tmp)
            tmp = 0.0792037*weight
            call this%stages(12)%init(nodes, usedcolors, tmp)
    end select
end subroutine FullExp_init

subroutine FullExp_dealloc(this)
    class(FullExp) :: this
    integer :: i
    do i = 1, this%evals
       call this%stages(i)%dealloc()
    enddo
    deallocate(this%stages)
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
    do i = 1, this%evals
       call this%stages(i)%vecmult(vec)
    enddo
end subroutine FullExp_vecmult

subroutine FullExp_vecmult_T(this, vec)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i
    do i = this%evals, 1, -1
       call this%stages(i)%vecmult_T(vec)
    enddo
end subroutine FullExp_vecmult_T

subroutine FullExp_lmult(this, mat)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), intent(inout) :: mat(:,:)
    integer :: i
    do i = this%evals-1, 1, -2
       call this%stages(i+1)%lmult_T(mat)
       call this%stages(i)%lmult(mat)
    enddo
end subroutine FullExp_lmult

subroutine FullExp_lmultinv(this, mat)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), intent(inout) :: mat(:,:)
    integer :: i
    do i = 1, this%evals, 2
       call this%stages(i)%lmultinv(mat)
       call this%stages(i+1)%lmultinv_T(mat)
    enddo
end subroutine FullExp_lmultinv

subroutine FullExp_rmult(this, mat)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = 1, this%evals,2
       call this%stages(i)%rmult(mat)
       call this%stages(i+1)%rmult_T(mat)
    enddo
end subroutine FullExp_rmult

subroutine FullExp_rmultinv(this, mat)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), intent(inout) :: mat(:,:)
    integer :: i
    do i = this%evals-1, 1, -2
       call this%stages(i+1)%rmultinv_T(mat)
       call this%stages(i)%rmultinv(mat)
    enddo
end subroutine FullExp_rmultinv

subroutine FullExp_lmult_T(this, mat)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = 1, this%evals, 2
       call this%stages(i)%lmult_T(mat)
       call this%stages(i+1)%lmult(mat)
    enddo
end subroutine FullExp_lmult_T


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
!> Perform the multiplication of this exponential with a matrix: out = this*mat
!
!> @param[in] this The exponential that we consider
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
subroutine SingleColExp_lmult(this, mat)
    class(SingleColExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2
    complex(kind=kind(0.D0)) :: t1(step), t2(step)
    
    ndim = size(mat,1)
    loopend = (ndim/step)*step
    do j = 1, loopend, step
        do i = 1, this%nrofentries! for every matrix
            do k = 1,step
                t1(k) = mat(this%nodes(i)%x, j+k-1)
                t2(k) = mat(this%nodes(i)%y, j+k-1)
            enddo
            do k = 1, step
                mat(this%nodes(i)%x, j+k-1) = this%nodes(i)%c * t1(k) + this%nodes(i)%s* t2(k)
                mat(this%nodes(i)%y, j+k-1) = this%nodes(i)%c * t2(k) + this%nodes(i)%s* t1(k)
            enddo
        enddo
    enddo
end subroutine SingleColExp_lmult

subroutine SingleColExp_lmultinv(this, mat)
    class(SingleColExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2
    complex(kind=kind(0.D0)) :: t1(step), t2(step)
    
    ndim = size(mat,1)
    loopend = (ndim/step)*step
    do j = 1, loopend, step
        do i = 1, this%nrofentries! for every matrix
            do k = 1,step
                t1(k) = mat(this%nodes(i)%x, j+k-1)
                t2(k) = mat(this%nodes(i)%y, j+k-1)
            enddo
            do k = 1, step
                mat(this%nodes(i)%x, j+k-1) = this%nodes(i)%c * t1(k) - this%nodes(i)%s* t2(k)
                mat(this%nodes(i)%y, j+k-1) = this%nodes(i)%c * t2(k) - this%nodes(i)%s* t1(k)
            enddo
        enddo
    enddo
end subroutine SingleColExp_lmultinv

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Perform the multiplication of this exponential with a matrix: out = mat*this
!
!> @param[in] this The exponential that we consider
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
subroutine SingleColExp_rmult(this, mat)
    class(SingleColExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i, j, k, ndim
    complex(kind=kind(0.D0)) :: t1, t2
    
    ndim = size(mat,1)
    do i = 1, this%nrofentries! for every matrix
        do j = 1, ndim
        t1 = mat(j, this%nodes(i)%x)
        t2 = mat(j, this%nodes(i)%y)
        mat(j, this%nodes(i)%x) = this%nodes(i)%c * t1 + this%nodes(i)%s* t2
        mat(j, this%nodes(i)%y) = this%nodes(i)%c * t2 + this%nodes(i)%s* t1
        enddo
    enddo
end subroutine SingleColExp_rmult

subroutine SingleColExp_rmultinv(this, mat)
    class(SingleColExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i, j, k, ndim
    complex(kind=kind(0.D0)) :: t1, t2
    
    ndim = size(mat,1)
    do i = 1, this%nrofentries! for every matrix
        do j = 1, ndim
        t1 = mat(j, this%nodes(i)%x)
        t2 = mat(j, this%nodes(i)%y)
        mat(j, this%nodes(i)%x) = this%nodes(i)%c * t1 - this%nodes(i)%s* t2
        mat(j, this%nodes(i)%y) = this%nodes(i)%c * t2 - this%nodes(i)%s* t1
        enddo
    enddo
end subroutine SingleColExp_rmultinv

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
!> @param[in] weight a prefactor for the exponent.
!--------------------------------------------------------------------
subroutine SingleColExp_init(this, nodes, nredges, weight)
    class(SingleColExp) :: this
    type(node), dimension(:), intent(in) :: nodes
    integer, intent(in) :: nredges
    complex (kind=kind(0.d0)), intent(in) :: weight
    integer :: i
    allocate(this%nodes(nredges))
    this%nrofentries = nredges
#ifndef NDEBUG
    write(*,*) "Setting up strict. sparse matrix with ", nredges, "edges"
#endif
    do i = 1, nredges
        this%nodes(i)%x = nodes(i)%x
        this%nodes(i)%y = nodes(i)%y
        this%nodes(i)%c = cosh(weight*nodes(i)%axy)
        this%nodes(i)%s = sinh(weight*nodes(i)%axy)
    enddo
! All nodes that we have been passed are now from a single color.
! They constitute now a strictly sparse matrix.
! Further processing of the entries could be done here.
end subroutine SingleColExp_init

subroutine SingleColExp_dealloc(this)
    class(SingleColExp) :: this
    deallocate(this%nodes)
end subroutine SingleColExp_dealloc

subroutine EulerExp_dealloc(this)
    class(EulerExp) :: this
    deallocate(this%singleexps)
end subroutine EulerExp_dealloc

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
subroutine EulerExp_vecmult(this, vec)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i
    do i = 1, this%nrofcols
       call this%singleexps(i)%vecmult(vec)
    enddo
end subroutine EulerExp_vecmult

subroutine EulerExp_vecmult_T(this, vec)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i
    do i = this%nrofcols, 1, -1
       call this%singleexps(i)%vecmult(vec)
    enddo
end subroutine EulerExp_vecmult_T

subroutine EulerExp_lmultinv(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%lmultinv(mat)
    enddo
end subroutine EulerExp_lmultinv

subroutine EulerExp_lmult(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = this%nrofcols, 1, -1
        call this%singleexps(i)%lmult(mat)
    enddo
end subroutine EulerExp_lmult

subroutine EulerExp_rmult(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%rmult(mat)
    enddo
end subroutine EulerExp_rmult

subroutine EulerExp_rmultinv(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = this%nrofcols, 1, -1
        call this%singleexps(i)%rmultinv(mat)
    enddo
end subroutine EulerExp_rmultinv

subroutine EulerExp_rmult_T(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = this%nrofcols, 1, -1
        call this%singleexps(i)%rmult(mat)
    enddo
end subroutine EulerExp_rmult_T

subroutine EulerExp_rmultinv_T(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%rmultinv(mat)
    enddo
end subroutine EulerExp_rmultinv_T

subroutine EulerExp_lmult_T(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%lmult(mat)
    enddo
end subroutine EulerExp_lmult_T

subroutine EulerExp_lmultinv_T(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = this%nrofcols, 1, -1
        call this%singleexps(i)%lmultinv(mat)
    enddo
end subroutine EulerExp_lmultinv_T
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
!> @param[in] weight a prefactor of the exponent
!--------------------------------------------------------------------
subroutine EulerExp_init(this, nodes, usedcolors, weight)
    class(EulerExp) :: this
    type(node), dimension(:), intent(in) :: nodes
    integer, intent(in) :: usedcolors
    complex (kind=kind(0.d0)), intent(in) :: weight
    integer, dimension(:), allocatable :: nredges, edgectr
    integer :: i, maxedges, k
    type(node), dimension(:, :), allocatable :: colsepnodes! An array of nodes separated by color
    character(len=64) :: filename
#ifndef NDEBUG
    write(*,*) "Setting up Euler Checkerboard exponential."
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
    do i = 1, usedcolors
    write (filename, "(A6,I3)") "matrix", i
    open(unit=5,file=filename)
    do k = 1, nredges(i)
    write (5, *) colsepnodes(i, k)%x, colsepnodes(i, k)%y, dble(colsepnodes(i, k)%axy)
    enddo
    enddo
!     do i = 1, usedcolors
!     write (*,*) edgectr(i), nredges(i)
!     enddo
    ! Now that we have properly separated which entry of a matrix belongs to
    ! which color we can create an exponential for each color that exploits
    ! the structure that the color decomposition creates strictly sparse matrices.
    allocate(this%singleexps(usedcolors))
    do i = 1, usedcolors
        call this%singleexps(i)%init(colsepnodes(i, :), nredges(i), weight)
    enddo
    deallocate(colsepnodes)
    deallocate(nredges, edgectr)
end subroutine EulerExp_init

end module Exponentials_mod
