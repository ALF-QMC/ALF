! MIT License
! 
! Copyright (c) 2018-2021 Florian Goth
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
    Use ZeroDiagSingleColExp_mod
    Use HomogeneousSingleColExp_mod
    Use TraceLessSingleColExp_mod
    Use GeneralSingleColExp_mod
    implicit none

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> A wrapper type such that Fortran can hold an array of 
!> base class pointers.
!--------------------------------------------------------------------    
    type :: SingleColeExpBaseWrapper
        class(SingleColExpBase), pointer :: dat => null()
    end type

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This holds together a set of exponentials that, if applied in the
!> correct order, approximate e^A to first order.
!> It provides functions for matrix-matrix and matrix-vector
!> multiplications in transposed and non-transposed manner.
!--------------------------------------------------------------------
    type :: EulerExp
        integer :: nrofcols ! The number of colors
        Type(SingleColeExpBaseWrapper), allocatable :: singleexps(:) ! an array of pointers to the actual polymorphic classes that do the work
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
        procedure :: adjointaction => EulerExp_adjointaction
        procedure :: adjoint_over_two => EulerExp_adjoint_over_two
        procedure :: adjoint_over_two_T => EulerExp_adjoint_over_two_T
        procedure :: rmultinv_T => EulerExp_rmultinv_T
        procedure :: lmultinv_T => EulerExp_lmultinv_T
    end type EulerExp

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This holds together a set of Euler exponentials
!> and applies them in the correct order to obtain higher order
!> approximations.
!> Note that due to the use of higher-order approximations we have successive
!> Strang pairings. Hence the Euler type method needs to be treated separately.
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
        procedure :: adjoint_over_two => FullExp_adjoint_over_two
    end type FullExp

contains

subroutine FullExp_init(this, nodes, usedcolors, dcols, method, weight)
    class(FullExp) :: this
    type(node), dimension(:), intent(in) :: nodes
    integer, intent(in) :: usedcolors, method
    real(kind=kind(0.D0)), intent(in), allocatable, dimension(:, :) :: dcols
    real (kind=kind(0.d0)), intent(in) :: weight
    real (kind=kind(0.d0)) :: tmp

    this%method = method
#ifndef NDEBUG
    write(*,*) "Setting up Full Checkerboard exponential."
#endif
    select case (method)
        case (2)! Strang
            this%evals = 2
            allocate(this%stages(this%evals))
            tmp = 1.D0/2.D0*weight
            call this%stages(1)%init(nodes, usedcolors, dcols, tmp)
            call this%stages(2)%init(nodes, usedcolors, dcols, tmp)
        case (3)! SE_2 2, Blanes
            this%evals = 4
            allocate(this%stages(this%evals))
            tmp = 0.21178*weight
            call this%stages(1)%init(nodes, usedcolors, dcols, tmp)
            tmp = 0.28822*weight
            call this%stages(2)%init(nodes, usedcolors, dcols, tmp)
            tmp = 0.28822*weight
            call this%stages(3)%init(nodes, usedcolors, dcols, tmp)
            tmp = 0.21178*weight
            call this%stages(4)%init(nodes, usedcolors, dcols, tmp)
        case (4)! S_3 4, Yoshida, Neri, Suzuki, 1990
            this%evals = 6
            allocate(this%stages(this%evals))
            tmp =  0.6756035959798*weight
            call this%stages(1)%init(nodes, usedcolors, dcols, tmp)
            tmp =  0.6756035959798*weight
            call this%stages(2)%init(nodes, usedcolors, dcols, tmp)
            tmp = -0.8512071919597*weight
            call this%stages(3)%init(nodes, usedcolors, dcols, tmp)
            tmp = -0.8512071919597*weight
            call this%stages(4)%init(nodes, usedcolors, dcols, tmp)
            tmp =  0.6756035959798*weight
            call this%stages(5)%init(nodes, usedcolors, dcols, tmp)
            tmp =  0.6756035959798*weight
            call this%stages(6)%init(nodes, usedcolors, dcols, tmp)
        case (5)! SE_6 4, Blanes, Blanes and Moan 2002
            this%evals = 12
            allocate(this%stages(this%evals))
            tmp = 0.079203696431195694249716154899943D0 *weight
            call this%stages(1)%init(nodes, usedcolors, dcols, tmp)
            tmp = 0.13031141018216631233261892930386D0 *weight
            call this%stages(2)%init(nodes, usedcolors, dcols, tmp)
            tmp = 0.22286149586760771457161212083520D0 *weight
            call this%stages(3)%init(nodes, usedcolors, dcols, tmp)
            tmp = -0.36671326904742572450057735977680D0 *weight
            call this%stages(4)%init(nodes, usedcolors, dcols, tmp)
            tmp = 0.32464818868970622689484883949262D0 *weight
            call this%stages(5)%init(nodes, usedcolors, dcols, tmp)
            tmp = 0.1096884778767497764517813152452D0 *weight
            call this%stages(6)%init(nodes, usedcolors, dcols, tmp)
            tmp = 0.1096884778767497764517813152452D0 *weight
            call this%stages(7)%init(nodes, usedcolors, dcols, tmp)
            tmp = 0.3246481886897062268948488394926D0 *weight
            call this%stages(8)%init(nodes, usedcolors, dcols, tmp)
            tmp = -0.3667132690474257245005773597768D0 *weight
            call this%stages(9)%init(nodes, usedcolors, dcols, tmp)
            tmp = 0.2228614958676077145716121208352D0 *weight
            call this%stages(10)%init(nodes, usedcolors, dcols, tmp)
            tmp = 0.1303114101821663123326189293039D0 *weight
            call this%stages(11)%init(nodes, usedcolors, dcols, tmp)
            tmp = 0.0792036964311956942497161548999D0 *weight
            call this%stages(12)%init(nodes, usedcolors, dcols, tmp)
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
    complex(kind=kind(0.D0)), dimension(:), intent(inout) :: vec
    integer :: i
    do i = 1, this%evals
       call this%stages(i)%vecmult(vec)
    enddo
end subroutine FullExp_vecmult

subroutine FullExp_vecmult_T(this, vec)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), dimension(:), intent(inout) :: vec
    integer :: i
    do i = this%evals, 1, -1
       call this%stages(i)%vecmult_T(vec)
    enddo
end subroutine FullExp_vecmult_T

subroutine FullExp_lmult(this, mat)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), intent(inout), contiguous :: mat(:,:)
    integer :: i
    do i = this%evals-1, 1, -2
       call this%stages(i+1)%lmult_T(mat)
       call this%stages(i)%lmult(mat)
    enddo
end subroutine FullExp_lmult

subroutine FullExp_adjoint_over_two(this, mat)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), intent(inout) :: mat(:,:)
    integer :: i
    do i = this%evals-1, 1, -2
       call this%stages(i+1)%adjoint_over_two_T(mat)
       call this%stages(i)%adjoint_over_two(mat)
    enddo
end subroutine FullExp_adjoint_over_two

subroutine FullExp_lmultinv(this, mat)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), intent(inout), contiguous :: mat(:,:)
    integer :: i
    do i = 1, this%evals, 2
       call this%stages(i)%lmultinv(mat)
       call this%stages(i+1)%lmultinv_T(mat)
    enddo
end subroutine FullExp_lmultinv

subroutine FullExp_rmult(this, mat)
    class(FullExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
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
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i
    do i = 1, this%evals, 2
       call this%stages(i)%lmult_T(mat)
       call this%stages(i+1)%lmult(mat)
    enddo
end subroutine FullExp_lmult_T

subroutine EulerExp_dealloc(this)
    class(EulerExp) :: this
    integer :: i

    do i = 1, this%nrofcols
        call this%singleexps(i)%dat%dealloc()
        deallocate(this%singleexps(i)%dat)
    enddo
    deallocate(this%singleexps)
end subroutine EulerExp_dealloc

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function multiplies this Euler exponential with a vector.
!
!> @param[in] this The exponential opbject
!> @param[in] vec The vector that we multiply
!--------------------------------------------------------------------
subroutine EulerExp_vecmult(this, vec)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i
    do i = 1, this%nrofcols
       call this%singleexps(i)%dat%vecmult(vec)
    enddo
end subroutine EulerExp_vecmult

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function multiplies this transposed full exponential,
!> i.e reverses the order of application of exponentials, with a vector.
!
!> @param[in] this The exponential opbject
!> @param[in] vec The vector that we multiply
!--------------------------------------------------------------------
subroutine EulerExp_vecmult_T(this, vec)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i
    do i = this%nrofcols, 1, -1
       call this%singleexps(i)%dat%vecmult(vec)
    enddo
end subroutine EulerExp_vecmult_T

subroutine EulerExp_lmultinv(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :), contiguous :: mat
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%dat%lmultinv(mat)
    enddo
end subroutine EulerExp_lmultinv

subroutine EulerExp_lmult(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :), contiguous :: mat
    integer :: i
    do i = this%nrofcols, 1, -1
        call this%singleexps(i)%dat%lmult(mat)
    enddo
end subroutine EulerExp_lmult

subroutine EulerExp_adjointaction(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = this%nrofcols, 1, -1
        call this%singleexps(i)%dat%adjointaction(mat)
    enddo
end subroutine EulerExp_adjointaction

subroutine EulerExp_adjoint_over_two(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = this%nrofcols, 1, -1
        call this%singleexps(i)%dat%adjoint_over_two(mat)
    enddo
end subroutine EulerExp_adjoint_over_two

subroutine EulerExp_adjoint_over_two_T(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%dat%adjoint_over_two(mat)
    enddo
end subroutine EulerExp_adjoint_over_two_T

subroutine EulerExp_rmult(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%dat%rmult(mat)
    enddo
end subroutine EulerExp_rmult

subroutine EulerExp_rmultinv(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = this%nrofcols, 1, -1
        call this%singleexps(i)%dat%rmultinv(mat)
    enddo
end subroutine EulerExp_rmultinv

subroutine EulerExp_rmult_T(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = this%nrofcols, 1, -1
        call this%singleexps(i)%dat%rmult(mat)
    enddo
end subroutine EulerExp_rmult_T

subroutine EulerExp_rmultinv_T(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%dat%rmultinv(mat)
    enddo
end subroutine EulerExp_rmultinv_T

subroutine EulerExp_lmult_T(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = 1, this%nrofcols
        call this%singleexps(i)%dat%lmult(mat)
    enddo
end subroutine EulerExp_lmult_T

subroutine EulerExp_lmultinv_T(this, mat)
    class(EulerExp) :: this
    complex(kind=kind(0.D0)), dimension(:, :) :: mat
    integer :: i
    do i = this%nrofcols, 1, -1
        call this%singleexps(i)%dat%lmultinv(mat)
    enddo
end subroutine EulerExp_lmultinv_T

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> A function to determine the equality of twor floating point numbers.

!> @param[in] a first number
!> @param[in] b second number
!> @return true if the numbers are equal up to 10^-15 of relative deviation.
!--------------------------------------------------------------------
function fpequal(a, b) result(isequal)
    real(kind=kind(0.D0)), intent(in) :: a,b
    logical :: isequal
    
    isequal = .true.
    if(abs(a - b) > max(abs(a), abs(b)) * 1E-15) then
        isequal = .false.
    endif
end function

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function classifies the diagonal.

!> @param[in] mys a vector containing the chemical potentials.
!> @return 0: ZeroDiag, 1:HomogeneousSingleColExp, 2:Traceless: 3: General
!--------------------------------------------------------------------
function determinediagtype(nodes, nrnodes, mys) result(diagtype)
    real(kind=kind(0.D0)), intent(in), dimension(:) :: mys
    type(node), dimension(:), intent(in) :: nodes
    integer, intent(in) :: nrnodes
    integer :: diagtype
    integer :: i
    real(kind=kind(0.D0)) :: localzero
    logical :: isequal, iszero, istraceless

    ! check for pairwise equality
    isequal = .true.
    do i = 1, nrnodes
        if ( .not. fpequal(mys(nodes(i)%x), mys(nodes(i)%y ) )) then
            isequal = .false.
        endif
    enddo
    
    if(isequal) then
    ! check whether they are as good as zero.
        iszero = .true.
        
        do i = 1, nrnodes
            localzero = 1E-15*abs(nodes(i)%axy)*sqrt(2.D0*(mys(nodes(i)%x)/abs(nodes(i)%axy) )**2 + 1.D0)
            if (abs(mys(nodes(i)%x)) > localzero) iszero = .false.
        enddo
        if (iszero) then
            diagtype = 0
        else
            diagtype = 1
        endif
    else
    ! check whether all blocks are traceless
        istraceless = .true.
        
        do i = 1, nrnodes
            localzero = 1E-15*frobnorm(mys(nodes(i)%x), mys(nodes(i)%y), nodes(i)%axy)
            if (abs(mys(nodes(i)%x) + mys(nodes(i)%y)) > localzero) istraceless = .false.
        enddo
        if (istraceless) then
            diagtype = 2
        else
            diagtype = 3
        endif
    endif
end function

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function creates an exponential object from an array of nodes.
!
!> @param this The exponential object
!> @param[in] nodes The array of nodes
!> @param[in] usedcolors the number of used colors/terms in 
!>                       the decomposition.
!> @param[in] dcols a usedcolors x ndim array that contains for each color the diagonal entries to be used.
!> @param[in] weight a prefactor of the exponent.
!--------------------------------------------------------------------
subroutine EulerExp_init(this, nodes, usedcolors, dcols, weight)
    class(EulerExp) :: this
    type(node), dimension(:), intent(in) :: nodes
    integer, intent(in) :: usedcolors
    real(kind=kind(0.D0)), intent(in), allocatable, dimension(:, :) :: dcols
    real(kind=kind(0.d0)), intent(in) :: weight
    integer, dimension(:), allocatable :: nredges ! An array for determining how many edges are there for each color
    integer, dimension(:), allocatable :: edgectr ! A helper array for counting.
    integer :: i, maxedges
! !     integer :: k
! !     character(64) :: filename
    type(node), dimension(:, :), allocatable :: colsepnodes! An array of nodes separated by color
    class(ZeroDiagSingleColExp), pointer :: zerodiagexp => null()
    class(HomogeneousSingleColExp), pointer :: homexp => null()
    class(TraceLessSingleColExp), pointer :: tracelessexp => null()
    class(GeneralSingleColExp), pointer :: generalexp => null()
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

! ! Useful for generating input for mathematica
! !      do i = 1, usedcolors
! !      write (filename, "(A6,I3)") "matrix", i
! !      open(unit=5,file=filename)
! !      do k = 1, nredges(i)
! !      write (5, *) "{{", colsepnodes(i, k)%x, ",",colsepnodes(i, k)%y,"} -> ", dble(colsepnodes(i, k)%axy), "}"
! !      enddo
! !      close(unit=5)
! !      enddo

    ! Now that we have properly separated which entry of a matrix belongs to
    ! which color we can create an exponential for each color that exploits
    ! the structure, that the color decomposition creates strictly sparse matrices.


    allocate(this%singleexps(usedcolors))
    do i = 1, usedcolors

        ! In each color we have to determine which optimizations are possible
        select case(determinediagtype( colsepnodes(i, :), nredges(i), dcols(i,:) ))
        case(0)
            allocate(zerodiagexp)
            this%singleexps(i)%dat => zerodiagexp
        case(1)
            allocate(homexp)
            this%singleexps(i)%dat => homexp
        case(2)
            allocate(tracelessexp)
            this%singleexps(i)%dat => tracelessexp
        case(3)
            allocate(generalexp)
            this%singleexps(i)%dat => generalexp
        end select
        
        call this%singleexps(i)%dat%init(colsepnodes(i, :), nredges(i), dcols(i, :), weight)
    enddo
    deallocate(nredges, edgectr, colsepnodes)
end subroutine EulerExp_init

end module Exponentials_mod
