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

module NonHomogeneousSingleColExp_mod
    Use Node_mod
    Use SingleColExpBase_mod
    implicit none

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This holds together all the low-level routines for performing the
!> multiplications.
!> This particular class allows for non-identical chemical potentials,
!> in contrast to HomogeneousSingleColExp .
!--------------------------------------------------------------------
    type, extends(SingleColExpBase) :: NonHomogeneousSingleColExp
        integer :: nrofentries
        integer, allocatable :: x(:), y(:)
        complex (kind=kind(0.d0)), allocatable :: s(:), s2(:), p(:)
        real (kind=kind(0.d0)), allocatable :: c(:), c2(:) ! the cosh arrays are twice as big since we need two real values
    contains
        procedure :: init => NonHomogeneousSingleColExp_init
        procedure :: dealloc => NonHomogeneousSingleColExp_dealloc
        procedure :: vecmult => NonHomogeneousSingleColExp_vecmult
        procedure :: lmult => NonHomogeneousSingleColExp_lmult
        procedure :: lmultinv => NonHomogeneousSingleColExp_lmultinv
        procedure :: rmult => NonHomogeneousSingleColExp_rmult
        procedure :: rmultinv => NonHomogeneousSingleColExp_rmultinv
        procedure :: adjoint_over_two => NonHomogeneousSingleColExp_adjoint_over_two
        procedure :: adjointaction => NonHomogeneousSingleColExp_adjointaction
    end type

contains

subroutine NonHomogeneousSingleColExp_vecmult(this, vec)
    class(NonHomogeneousSingleColExp) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i
    complex(kind=kind(0.D0)) :: t1,t2
    do i = 1, this%nrofentries! for every matrix
        t1 = vec(this%x(i))
        t2 = vec(this%y(i))
        vec(this%x(i)) = this%c(2*i) * t1 + this%s(i) * t2
        vec(this%y(i)) = this%c(2*i+1) * t2 + conjg(this%s(i)) * t1
    enddo
end subroutine NonHomogeneousSingleColExp_vecmult

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
!
!> @param[in] this The exponential that we consider
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
subroutine NonHomogeneousSingleColExp_lmult(this, mat)
    class(NonHomogeneousSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout), contiguous :: mat
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2 ! determined to be fastest on 6x6 hubbard
    complex(kind=kind(0.D0)) :: t1(step), t2(step)
    integer, allocatable, dimension(:) :: xyarray
    complex(kind=kind(0.D0)), allocatable, dimension(:) :: snh
    real(kind=kind(0.D0)), allocatable, dimension(:) :: csh

! The intel compiler is really helped by using these temporary arrays
    allocate(xyarray(2*this%nrofentries), csh(2*this%nrofentries), snh(this%nrofentries) )
    xyarray = this%x
    csh = this%c
    snh = this%s

    ndim = size(mat,1)
    loopend = (ndim/step)*step

! ifort 2017
!DIR$ UNROLL_AND_JAM(4)
    do j = 1, loopend, step
        do i = 1, this%nrofentries! for every matrix
            do k = 1,step
                t1(k) = mat(xyarray(2*i-1), j+k-1)
                t2(k) = mat(xyarray(2*i), j+k-1)
            enddo
            do k = 1, step
                mat(xyarray(2*i-1), j+k-1) = csh(2*i) * t1(k) + snh(i) * t2(k)
                mat(xyarray(2*i), j+k-1) = csh(2*i+1) * t2(k) + conjg(snh(i)) * t1(k)
            enddo
        enddo
    enddo
    deallocate(xyarray, csh, snh)
end subroutine NonHomogeneousSingleColExp_lmult

subroutine NonHomogeneousSingleColExp_lmultinv(this, mat)
    class(NonHomogeneousSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2
    complex(kind=kind(0.D0)) :: t1(step), t2(step)
    
    ndim = size(mat,1)
    loopend = (ndim/step)*step
    do j = 1, loopend, step
        do i = 1, this%nrofentries! for every matrix
            do k = 1,step
                t1(k) = mat(this%x(2*i-1), j+k-1)
                t2(k) = mat(this%x(2*i), j+k-1)
            enddo
            do k = 1, step
                mat(this%x(2*i-1), j+k-1) = this%c(2*i) * t1(k) - this%s(i) * t2(k)
                mat(this%x(2*i), j+k-1) = this%c(2*i+1) * t2(k) - conjg(this%s(i)) * t1(k)
            enddo
        enddo
    enddo
end subroutine NonHomogeneousSingleColExp_lmultinv

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> The routines for moving to an adjoint representation : out = this.mat.this^(-1)
!> If needed we could instead calculate an eigendecomposition and use that.
!> We could really invest in a diagonal calculation at every multiplication
!> The stability of this topic has been discussed in 
!> Hargreaves, G. (2005). Topics in matrix computations: Stability and efficiency of algorithms (Doctoral dissertation, University of Manchester).
!> and "Unifying unitary and hyperbolic transformations Adam Bojanczyka, Sanzheng Qiaob;;1, Allan O. Steinhardt"
!> For the future we might want to look into fast hyperbolic rotations of Hargreaves, G. (2005).
!
!> @param[in] this The exponential that we consider
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
subroutine NonHomogeneousSingleColExp_adjointaction(this, mat)
    class(NonHomogeneousSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2
    complex(kind=kind(0.D0)) :: t1(step), t2(step)
    
    call this%lmult(mat)
    call this%rmultinv(mat)
end subroutine NonHomogeneousSingleColExp_adjointaction

subroutine NonHomogeneousSingleColExp_adjoint_over_two(this, mat)
    class(NonHomogeneousSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2
    complex(kind=kind(0.D0)) :: t1(step), t2(step), t1scal, t2scal, mys
    real(kind=kind(0.D0)) :: myc(2)
    
    ! lmult part
    ndim = size(mat,1)
    loopend = (ndim/step)*step
    do j = 1, loopend, step
        do i = 1, this%nrofentries! for every matrix
            mys = this%s2(i)
            myc(1) = this%c2(2*i)
            myc(2) = this%c2(2*i+1)
            do k = 1,step
                t1(k) = mat(this%x(2*i-1), j+k-1)
                t2(k) = mat(this%x(2*i), j+k-1)
            enddo
            do k = 1, step
                mat(this%x(2*i-1), j+k-1) = myc(1) * t1(k) + mys * t2(k)
                mat(this%x(2*i), j+k-1) = myc(2) * t2(k) + conjg(mys) * t1(k)
            enddo
        enddo
    enddo

    ! rmultinv part
    do i = 1, this%nrofentries! for every matrix
            myc(1) = this%c2(2*i)
            myc(2) = this%c2(2*i+1)
            mys = this%s2(i)
        do j = 1, ndim
            t1scal = mat(j, this%x(2*i-1))
            t2scal = mat(j, this%x(2*i))
            mat(j, this%x(2*i-1)) = myc(1) * t1scal - mys * t2scal
            mat(j, this%x(2*i)) = myc(2) * t2scal - conjg(mys) * t1scal
        enddo
    enddo
end subroutine NonHomogeneousSingleColExp_adjoint_over_two

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
subroutine NonHomogeneousSingleColExp_rmult(this, mat)
    class(NonHomogeneousSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, k, ndim
    complex(kind=kind(0.D0)) :: t1, t2
    
    ndim = size(mat,1)
    do i = 1, this%nrofentries! for every matrix
        do j = 1, ndim
        t1 = mat(j, this%x(2*i-1))
        t2 = mat(j, this%x(2*i))
        mat(j, this%x(2*i-1)) = this%c(2*i) * t1 + this%s(i)* t2
        mat(j, this%x(2*i)) = this%c(2*i+1) * t2 + conjg(this%s(i))* t1
        enddo
    enddo
end subroutine NonHomogeneousSingleColExp_rmult

subroutine NonHomogeneousSingleColExp_rmultinv(this, mat)
    class(NonHomogeneousSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, k, ndim
    complex(kind=kind(0.D0)) :: t1, t2
    
    ndim = size(mat,1)
    do i = 1, this%nrofentries! for every matrix
        do j = 1, ndim
        t1 = mat(j, this%x(2*i-1))
        t2 = mat(j, this%x(2*i))
        mat(j, this%x(2*i-1)) = this%c(2*i) * t1 - this%s(i) * t2
        mat(j, this%x(2*i)) = this%c(2*i+1) * t2 - conjg(this%s(i)) * t1
        enddo
    enddo
end subroutine NonHomogeneousSingleColExp_rmultinv

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This sets up the data to perform the exponentiation of a 
!> strictly sparse matrix.
!> The internal layout is that the non-zero element a_xy stored at (x,y) in the matrix
!> has x(i) at x(2i-1) and y(i) at x(2i).
!> the two values for the diagonals are stored in c(2i) and c(2i+1)
!
!> @param[inout] this the NonHomogeneousSingleColExp object.
!> @param[in] nodes The nodes that belng to this color.
!> @param[in] nredges how many nodes of this color.
!> @param[in] weight a prefactor for the exponent.
!--------------------------------------------------------------------
subroutine NonHomogeneousSingleColExp_init(this, nodes, nredges, mys, weight)
    class(NonHomogeneousSingleColExp), intent(inout) :: this
    type(node), dimension(:), intent(in) :: nodes
    real(kind=kind(0.D0)), intent(in), allocatable, dimension(:) :: mys
    integer, intent(in) :: nredges
    real (kind=kind(0.d0)), intent(in) :: weight
    integer :: i
    real (kind=kind(0.d0)) :: nf, my1, my2, localzero
    allocate(this%x(2*nredges), this%y(nredges), this%c(2*nredges), this%s(nredges))
    allocate(this%c2(2*nredges), this%s2(nredges), this%p(nredges))
    this%nrofentries = nredges
#ifndef NDEBUG
    write(*,*) "Setting up strict. sparse matrix with ", nredges, "edges"
#endif
    do i = 1, nredges
        this%x(2*i-1) = nodes(i)%x
        this%x(2*i) = nodes(i)%y
        this%y(i) = nodes(i)%y
        this%p(i) = weight*nodes(i)%axy
        !calculate Frobenius norm
        my1 = mys(nodes(i)%x)
        my2 = mys(nodes(i)%y)
        nf = sqrt(my1*my1+my2*my2 + 2*dble(this%p(i) * conjg(this%p(i))))
        localzero = 1E-15*nf ! definition of my local scale that defines zero
        if (abs(my1-my2) < localzero) then
            write(*,*) "[NonHomogeneousSingleColExp_init]: Identical diagonals found! This should go in another class!"
        endif
        ! This is the order of operations that yields stable matrix inversions
        ! We assume that the matrix that we have decomposed is hermitian:
        ! M=(d  , b)
        !   (b^*, -d) then the below entries follow for the exponential and cosh is real.
        ! The fixup for generic chemical potentials happens at the end.
        this%c(i) = cosh(abs(weight*nodes(i)%axy))
        this%c2(i) = cosh(abs(weight*nodes(i)%axy)/2)
        ! I got the most reliable results if the hyperbolic pythagoras is best fulfilled.
        ! If we generalize this, to non-zero diagonals, this means 
        this%s(i) = sqrt(this%c(i)**2-1.0)*weight*nodes(i)%axy/abs(weight*nodes(i)%axy)
        this%s2(i) = sqrt(this%c2(i)**2-1.0)*weight*nodes(i)%axy/abs(weight*nodes(i)%axy)
    enddo
! All nodes that we have been passed are now from a single color.
! They constitute now a strictly sparse matrix.
! Further processing of the entries could be done here.
end subroutine NonHomogeneousSingleColExp_init

subroutine NonHomogeneousSingleColExp_dealloc(this)
    class(NonHomogeneousSingleColExp), intent(inout) :: this
    deallocate(this%x, this%y, this%c, this%s, this%c2, this%s2)
end subroutine NonHomogeneousSingleColExp_dealloc

end module NonHomogeneousSingleColExp_mod
