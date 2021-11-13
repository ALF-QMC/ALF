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

module HomogeneousSingleColExp_mod
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
!> This particular class is specialized to the case that in each
!> checkerboard matrix, the chemical potentials are identical.

!--------------------------------------------------------------------
    type, extends(SingleColExpBase) :: HomogeneousSingleColExp
!         integer :: nrofentries
!         integer, allocatable :: x(:), y(:) ! the y array is still around but often(?) unused
        complex (kind=kind(0.d0)), allocatable :: s2(:), p(:)
        real (kind=kind(0.d0)), allocatable :: c2(:)
    contains
        procedure :: init => HomogeneousSingleColExp_init
        procedure :: dealloc => HomogeneousSingleColExp_dealloc
        procedure :: vecmult => HomogeneousSingleColExp_vecmult
        procedure :: lmult => HomogeneousSingleColExp_lmult
        procedure :: lmultinv => HomogeneousSingleColExp_lmultinv
        procedure :: rmult => HomogeneousSingleColExp_rmult
        procedure :: rmultinv => HomogeneousSingleColExp_rmultinv
        procedure :: adjoint_over_two => HomogeneousSingleColExp_adjoint_over_two
        procedure :: adjointaction => HomogeneousSingleColExp_adjointaction
    end type

contains

subroutine HomogeneousSingleColExp_vecmult(this, vec)
    class(HomogeneousSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:) :: vec
    integer :: i
    complex(kind=kind(0.D0)) :: t1,t2
    do i = 1, this%nrofentries! for every matrix
        t1 = vec(this%x(i))
        t2 = vec(this%y(i))
        vec(this%x(i)) = this%c(i) * t1 + this%s(i) * t2
        vec(this%y(i)) = this%c(i) * t2 + conjg(this%s(i)) * t1
    enddo
end subroutine HomogeneousSingleColExp_vecmult

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
subroutine HomogeneousSingleColExp_lmult(this, mat)
    class(HomogeneousSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout), contiguous :: mat
    
    call lmultbase(this%c, this%s, this%x, this%nrofentries, mat)
end subroutine HomogeneousSingleColExp_lmult

subroutine HomogeneousSingleColExp_lmultinv(this, mat)
    class(HomogeneousSingleColExp), intent(in) :: this
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
                mat(this%x(2*i-1), j+k-1) = this%c(i) * t1(k) - this%s(i) * t2(k)
                mat(this%x(2*i), j+k-1) = this%c(i) * t2(k) - conjg(this%s(i)) * t1(k)
            enddo
        enddo
    enddo
end subroutine HomogeneousSingleColExp_lmultinv

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
subroutine HomogeneousSingleColExp_adjointaction(this, mat)
    class(HomogeneousSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    
    call this%lmult(mat)
    call this%rmultinv(mat)
end subroutine HomogeneousSingleColExp_adjointaction

subroutine HomogeneousSingleColExp_adjoint_over_two(this, mat)
    class(HomogeneousSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2
    complex(kind=kind(0.D0)) :: t1(step), t2(step), t1scal, t2scal, mys
    real(kind=kind(0.D0)) :: myc
    
    ! lmult part
    ndim = size(mat,1)
    loopend = (ndim/step)*step
    do j = 1, loopend, step
        do i = 1, this%nrofentries! for every matrix
            mys = this%s2(i)
            myc = this%c2(i)
            do k = 1,step
                t1(k) = mat(this%x(2*i-1), j+k-1)
                t2(k) = mat(this%x(2*i), j+k-1)
            enddo
            do k = 1, step
                mat(this%x(2*i-1), j+k-1) = myc * t1(k) + mys * t2(k)
                mat(this%x(2*i), j+k-1) = myc * t2(k) + conjg(mys) * t1(k)
            enddo
        enddo
    enddo

    ! rmultinv part
    do i = 1, this%nrofentries! for every matrix
            myc = this%c2(i)
            mys = this%s2(i)
        do j = 1, ndim
            t1scal = mat(j, this%x(2*i-1))
            t2scal = mat(j, this%x(2*i))
            mat(j, this%x(2*i-1)) = myc * t1scal - mys * t2scal
            mat(j, this%x(2*i)) = myc * t2scal - conjg(mys) * t1scal
        enddo
    enddo
end subroutine HomogeneousSingleColExp_adjoint_over_two

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
subroutine HomogeneousSingleColExp_rmult(this, mat)
    class(HomogeneousSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    
    call rmultbase(this%c, this%s, this%x, this%nrofentries, mat)
end subroutine HomogeneousSingleColExp_rmult

subroutine HomogeneousSingleColExp_rmultinv(this, mat)
    class(HomogeneousSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, ndim
    complex(kind=kind(0.D0)) :: t1, t2
    
    ndim = size(mat,1)
    do i = 1, this%nrofentries! for every matrix
        do j = 1, ndim
        t1 = mat(j, this%x(2*i-1))
        t2 = mat(j, this%x(2*i))
        mat(j, this%x(2*i-1)) = this%c(i) * t1 - this%s(i) * t2
        mat(j, this%x(2*i)) = this%c(i) * t2 - conjg(this%s(i)) * t1
        enddo
    enddo
end subroutine HomogeneousSingleColExp_rmultinv

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This sets up the data to perform the exponentiation of a 
!> strictly sparse matrix.
!> The internal layout is that the non-zero element a_xy stored at (x,y) in the matrix
!> has x(i) at x(2i-1) and y(i) at x(2i)
!
!> @param[inout] this the HomogeneousSingleColExp object.
!> @param[in] nodes The nodes that belng to this color.
!> @param[in] nredges how many nodes of this color.
!> @param[in] weight a prefactor for the exponent.
!--------------------------------------------------------------------
subroutine HomogeneousSingleColExp_init(this, nodes, nredges, mys, weight)
    class(HomogeneousSingleColExp), intent(inout) :: this
    type(node), dimension(:), intent(in) :: nodes
    real(kind=kind(0.D0)), intent(in), allocatable, dimension(:) :: mys
    integer, intent(in) :: nredges
    real (kind=kind(0.d0)), intent(in) :: weight
    integer :: i
    real (kind=kind(0.d0)) :: nf, my1, my2, localzero
    allocate(this%x(2*nredges), this%y(nredges), this%c(nredges), this%s(nredges))
    allocate(this%c2(nredges), this%s2(nredges), this%p(nredges))
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
        if (abs(my1-my2) > localzero) then
            write(*,*) "[HomogeneousSingleColExp_init]: Unequal diagonals found. This should not happen here."
            error stop 1
        endif
        ! This is the order of operations that yields stable matrix inversions
        ! We assume that the matrix that we have decomposed is hermitian:
        ! M=(0  , b)
        !   (b^*, 0) then the below entries follow for the exponential and cosh is real.
        ! The case of a uniform chemical potential is fixed up later.
        this%c(i) = cosh(abs(weight*nodes(i)%axy))
        this%c2(i) = cosh(abs(weight*nodes(i)%axy)/2)
        ! I got the most reliable results if the hyperbolic pythagoras is best fulfilled.
        ! If we generalize this, to non-zero diagonals, this means 
        this%s(i) = sqrt(this%c(i)**2-1.0)*weight*nodes(i)%axy/abs(weight*nodes(i)%axy)
        this%s2(i) = sqrt(this%c2(i)**2-1.0)*weight*nodes(i)%axy/abs(weight*nodes(i)%axy)
        if (abs(my1+my2) > 2*localzero) then ! chemical potential is actually different from zero
            this%c(i) = this%c(i) * exp(my1)
            this%c2(i) = this%c2(i) * exp(my1/2)
            this%s(i) = this%s(i) * exp(my1)
            this%s2(i) = this%s2(i) * exp(my1/2)
        endif
    enddo
! All nodes that we have been passed are now from a single color.
! They constitute now a strictly sparse matrix.
! Further processing of the entries could be done here.
end subroutine HomogeneousSingleColExp_init

subroutine HomogeneousSingleColExp_dealloc(this)
    class(HomogeneousSingleColExp), intent(inout) :: this
    deallocate(this%x, this%y, this%c, this%s, this%c2, this%s2)
end subroutine HomogeneousSingleColExp_dealloc

end module HomogeneousSingleColExp_mod
