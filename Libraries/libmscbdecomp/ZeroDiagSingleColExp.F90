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

module ZeroDiagSingleColExp_mod
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
!> This particular class is specialized to the case that in this
!> particular color all chemical potentials vanish.
!--------------------------------------------------------------------
    type, extends(SingleColExpBase) :: ZeroDiagSingleColExp
    contains
        procedure :: init => ZeroDiagSingleColExp_init
        procedure :: dealloc => ZeroDiagSingleColExp_dealloc
        procedure :: vecmult => ZeroDiagSingleColExp_vecmult
        procedure :: lmult => ZeroDiagSingleColExp_lmult
        procedure :: lmultinv => ZeroDiagSingleColExp_lmultinv
        procedure :: rmult => ZeroDiagSingleColExp_rmult
        procedure :: rmultinv => ZeroDiagSingleColExp_rmultinv
        procedure :: adjoint_over_two => ZeroDiagSingleColExp_adjoint_over_two
        procedure :: adjointaction => ZeroDiagSingleColExp_adjointaction
    end type

contains

subroutine ZeroDiagSingleColExp_vecmult(this, vec)
    class(ZeroDiagSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:),intent(inout) :: vec
    integer :: i
    complex(kind=kind(0.D0)) :: t1, t2
    do i = 1, this%nrofentries! for every matrix
        t1 = vec(this%xy(2*i-1))
        t2 = vec(this%xy(2*i))
        vec(this%xy(2*i-1)) = this%c(i) * t1 + this%s(i) * t2
        vec(this%xy(2*i)) = this%c(i) * t2 + conjg(this%s(i)) * t1
    enddo
end subroutine ZeroDiagSingleColExp_vecmult

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
subroutine ZeroDiagSingleColExp_lmult(this, mat)
    class(ZeroDiagSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout), contiguous :: mat

    call lmultbase(this%c, this%s, this%xy, this%nrofentries, mat)
end subroutine ZeroDiagSingleColExp_lmult

subroutine ZeroDiagSingleColExp_lmultinv(this, mat)
    class(ZeroDiagSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout), contiguous :: mat
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2
    complex(kind=kind(0.D0)) :: t1(step), t2(step)

    ndim = size(mat,1)
    loopend = (ndim/step)*step
    do j = 1, loopend, step
        do i = 1, this%nrofentries! for every matrix
            do k = 1,step
                t1(k) = mat(this%xy(2*i-1), j+k-1)
                t2(k) = mat(this%xy(2*i), j+k-1)
            enddo
            do k = 1, step
                mat(this%xy(2*i-1), j+k-1) = this%c(i) * t1(k) - this%s(i) * t2(k)
                mat(this%xy(2*i), j+k-1) = this%c(i) * t2(k) - conjg(this%s(i)) * t1(k)
            enddo
        enddo
    enddo

    ! remainder loop
    if ((ndim - loopend) .ne. 0) then
        do i = 1, this%nrofentries! for every matrix
            t1(1) = mat(this%xy(2*i-1), ndim)
            t2(1) = mat(this%xy(2*i), ndim)
            mat(this%xy(2*i-1), ndim) = this%c(i) * t1(1) - this%s(i) * t2(1)
            mat(this%xy(2*i), ndim) = this%c(i) * t2(1) - conjg(this%s(i)) * t1(1)
        enddo
    endif
end subroutine ZeroDiagSingleColExp_lmultinv

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
!> @param[in] this The exponential that we consider.
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
subroutine ZeroDiagSingleColExp_adjointaction(this, mat)
    class(ZeroDiagSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    
    call this%lmult(mat)
    call this%rmultinv(mat)
end subroutine ZeroDiagSingleColExp_adjointaction

subroutine ZeroDiagSingleColExp_adjoint_over_two(this, mat)
    class(ZeroDiagSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, ndim
    complex(kind=kind(0.D0)) :: t1scal, t2scal, mys
    real(kind=kind(0.D0)) :: myc
    
    ! lmult part
    ndim = size(mat,1)
    do i = 1, this%nrofentries! for every matrix
        ! get data
        mys = this%s2(i)
        myc = this%c2(i)
        do j = 1, ndim
                t1scal = mat(this%xy(2*i-1), j)
                t2scal = mat(this%xy(2*i), j)
                mat(this%xy(2*i-1), j) = myc * t1scal + mys * t2scal
                mat(this%xy(2*i), j) = myc * t2scal + conjg(mys) * t1scal
        enddo
!     enddo
! 
!     ! rmultinv part
!     do i = 1, this%nrofentries! for every matrix
        do j = 1, ndim
            t1scal = mat(j, this%xy(2*i-1))
            t2scal = mat(j, this%xy(2*i))
            mat(j, this%xy(2*i-1)) = myc * t1scal - mys * t2scal
            mat(j, this%xy(2*i)) = myc * t2scal - conjg(mys) * t1scal
        enddo
    enddo
end subroutine ZeroDiagSingleColExp_adjoint_over_two

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
subroutine ZeroDiagSingleColExp_rmult(this, mat)
    class(ZeroDiagSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat

    call rmultbase(this%c, this%s, this%xy, this%nrofentries, mat)
end subroutine ZeroDiagSingleColExp_rmult

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Perform the multiplication of this inverted exponential with a matrix: out = mat*(this)^-1
!
!> @param[in] this The exponential that we consider.
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
subroutine ZeroDiagSingleColExp_rmultinv(this, mat)
    class(ZeroDiagSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, ndim
    complex(kind=kind(0.D0)) :: t1, t2
    
    ndim = size(mat,1)
    do i = 1, this%nrofentries! for every matrix
        do j = 1, ndim
        t1 = mat(j, this%xy(2*i-1))
        t2 = mat(j, this%xy(2*i))
        mat(j, this%xy(2*i-1)) = this%c(i) * t1 - this%s(i) * t2
        mat(j, this%xy(2*i)  ) = this%c(i) * t2 - conjg(this%s(i)) * t1
        enddo
    enddo
end subroutine ZeroDiagSingleColExp_rmultinv

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
!> @param[inout] this the ZeroDiagSingleColExp object.
!> @param[in] nodes The nodes that belng to this color.
!> @param[in] nredges how many nodes of this color.
!> @param[in] weight a prefactor for the exponent.
!--------------------------------------------------------------------
subroutine ZeroDiagSingleColExp_init(this, nodes, nredges, mys, weight)
    implicit none
    class(ZeroDiagSingleColExp), intent(inout) :: this
    type(node), dimension(:), intent(in) :: nodes
    real(kind=kind(0.D0)), intent(in), dimension(:) :: mys
    integer, intent(in) :: nredges
    real (kind=kind(0.d0)), intent(in) :: weight
    integer :: i
    real (kind=kind(0.d0)) :: my1, my2, localzero
    allocate(this%xy(2*nredges), this%c(nredges), this%s(nredges), this%c2(nredges), this%s2(nredges))
    this%nrofentries = nredges
#ifndef NDEBUG
    write(*,*) "[ZeroDiagSingleColExp] Setting up strict. sparse matrix with ", nredges, "edges"
#endif
    do i = 1, nredges
        this%xy(2*i-1) = nodes(i)%x
        this%xy(2*i) = nodes(i)%y

        !calculate Frobenius norm
        my1 = mys(nodes(i)%x)
        my2 = mys(nodes(i)%y)
        ! dependence on weight cancels in all comparisons.
        localzero = 1E-15*frobnorm(my1, my2, nodes(i)%axy) ! definition of my local scale that defines zero
        if ((abs(my1) > localzero) .or. (abs(my2) > localzero)) then
            write(*,*) "[ZeroDiagSingleColExp_init]: Diagonal NOT zero. This should not happen here."
            error stop 1
        endif
        ! This is the order of operations that yields stable matrix inversions
        ! We assume that the matrix that we have decomposed is hermitian:
        ! M=(0  , b)
        !   (b^*, 0) then the below entries follow for the exponential and cosh is real.
        ! chemical potentials are deferred to different classes
        
        this%c(i) = cosh(abs(weight*nodes(i)%axy))
        this%c2(i) = cosh(abs(weight*nodes(i)%axy)/2.D0)
        if (abs(weight*nodes(i)%axy) > 0.D0) then ! weight required here, since it could be zero
            ! I got the most reliable results if the hyperbolic pythagoras is best fulfilled.
            this%s(i) = sqrt(this%c(i)**2-1.D0)*(weight*nodes(i)%axy/abs(weight*nodes(i)%axy))
            this%s2(i) = sqrt(this%c2(i)**2-1.D0)*(weight*nodes(i)%axy/abs(weight*nodes(i)%axy))
        else
            this%s(i) = 0
            this%s2(i) = 0
        endif
    enddo
! All nodes that we have been passed are now from a single color.
! They constitute now a strictly sparse matrix.
end subroutine ZeroDiagSingleColExp_init

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Deallocate our memory.
!
!> @param[in] this The exponential that we consider.
!--------------------------------------------------------------------
subroutine ZeroDiagSingleColExp_dealloc(this)
    class(ZeroDiagSingleColExp), intent(inout) :: this

    deallocate(this%xy, this%c, this%s, this%c2, this%s2)
end subroutine ZeroDiagSingleColExp_dealloc

end module ZeroDiagSingleColExp_mod
