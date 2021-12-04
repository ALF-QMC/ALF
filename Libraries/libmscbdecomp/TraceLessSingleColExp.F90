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

module TraceLessSingleColExp_mod
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
!> This particular class is specialized to the case that each
!> 2x2 block is a traceless matrix.
!--------------------------------------------------------------------
    type, extends(SingleColExpBase) :: TraceLessSingleColExp
    contains
        procedure :: init => TraceLessSingleColExp_init
        procedure :: dealloc => TraceLessSingleColExp_dealloc
        procedure :: vecmult => TraceLessSingleColExp_vecmult
        procedure :: lmult => TraceLessSingleColExp_lmult
        procedure :: lmultinv => TraceLessSingleColExp_lmultinv
        procedure :: rmult => TraceLessSingleColExp_rmult
        procedure :: rmultinv => TraceLessSingleColExp_rmultinv
        procedure :: adjoint_over_two => TraceLessSingleColExp_adjoint_over_two
        procedure :: adjointaction => TraceLessSingleColExp_adjointaction
    end type

contains

subroutine TraceLessSingleColExp_vecmult(this, vec)
    class(TraceLessSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:), intent(inout) :: vec
    integer :: i
    complex(kind=kind(0.D0)) :: t1,t2
    do i = 1, this%nrofentries! for every matrix
        t1 = vec(this%xy(2*i))
        t2 = vec(this%xy(2*i-1))
        vec(this%xy(2*i)) = this%c(2*i-1) * t1 + this%s(i) * t2
        vec(this%xy(2*i-1)) = this%c(2*i) * t2 + conjg(this%s(i)) * t1
    enddo
end subroutine TraceLessSingleColExp_vecmult

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
!> This is an internal helper function that finds reuse in multiple places.
!
!> @param[in] c the diagonal data
!> @param[in] s the off-diagonal data
!> @param[in] x the used matrix positions
!> @param[in] nrofentries how many vertices are in this family.
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------

pure subroutine lmultthreeelementbase(c, s, x, nrofentries, mat)
    real (kind=kind(0.d0)), allocatable, intent(in) :: c(:)
    complex (kind=kind(0.d0)), allocatable, intent(in) :: s(:)
    integer, allocatable, intent(in) :: x(:)
    integer, intent(in) ::nrofentries
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout), contiguous :: mat
    
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2 ! determined to be fastest on 6x6 hubbard
    complex(kind=kind(0.D0)) :: t1(step), t2(step)
    integer, allocatable, dimension(:) :: xyarray
    complex(kind=kind(0.D0)), allocatable, dimension(:) :: snh
    real(kind=kind(0.D0)), allocatable, dimension(:) :: csh

! The intel compiler is really helped by using these temporary arrays
    allocate(xyarray(nrofentries), csh(nrofentries), snh(nrofentries) )
    xyarray = x
    csh = c
    snh = s

    ndim = size(mat,1)
    loopend = (ndim/step)*step

! ifort 2017
!DIR$ UNROLL_AND_JAM(4)
    do j = 1, loopend, step
        do i = 1, nrofentries! for every matrix
            do k = 1,step
                t1(k) = mat(xyarray(2*i-1), j+k-1)
                t2(k) = mat(xyarray(2*i), j+k-1)
            enddo
            do k = 1, step
                mat(xyarray(2*i-1), j+k-1) = csh(2*i-1) * t1(k) + snh(i) * t2(k)
                mat(xyarray(2*i), j+k-1) = csh(2*i) * t2(k) + conjg(snh(i)) * t1(k)
            enddo
        enddo
    enddo
    
    ! remainder loop
    if ((ndim - loopend) .ne. 0) then
        do i = 1, nrofentries! for every matrix
            t1(1) = mat(xyarray(2*i-1), ndim)
            t2(1) = mat(xyarray(2*i), ndim)
            mat(xyarray(2*i-1), ndim) = csh(2*i-1) * t1(1) + snh(i) * t2(1)
            mat(xyarray(2*i), ndim) = csh(2*i) * t2(1) + conjg(snh(i)) * t1(1)
        enddo
    endif
    deallocate(xyarray, csh, snh)
end subroutine

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
subroutine TraceLessSingleColExp_lmult(this, mat)
    class(TraceLessSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout), contiguous :: mat

    call lmultthreeelementbase(this%c, this%s, this%xy, this%nrofentries, mat)
end subroutine TraceLessSingleColExp_lmult

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Perform the multiplication of this inverted exponential with a matrix:
!>  out = this*mat
!
!> @param[in] this The exponential that we consider
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
subroutine TraceLessSingleColExp_lmultinv(this, mat)
    class(TraceLessSingleColExp), intent(in) :: this
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
                mat(this%xy(2*i-1), j+k-1) = this%c(2*i) * t1(k) - this%s(i) * t2(k)
                mat(this%xy(2*i), j+k-1) = this%c(2*i-1) * t2(k) - conjg(this%s(i)) * t1(k)
            enddo
        enddo
    enddo
    
    ! remainder loop
    if ((ndim - loopend) .ne. 0) then
        do i = 1, this%nrofentries! for every matrix
            t1(1) = mat(this%xy(2*i-1), ndim)
            t2(1) = mat(this%xy(2*i), ndim)
            mat(this%xy(2*i-1), ndim) = this%c(2*i) * t1(1) - this%s(i) * t2(1)
            mat(this%xy(2*i), ndim) = this%c(2*i-1) * t2(1) - conjg(this%s(i)) * t1(1)
        enddo
    endif
end subroutine TraceLessSingleColExp_lmultinv

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
subroutine TraceLessSingleColExp_adjointaction(this, mat)
    class(TraceLessSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    
    call this%lmult(mat)
    call this%rmultinv(mat)
end subroutine TraceLessSingleColExp_adjointaction

subroutine TraceLessSingleColExp_adjoint_over_two(this, mat)
    class(TraceLessSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2
    complex(kind=kind(0.D0)) :: t1(step), t2(step), t1scal, t2scal, mys

    ! lmult part
    ndim = size(mat,1)
    loopend = (ndim/step)*step
    do j = 1, loopend, step
        do i = 1, this%nrofentries! for every matrix
            mys = this%s2(i)
            do k = 1,step
                t1(k) = mat(this%xy(2*i-1), j+k-1)
                t2(k) = mat(this%xy(2*i), j+k-1)
            enddo
            do k = 1, step
                mat(this%xy(2*i-1), j+k-1) = this%c2(2*i-1) * t1(k) + mys * t2(k)
                mat(this%xy(2*i), j+k-1) = this%c2(2*i) * t2(k) + conjg(mys) * t1(k)
            enddo
        enddo
    enddo

    ! remainder loop
    if ((ndim - loopend) .ne. 0) then
        do i = 1, this%nrofentries! for every matrix
            t1(1) = mat(this%xy(2*i-1), ndim)
            t2(1) = mat(this%xy(2*i), ndim)
            mat(this%xy(2*i-1), ndim) = this%c2(2*i-1) * t1(1) + this%s2(i) * t2(1)
            mat(this%xy(2*i), ndim) = this%c2(2*i) * t2(1) + conjg(this%s2(i)) * t1(1)
        enddo
    endif

    ! rmultinv part
    do i = 1, this%nrofentries! for every matrix
        mys = this%s2(i)
        do j = 1, ndim
            t1scal = mat(j, this%xy(2*i-1))
            t2scal = mat(j, this%xy(2*i))
            mat(j, this%xy(2*i-1)) = this%c2(2*i) * t1scal - mys * t2scal
            mat(j, this%xy(2*i)) = this%c2(2*i-1) * t2scal - conjg(mys) * t1scal
        enddo
    enddo
end subroutine TraceLessSingleColExp_adjoint_over_two

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Perform the multiplication with a matrix.
!> This version requires three data elements
!> This is an internal helper function that finds reuse in multiple places.
!
!> @param[in] c the diagonal data
!> @param[in] s the off-diagonal data
!> @param[in] x the used matrix positions
!> @param[in] nrofentries how many vertices are in this family.
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
pure subroutine rmultthreeelementbase(c, s, x, nrofentries, mat)
    real (kind=kind(0.d0)), allocatable, intent(in) :: c(:)
    complex (kind=kind(0.d0)), allocatable, intent(in) :: s(:)
    integer, allocatable, intent(in) :: x(:)
    integer, intent(in) ::nrofentries
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, ndim
    complex(kind=kind(0.D0)) :: t1, t2

    ndim = size(mat,1)
    do i = 1, nrofentries! for every matrix
        do j = 1, ndim
        t1 = mat(j, x(2*i-1))
        t2 = mat(j, x(2*i))
        mat(j, x(2*i-1)) = c(2*i-1) * t1 + s(i)* t2
        mat(j, x(2*i)) = c(2*i) * t2 + conjg(s(i))* t1
        enddo
    enddo
end subroutine

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
subroutine TraceLessSingleColExp_rmult(this, mat)
    class(TraceLessSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    
    call rmultthreeelementbase(this%c, this%s, this%xy, this%nrofentries, mat)
end subroutine TraceLessSingleColExp_rmult

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> multiplication with inverse.
!> We know that Inverse of
!> M= ( a, x  )   ( b, -x )
!>    (x^*, b) =  ( -x^*, a)
!> if det(M) = 1.
!
!> @param[in] this The exponential that we consider
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
subroutine TraceLessSingleColExp_rmultinv(this, mat)
    class(TraceLessSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    integer :: i, j, ndim
    complex(kind=kind(0.D0)) :: t1, t2
    
    ndim = size(mat,1)
    do i = 1, this%nrofentries! for every matrix
        do j = 1, ndim
        t1 = mat(j, this%xy(2*i-1))
        t2 = mat(j, this%xy(2*i))
        mat(j, this%xy(2*i-1)) = this%c(2*i) * t1 - this%s(i) * t2
        mat(j, this%xy(2*i)) = this%c(2*i-1) * t2 - conjg(this%s(i)) * t1
        enddo
    enddo
end subroutine TraceLessSingleColExp_rmultinv

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This calculates the input data of a checkerboard matrix, hence
!> the entries of C=exp({{d,o},{o^*,-d}})
!> While best preserving det(C) = 1.
!
!> @param [out] diag1 first diagonal
!> @param [out] diag2 second diagonal
!> @param [out] offout resulting off-diagonal
!> @param [in] diag diagonal d
!> @param [in] offinp The off-diagonal o
!> @param [in] weight a real prefactor that can be in front of the entire 2x2 sub-block.
!--------------------------------------------------------------------
subroutine expof2x2tracelesshermitianmatrix(diag1, diag2, offout, diag, offinp, weight)
    real (kind=kind(0.d0)), intent(in) :: weight, diag
    real (kind=kind(0.d0)), intent(out) :: diag1, diag2
    complex(kind=kind(0.D0)), intent(in) :: offinp
    complex(kind=kind(0.D0)), intent(out) :: offout
    real (kind=kind(0.d0)) :: sinhlocal, angle

    angle = sqrt(diag*diag + DBLE(offinp * conjg(offinp)))
    sinhlocal = sinh(weight*angle) ! weight does not seem to get an abs() here.
    diag1 = sqrt(1.0 + sinhlocal**2) ! solve for cosh
    diag2 = diag1
    diag1 = diag1 + diag*sinhlocal/angle
    diag2 = diag2 - diag*sinhlocal/angle
    offout = weight*offinp/abs(weight*offinp)*sqrt(diag1*diag2 - 1.0)! rescale offdiagonal such that det(exp(C)) == 1
end subroutine


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
!> @param[inout] this the TraceLessSingleColExp object.
!> @param[in] nodes The nodes that belng to this color.
!> @param[in] nredges how many nodes of this color.
!> @param[in] weight a prefactor for the exponent.
!--------------------------------------------------------------------
subroutine TraceLessSingleColExp_init(this, nodes, nredges, mys, weight)
    class(TraceLessSingleColExp), intent(inout) :: this
    type(node), dimension(:), intent(in) :: nodes
    real(kind=kind(0.D0)), intent(in), dimension(:) :: mys
    integer, intent(in) :: nredges
    real (kind=kind(0.d0)), intent(in) :: weight
    integer :: i
    real (kind=kind(0.d0)) :: my1, my2, localzero, tmp
    ! We need twice the amount of storage for the diagonal for this traceless case.
    allocate(this%xy(2*nredges), this%s(nredges), this%c(2*nredges), this%c2(2*nredges), this%s2(nredges))
    this%nrofentries = nredges
#ifndef NDEBUG
    write(*,*) "[TraceLessSingleColExp] Setting up strict. sparse matrix with ", nredges, "edges"
#endif
    do i = 1, nredges
        this%xy(2*i-1) = nodes(i)%x
        this%xy(2*i) = nodes(i)%y
        !calculate Frobenius norm
        my1 = mys(nodes(i)%x)
        my2 = mys(nodes(i)%y)
        ! dependence on weight cancels in all comparisons.
        localzero = 1E-15*frobnorm(my1, my2, nodes(i)%axy) ! definition of my local scale that defines zero
        if (abs(my1+my2) > localzero) then
            write(*,*) "[TraceLessSingleColExp_init]: Matrix not traceless. This should not happen here."
            error stop 1
        endif
        ! This is the order of operations that yields stable matrix inversions
        ! We assume that the matrix that we have decomposed is hermitian:
        ! M=(d  , b)
        !   (b^*, -d)
        ! with a real d.
        ! then the below entries follow for the exponential and cosh is real.
        ! more general chemical potentials are deferred to different classes
        
        call expof2x2tracelesshermitianmatrix(this%c(2*i-1), this%c(2*i), this%s(i), my1, nodes(i)%axy, weight)
        tmp = weight/2.D0
        call expof2x2tracelesshermitianmatrix(this%c2(2*i-1), this%c2(2*i), this%s2(i), my1, nodes(i)%axy, tmp)
    enddo
! All nodes that we have been passed are now from a single color.
! They constitute now a strictly sparse matrix.
end subroutine TraceLessSingleColExp_init

subroutine TraceLessSingleColExp_dealloc(this)
    class(TraceLessSingleColExp), intent(inout) :: this
    deallocate(this%xy, this%c, this%s, this%c2, this%s2)
end subroutine TraceLessSingleColExp_dealloc

end module TraceLessSingleColExp_mod
