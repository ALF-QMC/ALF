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

module GeneralSingleColExp_mod
    Use Node_mod
    Use TraceLessSingleColExp_mod
    implicit none

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @class GeneralSingleColExp
!> @brief This holds together all the low-level routines for performing the multiplications.
!>
!> This particular class allows for non-identical chemical potentials,
!> in contrast to HomogeneousSingleColExp and hence is the most general
!> implementation.
!--------------------------------------------------------------------
    type, extends(TraceLessSingleColExp) :: GeneralSingleColExp
        complex (kind=kind(0.d0)), allocatable :: sinv(:), s2inv(:) ! data for storing the sinh
        real (kind=kind(0.d0)), allocatable :: cinv(:), c2inv(:)! the cosh array is twice as big since we need two real values
    contains
        procedure :: init => GeneralSingleColExp_init
        procedure :: dealloc => GeneralSingleColExp_dealloc
        procedure :: lmultinv => GeneralSingleColExp_lmultinv
        procedure :: rmultinv => GeneralSingleColExp_rmultinv
        procedure :: adjoint_over_two => GeneralSingleColExp_adjoint_over_two
        procedure :: adjointaction => GeneralSingleColExp_adjointaction
    end type

contains

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Perform the multiplication of the inverse of this 
!> exponential with a matrix: out = this*mat
!
!> @param[in] this The exponential that we consider
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
subroutine GeneralSingleColExp_lmultinv(this, mat)
    class(GeneralSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout), contiguous :: mat
    
    call lmultthreeelementbase(this%cinv, this%sinv, this%xy, this%nrofentries, mat)
end subroutine GeneralSingleColExp_lmultinv

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
!> @param[inout] mat The matrix that we modify.
!--------------------------------------------------------------------
subroutine GeneralSingleColExp_adjointaction(this, mat)
    class(GeneralSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat
    
    call this%lmult(mat)
    call this%rmultinv(mat)
end subroutine GeneralSingleColExp_adjointaction

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> The routines for moving to an adjoint representation : out = this.mat.this^(-1),
!> but this time with a halved prefactor.
!
!> @param[in] this The exponential that we consider.
!> @param[inout] mat the matrix that we modify.
!--------------------------------------------------------------------
subroutine GeneralSingleColExp_adjoint_over_two(this, mat)
    class(GeneralSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat ! input matrix
    integer :: i, j, k, ndim, loopend
    integer, parameter :: step = 2
    complex(kind=kind(0.D0)) :: t1(step), t2(step), mys
    real(kind=kind(0.D0)) :: myc(2)
    
    ! lmult part
    ndim = size(mat,1)
    loopend = (ndim/step)*step
    do j = 1, loopend, step
        do i = 1, this%nrofentries! for every matrix
            mys = this%s2(i)
            myc(1) = this%c2(2*i-1)
            myc(2) = this%c2(2*i)
            do k = 1,step
                t1(k) = mat(this%xy(2*i-1), j+k-1)
                t2(k) = mat(this%xy(2*i), j+k-1)
            enddo
            do k = 1, step
                mat(this%xy(2*i-1), j+k-1) = myc(1) * t1(k) + mys * t2(k)
                mat(this%xy(2*i), j+k-1) = myc(2) * t2(k) + conjg(mys) * t1(k)
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

    ! rmultinv part with new data
    call rmultthreeelementbase(this%c2inv, this%s2inv, this%xy, this%nrofentries, mat)
end subroutine GeneralSingleColExp_adjoint_over_two

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Perform the multiplication of the inverse exponential with a matrix:
!> out = mat^-1*this.
!
!> @param[in] this The exponential that we consider.
!> @param[inout] mat The matrix that we modify.
!--------------------------------------------------------------------
subroutine GeneralSingleColExp_rmultinv(this, mat)
    class(GeneralSingleColExp), intent(in) :: this
    complex(kind=kind(0.D0)), dimension(:, :), intent(inout) :: mat

    call rmultthreeelementbase(this%cinv, this%sinv, this%xy, this%nrofentries, mat)
end subroutine GeneralSingleColExp_rmultinv

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
!> @param [in] weight a real prefactor
!> @param [in] mav the average chemical potential
!> @param [in] eps my definition of a local zero.
!--------------------------------------------------------------------
subroutine expof2x2hermitianmatrix(diag1, diag2, offout, diag, offinp, weight, mav, eps)
    real (kind=kind(0.d0)), intent(out) :: diag1, diag2
    complex(kind=kind(0.D0)), intent(out) :: offout
    real (kind=kind(0.d0)), intent(in) :: weight, diag, mav, eps
    complex(kind=kind(0.D0)), intent(in) :: offinp
    real (kind=kind(0.d0)) :: myexp
    
    call expof2x2tracelesshermitianmatrix(diag1, diag2, offout, diag, offinp, weight)
    
    if(abs(mav) > eps) then ! fixup chemical potential
        myexp = exp(weight * mav)
        diag1 = diag1 * myexp
        diag2 = diag2 * myexp
            
        offout = offout * myexp
    endif
end subroutine

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This sets up the data to perform the exponentiation of a 
!> strictly sparse matrix.
!> The internal layout is that the non-zero element a_xy stored at (x,y) in the matrix
!> has x(i) at x(2i-1) and y(i) at x(2i).
!> The two values for the diagonals are stored in c(2i) and c(2i+1).
!
!> @param[inout] this The GeneralSingleColExp object.
!> @param[in] nodes The nodes that belong to this color.
!> @param[in] nredges How many nodes of this color.
!> @param[in] mys A vector containing the diagonal of the matrix
!> @param[in] weight A prefactor for the exponent.
!--------------------------------------------------------------------
subroutine GeneralSingleColExp_init(this, nodes, nredges, mys, weight)
    class(GeneralSingleColExp), intent(inout) :: this
    type(node), dimension(:), intent(in) :: nodes
    real(kind=kind(0.D0)), intent(in), dimension(:) :: mys
    integer, intent(in) :: nredges
    real (kind=kind(0.d0)), intent(in) :: weight
    integer :: i
    real (kind=kind(0.d0)) :: my1, my2, localzero, md, mav, dweight

    allocate(this%xy(2*nredges), this%c(2*nredges), this%s(nredges), this%cinv(2*nredges), this%sinv(nredges))
    allocate(this%c2(2*nredges), this%c2inv(2*nredges), this%s2(nredges), this%s2inv(nredges))
    this%nrofentries = nredges
#ifndef NDEBUG
    write(*,*) " [GeneralSingleColExp_init]: Setting up strict. sparse matrix with ", nredges, "edges"
#endif
    do i = 1, nredges
        this%xy(2*i-1) = nodes(i)%x
        this%xy(2*i) = nodes(i)%y
        !calculate Frobenius norm
        my1 = mys(nodes(i)%x)
        my2 = mys(nodes(i)%y)
        ! dependence on weight cancels in all comps
        localzero = 1E-15*frobnorm(my1, my2, nodes(i)%axy) ! definition of my local scale that defines zero
        md = 0.5*(my1 - my2)
        mav = 0.5*(my1 + my2)
        ! This is the order of operations that yields stable matrix inversions
        ! We assume that the matrix that we have decomposed is hermitian:
        ! M=(d  , b)
        !   (b^*, -d)
        ! with d = (my1-m2)/2 and mav = (my1+m2)/2
        
        call expof2x2hermitianmatrix(this%c(2*i-1), this%c(2*i), this%s(i), md, nodes(i)%axy, weight, mav, localzero)
        
        dweight = -weight
        call expof2x2hermitianmatrix(this%cinv(2*i-1), this%cinv(2*i), this%sinv(i), md, nodes(i)%axy, dweight, mav, localzero)
        
        dweight = 0.5*weight
        call expof2x2hermitianmatrix(this%c2(2*i-1), this%c2(2*i), this%s2(i), md, nodes(i)%axy, dweight, mav, localzero)
        
        dweight = -0.5*weight
        call expof2x2hermitianmatrix(this%c2inv(2*i-1), this%c2inv(2*i), this%s2inv(i), md, nodes(i)%axy, dweight, mav, localzero)

    enddo
! All nodes that we have been passed are now from a single color.
! They constitute now a strictly sparse matrix adapted to a chemical potential..
end subroutine GeneralSingleColExp_init

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Tidy up internal data structures.
!
!> @param[inout] this the GeneralSingleColExp object.
!--------------------------------------------------------------------
subroutine GeneralSingleColExp_dealloc(this)
    class(GeneralSingleColExp), intent(inout) :: this
    deallocate(this%cinv, this%sinv, this%c2inv, this%s2inv)
    call this%TraceLessSingleColExp%dealloc()
end subroutine GeneralSingleColExp_dealloc

end module GeneralSingleColExp_mod
