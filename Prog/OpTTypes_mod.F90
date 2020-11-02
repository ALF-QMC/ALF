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


module OpTTypes_mod
    use ContainerElementBase_mod
    use Operator_mod
    implicit none
    
    !--------------------------------------------------------------------
    !> @author
    !> ALF-project
    !> @brief
    !> Encapsulates Operations for real OpTs
    !>
    !--------------------------------------------------------------------
    type, extends(ContainerElementBase) :: RealOpT
        Real(kind=kind(0.d0)), allocatable, dimension(:,:) :: mat, invmat, mat_1D2, invmat_1D2 !> We store the matrix here in the class
        Real(kind=kind(0.d0)) :: g, Zero
        integer, pointer :: P(:)
        Integer :: m, n, Ndim_hop
        
    contains
        procedure :: init => RealOpT_init ! initialize and allocate matrices
        procedure :: dealloc => RealOpT_dealloc ! dealloc matrices
        procedure :: rmult => RealOpT_rmult ! right multiplication with Op_T
        procedure :: lmult => RealOpT_lmult
        procedure :: rmultinv => RealOpT_rmultinv ! right multiplication with Op_T inverse
        procedure :: lmultinv => RealOpT_lmultinv
        procedure :: adjointaction => RealOpT_adjointaction
        procedure :: dump => RealOpT_dump ! dump matrices for debugging to screen
    end type RealOpT

    !--------------------------------------------------------------------
    !> @author
    !> ALF-project
    !> @brief
    !> Encapsulates Operations for Complex OpTs
    !>
    !--------------------------------------------------------------------
    type, extends(ContainerElementBase) :: CmplxOpT
        Complex(kind=kind(0.d0)),allocatable, dimension(:,:) :: mat, invmat, mat_1D2, invmat_1D2 !> We store the matrix here in the class
        Complex(kind=kind(0.d0)) :: g
        Real(kind=kind(0.d0)) :: Zero
        integer, pointer :: P(:)
        Integer :: m, n, Ndim_hop
    contains
        procedure :: init => CmplxOpT_init ! initialize and allocate matrices
        procedure :: dealloc => CmplxOpT_dealloc ! dealloc matrices
        procedure :: rmult => CmplxOpT_rmult ! right multiplication with Op_T
        procedure :: lmult => CmplxOpT_lmult
        procedure :: rmultinv => CmplxOpT_rmultinv ! right multiplication with Op_T inverse
        procedure :: lmultinv => CmplxOpT_lmultinv
        procedure :: adjointaction => CmplxOpT_adjointaction
        procedure :: dump => CmplxOpT_dump ! dump matrices for debugging to screen
    end type CmplxOpT

contains
    subroutine RealOpT_init(this, Op_T)
        class(RealOpT) :: this
        Type(Operator), intent(in) :: Op_T
        Complex(kind=kind(0.d0)),allocatable, dimension(:,:) :: cmat, cinvmat
        Complex(kind=kind(0.d0)) :: cg
        Integer :: i, j
        
        this%Zero = 1.E-12
        this%Ndim_hop = Op_T%N
        allocate(this%mat(this%Ndim_hop, this%Ndim_hop), this%invmat(this%Ndim_hop, this%Ndim_hop))
        allocate(this%mat_1D2(this%Ndim_hop, this%Ndim_hop), this%invmat_1D2(this%Ndim_hop, this%Ndim_hop))
        allocate(cmat(this%Ndim_hop, this%Ndim_hop), cinvmat(this%Ndim_hop, this%Ndim_hop))

        cg = -Op_T%g
        Call  Op_exp(cg, Op_T, cinvmat)
        cg = Op_T%g
        Call  Op_exp(cg, Op_T, cmat)
        ! copy over the data to the real storage
        this%mat = DBLE(cmat)
        this%invmat = DBLE(cinvmat)
        
        cg = -Op_T%g/2.0
        Call  Op_exp(cg, Op_T, cinvmat)
        cg = Op_T%g/2.0
        Call  Op_exp(cg, Op_T, cmat)
        ! copy over the data to the real storage
        this%mat_1D2 = DBLE(cmat)
        this%invmat_1D2 = DBLE(cinvmat)
        
        DO i = 1, this%Ndim_hop
            DO j = i, this%Ndim_hop
                this%mat(i, j) = (this%mat(i, j) + this%mat(j, i))/2.D0
                this%invmat(i, j) = (this%invmat(i, j) + this%invmat(j, i))/2.D0
                
                this%mat_1D2(i, j) = (this%mat_1D2(i, j) + this%mat_1D2(j, i))/2.D0
                this%invmat_1D2(i, j) = (this%invmat_1D2(i, j) + this%invmat_1D2(j, i))/2.D0
            ENDDO
        ENDDO
        this%P => Op_T%P
        this%g = DBLE(Op_T%g)
        deallocate(cmat, cinvmat)
    end subroutine
    
    subroutine RealOpT_adjointaction(this, arg)
        class(RealOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer :: n1, n2
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('L', 'U', this%Ndim_hop, n1, n2, this%mat_1D2, this%P, arg)
            call ZDSLSYMM('R', 'U', this%Ndim_hop, n1, n2, this%invmat_1D2, this%P, arg)
        Endif
    end subroutine
    
    subroutine RealOpT_rmult(this, arg)
        class(RealOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout),  dimension(:,:) :: arg
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: out
        Real(kind=kind(0.D0)), allocatable, dimension(:) :: rwork
        Integer :: n1, n2
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('R', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine
    
    subroutine RealOpT_rmultinv(this, arg)
        class(RealOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout),  dimension(:,:) :: arg
        Integer :: n1, n2
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('R', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine
    
    subroutine RealOpT_lmult(this, arg)
        class(RealOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout),  dimension(:,:) :: arg
        integer :: n1, n2
        
        ! taken from mmthr
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('L', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine
    
    subroutine RealOpT_lmultinv(this, arg)
        class(RealOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        integer :: n1, n2
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('L', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine
    
    subroutine CmplxOpT_init(this, Op_T)
        class(CmplxOpT) :: this
        Type(Operator), intent(in) :: Op_T
        Integer :: i, j
        
        this%Zero = 1.E-12
        this%Ndim_hop = Op_T%N
        
        allocate(this%mat(this%Ndim_hop, this%Ndim_hop), this%invmat(this%Ndim_hop, this%Ndim_hop))
        allocate(this%mat_1D2(this%Ndim_hop, this%Ndim_hop), this%invmat_1D2(this%Ndim_hop, this%Ndim_hop))
        
        this%g = -Op_T%g/2.0
        Call  Op_exp(this%g, Op_T, this%invmat_1D2)
        this%g = Op_T%g/2.0
        Call  Op_exp(this%g, Op_T, this%mat_1D2)
        
        this%g = -Op_T%g
        Call  Op_exp(this%g, Op_T, this%invmat)
        this%g = Op_T%g
        Call  Op_exp(this%g, Op_T, this%mat)
        
        DO i = 1, this%Ndim_hop
            DO j = i, this%Ndim_hop
                this%mat(i, j) = (this%mat(i, j) + Conjg(this%mat(j, i)))/2.D0
                this%invmat(i, j) = (this%invmat(i, j) + Conjg(this%invmat(j, i)))/2.D0
                
                this%mat_1D2(i, j) = (this%mat_1D2(i, j) + Conjg(this%mat_1D2(j, i)))/2.D0
                this%invmat_1D2(i, j) = (this%invmat_1D2(i, j) + Conjg(this%invmat_1D2(j, i)))/2.D0
            ENDDO
        ENDDO
        this%P => Op_T%P

    end subroutine

    subroutine CmplxOpT_adjointaction(this, arg)
        class(CmplxOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer :: n1, n2
        
!       Integer :: info, lwork
!       Integer, allocatable, dimension(:) :: ipiv
!       Complex(kind=kind(0.D0)), allocatable, dimension(:) :: work2
!       Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: work, work3
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            ! this%Ndim_hop
            ! n1 : arg = n1 x n2 matrix
            ! n2 : 
            ! mat and invmat: ndimhop x ndimhop matrices
            
            ! subroutine ZSLHEMM(side, uplo, N, M1, M2, A, P, Mat)
            ! NUMBLOCKS = 1
            ! IDXLIST(1) = 1
            ! DIMLIST(1) = 1
            
!            ALLOCATE (WORK(this%Ndim_hop, n2))
!            CALL ZLACPY('A', n1, n2, arg(this%P(1), 1), n1, work(1,1), this%Ndim_hop )
            ! mat^-1 * work = X <=>  mat * X = work
            ! call zhemm('L', 'U', ndimhop, n2, alpha, mat_1D2(1,1), ndimhop, work(1, 1), ndimhop, beta, arg(P(1), n1))
            ! arg = 0*arg + mat^-1 * work
            ! deallocate (work)
            
            ! the following should be equivalent.
!            allocate (ipiv(n1), work2(this%Ndim_hop*8) )
!            work3 = this%invmat_1D2
!            lwork = this%Ndim_hop*8
!            call zhesv('U', this%Ndim_hop, n2, work3, this%Ndim_hop, ipiv, work(1,1), this%Ndim_hop, work2, lwork, info)
!            call zlacpy('A', this%Ndim_hop, n2, work(1,1), this%Ndim_hop, arg(this%P(1), 1), n1)



          call ZSLHEMM('L', 'U', this%Ndim_hop, n1, n2, this%mat_1D2, this%P, arg)

            call ZSLHEMM('R', 'U', this%Ndim_hop, n1, n2, this%invmat_1D2, this%P, arg)
!           deallocate(work, work2, ipiv)
        Endif
        
    end subroutine
    
    subroutine CmplxOpT_rmult(this, arg)
        class(CmplxOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer :: n1, n2
        
        ! taken from mmthl
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call ZSLHEMM('R', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine
    
        subroutine CmplxOpT_rmultinv(this, arg)
        class(CmplxOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer :: n1, n2
        
        ! taken from mmthl_m1
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call ZSLHEMM('R', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine
    
    subroutine CmplxOpT_lmult(this, arg)
        class(CmplxOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        integer :: n1, n2
        
        ! taken from mmthr
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call ZSLHEMM('L', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine
    
    subroutine CmplxOpT_lmultinv(this, arg)
        class(CmplxOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        integer :: n1, n2
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call ZSLHEMM('L', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine

    subroutine CmplxOpT_dump(this)
        class(CmplxOpT), intent(in) :: this
        integer :: i,j

        do i = 1, size(this%mat, 1)
            write (*,*) (dble(this%mat(i,j)), j = 1,size(this%mat,2) )
        enddo
        write (*,*) "---------------"
        do i = 1, size(this%mat, 1)
            write (*,*) (dble(this%invmat(i,j)), j = 1,size(this%mat,2) )
        enddo
    end subroutine

    subroutine RealOpT_dump(this)
        class(RealOpT), intent(in) :: this
        integer :: i,j

        do i = 1, size(this%mat, 1)
            write (*,*) (dble(this%mat(i,j)), j = 1,size(this%mat,2) )
        enddo
        write (*,*) "---------------"
        do i = 1, size(this%mat, 1)
            write (*,*) (dble(this%invmat(i,j)), j = 1,size(this%mat,2) )
        enddo
    end subroutine
    
        subroutine CmplxOpT_dealloc(this)
        class(CmplxOpT), intent(inout) :: this
        
        deallocate(this%mat, this%invmat, this%mat_1D2, this%invmat_1D2)
    end subroutine

    subroutine RealOpT_dealloc(this)
        class(RealOpT), intent(inout) :: this
        
        deallocate(this%mat, this%invmat, this%mat_1D2, this%invmat_1D2)
    end subroutine
    
end module OpTTypes_mod
