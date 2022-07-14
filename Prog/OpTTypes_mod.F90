!  Copyright (C) 2020 - 2022 The ALF project
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
    use mat_subroutines
    implicit none

    private
    public :: RealExpOpT, CmplxExpOpT
    
    !--------------------------------------------------------------------
    !> @author
    !> ALF-project
    !> @brief
    !> Encapsulates operations for real exponentiated OpTs.
    !>
    !--------------------------------------------------------------------
    type, extends(ContainerElementBase) :: RealExpOpT
        Real(kind=kind(0.d0)), allocatable, dimension(:,:) :: mat, invmat, mat_1D2, invmat_1D2 !>We store the matrix in the class
        Real(kind=kind(0.d0)) :: g, Zero
        integer, allocatable :: P(:)
        Integer :: m, n, Ndim_hop
        
    contains
        procedure :: init => RealExpOpT_init ! initialize and allocate matrices
        procedure :: dealloc => RealExpOpT_dealloc ! dealloc matrices
        procedure :: rmult => RealExpOpT_rmult ! right multiplication with Op_T
        procedure :: lmult => RealExpOpT_lmult
        procedure :: rmultinv => RealExpOpT_rmultinv ! right multiplication with Op_T inverse
        procedure :: lmultinv => RealExpOpT_lmultinv
        procedure :: adjointaction => RealExpOpT_adjointaction
        procedure :: dump => RealExpOpT_dump ! dump matrices for debugging to screen
    end type RealExpOpT

    !--------------------------------------------------------------------
    !> @author
    !> ALF-project
    !> @brief
    !> Encapsulates operations for complex exponentiated OpTs.
    !>
    !--------------------------------------------------------------------
    type, extends(ContainerElementBase) :: CmplxExpOpT
        Complex(kind=kind(0.d0)), allocatable, dimension(:,:) :: mat, invmat, mat_1D2, invmat_1D2 !>We store the matrix inclass
        Complex(kind=kind(0.d0)) :: g
        Real(kind=kind(0.d0)) :: Zero
        integer, allocatable :: P(:)
        Integer :: m, n, Ndim_hop
    contains
        procedure :: init => CmplxExpOpT_init ! initialize and allocate matrices
        procedure :: dealloc => CmplxExpOpT_dealloc ! dealloc matrices
        procedure :: rmult => CmplxExpOpT_rmult ! right multiplication with Op_T
        procedure :: lmult => CmplxExpOpT_lmult
        procedure :: rmultinv => CmplxExpOpT_rmultinv ! right multiplication with Op_T inverse
        procedure :: lmultinv => CmplxExpOpT_lmultinv
        procedure :: adjointaction => CmplxExpOpT_adjointaction
        procedure :: dump => CmplxExpOpT_dump ! dump matrices for debugging to screen
    end type CmplxExpOpT

contains
    subroutine RealExpOpT_init(this, Op_T)
        class(RealExpOpT) :: this
        Type(Operator), intent(in) :: Op_T
        Complex(kind=kind(0.d0)), allocatable, dimension(:,:) :: cmat, cinvmat
        Complex(kind=kind(0.d0)) :: cg
        Integer :: i, j
        
        this%Zero = 1.E-12
        this%Ndim_hop = Op_T%N
        allocate(this%mat(this%Ndim_hop, this%Ndim_hop), this%invmat(this%Ndim_hop, this%Ndim_hop))
        allocate(this%mat_1D2(this%Ndim_hop, this%Ndim_hop), this%invmat_1D2(this%Ndim_hop, this%Ndim_hop))
        allocate(cmat(this%Ndim_hop, this%Ndim_hop), cinvmat(this%Ndim_hop, this%Ndim_hop))

        cg = -Op_T%g
        Call Op_exp(cg, Op_T, cinvmat)
        cg = Op_T%g
        Call Op_exp(cg, Op_T, cmat)
        ! copy over the data to the real storage
        this%mat = DBLE(cmat)
        this%invmat = DBLE(cinvmat)
        
        cg = -Op_T%g/2.0
        Call Op_exp(cg, Op_T, cinvmat)
        cg = Op_T%g/2.0
        Call Op_exp(cg, Op_T, cmat)
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
        this%P = Op_T%P ! copy all data locally to be consistent and less error prone
        this%g = DBLE(Op_T%g)
        deallocate(cmat, cinvmat)
    end subroutine
    
    subroutine RealExpOpT_adjointaction(this, arg, t)
        class(RealExpOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer :: n1, n2
        Integer, intent(in) :: t
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('L', 'U', this%Ndim_hop, n1, n2, this%mat_1D2, this%P, arg)
            call ZDSLSYMM('R', 'U', this%Ndim_hop, n1, n2, this%invmat_1D2, this%P, arg)
        Endif
    end subroutine
    
    subroutine RealExpOpT_rmult(this, arg, t)
        class(RealExpOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout),  dimension(:,:) :: arg
        Integer :: n1, n2
        Integer, intent(in) :: t
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('R', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine
    
    subroutine RealExpOpT_rmultinv(this, arg, t)
        class(RealExpOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout),  dimension(:,:) :: arg
        Integer :: n1, n2
        Integer, intent(in) :: t
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('R', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine
    
    subroutine RealExpOpT_lmult(this, arg, t)
        class(RealExpOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout),  dimension(:,:) :: arg
        Integer, intent(in) :: t
        integer :: n1, n2
        
        ! taken from mmthr
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('L', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine
    
    subroutine RealExpOpT_lmultinv(this, arg, t)
        class(RealExpOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer, intent(in) :: t
        integer :: n1, n2
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('L', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine
    
    subroutine CmplxExpOpT_init(this, Op_T)
        class(CmplxExpOpT) :: this
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
        this%P = Op_T%P ! copy all data locally to be consistent and less error prone
    end subroutine

    subroutine CmplxExpOpT_adjointaction(this, arg, t)
        class(CmplxExpOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer, intent(in) :: t
        Integer :: n1, n2

        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then            
            call ZSLHEMM('L', 'U', this%Ndim_hop, n1, n2, this%mat_1D2, this%P, arg)
            call ZSLHEMM('R', 'U', this%Ndim_hop, n1, n2, this%invmat_1D2, this%P, arg)
        Endif
        
    end subroutine
    
    subroutine CmplxExpOpT_rmult(this, arg, t)
        class(CmplxExpOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer, intent(in) :: t
        Integer :: n1, n2
        
        ! taken from mmthl
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call ZSLHEMM('R', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine
    
        subroutine CmplxExpOpT_rmultinv(this, arg, t)
        class(CmplxExpOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer, intent(in) :: t
        Integer :: n1, n2
        
        ! taken from mmthl_m1
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call ZSLHEMM('R', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine
    
    subroutine CmplxExpOpT_lmult(this, arg, t)
        class(CmplxExpOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer, intent(in) :: t
        integer :: n1, n2
        
        ! taken from mmthr
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call ZSLHEMM('L', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine
    
    subroutine CmplxExpOpT_lmultinv(this, arg, t)
        class(CmplxExpOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer, intent(in) :: t
        integer :: n1, n2
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call ZSLHEMM('L', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine

    subroutine CmplxExpOpT_dump(this)
        class(CmplxExpOpT), intent(in) :: this
        integer :: i,j

        do i = 1, size(this%mat, 1)
            write (*,*) (dble(this%mat(i,j)), j = 1,size(this%mat,2) )
        enddo
        write (*,*) "---------------"
        do i = 1, size(this%mat, 1)
            write (*,*) (dble(this%invmat(i,j)), j = 1,size(this%mat,2) )
        enddo
    end subroutine

    subroutine RealExpOpT_dump(this)
        class(RealExpOpT), intent(in) :: this
        integer :: i,j

        do i = 1, size(this%mat, 1)
            write (*,*) (dble(this%mat(i,j)), j = 1,size(this%mat,2) )
        enddo
        write (*,*) "---------------"
        do i = 1, size(this%mat, 1)
            write (*,*) (dble(this%invmat(i,j)), j = 1,size(this%mat,2) )
        enddo
    end subroutine

    subroutine CmplxExpOpT_dealloc(this)
        class(CmplxExpOpT), intent(inout) :: this

        deallocate(this%mat, this%invmat, this%mat_1D2, this%invmat_1D2, this%P)
    end subroutine

    subroutine RealExpOpT_dealloc(this)
        class(RealExpOpT), intent(inout) :: this

        deallocate(this%mat, this%invmat, this%mat_1D2, this%invmat_1D2, this%P)
    end subroutine

end module OpTTypes_mod
