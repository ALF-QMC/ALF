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


module mscbOpT_mod
    use ContainerElementBase_mod
    use Operator_mod
    use MvG_mod
    use Exponentials_mod
    implicit none

    private
    public :: RealmscbOpT, CmplxmscbOpT
    
    !--------------------------------------------------------------------
    !> @author
    !> ALF-project
    !> @brief
    !> Encapsulates operations for real exponentiated OpTs.
    !>
    !--------------------------------------------------------------------
    type, extends(ContainerElementBase) :: RealmscbOpT
        Real(kind=kind(0.d0)), allocatable, dimension(:,:) :: mat, invmat, mat_1D2, invmat_1D2 !>We store the matrix in the class
        Real(kind=kind(0.d0)) :: g, Zero
        integer, allocatable :: P(:)
        type(FullExp) :: fe
        Integer :: m, n, Ndim_hop
        
    contains
        procedure :: init => RealmscbOpT_init ! initialize and allocate matrices
        procedure :: dealloc => RealmscbOpT_dealloc ! dealloc matrices
        procedure :: rmult => RealmscbOpT_rmult ! right multiplication with Op_T
        procedure :: lmult => RealmscbOpT_lmult
        procedure :: rmultinv => RealmscbOpT_rmultinv ! right multiplication with Op_T inverse
        procedure :: lmultinv => RealmscbOpT_lmultinv
        procedure :: adjointaction => RealmscbOpT_adjointaction
        procedure :: dump => RealmscbOpT_dump ! dump matrices for debugging to screen
    end type RealmscbOpT

    !--------------------------------------------------------------------
    !> @author
    !> ALF-project
    !> @brief
    !> Encapsulates operations for complex exponentiated OpTs.
    !>
    !--------------------------------------------------------------------
    type, extends(ContainerElementBase) :: CmplxmscbOpT
        Complex(kind=kind(0.d0)) :: g
        Real(kind=kind(0.d0)) :: Zero
        integer, allocatable :: P(:)
        type(FullExp) :: fe
        Integer :: m, n, Ndim_hop
    contains
        procedure :: init => CmplxmscbOpT_init ! initialize and allocate matrices
        procedure :: dealloc => CmplxmscbOpT_dealloc ! dealloc matrices
        procedure :: rmult => CmplxmscbOpT_rmult ! right multiplication with Op_T
        procedure :: lmult => CmplxmscbOpT_lmult
        procedure :: rmultinv => CmplxmscbOpT_rmultinv ! right multiplication with Op_T inverse
        procedure :: lmultinv => CmplxmscbOpT_lmultinv
        procedure :: adjointaction => CmplxmscbOpT_adjointaction
        procedure :: dump => CmplxmscbOpT_dump ! dump matrices for debugging to screen
    end type CmplxmscbOpT

contains
    subroutine RealmscbOpT_init(this, Op_T)
        class(RealmscbOpT) :: this
        Type(Operator), intent(in) :: Op_T
        Complex(kind=kind(0.d0)), allocatable, dimension(:,:) :: cmat, cinvmat
        Complex(kind=kind(0.d0)) :: cg
        Integer :: i, j
        type(GraphData) :: gd
        
        gd = mat2verts(Op_T%O) ! convert to graphdata structure
        call MvG_decomp(gd%verts) ! perform the decomposition
        
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
        this%P = Op_T%P
        this%g = DBLE(Op_T%g)
        deallocate(cmat, cinvmat)
    end subroutine
    
    subroutine RealmscbOpT_adjointaction(this, arg)
        class(RealmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        Integer :: n1, n2
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('L', 'U', this%Ndim_hop, n1, n2, this%mat_1D2, this%P, arg)
            call ZDSLSYMM('R', 'U', this%Ndim_hop, n1, n2, this%invmat_1D2, this%P, arg)
        Endif
    end subroutine
    
    subroutine RealmscbOpT_rmult(this, arg)
        class(RealmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout),  dimension(:,:) :: arg
        Integer :: n1, n2
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('R', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine
    
    subroutine RealmscbOpT_rmultinv(this, arg)
        class(RealmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout),  dimension(:,:) :: arg
        Integer :: n1, n2
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('R', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine
    
    subroutine RealmscbOpT_lmult(this, arg)
        class(RealmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout),  dimension(:,:) :: arg
        integer :: n1, n2
        
        ! taken from mmthr
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('L', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine
    
    subroutine RealmscbOpT_lmultinv(this, arg)
        class(RealmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        integer :: n1, n2
        
        n1 = size(arg,1)
        n2 = size(arg,2)
        If ( this%g*this%g > this%Zero ) then
            call ZDSLSYMM('L', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine
    
    subroutine CmplxmscbOpT_init(this, Op_T, method)
        class(CmplxmscbOpT) :: this
        Type(Operator), intent(in) :: Op_T
        integer, intent(in) :: method
        Integer :: i, k
        type(GraphData) :: gd
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: tmp
        
        this%Zero = 1.E-12
        this%Ndim_hop = Op_T%N
        this%g = Op_T%g
        allocate(tmp(Op_T%N, Op_T%N))
        tmp = this%g* Op_T%O
        
        gd = mat2verts(Op_T%O) ! convert to graphdata structure
        deallocate(tmp)
        call MvG_decomp(gd%verts) ! perform the decomposition
        
        ! some sanity checks and status informations
        gd%usedcolors = 0
        gd%nredges = 0
        do i = 1, gd%ndim
            gd%deltag = max(gd%deltag, gd%verts(i)%degree)
            do k = 1, gd%verts(i)%degree
                if (gd%verts(i)%nbrs(k) > i) gd%nredges = gd%nredges + 1
                if (gd%verts(i)%nbrs(k) > gd%ndim) then
                    write(*,*) "invalid nbr!!!"
                    STOP
                endif
                gd%usedcolors = max(gd%usedcolors, gd%verts(i)%cols(k))
            enddo
        enddo
        write (*,*) "Nr edges: ", gd%nredges
        if (gd%usedcolors == gd%deltag) then
            write(*,*) "Maximum Degree", gd%deltag, ". Found", gd%usedcolors," Families -> optimal decomposition"
        else
            write(*,*) "Maximum Degree", gd%deltag, ". Found", gd%usedcolors," Families"
        endif

        this%fe = createFullExponentialfromGraphData(gd, method)
        
        ! check wether it is supported behaviour
        do i = 1, size(Op_T%P)
        if (Op_T%P(i) /= i) then
        write (*,*) "unsupported case."
        endif
        enddo
        this%P = Op_T%P

    end subroutine

    subroutine CmplxmscbOpT_adjointaction(this, arg)
        class(CmplxmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg

        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%fe%adjoint_over_two(arg)
        Endif

    end subroutine
    
    subroutine CmplxmscbOpT_rmult(this, arg)
        class(CmplxmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        
        !FIXME: P
        
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%fe%rmult(arg)
!            call ZSLHEMM('R', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine

    subroutine CmplxmscbOpT_rmultinv(this, arg)
        class(CmplxmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        
        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%fe%rmultinv(arg)
!            call ZSLHEMM('R', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine
    
    subroutine CmplxmscbOpT_lmult(this, arg)
        class(CmplxmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        
        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%fe%lmult(arg)
!             call ZSLHEMM('L', 'U', this%Ndim_hop, n1, n2, this%mat, this%P, arg)
        Endif
    end subroutine

    subroutine CmplxmscbOpT_lmultinv(this, arg)
        class(CmplxmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg

        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%fe%lmultinv(arg)
!             call ZSLHEMM('L', 'U', this%Ndim_hop, n1, n2, this%invmat, this%P, arg)
        Endif
    end subroutine

    subroutine CmplxmscbOpT_dump(this)
        class(CmplxmscbOpT), intent(in) :: this
        integer :: i, j
        
        write (*,*) "method: ", this%fe%method
        write (*,*) "colors: ", this%fe%stages(1)%nrofcols

    end subroutine

    subroutine RealmscbOpT_dump(this)
        class(RealmscbOpT), intent(in) :: this
        integer :: i,j

    end subroutine
    
    subroutine CmplxmscbOpT_dealloc(this)
        class(CmplxmscbOpT), intent(inout) :: this
        
        call this%fe%dealloc()
    end subroutine

    subroutine RealmscbOpT_dealloc(this)
        class(RealmscbOpT), intent(inout) :: this
        
    end subroutine
    
end module mscbOpT_mod
