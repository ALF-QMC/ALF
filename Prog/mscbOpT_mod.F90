!  Copyright (C) 2020 -2021 The ALF project
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
    use Exponentials_mod, only: EulerExp, FullExp
    implicit none

    private
    public :: CmplxmscbOpT, CmplxEulermscbOpT

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
        procedure :: lmult => CmplxmscbOpT_lmult  ! left multiplication with Op_T
        procedure :: rmultinv => CmplxmscbOpT_rmultinv ! right multiplication with Op_T inverse
        procedure :: lmultinv => CmplxmscbOpT_lmultinv  ! left multiplication with Op_T inverse
        procedure :: adjoint => CmplxmscbOpT_adjoint
        procedure :: reverseadjoint => CmplxmscbOpT_reverseadjoint
        procedure :: adjoint_over_two => CmplxmscbOpT_adjoint_over_two
        procedure :: dump => CmplxmscbOpT_dump ! dump matrices for debugging to screen
    end type CmplxmscbOpT

    !--------------------------------------------------------------------
    !> @author
    !> ALF-project
    !> @brief
    !> A specialized class for encapsulating the Euler approximation.
    !>
    !--------------------------------------------------------------------
    type, extends(ContainerElementBase) :: CmplxEulermscbOpT
        Complex(kind=kind(0.d0)) :: g
        Real(kind=kind(0.d0)) :: Zero
        integer, allocatable :: P(:)
        type(EulerExp) :: ee
        Integer :: m, n, Ndim_hop
    contains
        procedure :: init => CmplxEulermscbOpT_init ! initialize and allocate matrices
        procedure :: dealloc => CmplxEulermscbOpT_dealloc ! dealloc matrices
        procedure :: rmult => CmplxEulermscbOpT_rmult ! right multiplication with Op_T
        procedure :: lmult => CmplxEulermscbOpT_lmult
        procedure :: rmultinv => CmplxEulermscbOpT_rmultinv ! right multiplication with Op_T inverse
        procedure :: lmultinv => CmplxEulermscbOpT_lmultinv
        procedure :: adjoint => CmplxEulermscbOpT_adjoint
        procedure :: reverseadjoint => CmplxEulermscbOpT_reverseadjoint
        procedure :: adjoint_over_two => CmplxEulermscbOpT_adjoint_over_two
        procedure :: dump => CmplxEulermscbOpT_dump ! dump matrices for debugging to screen
    end type CmplxEulermscbOpT

contains
    !--------------------------------------------------------------------
    !> @author
    !> ALF-project
    !> @brief
    !> Transform an Op_T object to a GraphData object and a vector containing
    !> the matrix diagonal.
    !
    !> @param Op_T the operator.
    !> @param gd the Graphdata object. We take care to create one.
    !> @param diags the vector of the diagonal. We allocate it.
    !>
    !--------------------------------------------------------------------
    subroutine Op_T_to_graphdata(Op_T, gd, diags)
        use Operator_mod
        use MvG_mod
        use graphdata_mod, only: GraphData, mat2verts
        implicit none

        Type(Operator), intent(in) :: Op_T
        type(GraphData), intent(inout) :: gd
        Real(kind=kind(0.D0)), allocatable, dimension(:), intent(inout) :: diags
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: tmp
        integer :: i
        
        allocate(tmp(Op_T%N, Op_T%N), diags(Op_T%N))
        tmp = Op_T%g* Op_T%O
        
        ! Let's retrieve the diagonal of the matrix, hence the chemical potential
        do i = 1, Op_T%N
            diags(i) = DBLE(tmp(i, i))!If the input matrix is hermitian, this has to be real.
            tmp(i, i) = 0 ! Set to zero afterwards to disable any checks in the MvG decomposition
        enddo
        gd = mat2verts(tmp) ! convert to graphdata structure
        deallocate(tmp)
        call MvG_decomp(gd%verts) ! perform the decomposition
    end subroutine
    
    subroutine CmplxmscbOpT_init(this, Op_T, method)
        use Operator_mod
        use graphdata_mod
        implicit none

        class(CmplxmscbOpT) :: this
        type(Operator), intent(in) :: Op_T
        integer, intent(in) :: method
        integer :: i
        type(GraphData) :: gd
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: tmp
        Real(kind=kind(0.D0)), allocatable, dimension(:) :: diags

        this%Zero = 1.E-12
        this%Ndim_hop = Op_T%N
        this%g = Op_T%g

        call Op_T_to_graphdata(Op_T, gd, diags)
        ! some sanity checks and status informations
        call determine_used_colors_of_graph(gd)
        write (*,*) "Nr edges: ", gd%nredges
        if (gd%usedcolors == gd%deltag) then
            write(*,*) "Maximum Degree", gd%deltag, ". Found", gd%usedcolors," Families -> optimal decomposition"
        else
            write(*,*) "Maximum Degree", gd%deltag, ". Found", gd%usedcolors," Families"
        endif

        this%fe = createFullExponentialfromGraphData(gd, diags, method)

        ! check wether it is supported behaviour
        do i = 1, size(Op_T%P)
            if (Op_T%P(i) /= i) then
                write (*,*) "P unsupported"
            endif
        enddo
        this%P = Op_T%P
        call dealloc_graphdata(gd)
        deallocate(diags)
    end subroutine

    subroutine CmplxmscbOpT_reverseadjoint(this, arg)
        class(CmplxmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg

        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%fe%reverseadjoint(arg)
        Endif
    end subroutine

    subroutine CmplxmscbOpT_adjoint(this, arg)
        class(CmplxmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg

        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%fe%adjoint(arg)
        Endif
    end subroutine
    
    subroutine CmplxmscbOpT_adjoint_over_two(this, arg)
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
        Endif
    end subroutine

    subroutine CmplxmscbOpT_rmultinv(this, arg)
        class(CmplxmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        
        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%fe%rmultinv(arg)
        Endif
    end subroutine

    subroutine CmplxmscbOpT_lmult(this, arg)
        class(CmplxmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        
        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%fe%lmult(arg)
        Endif
    end subroutine

    subroutine CmplxmscbOpT_lmultinv(this, arg)
        class(CmplxmscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg

        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%fe%lmultinv(arg)
        Endif
    end subroutine

    subroutine CmplxmscbOpT_dump(this)
        class(CmplxmscbOpT), intent(in) :: this
        integer :: i, j
        
        write (*,*) "method: ", this%fe%method
        write (*,*) "colors: ", this%fe%stages(1)%nrofcols

    end subroutine
    
    subroutine CmplxmscbOpT_dealloc(this)
        class(CmplxmscbOpT), intent(inout) :: this
        
        call this%fe%dealloc()
    end subroutine
    
    subroutine CmplxEulermscbOpT_init(this, Op_T)
        use Operator_mod
        use graphdata_mod
        implicit none

        class(CmplxEulermscbOpT) :: this
        Type(Operator), intent(in) :: Op_T
        Integer :: i, k
        type(GraphData) :: gd
        real(kind=kind(0.D0)), allocatable, dimension(:) :: diags
        Complex(kind=kind(0.D0)), allocatable, dimension(:,:) :: tmp

        this%Zero = 1.E-12
        this%Ndim_hop = Op_T%N
        this%g = Op_T%g

        call Op_T_to_graphdata(Op_T, gd, diags)

        ! some sanity checks and status informations
        call determine_used_colors_of_graph(gd)
        write (*,*) "Nr edges: ", gd%nredges
        if (gd%usedcolors == gd%deltag) then
            write(*,*) "Maximum Degree", gd%deltag, ". Found", gd%usedcolors," Families -> optimal decomposition"
        else
            write(*,*) "Maximum Degree", gd%deltag, ". Found", gd%usedcolors," Families"
        endif

        this%ee = createEulerExponentialfromGraphData(gd, diags)

        ! check wether it is supported behaviour
        do i = 1, size(Op_T%P)
        if (Op_T%P(i) /= i) then
            write (*,*) "P unsupported."
        endif
        enddo
        this%P = Op_T%P
        call dealloc_graphdata(gd)
        deallocate(diags)
    end subroutine

    subroutine CmplxEulermscbOpT_reverseadjoint(this, arg)
        class(CmplxEulermscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg

        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%ee%reverseadjoint(arg)
        Endif

    end subroutine

    subroutine CmplxEulermscbOpT_adjoint(this, arg)
        class(CmplxEulermscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg

        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%ee%adjoint(arg)
        Endif

    end subroutine

    subroutine CmplxEulermscbOpT_adjoint_over_two(this, arg)
        class(CmplxEulermscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg

        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%ee%adjoint_over_two(arg)
        Endif

    end subroutine
    
    subroutine CmplxEulermscbOpT_rmult(this, arg)
        class(CmplxEulermscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        
        !FIXME: P
        
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%ee%rmult(arg)
        Endif
    end subroutine

    subroutine CmplxEulermscbOpT_rmultinv(this, arg)
        class(CmplxEulermscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        
        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%ee%rmultinv(arg)
        Endif
    end subroutine
    
    subroutine CmplxEulermscbOpT_lmult(this, arg)
        class(CmplxEulermscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg
        
        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%ee%lmult(arg)
        Endif
    end subroutine

    subroutine CmplxEulermscbOpT_lmultinv(this, arg)
        class(CmplxEulermscbOpT), intent(in) :: this
        Complex(kind=kind(0.D0)), intent(inout), dimension(:,:) :: arg

        !FIXME: P
        If ( dble(this%g*conjg(this%g)) > this%Zero ) then
            call this%ee%lmultinv(arg)
        Endif
    end subroutine

    subroutine CmplxEulermscbOpT_dump(this)
        class(CmplxEulermscbOpT), intent(in) :: this
        integer :: i, j
        
        write (*,*) "colors: ", this%ee%nrofcols

    end subroutine
    
    subroutine CmplxEulermscbOpT_dealloc(this)
        class(CmplxEulermscbOpT), intent(inout) :: this
        
        call this%ee%dealloc()
    end subroutine
    
end module mscbOpT_mod
