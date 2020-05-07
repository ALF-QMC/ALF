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


!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Defines a basis spanning the Fock space, for exact diagonalisation
!>
!>
!--------------------------------------------------------------------


MODULE ed_mod
    use Operator_mod
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ed_ham


    TYPE ed_state
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>   
!> @brief 
!> Defines a state in the Fock space, for exact diagonalisation. 
        complex(kind=kind(0.d0)) :: factor
        integer :: i
        !private
        !logical, allocatable :: v(:)
        
        CONTAINS
            PROCEDURE :: N_fermions => ed_state_N_fermions
            PROCEDURE :: annihil_e => ed_state_annihil_e
            PROCEDURE :: create_e  => ed_state_create_e
            
            !GENERIC :: ASSIGNMENT(=) => assign
    END TYPE ed_state


    TYPE ed_ham
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>   
!> @brief 
!> Defines a basis spanning the Fock space, for exact diagonalisation. 
        private
        INTEGER   :: N_orbitals, N_states
        complex(kind=kind(0.d0)), allocatable :: H(:,:)
        CONTAINS
            PROCEDURE :: init => ed_ham_init
            PROCEDURE :: add_op_t => ed_ham_add_op_t
            PROCEDURE :: add_op_v => ed_ham_add_op_v
            PROCEDURE :: build_h => ed_ham_build_h
    END TYPE ed_ham

CONTAINS

    function N_fermions(i)
        IMPLICIT NONE
        integer, intent(in) :: i
        
        integer :: N_fermions
        integer :: n
        
        N_fermions = 0
        do n = 0, 31
            if( btest(i, n) ) N_fermions = N_fermions + 1
        enddo
    
    end function N_fermions
    

    function ed_state_N_fermions(this)
        IMPLICIT NONE
        class(ed_state), INTENT(INOUT) :: this
        
        integer :: ed_state_N_fermions
        
        ed_state_N_fermions = N_fermions(this%i)
    
    end function ed_state_N_fermions
    

    subroutine ed_state_annihil_e(this, i_e)
        IMPLICIT NONE
        class(ed_state), INTENT(INOUT) :: this
        integer, intent(in) :: i_e
        
        if( .not. btest(this%i,i_e) ) then
            this%factor = cmplx(0.d0, 0.d0, kind=kind(0.d0))
            return
        endif
        
        if( mod( N_fermions(ibits(this%i,0,i_e-1)), 2) == 1 ) then
            this%factor = this%factor * cmplx(-1.d0, 0.d0, kind=kind(0.d0))
        endif
        
        this%i = IBCLR(this%i,i_e)
    
    end subroutine ed_state_annihil_e
    

    subroutine ed_state_create_e(this, i_e)
        IMPLICIT NONE
        class(ed_state), INTENT(INOUT) :: this
        integer, intent(in) :: i_e
        
        if( btest(this%i,i_e) ) then
            this%factor = cmplx(0.d0, 0.d0, kind=kind(0.d0))
            return
        endif
        
        if( mod( N_fermions(ibits(this%i,0,i_e-1)), 2) == 1 ) then
            this%factor = this%factor * cmplx(-1.d0, 0.d0, kind=kind(0.d0))
        endif
        
        this%i = IBSET(this%i,i_e)
    
    end subroutine ed_state_create_e
    

    subroutine ed_ham_init(this, N_orbitals)
        IMPLICIT NONE
        class(ed_ham), INTENT(INOUT) :: this
        integer, intent(in) :: N_orbitals
        
        this%N_orbitals = N_orbitals
        this%N_states   = 2**N_orbitals
        allocate( this%H(N_orbitals, N_orbitals) )
        this%H(:,:) = cmplx(0.d0, 0.d0, kind=kind(0.d0))
    
    end subroutine ed_ham_init
    

    subroutine ed_ham_add_op_t(this, op_t, dtau)
        IMPLICIT NONE
        class(ed_ham)  , intent(inout) :: this
        type(Operator), intent(in)    :: Op_T
        real(Kind=Kind(0.d0)), intent(in) :: dtau
        
        integer :: i, n1, n2
        type(ed_state) :: state
        
        do i=0, this%N_states-1
            do n1=1, op_t%N
                do n2=1, op_t%N
                    state%i = i
                    state%factor = op_t%O(n2,n1) * op_t%g / (-dtau)
                    call state%annihil_e( op_t%p(n1)-1 )
                    call state%create_e ( op_t%p(n2)-1 )
                    
                    this%H(state%i+1, i+1) = this%H(state%i+1, i+1) + state%factor
                enddo
            enddo
        enddo
    
    end subroutine ed_ham_add_op_t
    

    subroutine ed_ham_add_op_v(this, op_v, dtau)
        IMPLICIT NONE
        class(ed_ham)  , intent(inout) :: this
        type(Operator), intent(in)    :: op_v
        real(Kind=Kind(0.d0)), intent(in) :: dtau
        
        integer :: i, n1, n2, m1, m2
        type(ed_state) :: state
        
        do i=0, this%N_states-1
            this%H(i+1, i+1) = this%H(i+1, i+1) + op_v%g**2 * op_v%alpha**2 / (-dtau)
            do n1=1, op_v%N
                do n2=1, op_v%N
                    state%i = i
                    state%factor = op_v%O(n2,n1) * op_v%g**2 * op_v%alpha * 2.d0 / (-dtau)
                    call state%annihil_e( op_v%p(n1)-1 )
                    call state%create_e ( op_v%p(n2)-1 )
                    
                    this%H(state%i+1, i+1) = this%H(state%i+1, i+1) + state%factor
                    
                    do m1=1, op_v%N
                        do m2=1, op_v%N
                            state%i = i
                            state%factor = op_v%O(n2,n1) * op_v%O(m2,m1) * op_v%g**2 / (-dtau)
                            call state%annihil_e( op_v%p(m1)-1 )
                            call state%create_e ( op_v%p(m2)-1 )
                            call state%annihil_e( op_v%p(n1)-1 )
                            call state%create_e ( op_v%p(n2)-1 )
                            
                            this%H(state%i+1, i+1) = this%H(state%i+1, i+1) + state%factor
                        enddo
                    enddo
                enddo
            enddo
        enddo
    
    end subroutine ed_ham_add_op_v
    
    
          
    subroutine ed_ham_build_h(this, ndim, OP_T, OP_V, dtau)
        IMPLICIT NONE
        class(ed_ham)  , intent(inout) :: this
        integer, intent(in) :: ndim
        Type (Operator), intent(in) :: Op_V(:,:)
        Type (Operator), intent(in) :: Op_T(:,:)
        real(Kind=Kind(0.d0)), intent(in) :: dtau
        
        
        INTEGER :: nf, n
        print*, "Building ED-Hamiltonian"
        call this%init(2**ndim)
        do nf=1,1
            do n=1, Size(OP_T,1)
                call this%add_op_t( OP_T(n,nf), dtau )
            enddo
            do n=1, Size(OP_V,1)
                call this%add_op_v( OP_V(n,nf), dtau )
            enddo
        enddo
    
    end subroutine ed_ham_build_h
          
     
   END MODULE ed_mod
   
! program test
!     use ed_mod
!     IMPLICIT NONE
!     integer :: i
!     integer :: N_orbitals, N_fermions, N_states
!     type(ed_state) :: state
!     
!     N_orbitals = 20
!     N_fermions = 5
!     
!     N_states = 0
!     
!     do i=0, 2**N_orbitals-1
!         state%i = i
!         !print*, i, state%N_fermions()
!         if (state%N_fermions() == N_fermions) N_states = N_states + 1
! !         do n=0, 31
! !             print*,n, BTEST(i,n)
! !         enddo
!     enddo
!     print*, N_states
! 
! end program
