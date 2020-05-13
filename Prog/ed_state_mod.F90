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


MODULE ed_state_mod
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ed_state


    TYPE ed_state
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>   
!> @brief 
!> Defines a state in the Fock space, for exact diagonalisation.  
        private
        complex(kind=kind(0.d0)) :: factor
        integer :: N_orbitals, i
        
        CONTAINS
            PROCEDURE :: init => ed_state_init
            PROCEDURE :: set => ed_state_set
            PROCEDURE :: get_i => ed_state_get_i
            PROCEDURE :: get_factor => ed_state_get_factor
            PROCEDURE :: N_fermions => ed_state_N_fermions
            PROCEDURE :: print => ed_state_print
            PROCEDURE :: annihil_e => ed_state_annihil_e
            PROCEDURE :: create_e  => ed_state_create_e
            
            !GENERIC :: ASSIGNMENT(=) => assign
    END TYPE ed_state

CONTAINS

    subroutine ed_state_init(this, N_orbitals)
        IMPLICIT NONE
        class(ed_state), INTENT(INOUT) :: this
        integer, intent(in) :: N_orbitals
        
        this%N_orbitals = N_orbitals
        this%i = 0
    
    end subroutine ed_state_init
    

    subroutine ed_state_set(this, i0, s)
        IMPLICIT NONE
        class(ed_state), INTENT(INOUT) :: this
        integer, intent(in) :: i0
        integer, intent(in), optional :: s
        
        if( present(s) ) then
            this%i = 2**(this%N_orbitals * s) + i0
        else
            this%i = i0
        endif
        this%factor = cmplx( 1.d0, 0.d0, kind=kind(0.d0) )
    
    end subroutine ed_state_set
    

    function ed_state_get_i(this)
        IMPLICIT NONE
        class(ed_state), INTENT(INOUT) :: this
        
        integer :: ed_state_get_i
        
        ed_state_get_i = this%i
    
    end function ed_state_get_i
    

    function ed_state_get_factor(this)
        IMPLICIT NONE
        class(ed_state), INTENT(INOUT) :: this
        
        complex(kind=kind(0.d0)) :: ed_state_get_factor
        
        ed_state_get_factor = this%factor
    
    end function ed_state_get_factor
    

!     function ed_state_test(this, i0, s)
!         IMPLICIT NONE
!         class(ed_state), INTENT(INOUT) :: this
!         integer, intent(in) :: i0, s
!         logical :: ed_state_test
!         integer :: i_e
!         
!         ed_state_test = btest(this%i, this%N_orbitals*s + i0)
!     
!     end function ed_state_test

    
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
    

    subroutine ed_state_print(this)
        IMPLICIT NONE
        class(ed_state), INTENT(INOUT) :: this
        integer :: i
        character (len=32) :: str
        
        str = ''
        
        print*, this%factor
        do i=0, 31
            if ( btest(this%i, i) ) then
                write(str,'(A,A)') trim(str), '1'
            else
                write(str,'(A,A)') trim(str), '0'
            endif
        enddo
        print *, str
    
    end subroutine ed_state_print
    

    subroutine ed_state_annihil_e(this, i0, s)
        IMPLICIT NONE
        class(ed_state), INTENT(INOUT) :: this
        integer, intent(in) :: i0, s
        
        integer :: i_e
        
        i_e = this%N_orbitals*s + i0
        
        if( .not. btest(this%i, i_e) ) then
            this%factor = cmplx(0.d0, 0.d0, kind=kind(0.d0))
            return
        endif
        
        if( mod( N_fermions(ibits(this%i,0,i_e)), 2) == 1 ) then
            this%factor = this%factor * cmplx(-1.d0, 0.d0, kind=kind(0.d0))
        endif
        
        this%i = IBCLR(this%i,i_e)
    
    end subroutine ed_state_annihil_e
    

    subroutine ed_state_create_e(this, i0, s)
        IMPLICIT NONE
        class(ed_state), INTENT(INOUT) :: this
        integer, intent(in) :: i0, s
        
        integer :: i_e
        
        i_e = this%N_orbitals*s + i0
        
        if( btest(this%i,i_e) ) then
            this%factor = cmplx(0.d0, 0.d0, kind=kind(0.d0))
            return
        endif
        
        if( mod( N_fermions(ibits(this%i,0,i_e)), 2) == 1 ) then
            this%factor = this%factor * cmplx(-1.d0, 0.d0, kind=kind(0.d0))
        endif
        
        this%i = IBSET(this%i,i_e)
    
    end subroutine ed_state_create_e
     
   END MODULE ed_state_mod
