!  Copyright (C) 2017 The ALF project
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
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
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

MODULE MarkovPredictor_Mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: MarkovPredictor
    TYPE MarkovPredictor
        REAL(Kind=Kind(0.d0)), allocatable, Dimension(:, :) :: P
        INTEGER(Kind=Kind(0.d0)), allocatable, Dimension(:, :) :: Pint
        INTEGER(Kind=Kind(0.d0)), allocatable, Dimension(:) :: sums
        INTEGER, Dimension(:), allocatable :: previousMeasurements
        Integer, allocatable, Dimension(:) :: states
        Integer :: order, nrstates, nrofpreviousmeasurements, mcstates
        INTEGER :: nrofupdates

        CONTAINS
!            PROCEDURE :: alloc => alloc_UDV_state
            PROCEDURE :: init => init_MarkovPredictor
            PROCEDURE :: update => update_MarkovPredictor
            PROCEDURE :: predict => predict_MarkovPredictor
            PROCEDURE :: mapstatetoindex => mapstatetoindex_MarkovPredictor
!             PROCEDURE :: dealloc => dealloc_UDV_state
!             PROCEDURE :: reset => reset_UDV_state
!             PROCEDURE :: assign => assign_UDV_state
            PROCEDURE :: print => print_MarkovPredictor
!            GENERIC :: ASSIGNMENT(=) => assign
    END TYPE MarkovPredictor

CONTAINS
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function initializes the memory of an object.
!>
!> @param [inout] this The object to be modified.
!> @param [in] order The order of the matrix chain
!> @param [in] nrofstates The number of states
!> @param [in] states the set of used states
!-------------------------------------------------------------------
SUBROUTINE init_MarkovPredictor(this, order, nrstates, states)
    IMPLICIT NONE
    CLASS(MarkovPredictor), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: order, nrstates
    INTEGER, INTENT(IN), Dimension(:) :: states
    INTEGER :: I
    
    this%order = order
    this%nrstates = nrstates
    this%nrofpreviousmeasurements = order
    this%states = states
    this%mcstates = nrstates**order
!     DO I = 1, nrstates
!         this%statestoindex()
!     ENDDO
    ALLOCATE(this%P(nrstates, nrstates), this%previousMeasurements(order), this%Pint(nrstates, nrstates))
    ALLOCATE(this%sums(mcstates))
    this%previousMeasurements(1) = 1
    this%nrofupdates = nrstates
    this%sums = nrstates
    
    this%Pint = 1 ! All states are equal
END SUBROUTINE init_MarkovPredictor

INTEGER FUNCTION mapstatetoindex_MarkovPredictor(this, state)
IMPLICIT NONE
    CLASS(MarkovPredictor), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: state
INTEGER :: I, stateidx
DO I = 1, this%nrstates
IF(this%states(I) == state) stateidx = I
ENDDO
mapstatetoindex_MarkovPredictor = stateidx
END FUNCTION mapstatetoindex_MarkovPredictor

! update the transition matrix with new data
SUBROUTINE update_MarkovPredictor(this, state)
IMPLICIT NONE
    CLASS(MarkovPredictor), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: state
INTEGER :: I
INTEGER :: stateidx
this%nrofupdates = this%nrofupdates + 1
! map state to index
! DO I = 1, this%nrstates
! IF(this%states(I) == state) stateidx = I
! ENDDO
stateidx = this%mapstatetoindex(state)
this%sums(this%previousMeasurements(1)) = this%sums(this%previousMeasurements(1)) + 1
this%Pint(this%previousMeasurements(1), stateidx) = this%Pint(this%previousMeasurements(1), stateidx) + 1
this%previousMeasurements(1) = stateidx
END SUBROUTINE update_MarkovPredictor


! Draw a realization from transition matrix given the current state
INTEGER FUNCTION predict_MarkovPredictor(this, curstate)
IMPLICIT NONE
    CLASS(MarkovPredictor), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: curstate
INTEGER :: myrand, stateidx, mysum, I
REAL :: mynum
!CALL RANDOM_NUMBER(mynum)
mynum = rand()
! DO I = 1, this%nrstates
! IF(this%states(I) == curstate) stateidx = I
! ENDDO
stateidx = this%mapstatetoindex(curstate)
!write (*,*) mynum
myrand = mynum * this%sums(stateidx)
! find out which state we will predict:
mysum = this%Pint(stateidx, 1)
!write (*,*) myrand, stateidx
I = 1
DO WHILE (mysum < myrand)
I = I + 1
mysum = mysum + this%Pint(stateidx, I)
ENDDO
predict_MarkovPredictor = this%states(I)
END FUNCTION predict_MarkovPredictor

SUBROUTINE print_MarkovPredictor(this)
IMPLICIT NONE
    CLASS(MarkovPredictor), INTENT(INOUT) :: this
INTEGER :: I
write (*,*) "order: ", this%order, "Nr of states: ", this%nrstates, "states", this%states
write (*, *) 0, this%states
DO I = 1, this%nrstates
write (*, *) this%states(I), this%Pint(I, :, 1)/DBLE(this%sums(I,1))
ENDDO

END SUBROUTINE print_MarkovPredictor

END MODULE MarkovPredictor_Mod
