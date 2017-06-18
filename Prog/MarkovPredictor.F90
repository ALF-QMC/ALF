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
        INTEGER(Kind=Kind(0.d0)), allocatable, Dimension(:, :) :: Pint ! Pint holds the number of times a transition occured
        INTEGER(Kind=Kind(0.d0)), allocatable, Dimension(:) :: sums ! this holds the row sums 
        INTEGER, Dimension(:), allocatable :: previousMeasurements ! This is a vector containing previous occurences of the symbols
        Integer, allocatable, Dimension(:) :: states ! This is a vector containing the actual states that the caller employs
        Integer :: order, nrstates, nrofpreviousmeasurements, mcstates
        INTEGER :: nrofupdates, historyidx

        CONTAINS
!            PROCEDURE :: alloc => alloc_UDV_state
            PROCEDURE :: init => init_MarkovPredictor
            PROCEDURE :: update => update_MarkovPredictor
            PROCEDURE :: predict => predict_MarkovPredictor
            PROCEDURE :: mapstatetoindex => mapstatetoindex_MarkovPredictor
            PROCEDURE :: updatehistory => updatehistory_MarkovPredictor
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
Allocate (this%P(nrstates, this%mcstates), this%previousMeasurements(order), this%Pint(nrstates, this%mcstates))
    ALLOCATE(this%sums(nrstates))
    this%previousMeasurements(1) = 1
    this%nrofupdates = nrstates
    this%sums = this%mcstates
    this%historyidx = 1! some intial value
    this%Pint = 1 ! All states are equal initially
END SUBROUTINE init_MarkovPredictor

! This function updates the history array and sets the associated idx
SUBROUTINE updatehistory_MarkovPredictor(this, state)
    IMPLICIT NONE
    CLASS(MarkovPredictor), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: state
    INTEGER :: I, stateidx
    this%historyidx = this%nrstates**(this%order-1) * this%mapstatetoindex(state)
    ! now let's update the vector of previously drawn samples
    do i = 1, this%order - 1
        this%previousMeasurements(i+1) = this%previousMeasurements(i)
        this%historyidx = this%historyidx + this%nrstates**(this%order - I - 1) * this%mapstatetoindex(this%previousMeasurements(I))
    enddo
    this%previousMeasurements(1) = state
    
!     ! mapping essentially coresponds to number conversion...
!     DO I = 1, this%order-1
!         stateidx = stateidx + this%nrstates**(this%order - I - 1) * this%mapstatetoindex(this%previousMeasurements(I))
!     ENDDO
!     gethistoryidx_MarkovPredictor = stateidx
END SUBROUTINE updatehistory_MarkovPredictor

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
    
    this%nrofupdates = this%nrofupdates + 1 ! update the amount of received data
    stateidx = this%mapstatetoindex(state) ! get the idx of the current state. Will correspond to column.
    ! lastidx contains the row that corresponds to the state that is stored in previous measurements.
    this%sums(this%historyidx) = this%sums(this%historyidx) + 1 ! normalization changes
    this%Pint(this%historyidx, stateidx) = this%Pint(this%historyidx, stateidx) + 1 ! a new transition has occured
    ! now let's update the vector of previously drawn samples
    CALL this%updatehistory(state)
!     do i = 1, this%order - 1
!         this%previousMeasurements(i+1) = this%previousMeasurements(i)
!     enddo
!     this%previousMeasurements(1) = state
END SUBROUTINE update_MarkovPredictor

! Draw a realization from transition matrix given the current state
INTEGER FUNCTION predict_MarkovPredictor(this, curstate)
IMPLICIT NONE
    CLASS(MarkovPredictor), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: curstate
    Integer(Kind=Kind(0.D0)) :: mysum, myrand
    INTEGER :: stateidx, I
    REAL :: mynum
    !CALL RANDOM_NUMBER(mynum)
    mynum = rand()
    ! DO I = 1, this%nrstates
    ! IF(this%states(I) == curstate) stateidx = I
    ! ENDDO
    ! obtain row that coresponds to that state
    ! The row that corresponds to the history is stored in historyidx
    !write (*,*) mynum
    myrand = mynum * this%sums(this%historyidx)
    ! find out which state we will predict:
    mysum = this%Pint(this%historyidx, 1)
    !write (*,*) myrand, stateidx
    I = 1
    DO WHILE (mysum < myrand)
    I = I + 1
    mysum = mysum + this%Pint(this%historyidx, I)
    ENDDO
    predict_MarkovPredictor = this%states(I)
END FUNCTION predict_MarkovPredictor

SUBROUTINE print_MarkovPredictor(this)
IMPLICIT NONE
    CLASS(MarkovPredictor), INTENT(INOUT) :: this
    INTEGER :: I
    write (*,*) "order: ", this%order, "Nr of states: ", this%nrstates, "states", this%states
    write (*, *) 0, this%states
    DO I = 1, this%mcstates
        write (*, *) this%states(MOD(I, this%nrstates)), this%Pint(I, :)/DBLE(this%sums(I))
    ENDDO
END SUBROUTINE print_MarkovPredictor

END MODULE MarkovPredictor_Mod
