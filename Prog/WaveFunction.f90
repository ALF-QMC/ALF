!  Copyright (C) 2016, 2017 The ALF project
! 
!  This file is part of the ALF project.
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
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

Module WaveFunction_mod

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Defines the Operator type, and provides a number of operations on this
!> type. 
!
!--------------------------------------------------------------------

  
  Implicit none

  
  Type WaveFunction
     complex (Kind=Kind(0.d0)), allocatable :: P(:,:)
     ! P is an N X Ndim matrix containing the trial wave function as a slater determinant
     ! P is an Ndim x N_part matrix representing N_part particles in a Ndim dimensional space
     ! The columns have to be linearly independent
  end type WaveFunction

  
Contains
  
!--------------------------------------------------------------------

  Pure subroutine WF_alloc(WF, Ndim, N_part)
    Implicit none
    Type (WaveFunction), intent(INOUT) :: WF
    Integer, Intent(IN) :: Ndim, N_part
    Allocate (WF%P(Ndim, N_part))
    WF%P = cmplx(0.d0, 0.d0, kind(0.D0))
  end subroutine WF_alloc

!--------------------------------------------------------------------

  Pure subroutine WF_clear(WF)
    Implicit none
    Type (WaveFunction), intent(INOUT) :: WF
    Deallocate (WF%P)
  end subroutine WF_clear 

end Module WaveFunction_mod
