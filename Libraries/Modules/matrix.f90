
!  Copyright (C) 2016 The ALF project
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
     MODULE Matrix

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Definition of a matrix type. 
!> 
!--------------------------------------------------------------------

       
       Type Mat_C
          complex (Kind=Kind(0.d0)), pointer :: el(:,:)
          Integer :: dim
       end Type Mat_C

       Type Mat_R
          Real (Kind=Kind(0.d0)), pointer :: el(:,:)
          Integer :: dim
       end Type Mat_R

       Interface Make_Mat
          module procedure constructor_C, constructor_R
       end Interface
       Interface Clear_Mat
          module procedure Destroy_C, Destroy_R
       end Interface

       Contains
         subroutine constructor_C(Mat,N)
           type (Mat_C) :: Mat
           Integer :: N
           allocate (Mat%el(N,N))
           Mat%el = cmplx(0.0,0.0)
           Mat%dim = N
         end subroutine constructor_C

         subroutine constructor_R(Mat,N)
           type (Mat_R) :: Mat
           Integer :: N
           allocate (Mat%el(N,N))
           Mat%el = 0.0
           Mat%dim = N
         end subroutine constructor_R

         subroutine Destroy_C(Mat)
           type (Mat_C) :: Mat
           deallocate (Mat%el)
         end subroutine Destroy_C

         subroutine Destroy_R(Mat)
           type (Mat_R) :: Mat
           deallocate (Mat%el)
         end subroutine Destroy_R
     end MODULE Matrix



