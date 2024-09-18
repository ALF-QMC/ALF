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

Module Files_mod
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Manipulation of character strings
!> 
!
!--------------------------------------------------------------------
   contains

     Character (len=64) function File_i( file, I)
        character (len=64) :: file
        integer            :: i
        write(File_i,'(A,"_",I0)') trim(file),i
      end function File_i

     Character (len=64) function File_add( file, file1)
        character (len=64) :: file, file1
        write(File_add,'(A,A)') trim(file),Trim(file1)
      end function File_add
      
      function str_to_upper(strIn) result(strOut)
           implicit none
      
           character(len=*), intent(in) :: strIn
           character(len=len(trim(strIn))) :: strOut
           integer :: i,j
      
           do i = 1, len(trim(strIn))
                j = iachar(strIn(i:i))
                if (j>= iachar("a") .and. j<=iachar("z") ) then
                     strOut(i:i) = achar(iachar(strIn(i:i))-32)
                else
                     strOut(i:i) = strIn(i:i)
                end if
           end do
      end function str_to_upper

end Module Files_mod
