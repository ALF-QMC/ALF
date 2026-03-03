!  Copyright (C) 2016-2026 The ALF project
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
!> @author ALF-project
!> @brief String manipulation utilities for file naming and text processing.
!> 
!> @details
!> Provides utility functions for:
!> - Appending integers to filenames (for numbered output files)
!> - Concatenating strings/filenames
!> - Case conversion (uppercase)
!>
!> Used throughout ALF for generating output file names with run indices,
!> parameter values, or other identifiers.
!
!--------------------------------------------------------------------
   contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Append integer to filename with underscore separator.
!
!> @details
!> Creates numbered filename: "basename_N" where N is formatted without
!> leading zeros. Commonly used for:
!> - Numbered output files (data_1, data_2, ...)
!> - Run indices
!> - Parameter sweeps
!>
!> Example: File_i("output", 42) → "output_42"
!
!> @param[in] file Base filename (up to 64 characters)
!> @param[in] I Integer to append
!> @return Concatenated string "file_I"
!
!> @note Uses I0 format (no leading zeros/blanks)
!--------------------------------------------------------------------
     Character (len=64) function File_i( file, I)
        character (len=64) :: file   ! Base filename
        integer            :: i       ! Index to append
        write(File_i,'(A,"_",I0)') trim(file),i
      end function File_i

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Concatenate two strings/filenames.
!
!> @details
!> Joins two strings by trimming whitespace and concatenating.
!> Useful for:
!> - Building file paths (directory + filename)
!> - Adding extensions
!> - Combining filename components
!>
!> Example: File_add("output", ".dat") → "output.dat"
!
!> @param[in] file First string component
!> @param[in] file1 Second string component
!> @return Concatenated string "file||file1" (trimmed)
!--------------------------------------------------------------------
     Character (len=64) function File_add( file, file1)
        character (len=64) :: file, file1   ! Strings to concatenate
        write(File_add,'(A,A)') trim(file),Trim(file1)
      end function File_add
      
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Convert string to uppercase.
!
!> @details
!> Converts all lowercase ASCII letters (a-z) to uppercase (A-Z).
!> Non-alphabetic characters (numbers, punctuation, etc.) remain unchanged.
!>
!> Uses ASCII character codes:
!> - Lowercase: 'a'(97) to 'z'(122)
!> - Uppercase: 'A'(65) to 'Z'(90)
!> - Difference: 32
!>
!> Typical use: normalizing command-line arguments or keywords for
!> case-insensitive comparison.
!>
!> Example: str_to_upper("QMC_test") → "QMC_TEST"
!
!> @param[in] strIn Input string (arbitrary length)
!> @return Uppercase version of input string (trimmed length)
!
!> @note Only handles ASCII characters; extended character sets unchanged
!--------------------------------------------------------------------
      function str_to_upper(strIn) result(strOut)
           implicit none
      
           character(len=*), intent(in) :: strIn
           character(len=len(trim(strIn))) :: strOut
           integer :: i,j
      
           ! Convert each character
           do i = 1, len(trim(strIn))
                j = iachar(strIn(i:i))           ! ASCII code of current character
                if (j>= iachar("a") .and. j<=iachar("z") ) then
                     ! Lowercase letter: convert to uppercase (subtract 32)
                     strOut(i:i) = achar(iachar(strIn(i:i))-32)
                else
                     ! Not lowercase: keep unchanged
                     strOut(i:i) = strIn(i:i)
                end if
           end do
      end function str_to_upper

end Module Files_mod
