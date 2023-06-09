!  Copyright (C) 2023 The ALF project
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
 
Module runtime_error_mod
    
  !--------------------------------------------------------------------
  !> @author 
  !> ALF-project
  !
  !> @brief 
  !> This module handles runtime error codes and terminates MPI-save.
  !
  !--------------------------------------------------------------------

#if defined(MPI)
          Use mpi
#endif
          use iso_fortran_env, only: output_unit, error_unit

    !========================================================================

          !--------------------------------------------------------------------
          !> @brief
          !> This type defines the error codes.
          !> The error codes are defined as an enum type.
          !> The error codes can be extended by adding new enumerators.
          !--------------------------------------------------------------------
          ENUM, BIND(C)
            ENUMERATOR :: ERROR_NONE = 0
            ENUMERATOR :: ERROR_GENERIC
            ENUMERATOR :: ERROR_FILE_NOT_FOUND
            ENUMERATOR :: ERROR_UNSTABLE_MATRIX
            ENUMERATOR :: ERROR_MISSING_OBS
            ENUMERATOR :: ERROR_FIELDS
            ENUMERATOR :: ERROR_HAMILTONIAN
            ENUMERATOR :: ERROR_GLOBAL_UPDATES
            ENUMERATOR :: ERROR_MAXENT
            ENUMERATOR :: ERROR_RUNNING_FILE_FOUND
            ! TO BE EXTENDED: I'm using the generic error code for now
          END ENUM
    
          contains

            !--------------------------------------------------------------------
            !> @brief
            !> This subroutine terminates the program with an error message.
            !> This termination is MPI-save. 
            !> This is useful for debugging.
            !>
            !> @param error_code
            !> The error code
            !--------------------------------------------------------------------
            Subroutine Terminate_on_error(error_code, filename, linenum)
              Implicit none
              Integer, INTENT(IN) :: error_code
              character(len=*), intent(in) :: filename
              integer,          intent(in) :: linenum
              
#if defined(MPI)
              Integer             :: ierr
#endif
              Character(len=100)  :: error_message
              
              Select Case(error_code)
                Case(ERROR_NONE)
                  error_message = "No error"
                Case(ERROR_GENERIC)
                  error_message = "Generic error"
                Case(ERROR_FILE_NOT_FOUND)
                  error_message = "File not found"
                Case(ERROR_UNSTABLE_MATRIX)
                  error_message = "Unstable matrix"
                Case(ERROR_MISSING_OBS)
                  error_message = "Missing observables"
                Case(ERROR_FIELDS)
                  error_message = "Error in field setup"
                Case(ERROR_GLOBAL_UPDATES)
                  error_message = "Error in global updates"
                Case(ERROR_HAMILTONIAN)
                  error_message = "Error in Hamiltonian setup"
                Case(ERROR_MAXENT)
                  error_message = "Error in maximum entropy method"
                Case(ERROR_RUNNING_FILE_FOUND)
                  error_message = "" ! The error message is printed in the calling routine
                Case DEFAULT
                  error_message = "Unknown error"
              End Select
              

              if (error_code /= ERROR_NONE) then
                write(error_unit,'(" ",A," (",I0,"):")') trim(filename), linenum
                Write(error_unit,*) error_message
                Write(error_unit,*) "Terminating program"
#if !defined(MPI)
                error stop error_code
#else
                call MPI_ABORT(MPI_COMM_WORLD,error_code,ierr)
#endif
              end if
              
            end Subroutine Terminate_on_error
            
end Module runtime_error_mod
          
    