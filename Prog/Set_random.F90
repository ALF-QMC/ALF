!  Copyright (C) 2016 - 2018 The ALF project
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


!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This  routine reads the seeds from the file File_seeds,  distributes them to mpi  processes and
!> sets the random number generator.
!
!--------------------------------------------------------------------

     Subroutine Set_Random_number_Generator(File_seeds,Seed_in)

       
#ifdef MPI
        Use mpi
#endif
        Use Random_Wrap

        Use iso_fortran_env, only: output_unit, error_unit

        Implicit none

        Character (LEN=64), Intent(IN) :: File_seeds
        Integer,  Intent(out) :: SEED_IN
        Integer :: I, IERR
        Integer, allocatable :: SEED_VEC(:)
        
#ifdef MPI
        INTEGER        :: STATUS(MPI_STATUS_SIZE), irank_g, isize_g, igroup, ISIZE, IRANK
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
        
#if defined(MPI)       
       IF (IRANK == 0) THEN
          OPEN(UNIT=5,FILE=File_seeds,STATUS='OLD',ACTION='READ',IOSTAT=IERR)
          IF (IERR /= 0) THEN
             WRITE(error_unit,*) 'Fields_in: unable to open <seeds>',IERR
             error stop 1
          END IF
          DO I = ISIZE-1,1,-1
             READ (5,*) SEED_IN
             CALL MPI_SEND(SEED_IN,1,MPI_INTEGER, I, I+1024, MPI_COMM_WORLD,IERR)
          ENDDO
          READ(5,*) SEED_IN
          CLOSE(5)
       ELSE
          CALL MPI_RECV(SEED_IN, 1, MPI_INTEGER,0,  IRANK + 1024,  MPI_COMM_WORLD,STATUS,IERR)
       ENDIF
       ALLOCATE (SEED_VEC(1))
       SEED_VEC(1) = SEED_IN
       CALL RANSET(SEED_VEC)
       DEALLOCATE (SEED_VEC)
#else
       OPEN(UNIT=5,FILE=FILE_seeds,STATUS='OLD',ACTION='READ',IOSTAT=IERR)
       IF (IERR /= 0) THEN
          WRITE(*,*) 'Fields_in: unable to open <seeds>',IERR
          error stop 1
       END IF
       READ (5,*) SEED_IN
       CLOSE(5)
       ALLOCATE(SEED_VEC(1))
       SEED_VEC(1) = SEED_IN
       CALL RANSET (SEED_VEC)
       DEALLOCATE  (SEED_VEC)
#endif
       
     end Subroutine Set_Random_number_Generator
