!  Copyright (C) 2016 The ALF project
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


       SUBROUTINE CONFOUT

!--------------------------------------------------------------------
!
!> @brief 
!> Prints out the field configuration and seeds so as to be able to 
!> pursue the run.
!
!--------------------------------------------------------------------

         USE HAMILTONIAN

         IMPLICIT NONE

         

#ifdef MPI
         INCLUDE 'mpif.h'
#endif
         ! LOCAL
         INTEGER        :: I, IERR, ISIZE, IRANK, SEED_IN, K, ISEED, NT, NR
         INTEGER, DIMENSION(:), ALLOCATABLE :: SEED_VEC
         REAL (Kind=Kind(0.d0))  :: X
         LOGICAL        :: LCONF 
         CHARACTER (LEN=64) :: FILE_SR, FILE_TG

#if defined(MPI)
         INTEGER        :: STATUS(MPI_STATUS_SIZE)
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif 
         
#if defined(MPI) && !defined(TEMPERING)
         CALL GET_SEED_LEN(K)
         ALLOCATE(SEED_VEC(K))
         CALL RANGET(SEED_VEC)
         FILE_SR = "confout"
         FILE_TG = FILE_I(FILE_SR,IRANK)
         OPEN (UNIT = 10, FILE=FILE_TG, STATUS='UNKNOWN', ACTION='WRITE')
         WRITE(10,*) SEED_VEC
         DO NT = 1,LTROT
            DO I = 1,SIZE(NSIGMA,1)
               WRITE(10,*) NSIGMA(I,NT) 
            ENDDO
         ENDDO
         CLOSE(10)
         DEALLOCATE(SEED_VEC)
            
#else   
         CALL GET_SEED_LEN(K)
         ALLOCATE(SEED_VEC(K))
         CALL RANGET(SEED_VEC)
         FILE_TG = "confout_0"
#if defined(TEMPERING)
         write(File_TG,'(A,I0,A)') "Temp_",Irank,"/confout_0"
#else
         FILE_TG = "confout_0"
#endif
         OPEN (UNIT = 10, FILE=FILE_TG, STATUS='UNKNOWN', ACTION='WRITE')
         WRITE(10,*) SEED_VEC
         DO NT = 1,LTROT
            DO I = 1,SIZE(NSIGMA,1)
               WRITE(10,*) NSIGMA(I,NT) 
            ENDDO
         ENDDO
         CLOSE(10)
         DEALLOCATE(SEED_VEC)
#endif
         
       END SUBROUTINE CONFOUT
