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


       SUBROUTINE CONFOUT

!--------------------------------------------------------------------
!
!> @brief 
!> Prints out the field configuration and seeds so as to be able to 
!> pursue the run.
!
!--------------------------------------------------------------------

         USE HAMILTONIAN
         USE F95ZLIB
         USE IOPORTS
         USE ISO_C_BINDING
#ifdef MPI
         Use mpi
#endif
         IMPLICIT NONE

         ! LOCAL
         INTEGER        :: I, IERR, ISIZE, IRANK, K, NT
         INTEGER, DIMENSION(:), ALLOCATABLE :: SEED_VEC
         CHARACTER (LEN=64) :: FILE_TG
         TYPE(IOPORT) :: fd
         CHARACTER(LEN=255), TARGET :: LINE
         INTEGER        :: IOS

#if defined(MPI)
         INTEGER        :: irank_g, isize_g, igroup
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
         call MPI_Comm_rank(Group_Comm, irank_g, ierr)
         call MPI_Comm_size(Group_Comm, isize_g, ierr)
         igroup           = irank/isize_g
         !Write(6,*) "Group, rank :", igroup, irank_g

         CALL GET_SEED_LEN(K)
         ALLOCATE(SEED_VEC(K))
         CALL RANGET(SEED_VEC)
          
       
#if defined(MPI)
#if defined(TEMPERING) 
         write(FILE_TG,'(A,I0,A,I0,A)') "Temp_",igroup,"/confout_",irank_g,".gz"
#else
         write(FILE_TG,'(A,I0,A)') "confout_",irank_g,".gz"
#endif
#else
         FILE_TG = "confout_0.gz"
#endif

         CALL FGZ_OPEN(TRIM(ADJUSTL(FILE_TG)),'w6',fd,ios)
         WRITE(line,*) SEED_VEC
         CALL FGZ_WRITE(fd,TRIM(LINE),'yes',IOS)
         DO NT = 1,LTROT
            DO I = 1,SIZE(NSIGMA,1)
               WRITE(line,*) NSIGMA(I,NT)
               CALL FGZ_WRITE(fd,TRIM(LINE),'yes',IOS)
            ENDDO
         ENDDO
         CALL FGZ_CLOSE(fd,IOS)
         
         DEALLOCATE(SEED_VEC)

       END SUBROUTINE CONFOUT
