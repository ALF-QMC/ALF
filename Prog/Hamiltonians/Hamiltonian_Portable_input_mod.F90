!  Copyright (C) 2016 - 2022 The ALF project
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
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
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
!> Reads in the parameters for the construction of the Hamiltonian from the files geometry.txt (Ham_latt), hoppings.txt (Ham_hop) and potentials.txt (Ham_V) and for the observables.
!>
!--------------------------------------------------------------------

module Hamiltonian_Portable_input_mod

   Use runtime_error_mod

   implicit none

   contains

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reads in  the  Lattice parameters
!--------------------------------------------------------------------

        Subroutine read_latt(L1, l2, Norb, a_p, Orb_pos, Group_Comm)

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif

          implicit none

          Real (Kind=Kind(0.d0)), intent(out)  :: a_p(3,3)
          real (kind=kind(0.d0)), allocatable, intent(out) :: Orb_pos(:,:)
          integer, intent(out) :: L1, L2, Norb
          integer, intent(in)  :: Group_Comm
 
          integer                :: ierr, unit_latt, no, i
          Character (len=64)     :: file_latt
          real (kind=kind(0.d0)) :: x, y, z

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g

          If (Irank_g == 0) then
#endif
             File_latt = "geometry.txt"
#if defined(TEMPERING)
             write(File_latt,'(A,I0,A)') "Temp_",igroup,"/geometry.txt"
#endif
             Open(newunit=unit_latt, file=file_latt, status="old", action="read", iostat=ierr)
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'unable to open <geometry.txt>', ierr
                Call Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
             END IF

             read(unit_latt,*) L1, L2
             do i = 1, size(a_p,1)
                read(unit_latt,*) x, y, z
                a_p(i,1) = x; a_p(i,2) = y; a_p(i,3) = z
             enddo
             read(unit_latt,*) norb
             allocate(Orb_pos(Norb,3))
             do no = 1, Norb
                read(unit_latt,*) i, x, y, z
                Orb_pos(no,1) = x; Orb_pos(no,2) = y; Orb_pos(no,3) = z
             enddo

             Close(unit_latt)
#ifdef MPI
          Endif

          CALL MPI_BCAST(L1       ,  1           ,MPI_INTEGER  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(L2       ,  1           ,MPI_INTEGER  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(a_p      ,  size(a_p)   ,MPI_REAL8    ,0,Group_Comm,ierr)
          CALL MPI_BCAST(Norb     ,  1           ,MPI_INTEGER  ,0,Group_Comm,ierr)
          if (irank_g /= 0) allocate(Orb_pos(Norb,3))
          CALL MPI_BCAST(Orb_pos,  size(Orb_pos) ,MPI_REAL8    ,0,Group_Comm,ierr)
#endif

        end subroutine read_latt

end module Hamiltonian_Portable_input_mod
