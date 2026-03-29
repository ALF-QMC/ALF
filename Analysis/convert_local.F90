!  Copyright (C) 2022 The ALF project
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

Program convert_local
  
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Program that reads in Local-type observable given by first command line argument
!> and writes it to HDF5 file given by second command line argument
!>
!
!--------------------------------------------------------------------
  
  Use Lattices_v3, only: Unit_cell, Lattice
  use ana_mod, only: read_local, write_obs_local_hdf5
#ifdef _OPENMP
  use check_omp_num_threads_mod
#endif
  Implicit none
  
  INTEGER :: i
  Character (len=64) :: File_in, File_h5
  
  Real    (Kind=Kind(0.d0)), allocatable :: sgn(:)
  Complex (Kind=Kind(0.d0)), pointer     :: bins(:,:,:,:)
  Type (Lattice)                         :: Latt
  Type (Unit_cell)                       :: Latt_unit
  Real    (Kind=Kind(0.d0))              :: dtau
  Character (len=:),         allocatable :: Channel
#ifdef _OPENMP
  call check_omp_num_threads()
#endif
  
  i = 1
  CALL GET_COMMAND_ARGUMENT(i, File_in)
  i = 2
  CALL GET_COMMAND_ARGUMENT(i, File_h5)
  
  write(*,*) "reading from ", File_in
  call read_local(File_in, sgn, bins, Latt, Latt_unit, dtau, Channel)
  
  write(*,*) "writing to ", File_h5
  call write_obs_local_hdf5(File_h5, File_in, sgn, bins, Latt, Latt_unit, dtau, Channel)
  deallocate( sgn, bins )
  
end Program
