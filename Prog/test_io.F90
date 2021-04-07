!  Copyright (C) 2016 - 2020 The ALF project
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
!> Reads in the VAR_QMC namelist from the file parameters, calls Ham_set and  carries out the sweeps. If the
!> program is compiled with the Tempering flag on, then the VAR_TEMP namelist will also be read in.
!>



program test_io
  Use Fields_mod
#ifdef MPI
  Use mpi
#endif
  Implicit none

  Interface
     Subroutine Set_Random_number_Generator(File_seeds, Seed_in)
       Character (LEN=64), Intent(IN) :: File_seeds
       Integer,  Intent(out) :: SEED_IN
     end Subroutine Set_Random_number_Generator
  end Interface
  
  Type (Fields) :: nsigma
  Integer       :: N_op
  Integer       :: Ltrot
  Integer       :: Group_Comm
  Character (len=64) :: File_seeds
  Integer :: Seed_in
  
  NAMELIST /VAR/ N_op, Ltrot

#ifdef MPI
  Integer :: ierr, Isize, Irank, Irank_g, Isize_g, color, key, igroup, mpi_per_parameter_set

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
  mpi_per_parameter_set = Isize
  color = irank/mpi_per_parameter_set
  key   =  0
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Group_comm, ierr)
  call MPI_Comm_rank(Group_Comm, Irank_g, ierr)
  call MPI_Comm_size(Group_Comm, Isize_g, ierr)
  igroup           = irank/isize_g
#endif
  
  OPEN(UNIT=5, FILE='parameters_io', STATUS='old', ACTION='read')
  READ(5, NML=VAR)
  
  call nsigma%make(N_op, Ltrot)
  
  File_seeds="seeds"
  Call Set_Random_number_Generator(File_seeds, Seed_in)
         
  Call nsigma%in(Group_Comm)
  
  Call nsigma%out(Group_Comm)
  
#ifdef MPI
  CALL MPI_FINALIZE(ierr)
#endif
        
end program