!  Copyright (C) 2019 The ALF project
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

Program repack_latt
  
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Program that reads in Lattice-type observable given by first command line argument 
!> and writes it to file given by second command line argument
!> 
!
!--------------------------------------------------------------------
  
  use ana_mod
#ifdef _OPENMP
  use check_omp_num_threads_mod
#endif
  Implicit none
  
  INTEGER :: i, nb, no, no1, nt
  Character (len=64) :: File_in, File_out
  
  Real    (Kind=Kind(0.d0)), allocatable :: phase(:)
  Complex (Kind=Kind(0.d0)), pointer     :: bins(:,:,:,:,:), bins0(:,:)
  Complex (Kind=Kind(0.d0)), pointer     :: bin(:,:,:,:), bin0(:)
  Integer                                :: Norb, Nunit, Ntau, Nbins
  Real    (Kind=Kind(0.d0))              :: dtau, X_p(2)
  Type (Lattice)    :: Latt

#ifdef _OPENMP
  call check_omp_num_threads()
#endif
  
  i = 1
  CALL GET_COMMAND_ARGUMENT(i, File_in)
  i = 2
  CALL GET_COMMAND_ARGUMENT(i, File_out)
  
  write(*,*) "reading from ", File_in
  call read_latt(File_in, phase, bins, bins0, Latt, dtau)
  
  Nunit  = size(bins, 1) 
  Ntau   = size(bins, 2)
  Norb   = size(bins, 3)
  Nbins  = size(bins, 5)
  
  write(*,*) "writing to ", File_out
  
  Open (Unit=10,File=File_out, status="unknown",  position="append")
  do nb = 1,Nbins
    Write(10,*) phase(nb), Norb, Latt%N, Ntau, dtau, &
                & size(Latt%L1_p), Latt%L1_p(:), Latt%L2_p(:), Latt%a1_p(:), Latt%a2_p(:)
    Do no = 1,Norb
      Write(10,*) bins0(no,nb)
    enddo
    do I = 1,Latt%N
      x_p = dble(Latt%listk(i,1))*Latt%b1_p + dble(Latt%listk(i,2))*Latt%b2_p  
      Write(10,*) X_p(1), X_p(2)
      Do nt = 1,Ntau
        do no = 1,Norb
          do no1 = 1,Norb
            Write(10,*) bins(I,nt,no,no1,nb)
          enddo
        enddo
      enddo
    enddo
  enddo
  close(10)
  
end Program
