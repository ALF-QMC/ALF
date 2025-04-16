!  Copyright (C) 2020-2025 The ALF project
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


   Program ana_hdf5
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Program for postprocessing Monte Carlo bins saved in HDF5 file "data.h5".
!> Analyzes all observables in "data.h5" if no command line arguments are supplied,
!> otherwise analyzes the observables given as command line arguments.
!
!--------------------------------------------------------------------
      use hdf5
      use ana_mod
#ifdef _OPENMP
      use check_omp_num_threads_mod
#endif
      implicit none
      Integer                         :: i, hdferr, nargs, storage_type, nlinks, max_corder
      Character (len=64)              :: File_out, File_in, name
      Character (len=64), allocatable :: names(:)
      INTEGER(HID_T)                  :: file_id, group_id
      INTEGER(HSIZE_T)                :: n

#ifdef _OPENMP
      call check_omp_num_threads()
#endif

      File_in = 'data.h5'
      
      CALL h5open_f(hdferr)
      CALL h5fopen_f (File_in, H5F_ACC_RDONLY_F, file_id, hdferr)
      group_id = file_id
      call h5gget_info_f(group_id, storage_type, nlinks, max_corder, hdferr)
      
      nargs = COMMAND_ARGUMENT_COUNT()
      if ( nargs > 0 ) then
         allocate( names(nargs) )
         do i = 1, nargs
            CALL GET_COMMAND_ARGUMENT(i, name)
            names(i) = name
         enddo
      else
         allocate( names(nlinks) )
         do n=0, nlinks-1
            call h5lget_name_by_idx_f(group_id, '/', H5_INDEX_NAME_F, H5_ITER_INC_F, n, name, hdferr)
            names(n+1) = name
         enddo
      endif
      
      do n=1, size(names)
         name = names(n)
         i = len(trim(name)) -4
         if ( i > 1 ) then
           if ( name(i:) == '_scal' ) then
              print *, ''
              print '(A,A)', "analyzing scalar observable ", name
              call Cov_vec(name, File_in)
           endif
         endif
      enddo
      
      do n=1, size(names)
         name = names(n)
         i = len(trim(name)) -2
         if ( i > 1 ) then
           if ( name(i:) == '_eq' ) then
              print *, ''
              print '(A,A)', "analyzing equal time correlation ", name
              call Cov_eq(name, File_in)
           endif
         endif
      enddo
      
      do n=1, size(names)
         name = names(n)
         i = len(trim(name)) -5
         if ( i > 1 ) then
            if ( name(i:) == '_local' ) then
               print *, ''
               print '(A,A)', "analyzing equal time local observables ", name
               call Cov_local(name,File_in)
            endif
         endif
      enddo

      do n=1, size(names)
         name = names(n)
         i = len(trim(name)) -3
         if ( i > 1 ) then
           if ( name(i:) == '_tau' ) then
              print *, ''
              print '(A,A)', "analyzing time displaced correlation ", name
              call Cov_tau(name, File_in)
           endif
         endif
      enddo
      
   end Program ana_hdf5
       
