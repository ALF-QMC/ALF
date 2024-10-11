!  Copyright (C) 2020-2021 The ALF project
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


#if defined(HDF5)
     Module alf_hdf5
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Helper subroutines for using ALF with HDF5.
!
!--------------------------------------------------------------------
       use runtime_error_mod
       use iso_fortran_env, only: output_unit, error_unit
       USE ISO_C_BINDING

       Use hdf5
       use h5lt
       
       Use Lattices_v3
       
       implicit none
       private
       public :: write_attribute, read_attribute, test_attribute, &
         init_dset, append_dat, write_latt, write_comment
     
       interface write_attribute
         MODULE PROCEDURE write_attribute_double, write_attribute_int, write_attribute_string, write_attribute_logical
       end interface write_attribute
       interface read_attribute
         MODULE PROCEDURE read_attribute_double, read_attribute_int, read_attribute_string, read_attribute_logical
       end interface read_attribute
       interface test_attribute
         MODULE PROCEDURE test_attribute_double, test_attribute_int, test_attribute_string, test_attribute_logical
       end interface test_attribute
  
       contains
  
         Subroutine init_dset(file_id, dsetname, dims, is_complex, chunklen)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This subroutine creates a new dataset in an opened HDF5 file for the
!> purpsose of filling it with Monte-Carlo bins.
!
!> @param [IN] file_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of the opened HDF5 file
!> \endverbatim
!> @param [IN] dsetname Character(len=64)
!> \verbatim
!>  Name of the new dataset
!> \endverbatim
!> @param [IN] dims(:) INTEGER(HSIZE_T)
!> \verbatim
!>  Shape of one bin. Whith size(dims) = size(bin)+1 and dims(size(dims)) = 0
!> \endverbatim
!> @param [IN] is_complex logical
!> \verbatim
!>  True if values to be stored can be complex.
!> \endverbatim
!> @param [IN] chunklen INTEGER(HSIZE_T), optional
!> \verbatim
!>  Size of data chunks in number of bins, default = 1
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           
           INTEGER(HID_T),     intent(in) :: file_id
           Character (len=64), intent(in) :: dsetname
           INTEGER(HSIZE_T),   intent(in) :: dims(:)
           logical,            intent(in) :: is_complex
           INTEGER(HSIZE_T),   intent(in), optional :: chunklen
           
           INTEGER                       :: rank, hdferr
           INTEGER(HSIZE_T), allocatable :: dimsc(:), maxdims(:)
           INTEGER(HID_T)                :: dset_id, dataspace, crp_list
           
           !CALL h5open_f(hdferr)
           
           !Define size of dataset and of chunks
           rank = size(dims)
           allocate( dimsc(rank), maxdims(rank) )
           dimsc         = dims
           dimsc(rank)   = 1
           if ( present(chunklen) ) dimsc(rank) = chunklen
           maxdims       = dims
           maxdims(rank) = H5S_UNLIMITED_F
           
           !Check for dims(rank) = 0
           if (dims(rank) /= 0) then
             write(error_unit,*) 'Error in init_dset: dims(rank) /= 0'
             Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif
           
           !Create Dataspace
           CALL h5screate_simple_f(rank, dims, dataspace, hdferr, maxdims)
           
           !Modify dataset creation properties, i.e. enable chunking
           CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
           CALL h5pset_chunk_f(crp_list, rank, dimsc, hdferr)
#ifdef HDF5_ZLIB
           ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
           CALL h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
#endif
           
           !Create a dataset using cparms creation properties.
           CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace, &
                           dset_id, hdferr, crp_list )
           
           CALL write_attribute_logical(dset_id, '.', 'is_complex', is_complex, hdferr)
           
           !Close objects
           CALL h5sclose_f(dataspace, hdferr)
           CALL h5pclose_f(crp_list,  hdferr)
           CALL h5dclose_f(dset_id,   hdferr)
           deallocate( dimsc, maxdims )
           
         end Subroutine init_dset

!--------------------------------------------------------------------
         
         Subroutine append_dat(file_id, dsetname, dat_ptr, Nbins_in)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This subroutine appends one bin to an existing dataset of an opened HDF5 file.
!
!> @param [IN] file_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of the opened HDF5 file
!> \endverbatim
!> @param [IN] dsetname Character(len=64)
!> \verbatim
!>  Name of the dataset
!> \endverbatim
!> @param [IN] data_ptr TYPE(C_PTR)
!> \verbatim
!>  C-pointer to the first element of the data to write.
!>  data should be all double precision.
!>  The length of the data is assumed from the existing dataset.
!> \endverbatim
!> @param [IN] Nbins_in Integer, optional
!> \verbatim
!>  Number of bins to be written, default = 1
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           
           INTEGER(HID_T),     intent(in) :: file_id
           Character (len=64), intent(in) :: dsetname
           TYPE(C_PTR),        intent(in) :: dat_ptr
           INTEGER, optional,  intent(in) :: Nbins_in
           
           INTEGER                       :: rank, hdferr, i, Nbins
           INTEGER(HSIZE_T)              :: mem_dims(1)
           INTEGER(HSIZE_T), allocatable :: dims(:), maxdims(:), offset(:), count(:)
           INTEGER(HID_T)                :: dset_id, dataspace, memspace
           
           !CALL h5open_f(hdferr)
           if( present(Nbins_in) ) then
             Nbins = Nbins_in
           else
             Nbins = 1
           endif
           
           !Open the  dataset.
           CALL h5dopen_f(file_id, dsetname, dset_id, hdferr)
           
           !Get dataset's dataspace handle.
           CALL h5dget_space_f(dset_id, dataspace, hdferr)
           
           !Get dataspace's rank.
           CALL h5sget_simple_extent_ndims_f(dataspace, rank, hdferr)
           allocate( dims(rank), maxdims(rank), offset(rank), count(rank) )
           
           !Get dataspace's dimensions.
           CALL h5sget_simple_extent_dims_f(dataspace, dims, maxdims, hdferr)
           
           !Extent dataset and define hyperslab to write on
           offset(:)    = 0
           offset(rank) = dims(rank)
           count(:)     = dims(:)
           count(rank)  = Nbins
           dims(rank)   = dims(rank)+Nbins
           CALL h5dset_extent_f(dset_id, dims, hdferr)
           CALL h5sclose_f(dataspace, hdferr)
           CALL h5dget_space_f(dset_id, dataspace, hdferr)
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, hdferr)
           
           !Define memory space of data
           mem_dims = Nbins
           do i=1, rank-1
             mem_dims = mem_dims*dims(i)
           enddo
           CALL h5screate_simple_f (1, mem_dims, memspace, hdferr)
           
           !Write data
           CALL H5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr, &
                           mem_space_id = memspace, file_space_id = dataspace)
           
           !Close objects
           CALL h5sclose_f(memspace, hdferr)
           CALL h5sclose_f(dataspace, hdferr)
           CALL h5dclose_f(dset_id,   hdferr)
           
           deallocate( dims, maxdims, offset, count )
           
         end Subroutine append_dat

!--------------------------------------------------------------------
         
         Subroutine write_latt(obj_id, Latt, Latt_Unit)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This subroutine writes the lattice in an opened HDF5 object.
!
!> @param [IN] obj_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of the opened HDF5 object
!> \endverbatim
!> @param [IN] Latt Type(Lattice)
!> \verbatim
!>  The Bravais lattice
!> \endverbatim
!> @param [IN] Latt_unit Type(Unit_cell)
!> \verbatim
!>  The unit cell
!> \endverbatim
!-------------------------------------------------------------------
            Implicit none
            INTEGER(HID_T),   intent(in) :: obj_id
            Type (Lattice),   intent(in) :: Latt
            Type (Unit_cell), intent(in) :: Latt_unit
                  
            Character (len=64) :: group_name, dset_name, attr_name
            INTEGER(HID_T)  :: group_id
            LOGICAL :: link_exists
            INTEGER            :: ierr, no
            INTEGER(HSIZE_T)   :: size_dat
            Real (Kind=Kind(0.d0)), allocatable :: temp(:)
            
            group_name = "lattice"
            CALL h5lexists_f(obj_id, group_name, link_exists, ierr)
            if ( link_exists ) return
            call h5gcreate_f(obj_id, group_name, group_id, ierr)
            
            size_dat = size(Latt%L1_p)
            dset_name = "."
            attr_name = "a1"
            call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%a1_p, size_dat, ierr )
            attr_name = "a2"
            call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%a2_p, size_dat, ierr )
            attr_name = "L1"
            call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%L1_p, size_dat, ierr )
            attr_name = "L2"
            call h5LTset_attribute_double_f(group_id, dset_name, attr_name, Latt%L2_p, size_dat, ierr )
            
            attr_name = "N_coord"
            call write_attribute(group_id, '.', attr_name, Latt_unit%N_coord, ierr)
            attr_name = "Norb"
            call write_attribute(group_id, '.', attr_name, Latt_unit%Norb, ierr)
            attr_name = "Ndim"
            call write_attribute(group_id, '.', attr_name, size(Latt_unit%Orb_pos_p, 2), ierr)
            
            size_dat = size(Latt_unit%Orb_pos_p, 2)
            allocate(temp(size_dat))
            do no = 1, Latt_unit%Norb
               temp(:) = Latt_unit%Orb_pos_p(no,:)
               write(attr_name, '("Orbital", I0)') no
               call h5LTset_attribute_double_f(group_id, dset_name, attr_name, temp, size_dat, ierr )
            enddo
            deallocate(temp)
            
            call h5gclose_f(group_id, ierr)
           
         end Subroutine write_latt

!--------------------------------------------------------------------

!--------------------------------------------------------------------

         Subroutine write_attribute_double(loc_id, obj_name, attr_name, attr_value, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Write a double as attribute to an HDF5 object.
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object to be written to in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [IN] attr_value double
!> \verbatim
!>  Value of attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: loc_id
           CHARACTER(LEN=*), INTENT(IN) :: obj_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           real(Kind=Kind(0.d0)), INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           INTEGER(HID_T) :: space_id, attr_id
           INTEGER(HSIZE_T), parameter :: dims(1) = 1
           
           CALL h5screate_f (H5S_SCALAR_F, space_id, ierr)
           call h5acreate_by_name_f(loc_id, obj_name, attr_name, H5T_NATIVE_DOUBLE, &
                                    space_id, attr_id, ierr)
           call h5awrite_f  (attr_id, H5T_NATIVE_DOUBLE, attr_value, dims, ierr)
           call h5aclose_f  (attr_id, ierr)
           call h5sclose_f  (space_id, ierr)
         end Subroutine write_attribute_double
        
!--------------------------------------------------------------------

         Subroutine write_attribute_int(loc_id, obj_name, attr_name, attr_value, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Write an integer as attribute to an HDF5 object.
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object to be written to in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [IN] attr_value integer
!> \verbatim
!>  Value of attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: loc_id
           CHARACTER(LEN=*), INTENT(IN) :: obj_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           INTEGER,          INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           INTEGER(HID_T) :: space_id, attr_id
           INTEGER(HSIZE_T), parameter :: dims(1) = 1
           
           CALL h5screate_f (H5S_SCALAR_F, space_id, ierr)
           call h5acreate_by_name_f(loc_id, obj_name, attr_name, H5T_NATIVE_INTEGER, &
                                    space_id, attr_id, ierr)
           call h5awrite_f  (attr_id, H5T_NATIVE_INTEGER, attr_value, dims, ierr)
           call h5aclose_f  (attr_id, ierr)
           call h5sclose_f  (space_id, ierr)
         end Subroutine write_attribute_int
        
!--------------------------------------------------------------------

         Subroutine write_attribute_string(loc_id, obj_name, attr_name, attr_value, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Write a string as attribute to an HDF5 object.
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object to be written to in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [IN] attr_value CHARACTER(LEN=*)
!> \verbatim
!>  Value of attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: loc_id
           CHARACTER(LEN=*), INTENT(IN) :: obj_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr

           call h5ltset_attribute_string_f(loc_id, obj_name, attr_name, &
                                           attr_value, ierr)
         end Subroutine write_attribute_string
        
!--------------------------------------------------------------------

         Subroutine write_attribute_logical(loc_id, obj_name, attr_name, attr_value, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Write a boolean as attribute to an HDF5 object (stored as integer).
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object to be written to in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [IN] attr_value logical
!> \verbatim
!>  Value of attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: loc_id
           CHARACTER(LEN=*), INTENT(IN) :: obj_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           LOGICAL,          INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           INTEGER        :: attr_value2
           INTEGER(HID_T) :: space_id, attr_id
           INTEGER(HSIZE_T), parameter :: dims(1) = 1
           
           attr_value2 = 0
           if ( attr_value ) attr_value2 = 1
           
           CALL h5screate_f (H5S_SCALAR_F, space_id, ierr)
           call h5acreate_by_name_f(loc_id, obj_name, attr_name, H5T_NATIVE_INTEGER, &
                                    space_id, attr_id, ierr)
           call h5awrite_f  (attr_id, H5T_NATIVE_INTEGER, attr_value2, dims, ierr)
           call h5aclose_f  (attr_id, ierr)
           call h5sclose_f  (space_id, ierr)
         end Subroutine write_attribute_logical
        
!--------------------------------------------------------------------

         Subroutine read_attribute_double(loc_id, obj_name, attr_name, attr_value, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Read a double from attribute of HDF5 object.
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object to be read from in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [OUT] attr_value double
!> \verbatim
!>  Value of attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),        INTENT(IN)  :: loc_id
           CHARACTER(LEN=*),      INTENT(IN)  :: obj_name
           CHARACTER(LEN=*),      INTENT(IN)  :: attr_name
           real(Kind=Kind(0.d0)), INTENT(OUT) :: attr_value
           INTEGER,               INTENT(OUT) :: ierr
           
           INTEGER(HID_T)              :: attr_id
           INTEGER(HSIZE_T), parameter :: dims(1) = 1
           
           call h5aopen_by_name_f(loc_id, obj_name, attr_name, attr_id, ierr)
           call h5aread_f  (attr_id, H5T_NATIVE_DOUBLE, attr_value, dims, ierr)
           call h5aclose_f (attr_id, ierr)
         end Subroutine read_attribute_double
        
!--------------------------------------------------------------------

         Subroutine read_attribute_int(loc_id, obj_name, attr_name, attr_value, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Read an integer from attribute of HDF5 object.
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object to be read from in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [OUT] attr_value integer
!> \verbatim
!>  Value of attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),   INTENT(IN)  :: loc_id
           CHARACTER(LEN=*), INTENT(IN)  :: obj_name
           CHARACTER(LEN=*), INTENT(IN)  :: attr_name
           INTEGER,          INTENT(OUT) :: attr_value
           INTEGER,          INTENT(OUT) :: ierr
           
           INTEGER(HID_T)              :: attr_id
           INTEGER(HSIZE_T), parameter :: dims(1) = 1
           
           call h5aopen_by_name_f(loc_id, obj_name, attr_name, attr_id, ierr)
           call h5aread_f  (attr_id, H5T_NATIVE_INTEGER, attr_value, dims, ierr)
           call h5aclose_f (attr_id, ierr)
         end Subroutine read_attribute_int
        
!--------------------------------------------------------------------

         Subroutine read_attribute_string(loc_id, obj_name, attr_name, attr_value, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Read a string from attribute of HDF5 object.
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object to be read from in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [OUT] attr_value CHARACTER(LEN=*)
!> \verbatim
!>  Value of attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),   INTENT(IN)  :: loc_id
           CHARACTER(LEN=*), INTENT(IN)  :: obj_name
           CHARACTER(LEN=*), INTENT(IN)  :: attr_name
           CHARACTER(LEN=*), INTENT(OUT) :: attr_value
           INTEGER,          INTENT(OUT) :: ierr

           call h5ltget_attribute_string_f(loc_id, obj_name, attr_name, &
                                           attr_value, ierr)
         end Subroutine read_attribute_string
        
!--------------------------------------------------------------------

         Subroutine read_attribute_logical(loc_id, obj_name, attr_name, attr_value, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Read a boolean from attribute of HDF5 object (stored as integer).
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object to be read from in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [OUT] attr_value logical
!> \verbatim
!>  Value of attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),   INTENT(IN)  :: loc_id
           CHARACTER(LEN=*), INTENT(IN)  :: obj_name
           CHARACTER(LEN=*), INTENT(IN)  :: attr_name
           LOGICAL,          INTENT(OUT) :: attr_value
           INTEGER,          INTENT(OUT) :: ierr
           
           INTEGER                     :: attr_value2
           INTEGER(HID_T)              :: attr_id
           INTEGER(HSIZE_T), parameter :: dims(1) = 1
           
           call h5aopen_by_name_f(loc_id, obj_name, attr_name, attr_id, ierr)
           call h5aread_f  (attr_id, H5T_NATIVE_INTEGER, attr_value2, dims, ierr)
           call h5aclose_f (attr_id, ierr)
           
           if ( attr_value2 == 0 ) then
             attr_value = .false.
           elseif ( attr_value2 == 1 ) then
             attr_value = .true.
           else
             write(error_unit,*) "Error in read_attribute_logical: attr_value2 is neither 0 or 1, but", attr_value2
             Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif
         end Subroutine read_attribute_logical
        
!--------------------------------------------------------------------

         Subroutine test_attribute_double(loc_id, obj_name, attr_name, attr_value, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Test whether supplied double is identical to attribute stored in HDF5 file.
!> If not, triggers error stop.
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object attribute is attached to in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [IN] attr_value double
!> \verbatim
!>  Value of supplied attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: loc_id
           CHARACTER(LEN=*), INTENT(IN) :: obj_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           real(Kind=Kind(0.d0)), INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           LOGICAL :: attr_exists
           real(Kind=Kind(0.d0)), parameter :: ZERO = 10D-8
           real(Kind=Kind(0.d0)) :: test_double, diff
           
           call h5aexists_by_name_f(loc_id, obj_name, attr_name, attr_exists, ierr)
           
           if ( .not. attr_exists ) then
             call write_attribute_double(loc_id, obj_name, attr_name, attr_value, ierr)
           else
             call read_attribute_double(loc_id, obj_name, attr_name, test_double, ierr)
             diff = abs(attr_value - test_double)
             if (diff > ZERO) then
               write(error_unit,*) 'Error in test_attribute_double:', attr_name, ' = ', attr_value, '/=', test_double
               Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
             endif
           endif
         end Subroutine test_attribute_double
        
!--------------------------------------------------------------------

         Subroutine test_attribute_int(loc_id, obj_name, attr_name, attr_value, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Test whether supplied integer is identical to attribute stored in HDF5 file.
!> If not, triggers error stop.
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object attribute is attached to in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [IN] attr_value integer
!> \verbatim
!>  Value of supplied attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: loc_id
           CHARACTER(LEN=*), INTENT(IN) :: obj_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           INTEGER,          INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           LOGICAL :: attr_exists
           INTEGER :: test_int
           
           call h5aexists_by_name_f(loc_id, obj_name, attr_name, attr_exists, ierr)
           
           if ( .not. attr_exists ) then
             call write_attribute_int(loc_id, obj_name, attr_name, attr_value, ierr)
           else
             call read_attribute_int(loc_id, obj_name, attr_name, test_int, ierr)
             if (attr_value /= test_int) then
               write(error_unit,*) 'Error in test_attribute_int:', attr_name, ' = ', attr_value, '/=', test_int
               Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
             endif
           endif
         end Subroutine test_attribute_int
        
!--------------------------------------------------------------------

         Subroutine test_attribute_string(loc_id, obj_name, attr_name, attr_value, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Test whether supplied string is identical to attribute stored in HDF5 file.
!> If not, triggers error stop.
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object attribute is attached to in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [IN] attr_value CHARACTER(LEN=*)
!> \verbatim
!>  Value of supplied attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: loc_id
           CHARACTER(LEN=*), INTENT(IN) :: obj_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           LOGICAL :: attr_exists
           CHARACTER(LEN=64) :: test_string
           
           call h5aexists_by_name_f(loc_id, obj_name, attr_name, attr_exists, ierr)
           
           if ( .not. attr_exists ) then
             call write_attribute_string(loc_id, obj_name, attr_name, attr_value, ierr)
           else
             call read_attribute_string(loc_id, obj_name, attr_name, test_string, ierr)
             if (trim(attr_value) /= trim(test_string)) then
               write(error_unit,*) 'Error in test_attribute_string:', attr_name, ' = ', attr_value, '/=', test_string
               Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
             endif
           endif
         end Subroutine test_attribute_string
        
!--------------------------------------------------------------------

         Subroutine test_attribute_logical(loc_id, obj_name, attr_name, attr_value, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Test whether supplied boolean is identical to attribute stored in HDF5 file.
!> If not, triggers error stop.
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object attribute is attached to in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [IN] attr_value logical
!> \verbatim
!>  Value of supplied attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: loc_id
           CHARACTER(LEN=*), INTENT(IN) :: obj_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           LOGICAL,          INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           LOGICAL :: attr_exists
           LOGICAL :: test_bool
           
           call h5aexists_by_name_f(loc_id, obj_name, attr_name, attr_exists, ierr)
           
           if ( .not. attr_exists ) then
             call write_attribute_logical(loc_id, obj_name, attr_name, attr_value, ierr)
           else
             call read_attribute_logical(loc_id, obj_name, attr_name, test_bool, ierr)
             if (attr_value .neqv. test_bool) then
               write(error_unit,*) 'Error in test_attribute_logical:', attr_name, ' = ', attr_value, '/=', test_bool
               Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
             endif
           endif
         end Subroutine test_attribute_logical



         Subroutine write_comment(loc_id, obj_name, attr_name, comment, ierr)
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Write a comment (array of strings, each 64 characters in length) 
!> as attribute to an HDF5 object.
!
!> @param [IN] loc_id INTEGER(HID_T)
!> \verbatim
!>  Idendifier of opened HDF5 object
!> \endverbatim
!> @param [IN] obj_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of object to be written to in relation to loc_id
!> \endverbatim
!> @param [IN] attr_name CHARACTER(LEN=*)
!> \verbatim
!>  Name of attribute
!> \endverbatim
!> @param [IN] comment(:) CHARACTER(LEN=64)
!> \verbatim
!>  Value of attribute
!> \endverbatim
!> @param [OUT] ierr integer
!> \verbatim
!>  Error code
!> \endverbatim
!-------------------------------------------------------------------
           
           IMPLICIT NONE
           INTEGER(HID_T),    INTENT(IN) :: loc_id
           CHARACTER(LEN=*),  INTENT(IN) :: obj_name
           CHARACTER(LEN=*),  INTENT(IN) :: attr_name
           CHARACTER(len=64), INTENT(IN) :: comment(:)
           INTEGER, INTENT(OUT) ::   ierr
           
           INTEGER(HID_T)   :: attr_id       ! Attribute identifier
           INTEGER(HID_T)   :: space_id      ! Attribute Dataspace identifier
           INTEGER(HID_T)   :: type_id       ! Attribute datatype identifier
           INTEGER          :: rank = 1      ! Attribure rank
           INTEGER(HSIZE_T) :: dims(1)       ! Attribute dimensions
           INTEGER(SIZE_T)  :: attrlen = 64  ! Length of the attribute string
           
           ! Create scalar data space for the attribute.
           dims(1) = size(comment)
           CALL h5screate_simple_f(rank, dims, space_id, ierr)
           
           ! Create datatype for the attribute.
           CALL h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, ierr)
           CALL h5tset_size_f(type_id, attrlen, ierr)
           
           ! Create dataset attribute.
           call h5acreate_by_name_f(loc_id, obj_name, attr_name, type_id, &
                                    space_id, attr_id, ierr)
          
           ! Write the attribute data.
           CALL h5awrite_f(attr_id, type_id, comment, dims, ierr)
           
           ! Close the attribute, datatype and data space.
           CALL h5aclose_f(attr_id, ierr)
           CALL h5tclose_f(type_id, ierr)
           CALL h5sclose_f(space_id, ierr)
           
        end Subroutine write_comment
         
     end Module alf_hdf5
#endif
