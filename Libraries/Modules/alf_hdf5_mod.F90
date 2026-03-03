!  Copyright (C) 2020-2026 The ALF project
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
!> @author ALF-project
!> @brief Helper subroutines for using ALF with HDF5.
!
!> @details
!> This module provides a Fortran interface for HDF5 operations within the
!> ALF (Algorithms for Lattice Fermions) project. It includes functionality for:
!> - Creating and managing HDF5 datasets for Monte Carlo bin data
!> - Writing and reading attributes of various types (double, integer, string, logical)
!> - Storing lattice structure information in HDF5 format
!> - Validating consistency between supplied values and stored HDF5 attributes
!>
!> The module defines three generic interfaces:
!> - write_attribute: Write attributes to HDF5 objects (supports double, integer, string, logical)
!> - read_attribute: Read attributes from HDF5 objects (supports double, integer, string, logical)
!> - test_attribute: Verify attribute values or write them if they don't exist
!
!> @note
!> This module is only available when ALF is compiled with HDF5 support (HDF5 preprocessor flag).
!> Logical values are stored as integers (0 = false, 1 = true) in HDF5.
!
!> @see
!> For HDF5 library documentation, visit: https://portal.hdfgroup.org/
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
     
       !> Generic interface for writing attributes of various types to HDF5 objects
       interface write_attribute
         MODULE PROCEDURE write_attribute_double, write_attribute_int, write_attribute_string, write_attribute_logical
       end interface write_attribute
       
       !> Generic interface for reading attributes of various types from HDF5 objects
       interface read_attribute
         MODULE PROCEDURE read_attribute_double, read_attribute_int, read_attribute_string, read_attribute_logical
       end interface read_attribute
       
       !> Generic interface for testing/validating attributes or creating them if they don't exist
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
!> Creates a new extensible dataset in an opened HDF5 file for storing Monte Carlo bins.
!
!> @details
!> This subroutine creates an HDF5 dataset with unlimited size in the last dimension,
!> allowing bins to be appended incrementally via append_dat(). The dataset uses
!> chunking to enable efficient extensibility. All data is stored as double precision.
!>
!> The dims array defines the shape of a single bin, where:
!> - rank = size(dims)
!> - dims(1:rank-1) define the bin shape
!> - dims(rank) must be 0 (will be extended as bins are added)
!
!> @param[in] file_id HDF5 file identifier
!> @param[in] dsetname Name of the new dataset (max 64 characters)
!> @param[in] dims Shape of one bin with rank+1 dimensions where dims(rank) = 0
!> @param[in] is_complex If .true., indicates values can be complex (stored as metadata)
!> @param[in] chunklen Size of data chunks in number of bins (default = 1)
!
!> @note
!> The last dimension (dims(rank)) must be initialized to 0. It will be extended
!> automatically when data is appended via append_dat().
!
!> @warning
!> If HDF5_ZLIB is defined at compile time, DEFLATE compression will be applied
!> to the dataset using the specified compression level.
!
!> @see append_dat
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
           
           ! Define size of dataset and of chunks
           ! The last dimension is unlimited (for extensibility) and chunked
           rank = size(dims)
           allocate( dimsc(rank), maxdims(rank) )
           dimsc         = dims
           dimsc(rank)   = 1                      ! Chunk size: 1 bin by default
           if ( present(chunklen) ) dimsc(rank) = chunklen
           maxdims       = dims
           maxdims(rank) = H5S_UNLIMITED_F        ! Allow unlimited growth in last dimension
           
           ! Validate that dims(rank) = 0 (required for extensible datasets)
           if (dims(rank) /= 0) then
             write(error_unit,*) 'Error in init_dset: dims(rank) /= 0'
             Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif
           
           ! Create dataspace with extensible last dimension
           CALL h5screate_simple_f(rank, dims, dataspace, hdferr, maxdims)
           
           ! Enable chunking (required for extensible datasets)
           CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
           CALL h5pset_chunk_f(crp_list, rank, dimsc, hdferr)
#ifdef HDF5_ZLIB
           ! Apply ZLIB/DEFLATE compression if enabled at compile time
           ! HDF5_ZLIB macro defines the compression level (0-9)
           CALL h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
#endif
           
           ! Create the dataset with the specified properties
           CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dataspace, &
                           dset_id, hdferr, crp_list )
           
           ! Store whether data is complex as metadata attribute
           CALL write_attribute_logical(dset_id, '.', 'is_complex', is_complex, hdferr)
           
           ! Close objects
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
!> Appends Monte Carlo bin data to an existing extensible HDF5 dataset.
!
!> @details
!> This subroutine extends an HDF5 dataset created with init_dset() and writes
!> new bin data to it. It uses HDF5 hyperslab selection to append data to the
!> end of the dataset in the last dimension. The data layout is determined from
!> the existing dataset dimensions.
!>
!> Workflow:
!> 1. Open the dataset and get its current dimensions
!> 2. Extend the dataset size in the last dimension by Nbins
!> 3. Select a hyperslab at the end of the extended dataset
!> 4. Write the flattened data array to the selected region
!
!> @param[in] file_id HDF5 file identifier
!> @param[in] dsetname Name of the dataset (max 64 characters)
!> @param[in] dat_ptr C pointer to the first element of data (double precision array)
!> @param[in] Nbins_in Number of bins to append (default = 1)
!
!> @note
!> The data length is inferred from the dataset dimensions. The memory layout
!> is flattened into a 1D array = Nbins * product(dims(1:rank-1)).
!
!> @see init_dset
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
           
           ! Extend dataset in last dimension and select hyperslab for writing
           offset(:)    = 0
           offset(rank) = dims(rank)              ! Start writing after existing data
           count(:)     = dims(:)
           count(rank)  = Nbins                   ! Write Nbins slices
           dims(rank)   = dims(rank)+Nbins        ! New total size after extension
           CALL h5dset_extent_f(dset_id, dims, hdferr)
           CALL h5sclose_f(dataspace, hdferr)
           CALL h5dget_space_f(dset_id, dataspace, hdferr)
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, hdferr)
           
           ! Define memory space: flatten multi-dimensional bin data to 1D array
           ! Total elements = Nbins * product of bin dimensions
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
!> Writes lattice structure information to an HDF5 object.
!
!> @details
!> This subroutine creates a "lattice" group within the specified HDF5 object
!> and stores comprehensive lattice information as attributes:
!> - Bravais lattice vectors: a1, a2, L1, L2
!> - Unit cell properties: N_coord, Norb, Ndim
!> - Orbital positions: Orbital1, Orbital2, ... OrbitalN
!>
!> The group structure created is:
!> /[parent_object]/lattice/
!>   - Attributes: a1, a2, L1, L2, N_coord, Norb, Ndim, Orbital1, ...
!
!> @param[in] obj_id HDF5 object identifier (file or group)
!> @param[in] Latt Bravais lattice structure
!> @param[in] Latt_unit Unit cell structure
!
!> @note
!> If the "lattice" group already exists, the subroutine returns immediately
!> without modifying the file (idempotent behavior).
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
            ! Check if lattice group already exists (avoid duplicate writes)
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
!> Writes a double precision scalar attribute to an HDF5 object.
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of target object relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to create
!> @param[in] attr_value Double precision value to write
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @see read_attribute_double, test_attribute_double
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
!> Writes an integer scalar attribute to an HDF5 object.
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of target object relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to create
!> @param[in] attr_value Integer value to write
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @see read_attribute_int, test_attribute_int
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
!> Writes a string attribute to an HDF5 object.
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of target object relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to create
!> @param[in] attr_value String value to write
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @see read_attribute_string, test_attribute_string
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
!> Writes a logical (boolean) attribute to an HDF5 object.
!
!> @details
!> HDF5 does not have a native boolean type, so logical values are stored
!> as integers: .false. = 0, .true. = 1
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of target object relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to create
!> @param[in] attr_value Logical value to write
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @note
!> Uses integer representation: 0 for .false., 1 for .true.
!
!> @see read_attribute_logical, test_attribute_logical
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
           
           ! Convert logical to integer (HDF5 has no native boolean type)
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
!> Reads a double precision scalar attribute from an HDF5 object.
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of object to read from relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to read
!> @param[out] attr_value Double precision value read from attribute
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @see write_attribute_double, test_attribute_double
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
!> Reads an integer scalar attribute from an HDF5 object.
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of object to read from relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to read
!> @param[out] attr_value Integer value read from attribute
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @see write_attribute_int, test_attribute_int
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
!> Reads a string attribute from an HDF5 object.
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of object to read from relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to read
!> @param[out] attr_value String value read from attribute
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @see write_attribute_string, test_attribute_string
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
!> Reads a logical (boolean) attribute from an HDF5 object.
!
!> @details
!> Reads the integer representation stored by write_attribute_logical and
!> converts it back to a logical value: 0 = .false., 1 = .true.
!> Any other value triggers an error.
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of object to read from relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to read
!> @param[out] attr_value Logical value read from attribute
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @warning
!> The stored integer value must be 0 or 1. Other values will cause program termination.
!
!> @see write_attribute_logical, test_attribute_logical
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
           
           ! Convert integer back to logical (validate range)
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
!> Validates that a double attribute matches the supplied value, or creates it if absent.
!
!> @details
!> This subroutine checks if an attribute exists in the HDF5 file:
!> - If absent: writes the supplied value as a new attribute
!> - If present: reads the stored value and compares with supplied value
!>   - If values differ by more than tolerance (1e-8), terminates with error
!>   - If values match within tolerance, no action taken
!>
!> This is useful for ensuring consistency across multiple runs or restart files.
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of object relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to test
!> @param[in] attr_value Double precision value to test against
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @note
!> Comparison tolerance is 1e-8 (relative to absolute difference).
!
!> @warning
!> Terminates program if attribute exists but does not match supplied value.
!
!> @see write_attribute_double, read_attribute_double
!-------------------------------------------------------------------
           Implicit none
           INTEGER(HID_T),   INTENT(IN) :: loc_id
           CHARACTER(LEN=*), INTENT(IN) :: obj_name
           CHARACTER(LEN=*), INTENT(IN) :: attr_name
           real(Kind=Kind(0.d0)), INTENT(IN) :: attr_value
           INTEGER,   INTENT(OUT) :: ierr
           
           LOGICAL :: attr_exists
           real(Kind=Kind(0.d0)), parameter :: ZERO = 10D-8  ! Comparison tolerance
           real(Kind=Kind(0.d0)) :: test_double, diff
           
           call h5aexists_by_name_f(loc_id, obj_name, attr_name, attr_exists, ierr)
           
           if ( .not. attr_exists ) then
             call write_attribute_double(loc_id, obj_name, attr_name, attr_value, ierr)
           else
             call read_attribute_double(loc_id, obj_name, attr_name, test_double, ierr)
             diff = abs(attr_value - test_double)
             if (diff > ZERO) then
               write(error_unit,'(A)') 'Error in test_attribute_double:'
               write(error_unit,'(A)') '  Attribute name: "' // trim(attr_name) // '"'
               write(error_unit,'(A, F0.6)') '  Supplied value: ' , attr_value
                write(error_unit,'(A, F0.6)') '  Read value:     ' , test_double
               write(error_unit,'(A, F0.6)') '  Difference = ', diff
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
!> Validates that an integer attribute matches the supplied value, or creates it if absent.
!
!> @details
!> This subroutine checks if an attribute exists in the HDF5 file:
!> - If absent: writes the supplied value as a new attribute
!> - If present: reads the stored value and compares with supplied value
!>   - If values differ, terminates with error
!>   - If values match, no action taken
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of object relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to test
!> @param[in] attr_value Integer value to test against
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @warning
!> Terminates program if attribute exists but does not match supplied value.
!
!> @see write_attribute_int, read_attribute_int
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
               write(error_unit,'(A)') 'Error in test_attribute_int:'
               write(error_unit,'(A)') '  Attribute name: "' // trim(attr_name) // '"'
               write(error_unit,'(A, I0)') '  Supplied value: ' , attr_value
               write(error_unit,'(A, I0)') '  Read value:     ' , test_int
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
!> Validates that a string attribute matches the supplied value, or creates it if absent.
!
!> @details
!> This subroutine checks if an attribute exists in the HDF5 file:
!> - If absent: writes the supplied value as a new attribute
!> - If present: reads the stored value and compares with supplied value (trimmed)
!>   - If values differ, terminates with error
!>   - If values match, no action taken
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of object relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to test
!> @param[in] attr_value String value to test against
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @note
!> String comparison uses trim() to ignore trailing spaces.
!
!> @warning
!> Terminates program if attribute exists but does not match supplied value.
!
!> @see write_attribute_string, read_attribute_string
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
               write(error_unit,*) 'Error in test_attribute_string:'
               write(error_unit,*) '  Attribute name: "' // trim(attr_name) // '"'
               write(error_unit,*) '  Supplied value: "' // trim(attr_value) // '"'
               write(error_unit,*) '  Read value:     "' // trim(test_string) // '"'
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
!> Validates that a logical attribute matches the supplied value, or creates it if absent.
!
!> @details
!> This subroutine checks if an attribute exists in the HDF5 file:
!> - If absent: writes the supplied value as a new attribute
!> - If present: reads the stored value and compares with supplied value
!>   - If values differ, terminates with error
!>   - If values match, no action taken
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of object relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to test
!> @param[in] attr_value Logical value to test against
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @warning
!> Terminates program if attribute exists but does not match supplied value.
!
!> @see write_attribute_logical, read_attribute_logical
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
               write(error_unit,'(A)') 'Error in test_attribute_logical:'
               write(error_unit,'(A, A, A)') '  Attribute name: "' // trim(attr_name) // '"'
               write(error_unit,'(A,L1)') '  Supplied value: ', attr_value
               write(error_unit,'(A,L1)') '  Read value:     ', test_bool
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
!> Writes multi-line comment as an array of fixed-length strings to an HDF5 object.
!
!> @details
!> This subroutine stores documentation or metadata as an HDF5 attribute.
!> Each string in the array is stored with fixed length (64 characters).
!> This is useful for embedding parameter descriptions, citations, or
!> analysis notes directly in HDF5 output files.
!
!> @param[in] loc_id HDF5 object identifier (file, group, or dataset)
!> @param[in] obj_name Name of target object relative to loc_id (use "." for loc_id itself)
!> @param[in] attr_name Name of the attribute to create
!> @param[in] comment Array of strings, each exactly 64 characters long
!> @param[out] ierr HDF5 error code (0 on success)
!
!> @note
!> String length is fixed at 64 characters. Shorter strings are padded,
!> longer strings must be split across multiple array elements.
!
!> Example:
!> @code
!> character(len=64) :: comments(3)
!> comments(1) = "Simulation started at 2026-03-03"
!> comments(2) = "Parameters: U=4.0, beta=10.0"
!> comments(3) = "Lattice: 8x8 square"
!> call write_comment(file_id, ".", "simulation_info", comments, ierr)
!> @endcode
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
           INTEGER          :: rank = 1      ! Attribute rank
           INTEGER(HSIZE_T) :: dims(1)       ! Attribute dimensions
           INTEGER(SIZE_T)  :: attrlen = 64  ! Fixed string length (all strings 64 chars)
           
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
