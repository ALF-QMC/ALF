!  Copyright (C) 2021 The ALF project
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
!       https://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage https://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.


!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module provides an interface to allocate memory that is shared
!> between different MPI jobs from a single communicator on the same node
!
!--------------------------------------------------------------------


module mpi_shared_memory

    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER, C_SIZEOF
    use runtime_error_mod
#ifdef MPI
    use mpi
#endif
    use iso_fortran_env, only: output_unit, error_unit
    Implicit none
    private
    public :: mpi_shared_memory_init, allocate_shared_memory, deallocate_all_shared_memory, use_mpi_shm
  
    ! internal arrays that contain the shared memory chunk that can be distributed
    complex (Kind=Kind(0.d0)), POINTER, save :: shm_mem_chunk_cmplx(:)
    real    (Kind=Kind(0.d0)), POINTER, save :: shm_mem_chunk_real(:) 
    ! storing the MPI window ids that are use to release the memory at the end of the run / synchronization barriers
    integer, save, allocatable, dimension(:) :: mpi_wins_real, mpi_wins_cmplx
    ! internal members to manage mpi communication / memory distribution
    integer, save :: nodecomm, noderank, head_idx_cmplx, head_idx_real
    Real    (Kind=Kind(0.d0)), save :: chunk_size_gb
    integer, save :: num_chunks_real, num_chunks_cmplx
    integer(kind=8), save :: chunk_size_real(1), chunk_size_cmplx(1)
    logical, save :: initialized=.false.
    logical, save :: use_mpi_shm=.false. !> public variable to query if shared memory module is active
    
!--------------------------------------------------------------------
    !> @brief 
    !> interface to memory allocation routines.
    !> fortran array may be real or complex double; 1D to 4D;
    !> mpi_win_loc can be used for memory synchronization barrier;
    !> myrank to ensure that only one rank initializes the array
    !
    !> @param[out] fortran_array
    !> @param[out] mpi_win_loc
    !> @param[out] myrank
    !> @param[in]  arrayshape
    !--------------------------------------------------------------------
    INTERFACE allocate_shared_memory
      MODULE PROCEDURE allocate_shared_memory_1Dreal,  allocate_shared_memory_2Dreal, &
                    &  allocate_shared_memory_3Dreal,  allocate_shared_memory_4Dreal, &
                    &  allocate_shared_memory_1Dcmplx, allocate_shared_memory_2Dcmplx, &
                    &  allocate_shared_memory_3Dcmplx, allocate_shared_memory_4Dcmplx
    END INTERFACE
    
    Contains

      !--------------------------------------------------------------------
      !> @brief 
      !> initializes memory manager; can be called without MPI present.
      !> Does nothing if compiled without MPI
      !
      !> @param[in] mpi_communicator
      !> @param[in] chunk_size (in GB)
      !--------------------------------------------------------------------
      subroutine mpi_shared_memory_init(mpi_communicator, chunk_size)
        Implicit none
        integer, intent(in) :: mpi_communicator
        real(Kind=Kind(0.d0)), intent(in) :: chunk_size
        
#ifdef MPI
        integer :: ierr, tmp_int, status
        real(Kind=Kind(0.d0)) :: dummy_real_dp

        chunk_size_gb=chunk_size
        if (chunk_size_gb > 0) then
                use_mpi_shm=.true.
                CALL MPI_Comm_split_type(mpi_communicator, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, nodecomm,ierr)
                CALL MPI_Comm_rank(nodecomm, noderank,ierr)
                initialized=.true.
                num_chunks_real=0
                num_chunks_cmplx=0
                allocate(mpi_wins_real(10),mpi_wins_cmplx(10))
                if (noderank==0) write(*,*) "Chunk size for mpi shared memory is ", chunk_size_gb,"GB"
                if (MPI_ADDRESS_KIND < 8 .and. chunk_size_gb > 2) then
                        write(*,*) "Reducing chunk size to 2GB inorder to avoid integer overflow in MPI library."
                        chunk_size_gb = 2
                endif
        endif
#endif
  
      end subroutine mpi_shared_memory_init

      !--------------------------------------------------------------------
      !> @brief 
      !> internal helper routine that allocates chunks of shared MPI memory that can be carved out and distributed
      !
      !--------------------------------------------------------------------
      subroutine allocate_shm_chunk_real()
        Implicit none
        
#ifdef MPI
        INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
        INTEGER :: disp_unit, ierr, mpi_win_loc, tmp_int
        integer, allocatable :: tmp_int_array(:)
        real(Kind=Kind(0.d0)) :: dummy_real_dp
        TYPE(C_PTR) :: baseptr

        ! allocate GB(s) of memory as a 1D array of reals
        chunk_size_real(1) = nint(chunk_size_gb*1024.d0*1024.d0*1024.d0/dble(C_SIZEOF(dummy_real_dp)),8)
        
        if (.not. initialized) then
            WRITE(error_unit,*) 'Please initialize the mpi_shared_memory module before allocating the first array'
            Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        if (noderank == 0) then
            windowsize = int(chunk_size_real(1),MPI_ADDRESS_KIND)*int(C_SIZEOF(dummy_real_dp),MPI_ADDRESS_KIND) 
        else
            windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
        CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL, nodecomm, baseptr, mpi_win_loc, ierr)
        !count the number of alloacted chunks
        num_chunks_real=num_chunks_real+1
        !store mpi_win_loc, this is enough to deallocate all of them 
        tmp_int=size(mpi_wins_real)
        if (num_chunks_real >= tmp_int) then
            allocate(tmp_int_array(tmp_int))
            tmp_int_array=mpi_wins_real
            deallocate(mpi_wins_real)
            allocate(mpi_wins_real(tmp_int+10))
            mpi_wins_real(1:tmp_int)=tmp_int_array
            deallocate(tmp_int_array)
        endif
        mpi_wins_real(num_chunks_real)=mpi_win_loc

        ! Obtain the location of the memory segment
        if (noderank /= 0) then
            CALL MPI_Win_shared_query(mpi_win_loc, 0, windowsize, disp_unit, baseptr, ierr)
        end if

        ! baseptr can now be associated with a Fortran pointer
        ! and thus used to access the shared data
        CALL C_F_POINTER(baseptr, shm_mem_chunk_real, chunk_size_real)

        ! It's a fresh memory chunk, distribution will start at the beginning
        head_idx_real=1
#endif
      end subroutine allocate_shm_chunk_real
      
      !--------------------------------------------------------------------
      !> @brief 
      !> specific implementation of above interface for 1D real arrays
      !
      !> @param[out] fortran_array
      !> @param[out] mpi_win_loc
      !> @param[out] myrank
      !> @param[in]  arrayshape
      !--------------------------------------------------------------------
      subroutine allocate_shared_memory_1Dreal(fortran_array, mpi_win_loc, myrank, arrayshape)
        Implicit none
        real (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:)
        integer, intent(out) :: mpi_win_loc
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(1)
        
#ifdef MPI
        ! allocate a new chunk if the remaining space is not enough
        ! this routine takes no efford to retain accessability to the leftover, unused chunk
        ! but all mpi_windows are kept such that we can still properly deallocate them
        if (head_idx_real+arrayshape(1)>chunk_size_real(1)) then
            call allocate_shm_chunk_real
        endif

        ! pass on the shm rank (only one rank should write to the shared memory at the same time)
        myrank=noderank
        ! pass on the mpi_win (needed for syncronization) [possibly hide in an explicit sync routine]
        mpi_win_loc=mpi_wins_real(num_chunks_real)
        ! carve out a piece of the shm chunk and hand it out
        fortran_array (1:arrayshape(1)) => shm_mem_chunk_real(head_idx_real:head_idx_real+arrayshape(1)-1)
        ! mark it as distributed, i.e., shift the head accordingly
        head_idx_real=head_idx_real+arrayshape(1)
#else
        mpi_win_loc=-1
        myrank=-1
        !allocate plain old fortran array if run without MPI
!        allocate(fortran_array(arrayshape(1)))
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
  
      end subroutine allocate_shared_memory_1Dreal
      
      !--------------------------------------------------------------------
      !> @brief 
      !> specific implementation of above interface for 2D real arrays
      !
      !> @param[out] fortran_array
      !> @param[out] mpi_win_loc
      !> @param[out] myrank
      !> @param[in]  arrayshape
      !--------------------------------------------------------------------
      subroutine allocate_shared_memory_2Dreal(fortran_array, mpi_win_loc, myrank, arrayshape)
        Implicit none
        real (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:,:)
        integer, intent(out) :: mpi_win_loc
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(2)
        
#ifdef MPI
        ! allocate a new chunk if the remaining space is not enough
        ! this routine takes no efford to retain accessability to the leftover, unused chunk
        ! but all mpi_windows are kept such that we can still properly deallocate them
        if (head_idx_real+arrayshape(1)*arrayshape(2)>chunk_size_real(1)) then
            call allocate_shm_chunk_real
        endif

        ! pass on the shm rank (only one rank should write to the shared memory at the same time)
        myrank=noderank
        ! pass on the mpi_win (needed for syncronization) [possibly hide in an explicit sync routine]
        mpi_win_loc=mpi_wins_real(num_chunks_real)
        ! carve out a piece of the shm chunk and hand it out
        fortran_array (1:arrayshape(1),1:arrayshape(2)) => shm_mem_chunk_real(head_idx_real:head_idx_real+PRODUCT(arrayshape)-1)
        ! mark it as distributed, i.e., shift the head accordingly
        head_idx_real=head_idx_real+arrayshape(1)*arrayshape(2)
#else
        mpi_win_loc=-1
        myrank=-1
        !allocate plain old fortran array if run without MPI
!        allocate( fortran_array(arrayshape(1),1:arrayshape(2)) )
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
  
      end subroutine allocate_shared_memory_2Dreal
      
      !--------------------------------------------------------------------
      !> @brief 
      !> specific implementation of above interface for 3D real arrays
      !
      !> @param[out] fortran_array
      !> @param[out] mpi_win_loc
      !> @param[out] myrank
      !> @param[in]  arrayshape
      !--------------------------------------------------------------------
      subroutine allocate_shared_memory_3Dreal(fortran_array, mpi_win_loc, myrank, arrayshape)
        Implicit none
        real (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:,:,:)
        integer, intent(out) :: mpi_win_loc
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(3)
        
#ifdef MPI
        ! allocate a new chunk if the remaining space is not enough
        ! this routine takes no efford to retain accessability to the leftover, unused chunk
        ! but all mpi_windows are kept such that we can still properly deallocate them
        if (head_idx_real+PRODUCT(arrayshape)>chunk_size_real(1)) then
            call allocate_shm_chunk_real
        endif

        ! pass on the shm rank (only one rank should write to the shared memory at the same time)
        myrank=noderank
        ! pass on the mpi_win (needed for syncronization) [possibly hide in an explicit sync routine]
        mpi_win_loc=mpi_wins_real(num_chunks_real)
        ! carve out a piece of the shm chunk and hand it out
        fortran_array (1:arrayshape(1),1:arrayshape(2),1:arrayshape(3)) => &
            & shm_mem_chunk_real(head_idx_real:head_idx_real+PRODUCT(arrayshape)-1)
        ! mark it as distributed, i.e., shift the head accordingly
        head_idx_real=head_idx_real+PRODUCT(arrayshape)
#else
        mpi_win_loc=-1
        myrank=-1
        !allocate plain old fortran array if run without MPI
!        allocate( fortran_array(arrayshape(1),1:arrayshape(2),1:arrayshape(3)) )
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
  
      end subroutine allocate_shared_memory_3Dreal
      
      !--------------------------------------------------------------------
      !> @brief 
      !> specific implementation of above interface for 4D real arrays
      !
      !> @param[out] fortran_array
      !> @param[out] mpi_win_loc
      !> @param[out] myrank
      !> @param[in]  arrayshape
      !--------------------------------------------------------------------
      subroutine allocate_shared_memory_4Dreal(fortran_array, mpi_win_loc, myrank, arrayshape)
        Implicit none
        real (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:,:,:,:)
        integer, intent(out) :: mpi_win_loc
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(4)
     
#ifdef MPI
        ! allocate a new chunk if the remaining space is not enough
        ! this routine takes no efford to retain accessability to the leftover, unused chunk
        ! but all mpi_windows are kept such that we can still properly deallocate them
        if (head_idx_real+PRODUCT(arrayshape)>chunk_size_real(1)) then
            call allocate_shm_chunk_real
        endif

        ! pass on the shm rank (only one rank should write to the shared memory at the same time)
        myrank=noderank
        ! pass on the mpi_win (needed for syncronization) [possibly hide in an explicit sync routine]
        mpi_win_loc=mpi_wins_real(num_chunks_real)
        ! carve out a piece of the shm chunk and hand it out
        fortran_array (1:arrayshape(1),1:arrayshape(2),1:arrayshape(3),1:arrayshape(4)) => &
            & shm_mem_chunk_real(head_idx_real:head_idx_real+PRODUCT(arrayshape)-1)
        ! mark it as distributed, i.e., shift the head accordingly
        head_idx_real=head_idx_real+PRODUCT(arrayshape)
#else
        mpi_win_loc=-1
        myrank=-1
        !allocate plain old fortran array if run without MPI
!        allocate( fortran_array(arrayshape(1),1:arrayshape(2),1:arrayshape(3),1:arrayshape(4)) )
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
  
      end subroutine allocate_shared_memory_4Dreal



      !--------------------------------------------------------------------
      !> @brief 
      !> internal helper routine that allocates chunks of shared MPI complex memory that can be carved out and distributed
      !
      !--------------------------------------------------------------------
      subroutine allocate_shm_chunk_cmplx()
        Implicit none
        
#ifdef MPI
        INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
        INTEGER :: disp_unit, ierr, mpi_win_loc, tmp_int
        integer, allocatable :: tmp_int_array(:)
        complex(Kind=Kind(0.d0)) :: dummy_cmplx_dp
        TYPE(C_PTR) :: baseptr

        ! allocate GB(s) of memory as a 1D array of complex doubles
        chunk_size_cmplx(1) = nint(chunk_size_gb*1024.d0*1024.d0*1024.d0/dble(C_SIZEOF(dummy_cmplx_dp)),8)
        
        if (.not. initialized) then
            WRITE(error_unit,*) 'Please initialize the mpi_shared_memory module before allocating the first array'
            Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        if (noderank == 0) then
            windowsize = int(chunk_size_cmplx(1),MPI_ADDRESS_KIND)*int(C_SIZEOF(dummy_cmplx_dp),MPI_ADDRESS_KIND)
        else
            windowsize = 0_MPI_ADDRESS_KIND
        end if
        disp_unit = 1
        CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL, nodecomm, baseptr, mpi_win_loc, ierr)
        !count the number of alloacted chunks
        num_chunks_cmplx=num_chunks_cmplx+1
        !store mpi_win_loc, this is enough to deallocate all of them 
        tmp_int=size(mpi_wins_cmplx)
        if (num_chunks_cmplx >= tmp_int) then
            allocate(tmp_int_array(tmp_int))
            tmp_int_array=mpi_wins_cmplx
            deallocate(mpi_wins_cmplx)
            allocate(mpi_wins_cmplx(tmp_int+10))
            mpi_wins_cmplx(1:tmp_int)=tmp_int_array
            deallocate(tmp_int_array)
        endif
        mpi_wins_cmplx(num_chunks_cmplx)=mpi_win_loc

        ! Obtain the location of the memory segment
        if (noderank /= 0) then
            CALL MPI_Win_shared_query(mpi_win_loc, 0, windowsize, disp_unit, baseptr, ierr)
        end if

        ! baseptr can now be associated with a Fortran pointer
        ! and thus used to access the shared data
        CALL C_F_POINTER(baseptr, shm_mem_chunk_cmplx, chunk_size_cmplx)

        ! It's a fresh memory chunk, distribution will start at the beginning
        head_idx_cmplx=1
#endif
      end subroutine allocate_shm_chunk_cmplx
      
      
      !--------------------------------------------------------------------
      !> @brief 
      !> specific implementation of above interface for 1D complex arrays
      !
      !> @param[out] fortran_array
      !> @param[out] mpi_win_loc
      !> @param[out] myrank
      !> @param[in]  arrayshape
      !--------------------------------------------------------------------
      subroutine allocate_shared_memory_1Dcmplx(fortran_array, mpi_win_loc, myrank, arrayshape)
        Implicit none
        complex (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:)
        integer, intent(out) :: mpi_win_loc
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(1)
        
#ifdef MPI
        ! allocate a new chunk if the remaining space is not enough
        ! this routine takes no efford to retain accessability to the leftover, unused chunk
        ! but all mpi_windows are kept such that we can still properly deallocate them
        if (head_idx_cmplx+arrayshape(1)>chunk_size_cmplx(1)) then
            call allocate_shm_chunk_cmplx
        endif

        ! pass on the shm rank (only one rank should write to the shared memory at the same time)
        myrank=noderank
        ! pass on the mpi_win (needed for syncronization) [possibly hide in an explicit sync routine]
        mpi_win_loc=mpi_wins_cmplx(num_chunks_cmplx)
        ! carve out a piece of the shm chunk and hand it out
        fortran_array (1:arrayshape(1)) => shm_mem_chunk_cmplx(head_idx_cmplx:head_idx_cmplx+arrayshape(1)-1)
        ! mark it as distributed, i.e., shift the head accordingly
        head_idx_cmplx=head_idx_cmplx+arrayshape(1)
#else
        mpi_win_loc=-1
        myrank=-1
        !allocate plain old fortran array if run without MPI
!        allocate(fortran_array(arrayshape(1)))
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
  
      end subroutine allocate_shared_memory_1Dcmplx
      
      !--------------------------------------------------------------------
      !> @brief 
      !> specific implementation of above interface for 2D complex arrays
      !
      !> @param[out] fortran_array
      !> @param[out] mpi_win_loc
      !> @param[out] myrank
      !> @param[in]  arrayshape
      !--------------------------------------------------------------------
      subroutine allocate_shared_memory_2Dcmplx(fortran_array, mpi_win_loc, myrank, arrayshape)
        Implicit none
        complex (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:,:)
        integer, intent(out) :: mpi_win_loc
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(2)
                
#ifdef MPI
        ! allocate a new chunk if the remaining space is not enough
        ! this routine takes no efford to retain accessability to the leftover, unused chunk
        ! but all mpi_windows are kept such that we can still properly deallocate them
        if (head_idx_cmplx+PRODUCT(arrayshape)>chunk_size_cmplx(1)) then
            call allocate_shm_chunk_cmplx
        endif

        ! pass on the shm rank (only one rank should write to the shared memory at the same time)
        myrank=noderank
        ! pass on the mpi_win (needed for syncronization) [possibly hide in an explicit sync routine]
        mpi_win_loc=mpi_wins_cmplx(num_chunks_cmplx)
        ! carve out a piece of the shm chunk and hand it out
        fortran_array (1:arrayshape(1),1:arrayshape(2)) => &
                &shm_mem_chunk_cmplx(head_idx_cmplx:head_idx_cmplx+PRODUCT(arrayshape)-1)
        ! mark it as distributed, i.e., shift the head accordingly
        head_idx_cmplx=head_idx_cmplx+PRODUCT(arrayshape)
#else
        mpi_win_loc=-1
        myrank=-1
        !allocate plain old fortran array if run without MPI
!        allocate(fortran_array(arrayshape(1),arrayshape(2)))
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
  
      end subroutine allocate_shared_memory_2Dcmplx
      
      !--------------------------------------------------------------------
      !> @brief 
      !> specific implementation of above interface for 3D complex arrays
      !
      !> @param[out] fortran_array
      !> @param[out] mpi_win_loc
      !> @param[out] myrank
      !> @param[in]  arrayshape
      !--------------------------------------------------------------------
      subroutine allocate_shared_memory_3Dcmplx(fortran_array, mpi_win_loc, myrank, arrayshape)
        Implicit none
        complex (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:,:,:)
        integer, intent(out) :: mpi_win_loc
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(3)
                
#ifdef MPI
        ! allocate a new chunk if the remaining space is not enough
        ! this routine takes no efford to retain accessability to the leftover, unused chunk
        ! but all mpi_windows are kept such that we can still properly deallocate them
        if (head_idx_cmplx+PRODUCT(arrayshape)>chunk_size_cmplx(1)) then
            call allocate_shm_chunk_cmplx
        endif

        ! pass on the shm rank (only one rank should write to the shared memory at the same time)
        myrank=noderank
        ! pass on the mpi_win (needed for syncronization) [possibly hide in an explicit sync routine]
        mpi_win_loc=mpi_wins_cmplx(num_chunks_cmplx)
        ! carve out a piece of the shm chunk and hand it out
        fortran_array (1:arrayshape(1),1:arrayshape(2),1:arrayshape(3)) => &
                &shm_mem_chunk_cmplx(head_idx_cmplx:head_idx_cmplx+PRODUCT(arrayshape)-1)
        ! mark it as distributed, i.e., shift the head accordingly
        head_idx_cmplx=head_idx_cmplx+PRODUCT(arrayshape)
#else
        mpi_win_loc=-1
        myrank=-1
        !allocate plain old fortran array if run without MPI
!        allocate(fortran_array(arrayshape(1),arrayshape(2),arrayshape(3)))
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
  
      end subroutine allocate_shared_memory_3Dcmplx
      
      !--------------------------------------------------------------------
      !> @brief 
      !> specific implementation of above interface for 4D complex arrays
      !
      !> @param[out] fortran_array
      !> @param[out] mpi_win_loc
      !> @param[out] myrank
      !> @param[in]  arrayshape
      !--------------------------------------------------------------------
      subroutine allocate_shared_memory_4Dcmplx(fortran_array, mpi_win_loc, myrank, arrayshape)
        Implicit none
        complex (Kind=Kind(0.d0)), POINTER, intent(inout) :: fortran_array(:,:,:,:)
        integer, intent(out) :: mpi_win_loc
        integer, intent(out) :: myrank
        integer, intent(in)  :: arrayshape(4)
       
                
#ifdef MPI
        ! allocate a new chunk if the remaining space is not enough
        ! this routine takes no efford to retain accessability to the leftover, unused chunk
        ! but all mpi_windows are kept such that we can still properly deallocate them
        if (head_idx_cmplx+PRODUCT(arrayshape)>chunk_size_cmplx(1)) then
            call allocate_shm_chunk_cmplx
        endif

        ! pass on the shm rank (only one rank should write to the shared memory at the same time)
        myrank=noderank
        ! pass on the mpi_win (needed for syncronization) [possibly hide in an explicit sync routine]
        mpi_win_loc=mpi_wins_cmplx(num_chunks_cmplx)
        ! carve out a piece of the shm chunk and hand it out
        fortran_array (1:arrayshape(1),1:arrayshape(2),1:arrayshape(3),1:arrayshape(4)) => &
                &shm_mem_chunk_cmplx(head_idx_cmplx:head_idx_cmplx+PRODUCT(arrayshape)-1)
        ! mark it as distributed, i.e., shift the head accordingly
        head_idx_cmplx=head_idx_cmplx+PRODUCT(arrayshape)
#else
        mpi_win_loc=-1
        myrank=-1
        !allocate plain old fortran array if run without MPI
!        allocate(fortran_array(arrayshape(1),arrayshape(2),arrayshape(3),arrayshape))
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
  
      end subroutine allocate_shared_memory_4Dcmplx
      
      subroutine deallocate_all_shared_memory()
        Implicit none
        integer :: ierr, i, mpi_win_loc
        
#ifdef MPI
        external :: MPI_Win_free    ! This seems to be required by gfortran10 with OpenMPI on Fedora33 (should be part of MPI module)

        do i=1,num_chunks_real
            mpi_win_loc=mpi_wins_real(i)
            call MPI_WIN_FENCE(0, mpi_win_loc, ierr)
            call MPI_BARRIER(nodecomm,ierr)
            call MPI_Win_free(mpi_win_loc,ierr)
        enddo
        do i=1,num_chunks_cmplx
            mpi_win_loc=mpi_wins_cmplx(i)
            call MPI_WIN_FENCE(0, mpi_win_loc, ierr)
            call MPI_BARRIER(nodecomm,ierr)
            call MPI_Win_free(mpi_win_loc,ierr)
        enddo
#else
        WRITE(error_unit,*) 'This module requires MPI and should not be called without it'
        Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
  
      end subroutine deallocate_all_shared_memory
      
end module
