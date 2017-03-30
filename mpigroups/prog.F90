!***********************************************************************************
! Minimal program to test and learn the use of split communicators, groups
! Compilation: gfortran prog.F90
! optional configurations: -DMPI -DPTE
! where PTE = parallel tempering
! 
! Issues;
! 1. the number of groups is hard-coded 
! 2. the compiler options are either: -DMPI ,or: -DMPI -DPTE
!    There are no safe guards for the use of only -DPTE, so only -DPTE is probably not safe
!***********************************************************************************

program prog
#ifdef MPI
    include 'mpif.h'
#endif 
    real(kind=kind(0.d0)) :: p,j
    character(len=32)     :: filename
    integer               :: i,seed_in
    namelist /var_prog/ p
#ifdef MPI
    integer               :: world_size, world_rank,ierr, errcode
    real(kind=kind(0.d0)) :: x
#endif
#ifdef MPI
#ifdef PTE
    integer               :: part_size, part_rank,part_color,part_comm, number_of_groups
#endif
#endif


#ifdef MPI
    call mpi_init(ierr)
    if (ierr .ne. mpi_success) then
       write(*,*) 'error starting mpi program: terminating.'
       call mpi_abort(mpi_comm_world, errcode, ierr)
       stop
    end if
    call mpi_comm_size(mpi_comm_world,world_size,ierr)
    call mpi_comm_rank(mpi_comm_world,world_rank,ierr)
#endif

!***********************************************************************************
! set up the new groups. The group identifier is part_color
!***********************************************************************************
#ifdef MPI
#ifdef PTE
     number_of_groups = 4
     part_color = mod(world_rank,number_of_groups)
     call mpi_comm_split(mpi_comm_world,part_color, world_rank, part_comm,ierr)
     call mpi_comm_size(part_comm,part_size,ierr)
     call mpi_comm_rank(part_comm,part_rank,ierr)
#endif
#endif

!***********************************************************************************
! check if:
! a) plain MPI: parameter file exists
! b) MPI and split communicator parallel tempering flag)directories exist 
! and read parameters
!***********************************************************************************
#ifdef MPI
#ifndef PTE
    if ( world_rank .eq. 0 ) then 
       open(unit=5,file='parameters',status='old',action='read',iostat=ierr)
       if (ierr /= 0) then
           write(*,*) 'unable to open <parameters>',ierr
           call mpi_abort(mpi_comm_world, errcode, ierr)
           stop
       end if
       read(5,nml=var_prog)
       close(5)
    endif
#endif

#ifdef PTE
     if (part_rank .eq. 0) then
        write(filename,'(A,I1,A)') 'dir_',part_color,'/parameters'
        open(unit=5,file=filename,status='old',action='read',iostat=ierr)
        if (ierr /= 0) then
            write(*,*) 'unable to open: ',filename,ierr
            call mpi_abort(part_comm, errcode, ierr)
            call mpi_abort(mpi_comm_world, errcode, ierr)
            stop
        endif
        write(*,*) 'Read parameters from: ',filename
        read(5,nml=var_prog)
        close(5)
    endif
#endif

#endif
#ifndef MPI
open(unit=5,file='parameters',status='old',action='read',iostat=ierr)
    if (ierr /= 0) then
        write(*,*) 'unable to open: <parameters>',ierr
        call mpi_abort(mpi_comm_world, errcode, ierr)
        stop
    end if
    read(5,nml=var_prog)
    close(5)
#endif

!***********************************************************************************
!broadcast parameters
!***********************************************************************************
#ifdef MPI
#ifndef PTE
    call mpi_bcast(p ,1,mpi_real8,  0,mpi_comm_world,ierr)
#endif
#ifdef PTE
    call mpi_bcast(p ,1,mpi_real8,  0,part_comm,ierr)
#endif
#endif

#ifdef MPI
#ifdef PTE
    write(*,'(A,I4,I4,I4,I4,2x,F4.2)') 'rank, tasks, part_rank,part_size, p',world_rank, world_size, part_rank, part_size,p
!     write(*,*) 'Parameter is', p
#endif
#ifndef PTE
    write(*,'(A,I4,I4,2x,F4.2)') 'rank, tasks, p', world_rank, world_size,p
!     write(*,*) 'Parameter is', p
#endif
#endif


!***********************************************************************************
! Do some work, and in case of MPI, give each thread a slightly different task, depending on its rank 
!***********************************************************************************
    j = 0.d0
#ifdef MPI
#ifndef PTE
    j = 0.1d0*world_rank    
#endif
#ifdef PTE
    j = 0.1d0*part_rank
#endif
#endif
    do i = 1, 10
       j = j + p
    enddo
!***********************************************************************************
! Now I want to swap between groups
!***********************************************************************************
  
#ifdef MPI 
#ifdef PTE
! ideas where to look:
! mpi_comm_create, mpi_group_incl, mpi_comm_group, maybe mpi_car_create to set a topology, maybe MPI_Neighbor_alltoallw.

#endif
#endif
  
!***********************************************************************************
! collect the observable j
!***********************************************************************************
#ifdef MPI
#ifndef PTE
    x = 0.d0
    call mpi_reduce(j,x,1,mpi_real8,mpi_sum, 0,mpi_comm_world,ierr)
    j= x/dble(world_size)
#endif
#ifdef PTE
    x = 0.d0
    call mpi_reduce(j,x,1,mpi_real8,mpi_sum, 0,part_comm,ierr)
    j= x/dble(part_size)
#endif
#endif

!***********************************************************************************
! output the observable j
!***********************************************************************************
#ifndef PTE
#ifdef MPI
    if (world_rank .eq. 0) then
#endif
       open(unit=5,file='output',position='append',action='write')
       write(5,'(F9.4)') j
       close(5)
#ifdef MPI
    endif
#endif
#endif
       
#ifdef MPI
#ifdef PTE
    if (part_rank .eq. 0) then
        write(filename,'(A,I1,A)') 'dir_',part_color,'/output'
        open(unit=5,file=filename,position='append',action='write')
        write(5,'(F9.4)') j
        close(5)
    endif
#endif
#endif


!***********************************************************************************
!finalize the MPI
!***********************************************************************************
#ifdef MPI
#ifdef PTE
    call mpi_comm_free(part_comm);
#endif
    call mpi_finalize(err_mpi)
#endif

end program prog




