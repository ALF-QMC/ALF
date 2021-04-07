program test_io
  Use mpi
  Implicit none

  Integer       :: size_conf, i, ierr, Irank
  Real(Kind=Kind(0.d0)), allocatable :: conf(:)
  Real(Kind=Kind(0.d0)) :: x
  CHARACTER (LEN=64) :: file_conf
  
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
  
  size_conf = 64000
  allocate(conf(size_conf))
  
  do i=1, size_conf
    Call Random_Number(x)
    if (x > 0.5d0) then
      conf(i) = 1.d0
    else
      conf(i) = -1.d0
    endif
  enddo
  
  write(file_conf,'(A,I0)') "confout_", IRANK
  OPEN(UNIT=10, FILE=file_conf, STATUS='UNKNOWN', ACTION='WRITE')
  do i=1, size_conf
    write(10,*) nint(conf(i))
  enddo
  close(10)
  
  deallocate(conf)
  
  CALL MPI_FINALIZE(ierr)
        
end program