program test_io
  Use mpi
  Implicit none

  Integer       :: size_conf, i, ierr, Isize, Irank, SEED_IN, STATUS(MPI_STATUS_SIZE), k
  Real(Kind=Kind(0.d0)), allocatable :: conf(:)
  Real(Kind=Kind(0.d0)) :: x
  CHARACTER (LEN=64) :: file_conf
  Integer, allocatable :: SEED_VEC(:)
  
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
  
  CALL RANDOM_SEED (SIZE=K)
  ALLOCATE (SEED_VEC(K))
  
  IF (IRANK == 0) THEN
    OPEN(UNIT=5, FILE='seeds', STATUS='OLD', ACTION='READ')
    DO I = ISIZE-1,1,-1
      READ (5,*) SEED_IN
      CALL MPI_SEND(SEED_VEC, K, MPI_INTEGER, I, I+1024, MPI_COMM_WORLD,IERR)
    ENDDO
    READ(5,*) SEED_VEC
    CLOSE(5)
  ELSE
    CALL MPI_RECV(SEED_IN, K, MPI_INTEGER,0,  IRANK + 1024,  MPI_COMM_WORLD,STATUS,IERR)
  ENDIF
  CALL RANSET(SEED_VEC)
  CALL RANDOM_SEED (PUT=SEED_VEC)
  DEALLOCATE (SEED_VEC)
  
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


