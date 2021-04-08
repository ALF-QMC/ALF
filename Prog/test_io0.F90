program test_io
  Implicit none

  Integer       :: size_conf, i
  Real(Kind=Kind(0.d0)), allocatable :: conf(:)
  
  size_conf = 6400
  allocate(conf(size_conf))
  conf =  1.d0
  OPEN(UNIT=10, FILE="confout", STATUS='UNKNOWN', ACTION='WRITE')
  do i=1, size_conf
     write(10,*) nint(conf(i))
  enddo
  close(10)
  
  deallocate(conf)
  
        
end program
