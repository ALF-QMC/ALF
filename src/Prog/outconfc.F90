       SUBROUTINE confout

         Use Hamiltonian

         Implicit none

#include "machine"
         

#ifdef MPI
         INCLUDE 'mpif.h'
         ! Local
#endif
         
         Integer        :: I, IERR, ISIZE, IRANK, seed_in, K, iseed, Nt, nr
         Integer, dimension(:), allocatable :: Seed_vec
         Real (Kind=8)  :: X
         Logical        :: lconf 
         character (len=64) :: file_sr, File_tg

#ifdef MPI
         INTEGER        :: STATUS(MPI_STATUS_SIZE)
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
         
         Call Get_seed_Len(K)
         Allocate(Seed_vec(K))
         Call Ranget(Seed_vec)
         file_sr = "confout"
         file_tg = File_i(file_sr,IRANK)
         Open (Unit = 10, File=File_tg, status='unknown', ACTION='write')
         Write(10,*) Seed_vec
         do NT = 1,LTROT
            do I = 1,Size(Nsigma,1)
               write(10,*) NSIGMA(I,NT) 
            enddo
         enddo
         close(10)
         Deallocate(Seed_vec)
            
#else   
         Call Get_seed_Len(K)
         Allocate(Seed_vec(K))
         Call Ranget(Seed_vec)
         file_tg = "confout_0"
         Open (Unit = 10, File=File_tg, status='unknown', ACTION='write')
         Write(10,*) Seed_vec
         do NT = 1,LTROT
            do I = 1,Size(Nsigma,1)
               write(10,*) Nsigma(I,NT) 
            enddo
         enddo
         close(10)
         Deallocate(Seed_vec)
#endif
         
       END SUBROUTINE CONFOUT
