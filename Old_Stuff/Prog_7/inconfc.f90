        SUBROUTINE confin 

         Use Hamiltonian

         Implicit none

#include "machine"
         

         
#ifdef MPI
         INCLUDE 'mpif.h'
         ! Local
#endif   

         Integer        :: I, IERR, ISIZE, IRANK, seed_in, K, iseed, Nt
         Integer, dimension(:), allocatable :: Seed_vec
         Real (Kind=8)  :: X
         Logical ::   lconf 
         character (len=64) :: file_sr, File_tg

#ifdef MPI
         INTEGER        :: STATUS(MPI_STATUS_SIZE)
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif

         Allocate (Nsigma(Size(Op_V,1),Ltrot))
         
#ifdef MPI
         INQUIRE (FILE='confin_0', EXIST=lconf)
         If (lconf) Then
            file_sr = "confin"
            Call Get_seed_Len(K)
            Allocate(Seed_vec(K))
            file_tg = File_i(file_sr,IRANK)
            Open (Unit = 10, File=File_tg, status='old', ACTION='read')
            Read(10,*) Seed_vec
            Call Ranset(Seed_vec)
            do NT = 1,LTROT
               do I = 1,Size(Op_V,1)
                  Read(10,*) NSIGMA(I,NT) 
               enddo
            enddo
            close(10)
            Deallocate(Seed_vec)
         else
            If (Irank == 0) then
               Write(6,*) 'No initial configuration'
               OPEN(UNIT=5,FILE='seeds',STATUS='old',ACTION='read',IOSTAT=ierr)
               IF (ierr /= 0) THEN
                  WRITE(*,*) 'unable to open <seeds>',ierr
                  STOP
               END IF
               DO I = Isize-1,1,-1
                  Read (5,*) Seed_in
                  CALL MPI_SEND(Seed_in,1,MPI_INTEGER, I, I+1024,MPI_COMM_WORLD,IERR)
               enddo
               Read(5,*) Seed_in
               CLOSE(5)
            else
               CALL MPI_RECV(Seed_in, 1, MPI_INTEGER,0,  IRANK + 1024,  MPI_COMM_WORLD,STATUS,IERR)
            endif
            Call Get_seed_Len(K)
            !Write(6,*) K
            Allocate(Seed_vec(K))
            Do I = 1,K
               X =  Ranf_Imada(Seed_in)
               Seed_vec(I) = Seed_in
            enddo
            Call Ranset(Seed_vec)
            Deallocate(Seed_vec)
            do NT = 1,LTROT
               do I = 1,Size(Op_V,1)
                  X = RANF()
                  NSIGMA(I,NT) = 1
                  IF (X.GT.0.5) NSIGMA(I,NT) = -1
               enddo
            enddo
         endif
            
#else   
         INQUIRE (FILE='confin_0', EXIST=lconf)
         If (lconf) Then
            file_tg = "confin_0"
            Call Get_seed_Len(K)
            Allocate(Seed_vec(K))
            Open (Unit = 10, File=File_tg, status='old', ACTION='read')
            Read(10,*) Seed_vec
            Call Ranset(Seed_vec)
            do NT = 1,LTROT
               do I = 1,Size(Op_V,1)
                  Read(10,*) NSIGMA(I,NT) 
               enddo
            enddo
            close(10)
            Deallocate(Seed_vec)
         else
            Write(6,*) 'No initial configuration'
            OPEN(UNIT=5,FILE='seeds',STATUS='old',ACTION='read',IOSTAT=ierr)
            IF (ierr /= 0) THEN
               WRITE(*,*) 'unable to open <seeds>',ierr
               STOP
            END IF
            Read (5,*) Seed_in
            CLOSE(5)
            Call Get_seed_Len(K)
            !Write(6,*) K
            Allocate(Seed_vec(K))
            Do I = 1,K
               X =  Ranf_Imada(Seed_in)
               Seed_vec(I) = Seed_in
            enddo
            Call Ranset(Seed_vec)
            Deallocate(Seed_vec)
            do NT = 1,LTROT
               do I = 1,Size(Op_V,1)
                  X = RANF()
                  NSIGMA(I,NT) = 1
                  IF (X.GT.0.5) NSIGMA(I,NT) = -1
               enddo
            enddo
         endif
#endif
         
       END SUBROUTINE CONFIN
