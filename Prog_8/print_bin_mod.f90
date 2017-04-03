     Module Print_bin_mod

       Interface Print_bin
          module procedure Print_bin_C, Print_bin_R
       end Interface Print_bin


       Contains
         
         Subroutine  Print_bin_C(Dat_eq,Dat_eq0,Latt, Nobs, Phase_bin_tmp, file_pr)
           Use Lattices_v3
           Implicit none
#ifdef MPI
           include 'mpif.h'
#endif   

           Complex (Kind=Kind(0.d0)), Dimension(:,:,:), Intent(inout):: Dat_eq
           Complex (Kind=Kind(0.d0)), Dimension(:)    , Intent(inout):: Dat_eq0
           Type (Lattice),                     Intent(In)   :: Latt
           Complex (Kind=Kind(0.d0)),                   Intent(In)   :: Phase_bin_tmp
           Character (len=64),                 Intent(In)   :: File_pr
           Integer,                            Intent(In)   :: Nobs
          
           ! Local
           Character (len=64) :: File_tmp
           Integer :: Norb, I, no,no1
           Complex (Kind=Kind(0.d0)), allocatable :: Tmp(:,:,:), Tmp1(:)
           Real    (Kind=Kind(0.d0))              :: x_p(2) 
           Complex (Kind=Kind(0.d0))              :: Phase_bin
#ifdef MPI
           Complex (Kind=Kind(0.d0)):: Z
           Integer         :: Ierr, Isize, Irank
           INTEGER         :: STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif

           Phase_bin = Phase_bin_tmp
           Norb = size(Dat_eq,3)
           if ( .not. (Latt%N  == Size(Dat_eq,1) ) ) then 
              Write(6,*) 'Error in Print_bin' 
              Stop
           endif
           Allocate (Tmp(Latt%N,Norb,Norb), Tmp1(Norb) )
           Dat_eq = Dat_eq/dble(Nobs)
           Dat_eq0 = Dat_eq0/dble(Nobs*Latt%N)
           
#if defined(MPI) && !defined(TEMPERING)
           I = Latt%N*Norb*Norb
           Tmp = cmplx(0.d0, 0.d0, kind(0.D0))
           CALL MPI_REDUCE(Dat_eq,Tmp,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Dat_eq = Tmp/DBLE(ISIZE)
           I = 1
           Z = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Phase_bin,Z,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Phase_bin= Z/DBLE(ISIZE)

           I = Norb
           Tmp1 = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Dat_eq0,Tmp1,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Dat_eq0 = Tmp1/DBLE(ISIZE)

           If (Irank == 0 ) then
#endif
#if defined(TEMPERING)
           write(File_tmp,'(A,I0,A,A)') "Temp_",Irank,"/",trim(ADJUSTL(File_pr))
#else
           File_tmp=File_pr
#endif
              do no = 1,Norb
                 do no1 = 1,Norb
                    Call  Fourier_R_to_K(Dat_eq(:,no,no1), Tmp(:,no,no1), Latt)
                 enddo
              enddo
              Open (Unit=10,File=File_tmp, status="unknown",  position="append")
              Write(10,*) dble(Phase_bin),Norb,Latt%N
              do no = 1,Norb
                 Write(10,*) Dat_eq0(no)
              enddo
              do I = 1,Latt%N
                 x_p = dble(Latt%listk(i,1))*Latt%b1_p + dble(Latt%listk(i,2))*Latt%b2_p  
                 Write(10,*) X_p(1), X_p(2)
                 do no = 1,Norb
                    do no1 = 1,Norb
                       Write(10,*) tmp(I,no,no1)
                    enddo
                 enddo
              enddo
              close(10)
#if defined(MPI) && !defined(TEMPERING)
           Endif
#endif
              
           deallocate (Tmp, tmp1 )
          

         End Subroutine Print_bin_C
          

!=========================================================

         Subroutine  Print_bin_R(Dat_eq,Dat_eq0,Latt, Nobs, Phase_bin_tmp, file_pr)
           Use Lattices_v3
           Implicit none
#ifdef MPI
           include 'mpif.h'
#endif   
           
           Real    (Kind=Kind(0.d0)), Dimension(:,:,:), Intent(inout) :: Dat_eq
           Real    (Kind=Kind(0.d0)), Dimension(:)    , Intent(inout) :: Dat_eq0
           Type (Lattice),                     Intent(In)    :: Latt
           Complex (Kind=Kind(0.d0)),                   Intent(In)    :: Phase_bin_tmp
           Character (len=64),                 Intent(In)    :: File_pr
           Integer,                            Intent(In)    :: Nobs
           
           ! Local
           Character (len=64) :: File_tmp
           Integer :: Norb, I, no,no1
           Real    (Kind=Kind(0.d0)), allocatable :: Tmp(:,:,:), Tmp1(:)
           Real    (Kind=Kind(0.d0))              :: x_p(2) 
           Complex (Kind=Kind(0.d0))              :: Phase_bin
#ifdef MPI
           Integer        :: Ierr, Isize, Irank
           Complex (Kind=Kind(0.d0))              :: Z
           INTEGER        :: STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
           
           Phase_bin = Phase_bin_tmp
           Norb = size(Dat_eq,3)
           if ( .not. (Latt%N  == Size(Dat_eq,1) ) ) then 
              Write(6,*) 'Error in Print_bin' 
              Stop
           endif
           Allocate (Tmp(Latt%N,Norb,Norb), Tmp1(Norb) )
           Dat_eq  = Dat_eq/dble(Nobs)
           Dat_eq0 = Dat_eq0/(dble(Nobs)*dble(Latt%N))
#if defined(MPI) && !defined(TEMPERING)
           I = Latt%N*Norb*Norb
           Tmp = 0.d0
           CALL MPI_REDUCE(Dat_eq,Tmp,I,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Dat_eq = Tmp/DBLE(ISIZE)
           I = 1
           Z = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Phase_bin,Z,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Phase_bin= Z/DBLE(ISIZE)

           I = Norb 
           Tmp1 = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Dat_eq0,Tmp1,I,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Dat_eq0 = Tmp1/DBLE(ISIZE)

           If (Irank == 0 ) then
#endif
#if defined(TEMPERING)
           write(File_tmp,'(A,I0,A,A)') "Temp_",Irank,"/",trim(ADJUSTL(File_pr))
#else
           File_tmp=File_pr
#endif
           
              do no = 1,Norb
                 do no1 = 1,Norb
                    Call  Fourier_R_to_K(Dat_eq(:,no,no1), Tmp(:,no,no1), Latt)
                 enddo
              enddo
              Open (Unit=10,File=File_tmp, status="unknown",  position="append")
              Write(10,*) dble(Phase_bin),Norb,Latt%N
              do no = 1,Norb
                 Write(10,*) Dat_eq0(no)
              enddo
              do I = 1,Latt%N
                 x_p = dble(Latt%listk(i,1))*Latt%b1_p + dble(Latt%listk(i,2))*Latt%b2_p  
                 Write(10,*) X_p(1), X_p(2)
                 do no = 1,Norb
                    do no1 = 1,Norb
                       Write(10,*) tmp(I,no,no1)
                    enddo
                 enddo
              enddo
              close(10)
#if defined(MPI) && !defined(TEMPERING)
           endif
#endif
           deallocate (Tmp )
           
         End Subroutine Print_bin_R
!============================================================
         Subroutine  Print_scal(Obs, Nobs, file_pr)
           
           Implicit none
#ifdef MPI
           include 'mpif.h'
#endif   
           
           Complex   (Kind=Kind(0.d0)), Dimension(:), Intent(inout) :: Obs
           Character (len=64),               Intent(In)    :: File_pr
           Integer,                          Intent(In)    :: Nobs
           
           ! Local
           Character (len=64) :: File_tmp
           Integer :: Norb,I
           Complex  (Kind=Kind(0.d0)), allocatable :: Tmp(:)
#ifdef MPI
           Integer        :: Ierr, Isize, Irank
           INTEGER        :: STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
           
           Norb = size(Obs,1)
           Allocate ( Tmp(Norb) )
           Obs = Obs/dble(Nobs)
#if defined(MPI) && !defined(TEMPERING)
           Tmp = 0.d0
           CALL MPI_REDUCE(Obs,Tmp,Norb,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Obs = Tmp/DBLE(ISIZE)
           if (Irank == 0 ) then
#endif
#if defined(TEMPERING)
           write(File_tmp,'(A,I0,A,A)') "Temp_",Irank,"/",trim(ADJUSTL(File_pr))
#else
           File_tmp=File_pr
#endif
              Open (Unit=10,File=File_tmp, status="unknown",  position="append")
              WRITE(10,*) (Obs(I), I=1,size(Obs,1))
              close(10)
#if defined(MPI) && !defined(TEMPERING)
           endif
#endif
           deallocate (Tmp )
           
         End Subroutine Print_scal

!==============================================================
         Subroutine  Print_bin_tau(Dat_tau, Latt, Nobs, Phase_bin, file_pr, dtau, Dat0_tau)
           Use Lattices_v3
#ifdef ZLIB
           USE F95ZLIB
           USE IOPORTS
           USE ISO_C_BINDING
#endif  
           Implicit none
#ifdef MPI
           include 'mpif.h'
#endif   

           Complex (Kind=Kind(0.d0)), Dimension(:,:,:,:), Intent(inout):: Dat_tau   ! (Latt%N, Ltau,Norb, Norb)
           Complex (Kind=Kind(0.d0)), Dimension(:      ), Intent(inout), optional :: Dat0_tau  ! (Norb)
           Type (Lattice),                       Intent(In)   :: Latt
           Complex (Kind=Kind(0.d0)),                     Intent(In)   :: Phase_bin
           Character (len=64),                   Intent(In)   :: File_pr
           Integer,                              Intent(In)   :: Nobs
           Real (Kind=Kind(0.d0)),                        Intent(In)   :: dtau
          
           ! Local
           Character (len=64) :: File_tmp
           Integer :: Norb, I, no,no1, LT, nt,ios
           Complex (Kind=Kind(0.d0)), allocatable :: Tmp(:,:,:,:), Tmp0(:)
           Complex (Kind=Kind(0.d0)) :: Phase_mean 
           Real    (Kind=Kind(0.d0))              :: x_p(2) 
#ifdef ZLIB
           TYPE(IOPORT) :: fd
           CHARACTER(LEN=255), TARGET :: LINE
#endif  
#ifdef MPI
           Complex (Kind=Kind(0.d0)):: Z
           Integer         :: Ierr, Isize, Irank
           INTEGER         :: STATUS(MPI_STATUS_SIZE)
           CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif

           Phase_mean = Phase_bin
           Norb = size(Dat_tau,3)
           if ( .not. (Latt%N  == Size(Dat_tau,1) ) ) then 
              Write(6,*) 'Error in Print_bin' 
              Stop
           endif
           LT = Size(Dat_tau,2)
           Allocate (Tmp (Latt%N,LT,Norb,Norb) )
           Allocate (Tmp0(Norb) )

           Tmp0 = cmplx(0.d0, 0.d0, kind(0.D0))
           Dat_tau  = Dat_tau/dble(Nobs)
           If (Present(Dat0_tau) ) Dat0_tau = Dat0_tau/dble(Nobs*Latt%N*LT)
           
#if defined(MPI) && !defined(TEMPERING)
           I = Latt%N*Norb*Norb*LT
           Tmp = cmplx(0.d0, 0.d0, kind(0.D0))
           CALL MPI_REDUCE(Dat_tau,Tmp,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Dat_tau = Tmp/DBLE(ISIZE)
           
           If (Present(Dat0_tau) ) then
              I = Norb
              Tmp0 = cmplx(0.d0, 0.d0, kind(0.D0))
              CALL MPI_REDUCE(Dat0_tau,Tmp0,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
              Dat0_tau = Tmp0/DBLE(ISIZE)
           endif

           I = 1
           Z = cmplx(0.d0,0.d0,kind(0.d0))
           CALL MPI_REDUCE(Phase_mean,Z,I,MPI_COMPLEX16,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
           Phase_mean= Z/DBLE(ISIZE)
           If (Irank == 0 ) then
#endif
#if defined(TEMPERING)
           write(File_tmp,'(A,I0,A,A)') "Temp_",Irank,"/",trim(ADJUSTL(File_pr))
#else
           File_tmp=File_pr
#endif
              If (Present(Dat0_tau) ) Tmp0 = Dat0_tau
!$OMP parallel do default(shared) private(nt,no,no1)
              do nt = 1,LT
                 do no = 1,Norb
                    do no1 = 1,Norb
                       Call  Fourier_R_to_K(Dat_tau(:,nt,no,no1), Tmp(:,nt,no,no1), Latt)
                    enddo
                 enddo
              enddo
!$OMP end parallel do
#ifdef ZLIB
              write(File_tmp,*) TRIM(ADJUSTL(File_pr)),".gz"
              CALL FGZ_OPEN(TRIM(ADJUSTL(File_tmp)),'a6',fd,ios)
              Write(Line,*) dble(Phase_mean),Norb,Latt%N, LT, dtau
              CALL FGZ_WRITE(fd,TRIM(LINE),'yes',IOS)
#else
              Open (Unit=10,File=File_tmp, status="unknown",  position="append")
              Write(10,*) dble(Phase_mean),Norb,Latt%N, LT, dtau
#endif
              Do no = 1,Norb
#ifdef ZLIB
                 Write(LINE,*)  Tmp0(no)
                 CALL FGZ_WRITE(fd,TRIM(LINE),'yes',IOS)
#else
                 Write(10,*)  Tmp0(no)
#endif
              enddo
              do I = 1,Latt%N
                 x_p = dble(Latt%listk(i,1))*Latt%b1_p + dble(Latt%listk(i,2))*Latt%b2_p 
#ifdef ZLIB 
                 Write(Line,*) X_p(1), X_p(2)
                 CALL FGZ_WRITE(fd,TRIM(LINE),'yes',IOS)
#else
                 Write(10,*) X_p(1), X_p(2)
#endif
                 Do nt = 1,LT
                    do no = 1,Norb
                       do no1 = 1,Norb
#ifdef ZLIB 
                          Write(Line,*) tmp(I,nt,no,no1)
                          CALL FGZ_WRITE(fd,TRIM(LINE),'yes',IOS)
#else
                          Write(10,*) tmp(I,nt,no,no1)
#endif
                       enddo
                    enddo
                 enddo
              enddo
#ifdef ZLIB 
              CALL FGZ_CLOSE(fd,IOS)
#else
              close(10)
#endif
#if defined(MPI) && !defined(TEMPERING)
           Endif
#endif
              
           deallocate (Tmp, Tmp0 )

         End Subroutine Print_bin_tau


         
       end Module Print_bin_mod
