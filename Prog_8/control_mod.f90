  module Control

    Use MyMats
    Implicit none

    real (Kind=8)   , private, save :: XMEANG, XMAXG, XMAXP, CPU_time_st, CPU_time_en, Xmean_tau, Xmax_tau
    Integer         , private, save :: NCG, NCG_tau
    Integer (Kind=8), private, save :: NC_up, ACC_up
 
    
    Contains

      subroutine control_init
        Implicit none
        XMEANG    = 0.d0
        XMEAN_tau = 0.d0
        XMAXG     = 0.d0
        XMAX_tau  = 0.d0
        NCG       = 0
        NCG_tau   = 0
        NC_up     = 0
        ACC_up    = 0
        Call CPU_TIME(CPU_time_st)
      end subroutine control_init
      
      Subroutine Control_upgrade(Log) 
        Implicit none
        Logical :: Log
        NC_up = NC_up + 1
        if (Log) ACC_up = ACC_up + 1
      end Subroutine Control_upgrade

      Subroutine Control_PrecisionG(A,B,Ndim)
        Implicit none
        
        Integer :: Ndim
        Complex (Kind=8) :: A(Ndim,Ndim), B(Ndim,Ndim) 
        Real    (Kind=8) :: XMAX, XMEAN

        !Local 
        NCG = NCG + 1
        XMEAN = 0.d0
        XMAX  = 0.d0
        CALL COMPARE(A, B, XMAX, XMEAN)
        IF (XMAX  >  XMAXG) XMAXG = XMAX
        XMEANG = XMEANG + XMEAN
        !Write(6,*) 'Control', XMEAN, XMAX
      End Subroutine Control_PrecisionG

      Subroutine Control_Precision_tau(A,B,Ndim)
        Implicit none
        
        Integer :: Ndim
        Complex (Kind=8) :: A(Ndim,Ndim), B(Ndim,Ndim) 
        Real    (Kind=8) :: XMAX, XMEAN

        !Local 
        NCG_tau = NCG_tau + 1
        XMEAN = 0.d0
        XMAX  = 0.d0
        CALL COMPARE(A, B, XMAX, XMEAN)
        IF (XMAX  >  XMAX_tau) XMAX_tau = XMAX
        XMEAN_tau = XMEAN_tau + XMEAN
        !Write(6,*) 'Control_tau', XMEAN, XMAX
      End Subroutine Control_Precision_tau


      Subroutine Control_PrecisionP(Z,Z1)
        Implicit none
        Complex (Kind=8), INTENT(IN) :: Z,Z1
        Real    (Kind=8) :: X
        X = sqrt(dble((Z-Z1)*conjg(Z-Z1)))
        if ( X > XMAXP ) XMAXP = X
      End Subroutine Control_PrecisionP
      
      
      Subroutine control_Print
        Implicit none
#include "machine"
#ifdef MPI
        include 'mpif.h'
#endif
        Real (Kind=8) :: Time, Acc
#ifdef MPI
        REAL (KIND=8)  :: X
        Integer        :: Ierr, Isize, Irank
        INTEGER        :: STATUS(MPI_STATUS_SIZE)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
        
        ACC = 0.d0
        IF (NC_up > 0 )  ACC = dble(ACC_up)/dble(NC_up)
        Call CPU_TIME(CPU_time_en)
        Time = CPU_time_en -  CPU_time_st
#ifdef MPI
        X = 0.d0
        CALL MPI_REDUCE(XMEANG,X,1,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
        XMEANG = X/dble(Isize)
        X = 0.d0
        CALL MPI_REDUCE(XMEAN_tau,X,1,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
        XMEAN_tau = X/dble(Isize)
        X = 0.d0
        CALL MPI_REDUCE(ACC,X,1,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
        ACC = X/dble(Isize)

        X = 0.d0
        CALL MPI_REDUCE(Time,X,1,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
        Time = X/dble(Isize)

        
        CALL MPI_REDUCE(XMAXG,X,1,MPI_REAL8,MPI_MAX, 0,MPI_COMM_WORLD,IERR)
        XMAXG = X
        CALL MPI_REDUCE(XMAX_tau,X,1,MPI_REAL8,MPI_MAX, 0,MPI_COMM_WORLD,IERR)
        XMAX_tau= X


        CALL MPI_REDUCE(XMAXP,X,1,MPI_REAL8,MPI_MAX, 0,MPI_COMM_WORLD,IERR)
        XMAXP = X
        If (Irank == 0 ) then
#endif

           Open (Unit=50,file="info", status="unknown", position="append")
           If (NCG > 0 ) then 
              XMEANG = XMEANG/dble(NCG)
              Write(50,*) ' Precision Green  Mean, Max : ', XMEANG, XMAXG
              Write(50,*) ' Precision Phase, Max       : ', XMAXP
           endif
           If ( NCG_tau > 0 ) then
              XMEAN_tau = XMEAN_tau/dble(NCG_tau)
              Write(50,*) ' Precision tau    Mean, Max : ', XMEAN_tau, XMAX_tau
           endif
           Write(50,*) ' Acceptance                 : ', ACC
           Write(50,*) ' CPU Time                   : ', Time
           Close(50)
#ifdef MPI
        endif
#endif
      end Subroutine Control_Print
      
      subroutine make_truncation(prog_truncation,cpu_max,time_bin_start,time_bin_end)
      !!!!!!! Written by M. Bercx
      ! This subroutine checks if the conditions for a controlled termination of the program are met.
      ! The subroutine contains a hard-coded threshold (in unit of bins): 
      ! if time_remain/time_bin_duration < threshold the program terminates.

#include "machine"
#ifdef MPI
  include 'mpif.h'
#endif  
    
      logical, intent(out)     :: prog_truncation
      real(kind=8), intent(in) :: cpu_max, time_bin_start, time_bin_end
      real(kind=8)             :: time_bin_duration,time_alloc_end,time_remain,bins_remain,threshold
#ifdef MPI   
      real(kind=8)             :: bins_remain_mpi
      integer                  :: err_mpi,rank_mpi,tasks_mpi
#endif
      threshold = 1.5d0
      prog_truncation = .false.
      
#ifdef MPI
      call mpi_comm_size(mpi_comm_world, tasks_mpi, err_mpi)
      call mpi_comm_rank(mpi_comm_world, rank_mpi, err_mpi) 
#endif    
      time_alloc_end    = cpu_time_st + cpu_max*3600 
      time_bin_duration = time_bin_end-time_bin_start
      time_remain       = time_alloc_end - time_bin_end
      bins_remain       = time_remain/time_bin_duration
      
#ifdef MPI
      call mpi_reduce(bins_remain,bins_remain_mpi,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,err_mpi)
#endif
    
#ifdef MPI
      if (rank_mpi .eq. 0) bins_remain_mpi = bins_remain_mpi/tasks_mpi
      call mpi_bcast(bins_remain_mpi,1, mpi_double_precision,0, mpi_comm_world,err_mpi) 
      if (bins_remain_mpi .lt. threshold) prog_truncation = .true.
#else
      if (bins_remain .lt. threshold) prog_truncation = .true.
#endif
      end subroutine make_truncation
     

    end module control
  
  
