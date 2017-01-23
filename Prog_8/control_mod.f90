!  Copyright (C) 2016 The ALF project
! 
!  This file is part of the ALF project.
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

  module Control

!--------------------------------------------------------------------
!
!> @author 
!> ALF-project
!
!> @brief 
!> This module handles the  calculation of the acceptance ratio. It also  
!> monitors the precision of the code, as well as the timing.
!
!--------------------------------------------------------------------

    Use MyMats
    Implicit none

    real    (Kind=Kind(0.d0)) , private, save :: XMEANG, XMAXG, XMAXP, CPU_time_st, CPU_time_en, Xmean_tau, Xmax_tau
    Integer          , private, save :: NCG, NCG_tau
    Integer (Kind=Kind(0.d0)) , private, save :: NC_up, ACC_up
    Integer (Kind=kind(0.d0)),  private, save :: NC_Glob_up, ACC_Glob_up
    real    (Kind=Kind(0.d0)),  private, save :: XMAXP_Glob, XMEANP_Glob
    Integer (Kind=Kind(0.d0)),  private, save :: NC_Phase_GLob

    
    Contains

      subroutine control_init
        Implicit none
        XMEANG     = 0.d0
        XMEAN_tau  = 0.d0
        XMAXG      = 0.d0
        XMAX_tau   = 0.d0
        XMAXP      = 0.d0
        XMEANP_Glob= 0.d0 
        XMAXP_Glob = 0.d0

        NCG          = 0
        NCG_tau      = 0
        NC_up        = 0
        ACC_up       = 0
        NC_Glob_up   = 0
        ACC_Glob_up  = 0
        NC_Phase_GLob= 0
        Call CPU_TIME(CPU_time_st)
      end subroutine control_init
      
      Subroutine Control_upgrade(Log) 
        Implicit none
        Logical :: Log
        NC_up = NC_up + 1
        if (Log) ACC_up = ACC_up + 1
      end Subroutine Control_upgrade

      Subroutine Control_upgrade_Glob(Log) 
        Implicit none
        Logical :: Log
        NC_Glob_up = NC_Glob_up + 1
        if (Log) ACC_Glob_up = ACC_Glob_up + 1
      end Subroutine Control_upgrade_Glob

      Subroutine Control_PrecisionG(A,B,Ndim)
        Implicit none
        
        Integer :: Ndim
        Complex (Kind=Kind(0.d0)) :: A(Ndim,Ndim), B(Ndim,Ndim) 
        Real    (Kind=Kind(0.d0)) :: XMAX, XMEAN

        NCG = NCG + 1
        XMEAN = 0.d0
        XMAX  = 0.d0
        CALL COMPARE(A, B, XMAX, XMEAN)
        IF (XMAX  >  XMAXG) XMAXG = XMAX
        XMEANG = XMEANG + XMEAN
      End Subroutine Control_PrecisionG

      Subroutine Control_Precision_tau(A,B,Ndim)
        Implicit none
        
        Integer :: Ndim
        Complex (Kind=Kind(0.d0)) :: A(Ndim,Ndim), B(Ndim,Ndim) 
        Real    (Kind=Kind(0.d0)) :: XMAX, XMEAN

        NCG_tau = NCG_tau + 1
        XMEAN = 0.d0
        XMAX  = 0.d0
        CALL COMPARE(A, B, XMAX, XMEAN)
        IF (XMAX  >  XMAX_tau) XMAX_tau = XMAX
        XMEAN_tau = XMEAN_tau + XMEAN
      End Subroutine Control_Precision_tau


      Subroutine Control_PrecisionP(Z,Z1)
        Implicit none
        Complex (Kind=Kind(0.D0)), INTENT(IN) :: Z,Z1
        Real    (Kind=Kind(0.D0)) :: X
        X = ABS(Z-Z1)
        if ( X > XMAXP ) XMAXP = X
      End Subroutine Control_PrecisionP


      Subroutine Control_PrecisionP_Glob(Z,Z1)
        Implicit none
        Complex (Kind=Kind(0.D0)), INTENT(IN) :: Z,Z1
        Real    (Kind=Kind(0.D0)) :: X
        X = ABS(Z-Z1)
        if ( X > XMAXP_Glob ) XMAXP_Glob = X
        XMEANP_Glob = XMEANP_Glob + X
        NC_Phase_GLob = NC_Phase_GLob + 1
      End Subroutine Control_PrecisionP_Glob
      
      
      Subroutine control_Print
        Implicit none
#ifdef MPI
        include 'mpif.h'
#endif
        Real (Kind=Kind(0.d0)) :: Time, Acc, Acc_Glob
#ifdef MPI
        REAL (Kind=Kind(0.d0))  :: X
        Integer        :: Ierr, Isize, Irank
        INTEGER        :: STATUS(MPI_STATUS_SIZE)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
        
        ACC = 0.d0
        IF (NC_up > 0 )  ACC = dble(ACC_up)/dble(NC_up)
        ACC_Glob = 0.d0
        IF (NC_Glob_up    > 0 )  ACC_Glob    = dble(ACC_Glob_up)/dble(NC_Glob_up)
        IF (NC_Phase_GLob > 0 ) XMEANP_Glob  = XMEANP_Glob/dble(NC_Phase_GLob)

        Call CPU_TIME(CPU_time_en)
        Time = CPU_time_en -  CPU_time_st
#ifdef MPI


        X = 0.d0
        CALL MPI_REDUCE(ACC,X,1,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
        ACC = X/dble(Isize)
        X = 0.d0
        CALL MPI_REDUCE(ACC_Glob,X,1,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
        ACC_Glob = X/dble(Isize)


        X = 0.d0
        CALL MPI_REDUCE(XMEANG,X,1,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
        XMEANG = X/dble(Isize)
        X = 0.d0
        CALL MPI_REDUCE(XMEAN_tau,X,1,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
        XMEAN_tau = X/dble(Isize)

        X = 0.d0
        CALL MPI_REDUCE(XMEANP_Glob,X,1,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
        XMEANP_Glob = X/dble(Isize)

        X = 0.d0
        CALL MPI_REDUCE(Time,X,1,MPI_REAL8,MPI_SUM, 0,MPI_COMM_WORLD,IERR)
        Time = X/dble(Isize)

        
        CALL MPI_REDUCE(XMAXG,X,1,MPI_REAL8,MPI_MAX, 0,MPI_COMM_WORLD,IERR)
        XMAXG = X
        CALL MPI_REDUCE(XMAX_tau,X,1,MPI_REAL8,MPI_MAX, 0,MPI_COMM_WORLD,IERR)
        XMAX_tau= X

        CALL MPI_REDUCE(XMAXP,X,1,MPI_REAL8,MPI_MAX, 0,MPI_COMM_WORLD,IERR)
        XMAXP = X

        CALL MPI_REDUCE(XMAXP_GLOB,X,1,MPI_REAL8,MPI_MAX, 0,MPI_COMM_WORLD,IERR)
        XMAXP_GLOB = X
        

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
           !If (ACC_Glob > 1.D-200 ) then
              Write(50,*) ' Acceptance_Glob            : ', ACC_Glob
              Write(50,*) ' Mean Phase diff Glob       : ', XMEANP_Glob 
              Write(50,*) ' Max  Phase diff Glob       : ', XMAXP_Glob
           !endif
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

#ifdef MPI
  include 'mpif.h'
#endif  
    
      logical, intent(out)     :: prog_truncation
      real(Kind=Kind(0.d0)), intent(in) :: cpu_max, time_bin_start, time_bin_end
      real(Kind=Kind(0.d0))             :: time_bin_duration,time_alloc_end,time_remain,bins_remain,threshold
#ifdef MPI   
      real(Kind=Kind(0.d0))             :: bins_remain_mpi
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
  
  
