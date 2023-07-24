!  Copyright (C) 2016 - 2022 The ALF project
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
!     along with ALF. If not, see http://www.gnu.org/licenses/.
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


!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module handles the  calculation of the acceptance ratio.
!> It also monitors the precision of the code, as well as the timing.
!
!--------------------------------------------------------------------


module Control

    Use MyMats
    use iso_fortran_env, only: output_unit, error_unit
    Implicit none

    real    (Kind=Kind(0.d0)), private, save :: XMEANG, XMAXG, XMAXP,  Xmean_tau, Xmax_tau
    Integer (Kind=Kind(0.d0)), private, save :: count_CPU_start,count_CPU_end,count_rate,count_max
    Integer          , private, save :: NCG, NCG_tau
    Integer (Kind=Kind(0.d0)) , private, save :: NC_up, ACC_up
    Integer (Kind=Kind(0.d0)) , private, save :: NC_eff_up, ACC_eff_up
    Integer (Kind=kind(0.d0)),  private, save :: NC_Glob_up, ACC_Glob_up
    Integer (Kind=kind(0.d0)),  private, save :: NC_HMC_up, ACC_HMC_up
    Integer (Kind=kind(0.d0)),  private, save :: NC_Temp_up, ACC_Temp_up
    real    (Kind=Kind(0.d0)),  private, save :: XMAXP_Glob, XMEANP_Glob
    Integer (Kind=Kind(0.d0)),  private, save :: NC_Phase_GLob

    real    (Kind=Kind(0.d0)),  private, save :: XMAXP_HMC, XMEANP_HMC
    Integer (Kind=Kind(0.d0)),  private, save :: NC_Phase_HMC

    
    real    (Kind=Kind(0.d0)),  private, save :: size_clust_Glob_up, size_clust_Glob_ACC_up

    real    (Kind=Kind(0.d0)),  private, save :: Force_max, Force_mean
    Integer, private, save  :: Force_Count
#ifdef MPI
    Integer                  ,  private, save :: Ierr, Isize, Irank, irank_g, isize_g, igroup
#endif

    
    Contains

      subroutine control_init(Group_Comm)
#ifdef MPI
        Use mpi
#endif
        Implicit none
        Integer       , INTENT(IN)               :: Group_Comm

        XMEANG     = 0.d0
        XMEAN_tau  = 0.d0
        XMAXG      = 0.d0
        XMAX_tau   = 0.d0
        XMAXP      = 0.d0
        XMEANP_Glob= 0.d0
        XMAXP_Glob = 0.d0
        XMEANP_HMC = 0.d0
        XMAXP_HMC  = 0.d0

        
        NCG          = 0
        NCG_tau      = 0
        NC_up        = 0
        ACC_up       = 0
        NC_eff_up    = 0
        ACC_eff_up   = 0
        NC_Glob_up   = 0
        ACC_Glob_up  = 0
        NC_Phase_GLob= 0

        NC_Phase_HMC = 0
        NC_HMC_up    = 0
        ACC_HMC_up   = 0
        
        NC_Temp_up   = 0
        ACC_Temp_up  = 0

        size_clust_Glob_up    = 0.d0
        size_clust_Glob_ACC_up= 0.d0

        Force_max  = 0.d0
        Force_mean = 0.d0
        Force_count = 0
        
#ifdef MPI
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        call MPI_Comm_rank(Group_Comm, irank_g, ierr)
        call MPI_Comm_size(Group_Comm, isize_g, ierr)
        igroup           = irank/isize_g
#endif
        
        call system_clock(count_CPU_start,count_rate,count_max)
      end subroutine control_init


!-------------------------------------------------------------

      Subroutine Control_Langevin(Forces, Group_Comm)


        Implicit none
        
        Complex (Kind=Kind(0.d0)), Intent(In)  :: Forces(:,:)
        Integer, Intent(IN) :: Group_Comm
        
        Integer :: n1,n2, n, nt 
        Real (Kind = Kind(0.d0) ) :: X

        ! Test for not a  number
        n1 =  size(Forces,1)
        n2 =  size(Forces,2)
        Force_count =  Force_count  + 1

        X = 0.d0
        do  n = 1,n1
           do nt =1,n2
              If ( abs( Real(Forces(n,nt),kind(0.d0))) >=  Force_max  ) &
                   &  Force_max = abs( Real(Forces(n,nt),kind(0.d0)))
              X  = X + abs( Real(Forces(n,nt),kind(0.d0)) )
           enddo
        enddo
        Force_mean = Force_mean  +  X/Real(n1*n2,Kind(0.d0)) 
        
      end Subroutine Control_Langevin

      
      Subroutine Control_upgrade(toggle)
        Implicit none
        Logical :: toggle
        NC_up = NC_up + 1
        if (toggle) ACC_up = ACC_up + 1
      end Subroutine Control_upgrade

      Subroutine Control_upgrade_eff(toggle)
        Implicit none
        Logical :: toggle
        NC_eff_up = NC_eff_up + 1
        if (toggle) ACC_eff_up = ACC_eff_up + 1
      end Subroutine Control_upgrade_eff

      Subroutine Control_upgrade_Temp(toggle)
        Implicit none
        Logical :: toggle
        NC_Temp_up  =  NC_Temp_up    +  1
        if (toggle)  ACC_Temp_up  = ACC_Temp_up + 1
      end Subroutine Control_upgrade_Temp

      Subroutine Control_upgrade_Glob(toggle,size_clust)
        Implicit none
        Logical :: toggle
        Real (Kind=Kind(0.d0)), intent(in) :: size_clust
        NC_Glob_up = NC_Glob_up + 1
        size_clust_Glob_up = size_clust_Glob_up + size_clust
        if (toggle) then
          ACC_Glob_up = ACC_Glob_up + 1
          size_clust_Glob_ACC_up = size_clust_Glob_ACC_up + size_clust
        endif
      end Subroutine Control_upgrade_Glob

      Subroutine Control_upgrade_HMC(toggle)
        Implicit none
        Logical :: toggle
        NC_HMC_up = NC_HMC_up + 1
        if (toggle) then
           ACC_HMC_up = ACC_HMC_up + 1
        endif
      end Subroutine Control_upgrade_HMC


      Subroutine Control_PrecisionG(A,B,Ndim)
#ifdef MPI
        Use mpi
#endif
        use runtime_error_mod
        Implicit none

        Integer :: Ndim
        Complex (Kind=Kind(0.d0)) :: A(Ndim,Ndim), B(Ndim,Ndim)
        Real    (Kind=Kind(0.d0)) :: XMAX, XMEAN
        Character (len=64) :: file1 

        if (any(A /= A)) then
#if defined(TEMPERING) 
          write(File1,'(A,I0,A)') "Temp_",igroup,"/info"
#else
          File1 = "info"
#endif
          Open (Unit=50,file=file1, status="unknown", position="append")
          write(50,*)
#ifdef MPI
          write(50,*) "Task", Irank_g, "of group", igroup, "reports:"
#endif
          write(50,*) "Green function A contains NaN, calculation is being aborted!"
          write(50,*)
          close(50)
          write(error_unit,*)
#ifdef MPI
          write(error_unit,*) "Task", Irank_g, "of group", igroup, "reports:"
#endif
          write(error_unit,*) "Green function A contains NaN, calculation is being aborted!"
          write(error_unit,*)
          CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        if (any(B /= B)) then
#if defined(TEMPERING) 
          write(File1,'(A,I0,A)') "Temp_",igroup,"/info"
#else
          File1 = "info"
#endif
          Open (Unit=50,file=file1, status="unknown", position="append")
          write(50,*)
#ifdef MPI
          write(50,*) "Task", Irank_g, "of group", igroup, "reports:"
#endif
          write(50,*) "Green function B contains NaN, calculation is being aborted!"
          write(50,*)
          close(50)
          write(error_unit,*)
#ifdef MPI
          write(error_unit,*) "Task", Irank_g, "of group", igroup, "reports:"
#endif
          write(error_unit,*) "Green function B contains NaN, calculation is being aborted!"
          write(error_unit,*)
          CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        NCG = NCG + 1
        CALL COMPARE(A, B, XMAX, XMEAN)
        IF (XMAX  >  10.d0) then
#if defined(TEMPERING) 
          write(File1,'(A,I0,A)') "Temp_",igroup,"/info"
#else
          File1 = "info"
#endif
          Open (Unit=50,file=file1, status="unknown", position="append")
          write(50,*)
#ifdef MPI
          write(50,*) "Task", Irank_g, "of group", igroup, "reports:"
#endif
          write(50,*) XMAX, " is exceeding the threshold of 10 for G difference!"
          write(50,*) (XmeanG+Xmean)/ncg, " is the average deviation!"
          write(50,*) "This calculation is unstable and therefore aborted!!!"
          write(50,*)
          close(50)
          write(error_unit,*)
#ifdef MPI
          write(error_unit,*) "Task", Irank_g, "of group", igroup, "reports:"
#endif
          write(error_unit,*) XMAX, " is exceeding the threshold of 10 for G difference!"
          write(error_unit,*) (XmeanG+Xmean)/ncg, " is the average deviation!"
          write(error_unit,*) "This calculation is unstable and therefore aborted!!!"
          write(error_unit,*) 'Try with smaller Nwrap or dtau.'
          write(error_unit,*)
          
          CALL Terminate_on_error(ERROR_UNSTABLE_MATRIX,__FILE__,__LINE__)
          
        endif
        IF (XMAX  >  XMAXG) XMAXG = XMAX
        XMEANG = XMEANG + XMEAN
      End Subroutine Control_PrecisionG

      Subroutine Control_Precision_tau(A,B,Ndim)
        Implicit none

        Integer :: Ndim
        Complex (Kind=Kind(0.d0)) :: A(Ndim,Ndim), B(Ndim,Ndim)
        Real    (Kind=Kind(0.d0)) :: XMAX, XMEAN

        NCG_tau = NCG_tau + 1
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

      Subroutine Control_PrecisionP_HMC(Z,Z1)
        Implicit none
        Complex (Kind=Kind(0.D0)), INTENT(IN) :: Z,Z1
        Real    (Kind=Kind(0.D0)) :: X
        X = ABS(Z-Z1)
        if ( X > XMAXP_HMC ) XMAXP_HMC = X
        XMEANP_HMC = XMEANP_HMC + X
        NC_Phase_HMC = NC_Phase_HMC + 1
      End Subroutine Control_PrecisionP_HMC


      Subroutine Control_Print(Group_Comm, Global_update_scheme)
#ifdef MPI
        Use mpi
#endif
        Implicit none

        Integer, Intent(IN) :: Group_Comm
        Character (Len = 64), Intent(IN) :: Global_update_scheme
                

        Character (len=64) :: file1
        Real (Kind=Kind(0.d0)) :: Time, Acc, Acc_eff, Acc_Glob, Acc_Temp, size_clust_Glob, size_clust_Glob_ACC, Acc_HMC
#ifdef MPI
        REAL (Kind=Kind(0.d0))  :: X
        Integer        :: Ierr, Isize, Irank, irank_g, isize_g, igroup

        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        call MPI_Comm_rank(Group_Comm, irank_g, ierr)
        call MPI_Comm_size(Group_Comm, isize_g, ierr)
        igroup           = irank/isize_g
#endif

        ACC = 0.d0
        IF (NC_up > 0 )  ACC = dble(ACC_up)/dble(NC_up)
        ACC_eff = 0.d0
        IF (NC_eff_up > 0 )  ACC_eff = dble(ACC_eff_up)/dble(NC_eff_up)
        ACC_Glob = 0.d0
        size_clust_Glob = 0.d0
        size_clust_Glob_ACC = 0.d0
        IF (NC_Glob_up    > 0 )  then
          ACC_Glob    = dble(ACC_Glob_up)/dble(NC_Glob_up)
          size_clust_Glob     = size_clust_Glob_up     / dble( NC_Glob_up)
          size_clust_Glob_ACC = size_clust_Glob_ACC_up / dble(ACC_Glob_up)
        endif
        ACC_HMC = 0.d0
        IF (NC_HMC_up    > 0 )  then
           ACC_HMC    = dble(ACC_HMC_up)/dble(NC_HMC_up)
        endif

        
        ACC_TEMP = 0.d0
        IF (NC_Temp_up    > 0 )  ACC_Temp    = dble(ACC_Temp_up)/dble(NC_Temp_up)
        IF (NC_Phase_GLob > 0 ) XMEANP_Glob  = XMEANP_Glob/dble(NC_Phase_GLob)

        call system_clock(count_CPU_end)
        time = (count_CPU_end-count_CPU_start)/dble(count_rate)
        if (count_CPU_end .lt. count_CPU_start) time = (count_max+count_CPU_end-count_CPU_start)/dble(count_rate)
        If (trim(Global_update_scheme) == "Langevin") Force_mean =  Force_mean/real(Force_count,kind(0.d0)) 
        
#if defined(MPI)
        If (trim(Global_update_scheme) == "Langevin")  then
           X = 0.d0
           CALL MPI_REDUCE(Force_mean,X,1,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
           Force_mean= X/dble(Isize_g)
           CALL MPI_REDUCE(Force_max,X,1,MPI_REAL8,MPI_MAX, 0,Group_Comm,IERR)
           Force_max= X
        endif
        X = 0.d0
        CALL MPI_REDUCE(ACC,X,1,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
        ACC = X/dble(Isize_g)
        X = 0.d0
        CALL MPI_REDUCE(ACC_eff,X,1,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
        ACC_eff = X/dble(Isize_g)
        X = 0.d0
        CALL MPI_REDUCE(ACC_Glob,X,1,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
        ACC_Glob = X/dble(Isize_g)
        X = 0.d0
        CALL MPI_REDUCE(ACC_HMC,X,1,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
        ACC_HMC = X/dble(Isize_g)
        X = 0.d0
        CALL MPI_REDUCE(ACC_Temp ,X,1,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
        ACC_Temp  = X/dble(Isize_g)

        X = 0.d0
        CALL MPI_REDUCE(size_clust_Glob,X,1,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
        size_clust_Glob = X/dble(Isize_g)
        X = 0.d0
        CALL MPI_REDUCE(size_clust_Glob_ACC,X,1,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
        size_clust_Glob_ACC = X/dble(Isize_g)

        X = 0.d0
        CALL MPI_REDUCE(XMEANG,X,1,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
        XMEANG = X/dble(Isize_g)
        X = 0.d0
        CALL MPI_REDUCE(XMEAN_tau,X,1,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
        XMEAN_tau = X/dble(Isize_g)

        X = 0.d0
        CALL MPI_REDUCE(XMEANP_Glob,X,1,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
        XMEANP_Glob = X/dble(Isize_g)

        X = 0.d0
        CALL MPI_REDUCE(Time,X,1,MPI_REAL8,MPI_SUM, 0,Group_Comm,IERR)
        Time = X/dble(Isize_g)


        CALL MPI_REDUCE(XMAXG,X,1,MPI_REAL8,MPI_MAX, 0,Group_Comm,IERR)
        XMAXG = X
        CALL MPI_REDUCE(XMAX_tau,X,1,MPI_REAL8,MPI_MAX, 0,Group_Comm,IERR)
        XMAX_tau= X

        CALL MPI_REDUCE(XMAXP,X,1,MPI_REAL8,MPI_MAX, 0,Group_Comm,IERR)
        XMAXP = X

        CALL MPI_REDUCE(XMAXP_GLOB,X,1,MPI_REAL8,MPI_MAX, 0,Group_Comm,IERR)
        XMAXP_GLOB = X

        CALL MPI_REDUCE(XMAXP_HMC,X,1,MPI_REAL8,MPI_MAX, 0,Group_Comm,IERR)
        XMAXP_HMC = X

#endif

#if defined(TEMPERING)
        write(File1,'(A,I0,A)') "Temp_",igroup,"/info"
#else
        File1 = "info"
#endif

#if defined(MPI)
        If (Irank_g == 0 ) then
#endif

           Open (Unit=50,file=file1, status="unknown", position="append")
           If (NCG > 0 ) then
              XMEANG = XMEANG/dble(NCG)
              Write(50,*) ' Precision Green  Mean, Max : ', XMEANG, XMAXG
              Write(50,*) ' Precision Phase, Max       : ', XMAXP
           endif
           If ( NCG_tau > 0 ) then
              XMEAN_tau = XMEAN_tau/dble(NCG_tau)
              Write(50,*) ' Precision tau    Mean, Max : ', XMEAN_tau, XMAX_tau
           endif
           If ( NC_up > 0 ) then
              Write(50,*) ' Acceptance                 : ', ACC
           Endif
           If ( NC_eff_up > 0 ) then
              Write(50,*) ' Effective Acceptance       : ', ACC_eff
           Endif
#if defined(TEMPERING)
           Write(50,*) ' Acceptance Tempering       : ', ACC_Temp
#endif
           If (ACC_Glob > 1.D-200 ) then
              Write(50,*) ' Acceptance_Glob              : ', ACC_Glob
              Write(50,*) ' Mean Phase diff Glob         : ', XMEANP_Glob
              Write(50,*) ' Max  Phase diff Glob         : ', XMAXP_Glob
              Write(50,*) ' Average cluster size         : ', size_clust_Glob
              Write(50,*) ' Average accepted cluster size: ', size_clust_Glob_ACC
           endif
           if (trim(Global_update_scheme) == "Langevin") &
                &  Write(50,*) ' Langevin         Mean, Max : ', Force_mean,  Force_max
           
           if (trim(Global_update_scheme) == "HMC")   Then
              Write(50,*) ' Acceptance_HMC              : ', ACC_HMC
              Write(50,*) ' Mean Phase diff HMC         : ', XMEANP_HMC
              Write(50,*) ' Max  Phase diff HMC         : ', XMAXP_HMC
           Endif
           
           Write(50,*) ' CPU Time                   : ', Time
           Close(50)
#if defined(MPI)
        endif
#endif

      end Subroutine Control_Print


      subroutine make_truncation(prog_truncation,cpu_max,count_bin_start,count_bin_end, group_comm)
      !!!!!!! Written by M. Bercx, edited by J. Schwab
      ! This subroutine checks if the conditions for a controlled termination of the program are met.
      ! The subroutine contains a hard-coded threshold (in unit of bins):
      ! if time_remain/time_bin_duration < threshold the program terminates.

#ifdef MPI
      Use mpi
#endif

      logical, intent(out)                 :: prog_truncation
      real(kind=kind(0.d0)), intent(in)    :: cpu_max
      integer(kind=kind(0.d0)), intent(in) :: count_bin_start, count_bin_end
      integer, intent(in)                  :: group_comm
      real(kind=kind(0.d0))                :: count_alloc_end
      real(kind=kind(0.d0))                :: time_bin_duration,time_remain,bins_remain,threshold
#ifdef MPI
      real(kind=kind(0.d0))                :: bins_remain_mpi
      integer                              :: err_mpi, irank, isize, irank_g, isize_g
#endif
      threshold = 1.5d0
      prog_truncation = .false.

#ifdef MPI
      call mpi_comm_size(mpi_comm_world, isize, err_mpi)
      call mpi_comm_rank(mpi_comm_world, irank, err_mpi)
      call mpi_comm_size(group_comm, isize_g, err_mpi)
      call mpi_comm_rank(group_comm, irank_g, err_mpi)
#endif
      count_alloc_end   = count_CPU_start + cpu_max*3600*count_rate
      time_bin_duration = (count_bin_end-count_bin_start)/dble(count_rate)
      time_remain       = (count_alloc_end - count_bin_end)/dble(count_rate)
      if (count_bin_end .lt. count_bin_start) then ! the counter has wrapped around
         time_bin_duration = (count_max+count_bin_end-count_bin_start)/dble(count_rate)
         time_remain       = (count_alloc_end - count_bin_end-count_max)/dble(count_rate)
      endif
      bins_remain       = time_remain/time_bin_duration

#ifdef MPI
#ifdef PARALLEL_PARAMS
      call mpi_reduce(bins_remain,bins_remain_mpi,1,mpi_double_precision,mpi_sum,0,group_comm,err_mpi)
      if (irank_g .eq. 0) bins_remain_mpi = bins_remain_mpi/isize_g
      call mpi_bcast(bins_remain_mpi,1, mpi_double_precision,0, group_comm,err_mpi)
#else
      call mpi_reduce(bins_remain,bins_remain_mpi,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,err_mpi)
      if (irank .eq. 0) bins_remain_mpi = bins_remain_mpi/isize
      call mpi_bcast(bins_remain_mpi,1, mpi_double_precision,0, mpi_comm_world,err_mpi)
#endif
      bins_remain = bins_remain_mpi
#endif
      if (bins_remain .lt. threshold) prog_truncation = .true.
      end subroutine make_truncation

    end module control
