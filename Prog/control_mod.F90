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

   use files_mod
   use MyMats
   use iso_fortran_env, only: output_unit, error_unit
   implicit none

   real(Kind=kind(0.d0)), private, save :: XMEANG, XMAXG, XMAXP, Xmean_tau, Xmax_tau
   integer(Kind=kind(0.d0)), private, save :: count_CPU_start, count_CPU_end, count_rate, count_max
   integer, private, save :: NCG, NCG_tau
   integer(Kind=kind(0.d0)), private, save :: NC_up, ACC_up
   integer(Kind=kind(0.d0)), private, save :: NC_eff_up, ACC_eff_up
   integer(Kind=kind(0.d0)), private, save :: NC_Glob_up, ACC_Glob_up
   integer(Kind=kind(0.d0)), private, save :: NC_HMC_up, ACC_HMC_up
   integer(Kind=kind(0.d0)), private, save :: NC_Temp_up, ACC_Temp_up
   real(Kind=kind(0.d0)), private, save :: XMAXP_Glob, XMEANP_Glob
   integer(Kind=kind(0.d0)), private, save :: NC_Phase_GLob

   real(Kind=kind(0.d0)), private, save :: XMAXP_HMC, XMEANP_HMC
   integer(Kind=kind(0.d0)), private, save :: NC_Phase_HMC

   real(Kind=kind(0.d0)), private, save :: size_clust_Glob_up, size_clust_Glob_ACC_up

   real(Kind=kind(0.d0)), private, save :: Force_max, Force_mean
   integer, private, save  :: Force_Count
#ifdef MPI
   integer, private, save :: Ierr, Isize, Irank, irank_g, isize_g, igroup
#endif

contains

   subroutine control_init(Group_Comm)
#ifdef MPI
      use mpi
#endif
      implicit none
      integer, intent(IN)               :: Group_Comm

      XMEANG = 0.d0
      XMEAN_tau = 0.d0
      XMAXG = 0.d0
      XMAX_tau = 0.d0
      XMAXP = 0.d0
      XMEANP_Glob = 0.d0
      XMAXP_Glob = 0.d0
      XMEANP_HMC = 0.d0
      XMAXP_HMC = 0.d0

      NCG = 0
      NCG_tau = 0
      NC_up = 0
      ACC_up = 0
      NC_eff_up = 0
      ACC_eff_up = 0
      NC_Glob_up = 0
      ACC_Glob_up = 0
      NC_Phase_GLob = 0

      NC_Phase_HMC = 0
      NC_HMC_up = 0
      ACC_HMC_up = 0

      NC_Temp_up = 0
      ACC_Temp_up = 0

      size_clust_Glob_up = 0.d0
      size_clust_Glob_ACC_up = 0.d0

      Force_max = 0.d0
      Force_mean = 0.d0
      Force_count = 0

#ifdef MPI
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif

      call system_clock(count_CPU_start, count_rate, count_max)
   end subroutine control_init

!-------------------------------------------------------------

   subroutine Control_Langevin(Forces, Group_Comm)

      implicit none

      complex(Kind=kind(0.d0)), intent(In)  :: Forces(:, :)
      integer, intent(IN) :: Group_Comm

      integer :: n1, n2, n, nt
      real(Kind=kind(0.d0)) :: X

      ! Test for not a  number
      n1 = size(Forces, 1)
      n2 = size(Forces, 2)
      Force_count = Force_count + 1

      X = 0.d0
      do n = 1, n1
         do nt = 1, n2
            if (abs(real(Forces(n, nt), kind(0.d0))) >= Force_max) &
                 &  Force_max = abs(real(Forces(n, nt), kind(0.d0)))
            X = X + abs(real(Forces(n, nt), kind(0.d0)))
         end do
      end do
      Force_mean = Force_mean + X/real(n1*n2, kind(0.d0))

   end subroutine Control_Langevin

   subroutine Control_upgrade(toggle)
      implicit none
      logical :: toggle
      NC_up = NC_up + 1
      if (toggle) ACC_up = ACC_up + 1
   end subroutine Control_upgrade

   subroutine Control_upgrade_eff(toggle)
      implicit none
      logical :: toggle
      NC_eff_up = NC_eff_up + 1
      if (toggle) ACC_eff_up = ACC_eff_up + 1
   end subroutine Control_upgrade_eff

   subroutine Control_upgrade_Temp(toggle)
      implicit none
      logical :: toggle
      NC_Temp_up = NC_Temp_up + 1
      if (toggle) ACC_Temp_up = ACC_Temp_up + 1
   end subroutine Control_upgrade_Temp

   subroutine Control_upgrade_Glob(toggle, size_clust)
      implicit none
      logical :: toggle
      real(Kind=kind(0.d0)), intent(in) :: size_clust
      NC_Glob_up = NC_Glob_up + 1
      size_clust_Glob_up = size_clust_Glob_up + size_clust
      if (toggle) then
         ACC_Glob_up = ACC_Glob_up + 1
         size_clust_Glob_ACC_up = size_clust_Glob_ACC_up + size_clust
      end if
   end subroutine Control_upgrade_Glob

   subroutine Control_upgrade_HMC(toggle)
      implicit none
      logical :: toggle
      NC_HMC_up = NC_HMC_up + 1
      if (toggle) then
         ACC_HMC_up = ACC_HMC_up + 1
      end if
   end subroutine Control_upgrade_HMC

   subroutine Control_PrecisionG(A, B, Ndim)
#ifdef MPI
      use mpi
#endif
      use runtime_error_mod
      implicit none

      integer :: Ndim
      complex(Kind=kind(0.d0)) :: A(Ndim, Ndim), B(Ndim, Ndim)
      real(Kind=kind(0.d0)) :: XMAX, XMEAN
      character(len=64) :: file1

      if (any(A /= A)) then
#if defined(TEMPERING)
         write (File1, '(A,I0,A)') "Temp_", igroup, "/info"
#else
         File1 = "info"
#endif
         open (Unit=50, file=file1, status="unknown", position="append")
         write (50, *)
#ifdef MPI
         write (50, *) "Task", Irank_g, "of group", igroup, "reports:"
#endif
         write (50, *) "Green function A contains NaN, calculation is being aborted!"
         write (50, *)
         close (50)
         write (error_unit, *)
#ifdef MPI
         write (error_unit, *) "Task", Irank_g, "of group", igroup, "reports:"
#endif
         write (error_unit, *) "Green function A contains NaN, calculation is being aborted!"
         write (error_unit, *)
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      if (any(B /= B)) then
#if defined(TEMPERING)
         write (File1, '(A,I0,A)') "Temp_", igroup, "/info"
#else
         File1 = "info"
#endif
         open (Unit=50, file=file1, status="unknown", position="append")
         write (50, *)
#ifdef MPI
         write (50, *) "Task", Irank_g, "of group", igroup, "reports:"
#endif
         write (50, *) "Green function B contains NaN, calculation is being aborted!"
         write (50, *)
         close (50)
         write (error_unit, *)
#ifdef MPI
         write (error_unit, *) "Task", Irank_g, "of group", igroup, "reports:"
#endif
         write (error_unit, *) "Green function B contains NaN, calculation is being aborted!"
         write (error_unit, *)
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      NCG = NCG + 1
      call COMPARE(A, B, XMAX, XMEAN)
      if (XMAX > 10.d0) then
#if defined(TEMPERING)
         write (File1, '(A,I0,A)') "Temp_", igroup, "/info"
#else
         File1 = "info"
#endif
         open (Unit=50, file=file1, status="unknown", position="append")
         write (50, *)
#ifdef MPI
         write (50, *) "Task", Irank_g, "of group", igroup, "reports:"
#endif
         write (50, *) XMAX, " is exceeding the threshold of 10 for G difference!"
         write (50, *) (XmeanG + Xmean)/ncg, " is the average deviation!"
         write (50, *) "This calculation is unstable and therefore aborted!!!"
         write (50, *)
         close (50)
         write (error_unit, *)
#ifdef MPI
         write (error_unit, *) "Task", Irank_g, "of group", igroup, "reports:"
#endif
         write (error_unit, *) XMAX, " is exceeding the threshold of 10 for G difference!"
         write (error_unit, *) (XmeanG + Xmean)/ncg, " is the average deviation!"
         write (error_unit, *) "This calculation is unstable and therefore aborted!!!"
         write (error_unit, *) 'Try with smaller Nwrap or dtau.'
         write (error_unit, *)

         call Terminate_on_error(ERROR_UNSTABLE_MATRIX, __FILE__, __LINE__)

      end if
      if (XMAX > XMAXG) XMAXG = XMAX
      XMEANG = XMEANG + XMEAN
   end subroutine Control_PrecisionG

   subroutine Control_Precision_tau(A, B, Ndim)
      implicit none

      integer :: Ndim
      complex(Kind=kind(0.d0)) :: A(Ndim, Ndim), B(Ndim, Ndim)
      real(Kind=kind(0.d0)) :: XMAX, XMEAN

      NCG_tau = NCG_tau + 1
      call COMPARE(A, B, XMAX, XMEAN)
      if (XMAX > XMAX_tau) XMAX_tau = XMAX
      XMEAN_tau = XMEAN_tau + XMEAN
   end subroutine Control_Precision_tau

   subroutine Control_PrecisionP(Z, Z1)
      implicit none
      complex(Kind=kind(0.d0)), intent(IN) :: Z, Z1
      real(Kind=kind(0.d0)) :: X
      X = abs(Z - Z1)
      if (X > XMAXP) XMAXP = X
   end subroutine Control_PrecisionP

   subroutine Control_PrecisionP_Glob(Z, Z1)
      implicit none
      complex(Kind=kind(0.d0)), intent(IN) :: Z, Z1
      real(Kind=kind(0.d0)) :: X
      X = abs(Z - Z1)
      if (X > XMAXP_Glob) XMAXP_Glob = X
      XMEANP_Glob = XMEANP_Glob + X
      NC_Phase_GLob = NC_Phase_GLob + 1
   end subroutine Control_PrecisionP_Glob

   subroutine Control_PrecisionP_HMC(Z, Z1)
      implicit none
      complex(Kind=kind(0.d0)), intent(IN) :: Z, Z1
      real(Kind=kind(0.d0)) :: X
      X = abs(Z - Z1)
      if (X > XMAXP_HMC) XMAXP_HMC = X
      XMEANP_HMC = XMEANP_HMC + X
      NC_Phase_HMC = NC_Phase_HMC + 1
   end subroutine Control_PrecisionP_HMC

   subroutine Control_Print(Group_Comm, Global_update_scheme)
#ifdef MPI
      use mpi
#endif
      implicit none

      integer, intent(IN) :: Group_Comm
      character(Len=64), intent(IN) :: Global_update_scheme

      character(len=64) :: file1
      real(Kind=kind(0.d0)) :: Time, Acc, Acc_eff, Acc_Glob, Acc_Temp, size_clust_Glob, size_clust_Glob_ACC, Acc_HMC
#ifdef MPI
      real(Kind=kind(0.d0))  :: X
      integer        :: Ierr, Isize, Irank, irank_g, isize_g, igroup

      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup = irank/isize_g
#endif

      ACC = 0.d0
      if (NC_up > 0) ACC = dble(ACC_up)/dble(NC_up)
      ACC_eff = 0.d0
      if (NC_eff_up > 0) ACC_eff = dble(ACC_eff_up)/dble(NC_eff_up)
      ACC_Glob = 0.d0
      size_clust_Glob = 0.d0
      size_clust_Glob_ACC = 0.d0
      if (NC_Glob_up > 0) then
         ACC_Glob = dble(ACC_Glob_up)/dble(NC_Glob_up)
         size_clust_Glob = size_clust_Glob_up/dble(NC_Glob_up)
         size_clust_Glob_ACC = size_clust_Glob_ACC_up/dble(ACC_Glob_up)
      end if
      ACC_HMC = 0.d0
      if (NC_HMC_up > 0) then
         ACC_HMC = dble(ACC_HMC_up)/dble(NC_HMC_up)
      end if

      ACC_TEMP = 0.d0
      if (NC_Temp_up > 0) ACC_Temp = dble(ACC_Temp_up)/dble(NC_Temp_up)
      if (NC_Phase_GLob > 0) XMEANP_Glob = XMEANP_Glob/dble(NC_Phase_GLob)

      call system_clock(count_CPU_end)
      time = (count_CPU_end - count_CPU_start)/dble(count_rate)
      if (count_CPU_end .lt. count_CPU_start) time = (count_max + count_CPU_end - count_CPU_start)/dble(count_rate)
      if (str_to_upper(Global_update_scheme) == "LANGEVIN") Force_mean = Force_mean/real(Force_count, kind(0.d0))

#if defined(MPI)
      if (str_to_upper(Global_update_scheme) == "LANGEVIN") then
         X = 0.d0
         call MPI_REDUCE(Force_mean, X, 1, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
         Force_mean = X/dble(Isize_g)
         call MPI_REDUCE(Force_max, X, 1, MPI_REAL8, MPI_MAX, 0, Group_Comm, IERR)
         Force_max = X
      end if
      X = 0.d0
      call MPI_REDUCE(ACC, X, 1, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
      ACC = X/dble(Isize_g)
      X = 0.d0
      call MPI_REDUCE(ACC_eff, X, 1, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
      ACC_eff = X/dble(Isize_g)
      X = 0.d0
      call MPI_REDUCE(ACC_Glob, X, 1, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
      ACC_Glob = X/dble(Isize_g)
      X = 0.d0
      call MPI_REDUCE(ACC_HMC, X, 1, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
      ACC_HMC = X/dble(Isize_g)
      X = 0.d0
      call MPI_REDUCE(ACC_Temp, X, 1, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
      ACC_Temp = X/dble(Isize_g)

      X = 0.d0
      call MPI_REDUCE(size_clust_Glob, X, 1, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
      size_clust_Glob = X/dble(Isize_g)
      X = 0.d0
      call MPI_REDUCE(size_clust_Glob_ACC, X, 1, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
      size_clust_Glob_ACC = X/dble(Isize_g)

      X = 0.d0
      call MPI_REDUCE(XMEANG, X, 1, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
      XMEANG = X/dble(Isize_g)
      X = 0.d0
      call MPI_REDUCE(XMEAN_tau, X, 1, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
      XMEAN_tau = X/dble(Isize_g)

      X = 0.d0
      call MPI_REDUCE(XMEANP_Glob, X, 1, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
      XMEANP_Glob = X/dble(Isize_g)

      X = 0.d0
      call MPI_REDUCE(Time, X, 1, MPI_REAL8, MPI_SUM, 0, Group_Comm, IERR)
      Time = X/dble(Isize_g)

      call MPI_REDUCE(XMAXG, X, 1, MPI_REAL8, MPI_MAX, 0, Group_Comm, IERR)
      XMAXG = X
      call MPI_REDUCE(XMAX_tau, X, 1, MPI_REAL8, MPI_MAX, 0, Group_Comm, IERR)
      XMAX_tau = X

      call MPI_REDUCE(XMAXP, X, 1, MPI_REAL8, MPI_MAX, 0, Group_Comm, IERR)
      XMAXP = X

      call MPI_REDUCE(XMAXP_GLOB, X, 1, MPI_REAL8, MPI_MAX, 0, Group_Comm, IERR)
      XMAXP_GLOB = X

      call MPI_REDUCE(XMAXP_HMC, X, 1, MPI_REAL8, MPI_MAX, 0, Group_Comm, IERR)
      XMAXP_HMC = X

#endif

#if defined(TEMPERING)
      write (File1, '(A,I0,A)') "Temp_", igroup, "/info"
#else
      File1 = "info"
#endif

#if defined(MPI)
      if (Irank_g == 0) then
#endif

         open (Unit=50, file=file1, status="unknown", position="append")
         if (NCG > 0) then
            XMEANG = XMEANG/dble(NCG)
            write (50, *) ' Precision Green  Mean, Max : ', XMEANG, XMAXG
            write (50, *) ' Precision Phase, Max       : ', XMAXP
         end if
         if (NCG_tau > 0) then
            XMEAN_tau = XMEAN_tau/dble(NCG_tau)
            write (50, *) ' Precision tau    Mean, Max : ', XMEAN_tau, XMAX_tau
         end if
         if (NC_up > 0) then
            write (50, *) ' Acceptance                 : ', ACC
         end if
         if (NC_eff_up > 0) then
            write (50, *) ' Effective Acceptance       : ', ACC_eff
         end if
#if defined(TEMPERING)
         write (50, *) ' Acceptance Tempering       : ', ACC_Temp
#endif
         if (ACC_Glob > 1.d-200) then
            write (50, *) ' Acceptance_Glob              : ', ACC_Glob
            write (50, *) ' Mean Phase diff Glob         : ', XMEANP_Glob
            write (50, *) ' Max  Phase diff Glob         : ', XMAXP_Glob
            write (50, *) ' Average cluster size         : ', size_clust_Glob
            write (50, *) ' Average accepted cluster size: ', size_clust_Glob_ACC
         end if
         if (str_to_upper(Global_update_scheme) == "LANGEVIN") &
              &  write (50, *) ' Langevin         Mean, Max : ', Force_mean, Force_max

         if (str_to_upper(Global_update_scheme) == "HMC") then
            write (50, *) ' Acceptance_HMC              : ', ACC_HMC
            write (50, *) ' Mean Phase diff HMC         : ', XMEANP_HMC
            write (50, *) ' Max  Phase diff HMC         : ', XMAXP_HMC
         end if

         write (50, *) ' CPU Time                   : ', Time
         close (50)
#if defined(MPI)
      end if
#endif

   end subroutine Control_Print

   subroutine make_truncation(prog_truncation, cpu_max, count_bin_start, count_bin_end, group_comm)
      !!!!!!! Written by M. Bercx, edited by J. Schwab
      ! This subroutine checks if the conditions for a controlled termination of the program are met.
      ! The subroutine contains a hard-coded threshold (in unit of bins):
      ! if time_remain/time_bin_duration < threshold the program terminates.

#ifdef MPI
      use mpi
#endif

      logical, intent(out)                 :: prog_truncation
      real(kind=kind(0.d0)), intent(in)    :: cpu_max
      integer(kind=kind(0.d0)), intent(in) :: count_bin_start, count_bin_end
      integer, intent(in)                  :: group_comm
      real(kind=kind(0.d0))                :: count_alloc_end
      real(kind=kind(0.d0))                :: time_bin_duration, time_remain, bins_remain, threshold
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
      count_alloc_end = count_CPU_start + cpu_max*3600*count_rate
      time_bin_duration = (count_bin_end - count_bin_start)/dble(count_rate)
      time_remain = (count_alloc_end - count_bin_end)/dble(count_rate)
      if (count_bin_end .lt. count_bin_start) then ! the counter has wrapped around
         time_bin_duration = (count_max + count_bin_end - count_bin_start)/dble(count_rate)
         time_remain = (count_alloc_end - count_bin_end - count_max)/dble(count_rate)
      end if
      bins_remain = time_remain/time_bin_duration

#ifdef MPI
#ifdef PARALLEL_PARAMS
      call mpi_reduce(bins_remain, bins_remain_mpi, 1, mpi_double_precision, mpi_sum, 0, group_comm, err_mpi)
      if (irank_g .eq. 0) bins_remain_mpi = bins_remain_mpi/isize_g
      call mpi_bcast(bins_remain_mpi, 1, mpi_double_precision, 0, group_comm, err_mpi)
#else
      call mpi_reduce(bins_remain, bins_remain_mpi, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, err_mpi)
      if (irank .eq. 0) bins_remain_mpi = bins_remain_mpi/isize
      call mpi_bcast(bins_remain_mpi, 1, mpi_double_precision, 0, mpi_comm_world, err_mpi)
#endif
      bins_remain = bins_remain_mpi
#endif
      if (bins_remain .lt. threshold) prog_truncation = .true.
   end subroutine make_truncation

end module control
