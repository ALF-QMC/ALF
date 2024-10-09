!  Copyright (C) 2016 - 2022 The ALF project
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
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
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
!
!> @brief
!> Reads in the VAR_QMC and VAR_HAM_NAME namelists from the file parameters, calls Ham_set and carries out the sweeps.
!> If the program is compiled with the Tempering flag on, then the VAR_TEMP namelist will also be read in.
!> If the PARALLEL_PARAMS is set, VAR_QMC and VAR_HAM_NAME are read from Temp_{igroup}/parameters instead.
!>
!> @details
!> \verbatim
!>  The parameters in the VAR_QMC namelist read
!> \endverbatim
!> @param Nwrap Integer
!> \verbatim
!>  Number of time slices between stabilization (QR)
!>  Has to be specified.
!> \endverbatim
!> @param Nsweep Integer
!> \verbatim
!>  Number of sweeps per bin
!>  Has to be specified.
!> \endverbatim
!> @param Nbin Integer
!> \verbatim
!>  Number of bins
!>  Has to be specified.
!> \endverbatim
!> @param Ltau Integer
!> \verbatim
!>  If Ltau=1 time displaced correlations will be measured.
!>  Has to be specified.
!> \endverbatim
!> @param LOBS_ST LOBS_EN Integer
!> \verbatim
!>  Time slice interval for measurements
!>  Default values:  LOBS_ST = Thtrot +1,  LOBS_ST = Ltrot - Thtrot
!>  Note that Thtrot corresponds to the projection time in units of
!>  the time step  and is equal to zero for the finite temperature code.
!> \endverbatim
!> @param CPU_MAX Real
!> \verbatim
!>  Available Wallclock time. The program will carry as many bins as
!>  possible during this time
!>  If not specified the program will stop after NBIN bins are calculated
!> \endverbatim
!> @param Propose_S0 Logical
!> \verbatim
!>  If true, spin flips are proposed with probability exp(-S_0(C')). See documentation.
!>  Default:  Propose_S0=.false.
!> \endverbatim
!> @param Global_moves Logical
!> \verbatim
!>  If true, global moves will be carried out.
!>  Default: Global_moves=.false.
!> \endverbatim
!> @param N_Global Integer
!> \verbatim
!>  Number of global moves per  sequential sweep.
!>  Default: N_Global=0
!> \endverbatim
!> @param Global_tau_moves Logical
!> \verbatim
!>  If true, global moves on a given time slice will be carried out
!>  Default: Global_tau_moves=.false.
!> \endverbatim
!> @param N_Global_tau Integer
!> \verbatim
!>  Number of global_tau moves that will be carried out per time-slice.
!>  Default: N_Global_tau=0
!> \endverbatim
!> @param Nt_sequential_start  Integer
!> @param Nt_sequential_end  Integer
!> \verbatim
!> Interval over which one will carry out sequential updating on a single time slice.
!> Default: Nt_sequential_start = 1  Nt_sequential_end=size(OP_V,1)). This default is
!> automatically if Global_tau_moves=.false.
!> \endverbatim

!--------------------------------------------------------------------

program Main

   use runtime_error_mod
   use Operator_mod
   use Lattices_v3
   use MyMats
   use Hamiltonian_main
   use Control
   use Tau_m_mod
   use Tau_p_mod
   use Hop_mod
   use Global_mod
   use UDV_State_mod
   use Wrapgr_mod
   use Fields_mod
   use WaveFunction_mod
   use entanglement_mod
   use iso_fortran_env, only: output_unit, error_unit
   use Langevin_HMC_mod
   use wrapur_mod
   use wrapul_mod
   use cgr1_mod
   use set_random

#ifdef MPI
   use mpi
#endif
#ifdef HDF5
   use hdf5
   use h5lt
#endif
   implicit none

#include "git.h"

   complex(Kind=kind(0.d0)), dimension(:, :), allocatable   ::  TEST
   complex(Kind=kind(0.d0)), dimension(:, :, :), allocatable    :: GR, GR_Tilde
   class(UDV_State), dimension(:), allocatable :: udvl, udvr
   complex(Kind=kind(0.d0)), dimension(:), allocatable   :: Phase_array

   integer :: Nwrap, NSweep, NBin, NBin_eff, Ltau, NSTM, NT, NT1, NVAR, LOBS_EN, LOBS_ST, NBC, NSW
   integer :: NTAU, NTAU1
   real(Kind=kind(0.d0)) :: CPU_MAX
   character(len=64) :: file_seeds, file_para, file_dat, file_info, ham_name
   integer :: Seed_in
   complex(Kind=kind(0.d0)), allocatable, dimension(:, :) :: Initial_field

   ! Space for choosing sampling scheme
   logical :: Propose_S0, Tempering_calc_det
   logical :: Global_moves, Global_tau_moves
   integer :: N_Global
   integer :: Nt_sequential_start, Nt_sequential_end, mpi_per_parameter_set
   integer :: N_Global_tau
   logical :: Sequential
   real(Kind=kind(0.d0)) ::  Amplitude  !    Needed for  update of  type  3  and  4  fields.

#ifdef HDF5
   integer(HID_T) :: file_id
   logical :: file_exists
#endif
   !  Space for reading in Langevin & HMC  parameters
   logical                      :: Langevin, HMC
   integer                      :: Leapfrog_Steps, N_HMC_sweeps
   real(Kind=kind(0.d0))      :: Delta_t_Langevin_HMC, Max_Force

#if defined(TEMPERING)
   integer :: N_exchange_steps, N_Tempering_frequency
   namelist /VAR_TEMP/ N_exchange_steps, N_Tempering_frequency, mpi_per_parameter_set, Tempering_calc_det
#endif

   namelist /VAR_QMC/ Nwrap, NSweep, NBin, Ltau, LOBS_EN, LOBS_ST, CPU_MAX, &
        &               Propose_S0, Global_moves, N_Global, Global_tau_moves, &
        &               Nt_sequential_start, Nt_sequential_end, N_Global_tau, &
        &               sequential, Langevin, HMC, Delta_t_Langevin_HMC, &
        &               Max_Force, Leapfrog_steps, N_HMC_sweeps, Amplitude

   namelist /VAR_HAM_NAME/ ham_name

   !  General
   integer :: Ierr, I, nf, nf_eff, nst, n, n1, N_op
   logical :: Toggle, Toggle1
   complex(Kind=kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.d0)), Phase, Z, Z1
   real(Kind=kind(0.d0)) :: ZERO = 10d-8, X, X1
   real(Kind=kind(0.d0)) :: Mc_step_weight

   ! Storage for  stabilization steps
   integer, dimension(:), allocatable :: Stab_nt

   ! Space for storage.
   class(UDV_State), dimension(:, :), allocatable :: udvst

   ! For tests
   real(Kind=kind(0.d0)) :: Weight, Weight_tot

   ! For the truncation of the program:
   logical                   :: prog_truncation, run_file_exists
   integer(kind=kind(0.d0)) :: count_bin_start, count_bin_end

   ! For MPI shared memory
   character(64), parameter :: name = "ALF_SHM_CHUNK_SIZE_GB"
   character(64) :: chunk_size_str
   real(Kind=kind(0.d0)) :: chunk_size_gb

#ifdef MPI
   integer        :: Isize, Irank, Irank_g, Isize_g, color, key, igroup, MPI_COMM_i

   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
   call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)

   if (Irank == 0) then
#endif
      write (*, *) "ALF Copyright (C) 2016 - 2022 The ALF project contributors"
      write (*, *) "This Program comes with ABSOLUTELY NO WARRANTY; for details see license.GPL"
      write (*, *) "This is free software, and you are welcome to redistribute it under certain conditions."
#ifdef MPI
   end if
#endif

#if defined(TEMPERING) && defined(MPI)
   mpi_per_parameter_set = 1   ! Default value
   Tempering_calc_det = .true. ! Default value
   open (UNIT=5, FILE='parameters', STATUS='old', ACTION='read', IOSTAT=ierr)
   if (ierr /= 0) then
      write (error_unit, *) 'main: unable to open <parameters>', ierr
      call Terminate_on_error(ERROR_FILE_NOT_FOUND, __FILE__, __LINE__)
   end if
   read (5, NML=VAR_TEMP)
   close (5)
   call MPI_BCAST(N_exchange_steps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(N_Tempering_frequency, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(mpi_per_parameter_set, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(Tempering_calc_det, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   if (mod(ISIZE, mpi_per_parameter_set) .ne. 0) then
      write (error_unit, *) "mpi_per_parameter_set is not a multiple of total mpi processes"
      call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
   end if
   call Global_Tempering_setup
#elif !defined(TEMPERING)  && defined(MPI)
   mpi_per_parameter_set = Isize
#elif defined(TEMPERING)  && !defined(MPI)
   write (error_unit, *) 'Mpi has to be defined for tempering runs'
   call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
#endif
#if defined(PARALLEL_PARAMS) && !defined(TEMPERING)
   write (error_unit, *) 'TEMPERING has to be defined for PARALLEL_PARAMS'
   call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
#endif

#ifdef MPI
   color = irank/mpi_per_parameter_set
   key = 0
   call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, Group_comm, ierr)
   call MPI_Comm_rank(Group_Comm, Irank_g, ierr)
   call MPI_Comm_size(Group_Comm, Isize_g, ierr)
   igroup = irank/isize_g
   !Write(6,*) 'irank, Irank_g, Isize_g', irank, irank_g, isize_g
   !read environment variable called ALF_SHM_CHUNK_SIZE_GB
   !it should be a positive integer setting the chunk size of shared memory blocks in GB
   !if it is not set, or set to a non-positive (including 0) integer, the routine defaults back to the
   !usual Fortran allocation routines
   call get_environment_variable(Name, value=chunk_size_str, STATUS=ierr)
   if (ierr == 0) then
      read (chunk_size_str, *, IOSTAT=ierr) chunk_size_gb
   end if
   if (ierr /= 0 .or. chunk_size_gb < 0) then
      chunk_size_gb = 0
   end if
   call mpi_shared_memory_init(Group_Comm, chunk_size_gb)
#endif
   !Initialize entanglement pairs of MPI jobs
   !This routine can and should also be called if MPI is not activated
   !It will then deactivate the entanglement measurements, i.e., the user does not have to care about this
   call Init_Entanglement_replicas(Group_Comm)

#ifdef MPI
#ifdef PARALLEL_PARAMS
   MPI_COMM_i = Group_Comm
   if (irank_g == 0) then
      write (file_para, '(A,I0,A)') "Temp_", igroup, "/parameters"
#else
      MPI_COMM_i = MPI_COMM_WORLD
      if (Irank == 0) then
         file_para = "parameters"
#endif
#else
         file_para = "parameters"
#endif
         ! This is a set of variables that  identical for each simulation.
         Nwrap = 0; NSweep = 0; NBin = 0; Ltau = 0; LOBS_EN = 0; LOBS_ST = 0; CPU_MAX = 0.d0
         Propose_S0 = .false.; Global_moves = .false.; N_Global = 0
         Global_tau_moves = .false.; sequential = .true.; Langevin = .false.; HMC = .false.
         Delta_t_Langevin_HMC = 0.d0; Max_Force = 0.d0; Leapfrog_steps = 0; N_HMC_sweeps = 1
         Nt_sequential_start = 1; Nt_sequential_end = 0; N_Global_tau = 0; Amplitude = 1.d0
         open (UNIT=5, FILE=file_para, STATUS='old', ACTION='read', IOSTAT=ierr)
         if (ierr /= 0) then
            write (error_unit, *) 'main: unable to open <parameters>', file_para, ierr
            call Terminate_on_error(ERROR_FILE_NOT_FOUND, __FILE__, __LINE__)
         end if
         read (5, NML=VAR_QMC)
         rewind (5)
         read (5, NML=VAR_HAM_NAME)
         close (5)
         NBin_eff = NBin
#ifdef MPI
      end if
      call MPI_BCAST(Nwrap, 1, MPI_INTEGER, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(NSweep, 1, MPI_INTEGER, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(NBin, 1, MPI_INTEGER, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(Ltau, 1, MPI_INTEGER, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(LOBS_EN, 1, MPI_INTEGER, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(LOBS_ST, 1, MPI_INTEGER, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(CPU_MAX, 1, MPI_REAL8, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(Propose_S0, 1, MPI_LOGICAL, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(Global_moves, 1, MPI_LOGICAL, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(N_Global, 1, MPI_Integer, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(Global_tau_moves, 1, MPI_LOGICAL, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(Nt_sequential_start, 1, MPI_Integer, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(Nt_sequential_end, 1, MPI_Integer, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(N_Global_tau, 1, MPI_Integer, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(sequential, 1, MPI_LOGICAL, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(Langevin, 1, MPI_LOGICAL, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(HMC, 1, MPI_LOGICAL, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(Leapfrog_steps, 1, MPI_Integer, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(N_HMC_sweeps, 1, MPI_Integer, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(Max_Force, 1, MPI_REAL8, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(Delta_t_Langevin_HMC, 1, MPI_REAL8, 0, MPI_COMM_i, ierr)
      call MPI_BCAST(Amplitude, 1, MPI_REAL8, 0, MPI_COMM_i, ierr)

      call MPI_BCAST(ham_name, 64, MPI_CHARACTER, 0, MPI_COMM_i, ierr)
#endif
      call Fields_init(Amplitude)
      call Alloc_Ham(ham_name)
      leap_frog_bulk = .false.
      call ham%Ham_set()
      ! Test  if  user  has  specified  correct  array  size  for time dependent Hamiltonians
      N_op = size(OP_V, 1)
      do n = 1, N_op
         do i = 1, N_Fl
            if (Op_V(n, i)%get_g_t_alloc()) then
               if (size(Op_V(n, i)%g_t, 1) /= Ltrot) then
                  write (error_unit, *) "Array size of time-dependent coupling Op_V%g_t has to be Ltrot!"
                  call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
               end if
            end if
         end do
      end do
      do n = 1, size(Op_T, 1)
         do i = 1, N_Fl
            if (Op_T(n, i)%get_g_t_alloc()) then
               if (size(Op_T(n, i)%g_t, 1) /= Ltrot) then
                  write (error_unit, *) "Array size of time-dependent coupling Op_T%g_t has to be Ltrot!"
                  call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
               end if
            end if
         end do
      end do

      ! Test if user has initialized Calc_FL array
      if (.not. allocated(Calc_Fl)) then
         allocate (Calc_Fl(N_FL))
         Calc_Fl = .true.
      end if
      ! Count number of flavors to be calculated
      N_FL_eff = 0
      do I = 1, N_Fl
         if (Calc_Fl(I)) N_FL_eff = N_FL_eff + 1
      end do
      reconstruction_needed = .false.
      if (N_FL_eff /= N_FL) reconstruction_needed = .true.
      !initialize the flavor map
      allocate (Calc_Fl_map(N_FL_eff), Phase_array(N_FL))
      N_FL_eff = 0
      do I = 1, N_Fl
         if (Calc_Fl(I)) then
            N_FL_eff = N_FL_eff + 1
            Calc_Fl_map(N_FL_eff) = I
         end if
      end do

      if (Projector) then
         if (.not. allocated(WF_R) .or. .not. allocated(WF_L)) then
            write (error_unit, *) "Projector is selected but there are no trial wave functions!"
            call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
         end if
         do nf_eff = 1, N_fl_eff
            nf = Calc_Fl_map(nf_eff)
            if (.not. allocated(WF_R(nf)%P) .or. .not. allocated(WF_L(nf)%P)) then
               write (error_unit, *) "Projector is selected but there are no trial wave functions!"
               call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
            end if
         end do
      end if
      !  Default values of  measuring interval.
      if (Projector) then
         if (LOBS_ST == 0) then
            LOBS_ST = Thtrot + 1
         else
            if (LOBS_ST < Thtrot + 1) then
               write (error_unit, *) 'Measuring out of dedicating interval, LOBS_ST too small.'
               call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
            end if
         end if
         if (LOBS_EN == 0) then
            LOBS_EN = Ltrot - Thtrot
         else
            if (LOBS_EN > Ltrot - Thtrot) then
               write (error_unit, *) 'Measuring out of dedicating interval, LOBS_EN too big.'
               call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
            end if
         end if
      else
         if (LOBS_ST == 0) then
            LOBS_ST = 1
         end if
         if (LOBS_EN == 0) then
            LOBS_EN = Ltrot
         end if
      end if
      if (.not. Global_tau_moves) then
         ! This  corresponds to the default updating scheme
         Nt_sequential_start = 1
         Nt_sequential_end = size(OP_V, 1)
         N_Global_tau = 0
      else
         !  Gives the possibility to set parameters in the Hamiltonian file
         call ham%Overide_global_tau_sampling_parameters(Nt_sequential_start, Nt_sequential_end, N_Global_tau)
      end if

      call nsigma%make(N_op, Ltrot)
      do n = 1, N_op
         nsigma%t(n) = OP_V(n, 1)%type
         nsigma%flip_protocol(n) = OP_V(n, 1)%flip_protocol
      end do

      File_seeds = "seeds"
      call Set_Random_number_Generator(File_seeds, Seed_in)
      !Write(6,*) Seed_in

      call ham%Hamiltonian_set_nsigma(Initial_field)
      if (allocated(Initial_field)) then
         call nsigma%in(Group_Comm, Initial_field)
         deallocate (Initial_field)
      else
         call nsigma%in(Group_Comm)
      end if
      call Hop_mod_init

      if (abs(CPU_MAX) > Zero) NBIN = 10000000
      if (N_Global_tau > 0) then
         call Wrapgr_alloc
      end if

#if defined(HDF5)
#if defined(TEMPERING)
      write (file_dat, '(A,I0,A)') "Temp_", igroup, "/data.h5"
#else
      file_dat = "data.h5"
#endif
#if defined(MPI)
      if (Irank_g == 0) then
#endif
         call h5open_f(ierr)
         inquire (file=file_dat, exist=file_exists)
         if (.not. file_exists) then
            ! Create HDF5 file
            call h5fcreate_f(file_dat, H5F_ACC_TRUNC_F, file_id, ierr)
            call h5ltset_attribute_string_f(file_id, '/', 'program_name', 'ALF', ierr)
            call h5fclose_f(file_id, ierr)
         end if
         call ham%write_parameters_hdf5(file_dat)

#if defined(MPI)
      end if
#endif
#endif

      call control_init(Group_Comm)
      call ham%Alloc_obs(Ltau)

      if (mod(Ltrot, nwrap) == 0) then
         Nstm = Ltrot/nwrap
      else
         nstm = Ltrot/nwrap + 1
      end if
      allocate (Stab_nt(0:Nstm))
      Stab_nt(0) = 0
      do n = 1, Nstm - 1
         Stab_nt(n) = nwrap*n
      end do

      Stab_nt(Nstm) = Ltrot

      !   Sequential = .true.
      !TODO: check if sequential is done if some fields are discrete (Warning or error termination?)
      if (Langevin .or. HMC) then
         if (Langevin) then
#if defined(MPI)
            if (Irank_g == 0) then
#endif
               if (sequential) then
                  write (output_unit, *) "Langevin mode does not allow sequential updates."
                  write (output_unit, *) "Overriding Sequential=.True. from parameter files."
               end if
               if (HMC) then
                  write (output_unit, *) "Langevin mode does not allow HMC updates."
                  write (output_unit, *) "Overriding HMC=.True. from parameter files."
               end if
               if (Global_moves) then
                  write (output_unit, *) "Langevin mode does not allow global updates."
                  write (output_unit, *) "Overriding Global_moves=.True. from parameter files."
               end if
               if (Global_tau_moves) then
                  write (output_unit, *) "Langevin mode does not allow global tau updates."
                  write (output_unit, *) "Overriding Global_tau_moves=.True. from parameter files."
               end if
#if defined(TEMPERING)
               if (N_exchange_steps > 0) then
                  write (output_unit, *) "Langevin mode does not allow tempering updates."
                  write (output_unit, *) "Overwriting N_exchange_steps to 0."
               end if
#endif
#if defined(MPI)
            end if
#endif
            Sequential = .false.
            HMC = .false.
            Global_moves = .false.
            Global_tau_moves = .false.
#if defined(TEMPERING)
            N_exchange_steps = 0
#endif
         end if
         call Langevin_HMC%make(Langevin, HMC, Delta_t_Langevin_HMC, Max_Force, Leapfrog_steps)
      else
         call Langevin_HMC%set_Update_scheme(Langevin, HMC)
      end if

      if (.not. Sequential .and. Global_tau_moves) then
         write (output_unit, *) "Warning: Sequential = .False. and Global_tau_moves = .True."
         write (output_unit, *) "in the parameter file. Global tau updates will not occur if"
         write (output_unit, *) "Sequential is set to .False. ."
      end if

      if (.not. Sequential .and. .not. HMC .and. .not. Langevin .and. .not. Global_moves) then
         write (output_unit, *) "Warning: no updates will occur as Sequential, HMC, Langevin, and"
         write (output_unit, *) "Global_moves are all .False. in the parameter file."
      end if

      if (Sequential .and. Nt_sequential_end < Nt_sequential_start) then
         write (output_unit, *) "Warning: Nt_sequential_end is smaller than Nt_sequential_start"
      end if

#if defined(TEMPERING)
      write (file_info, '(A,I0,A)') "Temp_", igroup, "/info"
#else
      file_info = "info"
#endif

#if defined(MPI)
      if (Irank_g == 0) then
#endif
         open (Unit=50, file=file_info, status="unknown", position="append")
         write (50, *) 'Sweeps                              : ', Nsweep
         if (abs(CPU_MAX) < ZERO) then
            write (50, *) 'Bins                                : ', NBin
            write (50, *) 'No CPU-time limitation '
         else
            write (50, '(" Prog will stop after hours:",2x,F8.4)') CPU_MAX
         end if
         write (50, *) 'Measure Int.                        : ', LOBS_ST, LOBS_EN
         write (50, *) 'Stabilization,Wrap                  : ', Nwrap
         write (50, *) 'Nstm                                : ', NSTM
         write (50, *) 'Ltau                                : ', Ltau
         write (50, *) '# of interacting Ops per time slice : ', size(OP_V, 1)
         if (Propose_S0) &
              &  write (50, *) 'Propose Ising moves according to  bare Ising action'
         if (Global_moves) then
            write (50, *) 'Global moves are enabled   '
            write (50, *) '# of global moves / sweep :', N_Global
         end if
         if (sequential) then
            if (Global_tau_moves) then
               write (50, *) 'Nt_sequential_start: ', Nt_sequential_start
               write (50, *) 'Nt_sequential_end  : ', Nt_sequential_end
               write (50, *) 'N_Global_tau       : ', N_Global_tau
            else
               write (50, *) 'Default sequential updating '
            end if
         end if
         if (Langevin) then
            write (50, *) 'Langevin del_t: ', Delta_t_Langevin_HMC
            write (50, *) 'Max Force     : ', Max_Force
         end if
         if (HMC) then
            write (50, *) 'HMC del_t     : ', Delta_t_Langevin_HMC
            write (50, *) 'Leapfrog_Steps: ', Leapfrog_Steps
            write (50, *) 'HMC_Sweeps:     ', N_HMC_sweeps
         end if

         !Write out info  for  amplitude and flip_protocol
         Toggle = .false.
         do n = 1, N_op
            if (nsigma%t(n) == 3 .or. nsigma%t(n) == 4) Toggle = .true.
         end do
         if (Toggle) then
            write (50, *) 'Amplitude  for  t=3,4  vertices is  set to: ', Amplitude
         end if
         Toggle = .false.
         do n = 1, N_op
            if (nsigma%t(n) == 4) Toggle = .true.
         end do
         if (Toggle) then
            write (50, "('Flip protocal  for  t=4  vertices is  set to')", advance="no")
            do n1 = 1, 4
               Toggle1 = .false.
               do n = 1, N_op
                  if (nsigma%Flip_protocol(n) == n1 .and. nsigma%t(n) == 4) Toggle1 = .true.
               end do
               if (Toggle1) write (50, "(I2,2x)", advance="no") n1
            end do
            write (50, *)
         end if

#if defined(MPI)
         write (50, *) 'Number of mpi-processes : ', isize_g
         if (use_mpi_shm) write (50, *) 'Using mpi-shared memory in chunks of ', chunk_size_gb, 'GB.'
#endif
#if defined(GIT)
         write (50, *) 'This executable represents commit '&
              &      , GIT_COMMIT_HASH, ' of branch ', GIT_BRANCH, '.'
#endif
#if defined(STAB1)
         write (50, *) 'STAB1 is defined '
#endif
#if defined(STAB2)
         write (50, *) 'STAB2 is defined '
#endif
#if defined(STAB3)
         write (50, *) 'STAB3 is defined '
#endif
#if defined(STABLOG)
         write (50, *) 'LOG is defined '
#endif
#if defined(QRREF)
         write (50, *) 'QRREF is defined '
#endif
#if defined(TEMPERING) && !defined(PARALLEL_PARAMS)
         write (50, *) '# of exchange steps  ', N_exchange_steps
         write (50, *) 'Tempering frequency  ', N_Tempering_frequency
         write (50, *) 'Tempering Calc_det   ', Tempering_calc_det
#endif
         close (50)
#if defined(MPI)
      end if
#endif

      !Call Test_Hamiltonian
      allocate (Test(Ndim, Ndim), GR(NDIM, NDIM, N_FL), GR_Tilde(NDIM, NDIM, N_FL))
      allocate (udvl(N_FL_eff), udvr(N_FL_eff), udvst(NSTM, N_FL_eff))
      do nf_eff = 1, N_FL_eff
         nf = Calc_Fl_map(nf_eff)
         do n = 1, NSTM
            call udvst(n, nf_eff)%alloc(ndim)
         end do
         if (Projector) then
            call udvl(nf_eff)%init(ndim, 'l', WF_L(nf)%P)
            call udvr(nf_eff)%init(ndim, 'r', WF_R(nf)%P)
            call udvst(NSTM, nf_eff)%reset('l', WF_L(nf)%P)
         else
            call udvl(nf_eff)%init(ndim, 'l')
            call udvr(nf_eff)%init(ndim, 'r')
            call udvst(NSTM, nf_eff)%reset('l')
         end if
      end do

      do NST = NSTM - 1, 1, -1
         NT1 = Stab_nt(NST + 1)
         NT = Stab_nt(NST)
         !Write(6,*)'Hi', NT1,NT, NST
         call WRAPUL(NT1, NT, UDVL)
         do nf_eff = 1, N_FL_eff
            UDVST(NST, nf_eff) = UDVL(nf_eff)
         end do
      end do
      NT1 = stab_nt(1)
      call WRAPUL(NT1, 0, UDVL)

      NVAR = 1
      Phase_array = cmplx(1.d0, 0.d0, kind(0.d0))
      do nf_eff = 1, N_Fl_eff
         nf = Calc_Fl_map(nf_eff)
         call CGR(Z, NVAR, GR(:, :, nf), UDVR(nf_eff), UDVL(nf_eff))
         call Op_phase(Z, OP_V, Nsigma, nf)
         Phase_array(nf) = Z
      end do
      if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
      Phase = product(Phase_array)
      Phase = Phase**N_SUN
#ifdef MPI
      !WRITE(6,*) 'Phase is: ', Irank, PHASE, GR(1,1,1)
#else
      !WRITE(6,*) 'Phase is: ',  PHASE
#endif

      call Control_init(Group_Comm)

      do NBC = 1, NBIN
         ! Here, you have the green functions on time slice 1.
         ! Set bin observables to zero.

         call system_clock(count_bin_start)

         call ham%Init_obs(Ltau)
#if defined(TEMPERING)
         call Global_Tempering_init_obs
#endif

         do NSW = 1, NSWEEP

#if defined(TEMPERING) && !defined(PARALLEL_PARAMS)
            if (mod(NSW, N_Tempering_frequency) == 0) then
               !Write(6,*) "Irank, Call tempering", Irank, NSW, N_exchange_steps
               call Exchange_Step(Phase, GR, udvr, udvl, Stab_nt, udvst, N_exchange_steps, Tempering_calc_det)
            end if
#endif
            ! Global updates
            if (Global_moves) call Global_Updates(Phase, GR, udvr, udvl, Stab_nt, udvst, N_Global)

            if (str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN") then
               !  Carry out a Langevin update and calculate equal time observables.
               call Langevin_HMC%update(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, &
                    &                   LOBS_ST, LOBS_EN, LTAU)

               if (LTAU == 1) then
                  if (Projector) then
                     NST = 0
                     call Tau_p(udvl, udvr, udvst, GR, PHASE, NSTM, STAB_NT, NST, LOBS_ST, LOBS_EN)
                     call Langevin_HMC%set_L_Forces(.true.)
                  else
                     call Tau_m(udvst, GR, PHASE, NSTM, NWRAP, STAB_NT, LOBS_ST, LOBS_EN)
                     call Langevin_HMC%set_L_Forces(.true.)
                  end if
               end if
            end if

            if (str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC") then
               if (Sequential) call Langevin_HMC%set_L_Forces(.false.)
               do n = 1, N_HMC_sweeps
                  !  Carry out a Langevin update and calculate equal time observables.
                  call Langevin_HMC%update(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, &
                       &                   LOBS_ST, LOBS_EN, LTAU)
                  if (n /= N_HMC_sweeps) then
                     call Langevin_HMC%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
                          &  LOBS_ST, LOBS_EN, .true.)
                     call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
                     call Langevin_HMC%set_L_Forces(.true.)
                  end if
               end do

               !Do time-displaced measurements if needed, else set Calc_Obser_eq=.True. for the very first leapfrog ONLY
               if (.not. sequential) then
                  if (LTAU == 1) then
                     if (Projector) then
                        NST = 0
                        call Tau_p(udvl, udvr, udvst, GR, PHASE, NSTM, STAB_NT, NST, LOBS_ST, LOBS_EN)
                     else
                        call Tau_m(udvst, GR, PHASE, NSTM, NWRAP, STAB_NT, LOBS_ST, LOBS_EN)
                     end if
                  else
                     call Langevin_HMC%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
                     &  LOBS_ST, LOBS_EN, .true.)
                     call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
                  end if
                  call Langevin_HMC%set_L_Forces(.true.)
               end if
            end if

            if (Sequential) then
               ! Propagation from 1 to Ltrot
               ! Set the right storage to 1
               do nf_eff = 1, N_FL_eff
                  nf = Calc_Fl_map(nf_eff)
                  if (Projector) then
                     call udvr(nf_eff)%reset('r', WF_R(nf)%P)
                  else
                     call udvr(nf_eff)%reset('r')
                  end if
               end do

               NST = 1
               do NTAU = 0, LTROT - 1
                  NTAU1 = NTAU + 1
                  call WRAPGRUP(GR, NTAU, PHASE, Propose_S0, Nt_sequential_start, Nt_sequential_end, N_Global_tau)

                  if (NTAU1 == Stab_nt(NST)) then
                     NT1 = Stab_nt(NST - 1)
                     call WRAPUR(NT1, NTAU1, udvr)
                     Phase_array = cmplx(1.d0, 0.d0, kind(0.d0))
                     do nf_eff = 1, N_FL_eff
                        nf = Calc_Fl_map(nf_eff)
                        ! Read from storage left propagation from LTROT to  NTAU1
                        udvl(nf_eff) = udvst(NST, nf_eff)
                        ! Write in storage right prop from 1 to NTAU1
                        udvst(NST, nf_eff) = udvr(nf_eff)
                        NVAR = 1
                        if (NTAU1 .gt. LTROT/2) NVAR = 2
                        TEST(:, :) = GR(:, :, nf)
                        call CGR(Z1, NVAR, GR(:, :, nf), UDVR(nf_eff), UDVL(nf_eff))
                        call Control_PrecisionG(GR(:, :, nf), Test, Ndim)
                        call Op_phase(Z1, OP_V, Nsigma, nf)
                        Phase_array(nf) = Z1
                     end do
                     if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
                     Z = product(Phase_array)
                     Z = Z**N_SUN
                     call Control_PrecisionP(Z, Phase)
                     Phase = Z
                     NST = NST + 1
                  end if

                  if (NTAU1 .ge. LOBS_ST .and. NTAU1 .le. LOBS_EN) then
                     !Call  Global_tau_mod_Test(Gr,ntau1)
                     !Stop
                     !write(*,*) "GR before obser sum: ",sum(GR(:,:,1))
                     !write(*,*) "Phase before obser : ",phase
                     Mc_step_weight = 1.d0
                     if (Symm) then
                        call Hop_mod_Symm(GR_Tilde, GR, ntau1)
                        !reconstruction of NOT calculated block!!!
                        if (reconstruction_needed) call ham%GR_reconstruction(GR_Tilde)
                        call ham%Obser(GR_Tilde, PHASE, Ntau1, Mc_step_weight)
                     else
                        !reconstruction of NOT calculated block!!!
                        if (reconstruction_needed) call ham%GR_reconstruction(GR)
                        call ham%Obser(GR, PHASE, Ntau1, Mc_step_weight)
                     end if
                  end if
               end do

               do nf_eff = 1, N_FL_eff
                  nf = Calc_Fl_map(nf_eff)
                  if (Projector) then
                     call udvl(nf_eff)%reset('l', WF_L(nf)%P)
                  else
                     call udvl(nf_eff)%reset('l')
                  end if
               end do

               NST = NSTM - 1
               do NTAU = LTROT, 1, -1
                  NTAU1 = NTAU - 1
                  call WRAPGRDO(GR, NTAU, PHASE, Propose_S0, Nt_sequential_start, Nt_sequential_end, N_Global_tau)
                  if (NTAU1 .ge. LOBS_ST .and. NTAU1 .le. LOBS_EN) then
                     !write(*,*) "GR before obser sum: ",sum(GR(:,:,1))
                     !write(*,*) "Phase before obser : ",phase
                     Mc_step_weight = 1.d0
                     if (Symm) then
                        call Hop_mod_Symm(GR_Tilde, GR, ntau1)
                        !reconstruction of NOT calculated block!!!
                        if (reconstruction_needed) call ham%GR_reconstruction(GR_Tilde)
                        call ham%Obser(GR_Tilde, PHASE, Ntau1, Mc_step_weight)
                     else
                        !reconstruction of NOT calculated block!!!
                        if (reconstruction_needed) call ham%GR_reconstruction(GR)
                        call ham%Obser(GR, PHASE, Ntau1, Mc_step_weight)
                     end if
                  end if
                  if (Stab_nt(NST) == NTAU1 .and. NTAU1 .ne. 0) then
                     NT1 = Stab_nt(NST + 1)
                     !Write(6,*) 'Wrapul : ', NT1, NTAU1
                     call WRAPUL(NT1, NTAU1, udvl)
                     !Write(6,*)  'Write UL, read UR ', NTAU1, NST
                     Phase_array = cmplx(1.d0, 0.d0, kind(0.d0))
                     do nf_eff = 1, N_FL_eff
                        nf = Calc_Fl_map(nf_eff)
                        ! Read from store the right prop. from 1 to LTROT/NWRAP-1
                        udvr(nf_eff) = udvst(NST, nf_eff)
                        ! WRITE in store the left prop. from LTROT/NWRAP-1 to 1
                        udvst(NST, nf_eff) = udvl(nf_eff)
                        NVAR = 1
                        if (NTAU1 .gt. LTROT/2) NVAR = 2
                        TEST(:, :) = GR(:, :, nf)
                        call CGR(Z1, NVAR, GR(:, :, nf), UDVR(nf_eff), UDVL(nf_eff))
                        call Control_PrecisionG(GR(:, :, nf), Test, Ndim)
                        call Op_phase(Z1, OP_V, Nsigma, nf)
                        Phase_array(nf) = Z1
                     end do
                     if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
                     Z = product(Phase_array)
                     Z = Z**N_SUN
                     call Control_PrecisionP(Z, Phase)
                     Phase = Z
                     if (LTAU == 1 .and. Projector .and. Stab_nt(NST) <= THTROT + 1 .and. THTROT + 1 < Stab_nt(NST + 1)) then
                        call tau_p(udvl, udvr, udvst, GR, PHASE, NSTM, STAB_NT, NST, LOBS_ST, LOBS_EN)
                     end if
                     NST = NST - 1
                  end if
               end do

               !Calculate and compare green functions on time slice 0.
               NT1 = Stab_nt(0)
               NT = Stab_nt(1)
               call WRAPUL(NT, NT1, udvl)

               do nf_eff = 1, N_FL_eff
                  nf = Calc_Fl_map(nf_eff)
                  if (Projector) then
                     call udvr(nf_eff)%reset('r', WF_R(nf)%P)
                  else
                     call udvr(nf_eff)%reset('r')
                  end if
               end do
               Phase_array = cmplx(1.d0, 0.d0, kind(0.d0))
               do nf_eff = 1, N_FL_eff
                  nf = Calc_Fl_map(nf_eff)
                  TEST(:, :) = GR(:, :, nf)
                  NVAR = 1
                  call CGR(Z1, NVAR, GR(:, :, nf), UDVR(nf_eff), UDVL(nf_eff))
                  call Control_PrecisionG(GR(:, :, nf), Test, Ndim)
                  call Op_phase(Z1, OP_V, Nsigma, nf)
                  Phase_array(nf) = Z1
               end do
               if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
               Z = product(Phase_array)
               Z = Z**N_SUN
               call Control_PrecisionP(Z, Phase)
               Phase = Z
               NST = NSTM
               do nf_eff = 1, N_FL_eff
                  nf = Calc_Fl_map(nf_eff)
                  if (Projector) then
                     call udvst(NST, nf_eff)%reset('l', WF_L(nf)%P)
                  else
                     call udvst(NST, nf_eff)%reset('l')
                  end if
               end do

               if (LTAU == 1 .and. .not. Projector) then
                  call TAU_M(udvst, GR, PHASE, NSTM, NWRAP, STAB_NT, LOBS_ST, LOBS_EN)
               end if
            end if

         end do
         call ham%Pr_obs(Ltau)
#if defined(TEMPERING)
         call Global_Tempering_Pr
#endif

         call nsigma%out(Group_Comm)

         call system_clock(count_bin_end)
         prog_truncation = .false.
         if (abs(CPU_MAX) > Zero) then
            call make_truncation(prog_truncation, CPU_MAX, count_bin_start, count_bin_end, group_comm)
         end if
         if (prog_truncation) then
            Nbin_eff = nbc
            exit !exit the loop over the bin index, labelled NBC.
         end if
      end do

      ! Deallocate things
      do nf_eff = 1, N_FL_eff
         call udvl(nf_eff)%dealloc
         call udvr(nf_eff)%dealloc
         do n = 1, NSTM
            call udvst(n, nf_eff)%dealloc
         end do
      end do
      if (Projector) then
         do nf = 1, N_FL
            call WF_clear(WF_R(nf))
            call WF_clear(WF_L(nf))
         end do
      end if
      deallocate (udvl, udvr, udvst)
      deallocate (GR, TEST, Stab_nt, GR_Tilde)
      if (Projector) deallocate (WF_R, WF_L)
      if (N_Global_tau > 0) then
         call Wrapgr_dealloc
      end if
      do nf = 1, N_FL
         do n = 1, size(OP_V, 1)
            call Op_clear(Op_V(n, nf), Op_V(n, nf)%N)
         end do
         do n = 1, size(OP_T, 1)
            call Op_clear(Op_T(n, nf), Op_T(n, nf)%N)
         end do
      end do

#if defined(MPI)
      ! Gracefully deallocate all shared MPI memory (thw whole chunks)
      ! irrespective of where they actually have been used
      call deallocate_all_shared_memory
#endif

      call Control_Print(Group_Comm, Langevin_HMC%get_Update_scheme())

#if defined(MPI)
      if (Irank_g == 0) then
#endif
         if (abs(CPU_MAX) > Zero) then
#if defined(TEMPERING)
            write (file_info, '(A,I0,A)') "Temp_", igroup, "/info"
#else
            file_info = "info"
#endif
            open (Unit=50, file=file_info, status="unknown", position="append")
            write (50, *) ' Effective number of bins   : ', Nbin_eff
            close (50)
         end if
#if defined(MPI)
      end if
#endif

      call Langevin_HMC%clean()
      deallocate (Calc_Fl_map, Phase_array)

#ifdef MPI
      call MPI_FINALIZE(ierr)
#endif

      end program Main
