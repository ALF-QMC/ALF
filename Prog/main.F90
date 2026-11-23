!  Copyright (C) 2016 - 2026 The ALF project
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


Program Main

#ifdef MPI
        Use mpi
#endif
        Use runtime_error_mod, only: Terminate_on_error, ERROR_GENERIC, ERROR_HAMILTONIAN, &
             &                       ERROR_RUNNING_FILE_FOUND
        Use mpi_shared_memory, only: mpi_shared_memory_init, deallocate_all_shared_memory, use_mpi_shm
        Use QMC_runtime_var, only: check_compatibility_reconstruct_greens_function, &
             &                     check_langevin_schemes_and_variables, check_mala_variables, &
             &                     check_update_schemes_compatibility, get_amplitude, get_cpu_max, &
             &                     get_delta_t_langevin_hmc, get_delta_t_mala_global, &
             &                     get_delta_t_mala_global_tau, get_delta_t_mala_sequential, &
             &                     get_global_mala_moves, get_global_moves, get_global_tau_mala_moves, &
             &                     get_global_tau_moves, get_ham_name, get_hmc, get_langevin, &
             &                     get_leapfrog_steps, get_lobs_en, get_lobs_st, get_ltau, get_max_force, &
             &                     get_max_force_mala_global, get_max_force_mala_global_tau, &
             &                     get_max_force_mala_sequential, get_mpi_per_parameter_set, get_n_global, &
             &                     get_n_global_mala_sweeps, get_n_global_tau, get_n_global_tau_mala, &
             &                     get_n_hmc_sweeps, get_nbin, get_nsweep, get_nt_sequential_end, &
             &                     get_nt_sequential_start, get_nwrap, get_propose_s0, get_sequential, &
             &                     get_sequential_mala, get_tempering_calc_det, &
             &                     read_and_broadcast_qmc_var_and_ham_name, &
             &                     set_default_values_measuring_interval, apply_Langevin_overrides, &
             &                     set_default_global_tau_sampling, &
             &                     set_mpi_per_parameter_set, set_n_global_tau, &
             &                     set_n_global_tau_mala, apply_CPU_MAX_policy, set_nt_sequential_end, &
             &                     set_nt_sequential_start
        Use Hamiltonian_main, only: ham, Alloc_Ham, Validate_Ham_Variables, Setup_flavor_map, Clear_Ham_state, &
             &                      Set_Group_Comm, &
             &                      Ndim, N_FL, N_FL_eff, N_SUN, &
             &                      Ltrot, Thtrot, Op_T, Op_V, nsigma, WF_L, WF_R, Calc_Fl, Calc_Fl_map, &
             &                      Projector, Group_Comm, leap_frog_bulk
        Use files_mod, only: str_to_upper
        Use Natural_Constants, only: Eps_small
        Use Control, only: Control_init, Control_Print, &
             &             make_truncation
        Use Hop_mod, only: Hop_mod_init
        Use Global_mod, only: Global_Updates, Global_Tempering_setup, Global_Tempering_init_obs, &
             &                Global_Tempering_Pr
#if defined(TEMPERING)
        Use Global_mod, only: Exchange_Step
        Use QMC_runtime_var, only: read_and_broadcast_TEMPERING_var, get_N_exchange_steps, &
             &                     get_N_Tempering_frequency
#endif
        Use Wrapgr_mod, only: Wrapgr_alloc, Wrapgr_dealloc
        Use Fields_mod, only: Fields_init
        use entanglement_mod, only: Init_Entanglement_replicas
        use iso_fortran_env, only: output_unit, error_unit
        Use Langevin_HMC_mod, only: Langevin_HMC, Metropolis_Langevin
        use set_random, only: Set_Random_number_Generator
        use iso8601_datetime_mod, only: iso8601_datetime
        Use qmc_setup_mod, only: Setup_model, Setup_fields, Setup_update_schemes, Write_info_parameters, Write_info_start, Write_info_end
        Use comm_mod, only: is_root, comm_world
        Use qmc_driver_mod, only: Langevin_update, HMC_update, MALA_update
        Use qmc_sweep_mod, only: Sweep_up, Sweep_down
        Use qmc_walker_mod, only: qmc_walker_t, Walker_layout, Walker_alloc, Walker_dealloc, &
             &                    Walker_init_stack, Walker_init_GR
         
#ifdef HDF5
        use hdf5
        use h5lt
#endif
        Implicit none

#include "git.h"
        ! State of the Monte Carlo walker: Green function, phase, UDV propagations and storage
        Type (qmc_walker_t) :: walker
       
        Integer :: NBC, NSW
        Character (len=64) :: file_seeds, file_info, file_git_info
        Integer :: Seed_in
        Complex (Kind=Kind(0.d0)) , allocatable, dimension(:,:) :: Initial_field

#ifdef HDF5
        Character (len=64) :: file_dat
#endif

#ifdef HDF5
        INTEGER(HID_T) :: file_id
        Logical :: file_exists
#endif
         
        !General
        Integer :: I, nf, nf_eff, n, n1, N_op, NBin_eff
        Integer :: tmp_Nt_sequential_start, tmp_Nt_sequential_end, tmp_N_Global_tau, tmp_N_Global_tau_MALA
        Logical :: Toggle,  Toggle1
        Real    (Kind=Kind(0.d0)) :: chunk_size_gb   ! MPI shared memory chunk (GB); set with MPI only

#if defined(MPI) || defined(HDF5)
        Integer :: Ierr
#endif



        ! For the truncation of the program:
        logical                   :: prog_truncation, run_file_exists
        integer (kind=kind(0.d0)) :: count_bin_start, count_bin_end
        
        ! For MPI shared memory
#ifdef MPI
        character(64), parameter :: name="ALF_SHM_CHUNK_SIZE_GB"
        character(64) :: chunk_size_str
        Integer        :: Isize, Irank, Isize_g, color, key, igroup, comm_split

        CALL MPI_INIT(ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
        
        If (is_root(comm_world())) then
           write (*,*) "ALF Copyright (C) 2016 - 2026 The ALF project contributors"
           write (*,*) "This Program comes with ABSOLUTELY NO WARRANTY; for details see license.GPL"
           write (*,*) "This is free software, and you are welcome to redistribute it under certain conditions."

           ! Ensure that only one ALF is running at the same time, i.e. the file RUNNING is not present
           inquire (file='RUNNING', exist=run_file_exists)
           if (run_file_exists) then
             write (error_unit,*)
             write (error_unit,*) "ALF is already running or the previous run failed."
             write (error_unit,*) "Please ensure the following:"
             write (error_unit,*) " * Make sure no other simulation is currently running in this directory"
             write (error_unit,*) "   (Wait until the previous run is finished; it will automatically remove RUNNING)"
             write (error_unit,*) " * If the previous run crashed, make sure that"
             write (error_unit,*) "    1) the data files are not corrupted"
             write (error_unit,*) "       (run the analysis)"
             write (error_unit,*) "    2) the configuration files are not corrupted"
             write (error_unit,*) "       (e.g., h5dump confout_*.h5 or check number of lines in confout_*)"
             write (error_unit,*) "    3) If either data or configuration file are currupted (rare event), either"
             write (error_unit,*) "       * [PREFERED] remove them and start fresh (safe)"
             write (error_unit,*) "       * repair them (if you know what you are doing)"
             write (error_unit,*) "         (difficult or impossible; ensure data and configuration files synced)"
             write (error_unit,*) "    4) remove the file RUNNING manually before resubmition"
             write (error_unit,*) "Afterwards, you may rerun the simulation."
#ifdef MPI
             call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#else
   CALL Terminate_on_error(ERROR_RUNNING_FILE_FOUND,__FILE__,__LINE__)
#endif
           else
             open (unit=5, file='RUNNING', status='replace', action='write')
             write (5,*) "ALF is running"
             close (5)
           end if
        endif

#if defined(TEMPERING) && defined(MPI)
        call read_and_broadcast_TEMPERING_var()
        if ( mod(ISIZE,get_mpi_per_parameter_set()) .ne. 0 ) then
           Write (error_unit,*) "mpi_per_parameter_set is not a multiple of total mpi processes"
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        Call Global_Tempering_setup
#elif !defined(TEMPERING)  && defined(MPI)
        call set_mpi_per_parameter_set(Isize)
#elif defined(TEMPERING)  && !defined(MPI)
        Write(error_unit,*) 'Mpi has to be defined for tempering runs'
        CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
#if defined(PARALLEL_PARAMS) && !defined(TEMPERING)
        Write(error_unit,*) 'TEMPERING has to be defined for PARALLEL_PARAMS'
        CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif

#ifdef MPI
        color = irank/get_mpi_per_parameter_set()
        key   =  0
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,comm_split, ierr)
        call Set_Group_Comm(comm_split)
        call MPI_Comm_size(Group_Comm, Isize_g, ierr)
        igroup           = irank/isize_g
        !Write(6,*) 'irank, Irank_g, Isize_g', irank, irank_g, isize_g
        !read environment variable called ALF_SHM_CHUNK_SIZE_GB
        !it should be a positive integer setting the chunk size of shared memory blocks in GB
        !if it is not set, or set to a non-positive (including 0) integer, the routine defaults back to the
        !usual Fortran allocation routines
        CALL GET_ENVIRONMENT_VARIABLE(Name, VALUE=chunk_size_str, STATUS=ierr)
        if (ierr==0) then
           read(chunk_size_str,*,IOSTAT=ierr) chunk_size_gb
        endif
        if (ierr/=0 .or. chunk_size_gb<0) then
              chunk_size_gb=0
        endif
        CALL mpi_shared_memory_init(Group_Comm, chunk_size_gb)
#endif
        !Initialize entanglement pairs of MPI jobs
        !This routine can and should also be called if MPI is not activated
        !It will then deactivate the entanglement measurements, i.e., the user does not have to care about this
        call Init_Entanglement_replicas(Group_Comm)

#if defined(TEMPERING)
        write(file_info,'(A,I0,A)') "Temp_",igroup,"/info"
        write(file_git_info,'(A,I0,A)') "Temp_",igroup,"/info_git"
#else
        file_info = "info"
        file_git_info = "info_git"
#endif

        Call Write_info_start(file_info, file_git_info)
        call read_and_broadcast_QMC_var_and_ham_name(Group_Comm)
        NBin_eff = get_NBin()

        Call Setup_model(N_op)

        Call Setup_fields(N_op)
        
#if defined(HDF5)
#if defined(TEMPERING)
        write(file_dat,'(A,I0,A)') "Temp_",igroup,"/data.h5"
#else
        file_dat = "data.h5"
#endif
        if (is_root(Group_Comm)) then
          CALL h5open_f(ierr)
          inquire (file=file_dat, exist=file_exists)
          IF (.not. file_exists) THEN
            ! Create HDF5 file
            CALL h5fcreate_f(file_dat, H5F_ACC_TRUNC_F, file_id, ierr)
            call h5ltset_attribute_string_f(file_id, '/', 'program_name', 'ALF', ierr)
            call h5fclose_f(file_id, ierr)
          endif
          call ham%write_parameters_hdf5(file_dat)

        endif
#endif


        Call control_init(Group_Comm)
        Call ham%Alloc_obs(get_Ltau())

        Call Walker_layout(walker, Ltrot, get_Nwrap())

        Call Setup_update_schemes(N_op)

        Call Write_info_parameters(file_info, walker%NSTM, N_op, chunk_size_gb)
        
        !Call Test_Hamiltonian
        Call Walker_alloc(walker)

        ! Initial storage from Ltrot down to 0, then the Green function and phase on time slice 0
        Call Walker_init_stack(walker)
        Call Walker_init_GR(walker)
#ifdef MPI
        !WRITE(6,*) 'Phase is: ', Irank, PHASE, GR(1,1,1)
#else
        !WRITE(6,*) 'Phase is: ',  PHASE
#endif



        Call Control_init(Group_Comm)

        DO  NBC = 1, get_NBin()
           ! Here, you have the green functions on time slice 1.
           ! Set bin observables to zero.

           call system_clock(count_bin_start)

           Call ham%Init_obs(get_Ltau())
#if defined(TEMPERING)
           Call Global_Tempering_init_obs
#endif

           DO NSW = 1, get_NSweep()

#if defined(TEMPERING) && !defined(PARALLEL_PARAMS)
              IF (MOD(NSW,get_N_Tempering_frequency()) == 0) then
                 !Write(6,*) "Irank, Call tempering", Irank, NSW, get_N_exchange_steps()
                 CALL Exchange_Step(walker%Phase,walker%GR,walker%udvr, walker%udvl,walker%Stab_nt, walker%udvst, get_N_exchange_steps(), get_Tempering_calc_det())
              endif
#endif
              ! Global updates
              If ( get_Global_moves() ) Call Global_Updates(walker%Phase, walker%GR, walker%udvr, walker%udvl, walker%Stab_nt, walker%udvst, get_N_Global() )


              If ( str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" )  then
                 Call Langevin_update(walker)
              endif

              If (  str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC" )  then
                 Call HMC_update(walker)
              endif

              If (  str_to_upper(Metropolis_Langevin%get_Update_scheme()) == "MALA" )  then
                 Call MALA_update(walker, NSW)
              endif

              If (get_sequential())  then 
                 Call Sweep_up(walker)
                 
                 Call Sweep_down(walker)
              endif

           ENDDO
           Call ham%Pr_obs(get_Ltau())
#if defined(TEMPERING) && !defined(PARALLEL_PARAMS)
           Call Global_Tempering_Pr
#endif

           Call nsigma%out(Group_Comm)

           call system_clock(count_bin_end)
           prog_truncation = .false.
           if ( abs(get_CPU_MAX()) > Eps_small ) then
              Call make_truncation(prog_truncation,get_CPU_MAX(),count_bin_start,count_bin_end,group_comm)
           endif
           If (prog_truncation) then
              Nbin_eff = nbc
              exit !exit the loop over the bin index, labelled NBC.
           Endif
        Enddo

        ! Deallocate things
        Call Walker_dealloc(walker)
        If (get_N_Global_tau() > 0 .or. get_N_Global_tau_MALA() > 0 .or. get_Sequential_MALA()) then
           Call Wrapgr_dealloc
        endif
        ! Release the model state: trial wave functions, operators, flavour map
        Call Clear_Ham_state()

#if defined(MPI)  
        ! Gracefully deallocate all shared MPI memory (thw whole chunks)
        ! irrespective of where they actually have been used
        call deallocate_all_shared_memory
#endif

        Call Control_Print(Group_Comm, Langevin_HMC%get_Update_scheme(), get_Global_MALA_moves(), &
             &              get_Sequential_MALA(), get_Global_tau_MALA_moves())

        Call Write_info_end(file_info, NBin_eff)
        
        Call Langevin_HMC%clean()
        Call Metropolis_Langevin%clean()

         ! Delete the file RUNNING since the simulation finished successfully
        If (is_root(comm_world())) then
           open(unit=5, file='RUNNING', status='old')
           close(5, status='delete')
        endif

#ifdef MPI
        CALL MPI_FINALIZE(ierr)
#endif

      end Program Main
