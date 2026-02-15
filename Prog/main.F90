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


Program Main

#ifdef MPI
        Use mpi
#endif
        Use runtime_error_mod
        Use Operator_mod
        Use QMC_runtime_var
        Use Lattices_v3
        Use MyMats
        Use Hamiltonian_main
        Use Control
        Use Tau_m_mod
        Use Tau_p_mod
        Use Hop_mod
        Use Global_mod
        Use UDV_State_mod
        Use Wrapgr_mod
        Use Fields_mod
        Use WaveFunction_mod
        use entanglement_mod
        use iso_fortran_env, only: output_unit, error_unit
        Use Langevin_HMC_mod
        use wrapur_mod
        use wrapul_mod
        use cgr1_mod
        use set_random
         
#ifdef HDF5
        use hdf5
        use h5lt
#endif
        Implicit none

#include "git.h"
       COMPLEX (Kind=Kind(0.d0)), Dimension(:,:)  , Allocatable   ::  TEST
       COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable    :: GR, GR_Tilde
       CLASS(UDV_State), DIMENSION(:), ALLOCATABLE :: udvl, udvr
       COMPLEX (Kind=Kind(0.d0)), Dimension(:)  , Allocatable   :: Phase_array
       
       Integer :: NBC, NSW
       Integer :: NTAU, NTAU1
       Character (len=64) :: file_seeds, file_dat, file_info
       Integer :: Seed_in
       Complex (Kind=Kind(0.d0)) , allocatable, dimension(:,:) :: Initial_field


#ifdef HDF5
        INTEGER(HID_T) :: file_id
        Logical :: file_exists
#endif
         
        !General
        Integer :: NSTM, NT, NT1, NVAR
        Integer :: Ierr, I,nf, nf_eff, nst, n, n1, N_op, NBin_eff
        Integer :: tmp_Nt_sequential_start, tmp_Nt_sequential_end, tmp_N_Global_tau
        Logical :: Toggle,  Toggle1
        Complex (Kind=Kind(0.d0)) :: Phase, Z, Z1
        Real    (Kind=Kind(0.d0)) :: ZERO = 10D-8
        Real    (Kind=Kind(0.d0)) :: Mc_step_weight

        ! Storage for  stabilization steps
        Integer, dimension(:), allocatable :: Stab_nt 

        ! Space for storage.
        CLASS(UDV_State), Dimension(:,:), ALLOCATABLE :: udvst

        ! For the truncation of the program:
        logical                   :: prog_truncation, run_file_exists
        integer (kind=kind(0.d0)) :: count_bin_start, count_bin_end
        
        ! For MPI shared memory
        character(64), parameter :: name="ALF_SHM_CHUNK_SIZE_GB"
        character(64) :: chunk_size_str
        Real    (Kind=Kind(0.d0)) :: chunk_size_gb

#ifdef MPI
        Integer        :: Isize, Irank, Irank_g, Isize_g, color, key, igroup

        CALL MPI_INIT(ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        
        If (  Irank == 0 ) then
#endif
           write (*,*) "ALF Copyright (C) 2016 - 2022 The ALF project contributors"
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
#ifdef MPI
        endif
#endif

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
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Group_comm, ierr)
        call MPI_Comm_rank(Group_Comm, Irank_g, ierr)
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

        call read_and_broadcast_QMC_var_and_ham_name(Group_Comm)
        NBin_eff = get_NBin()

        Call Fields_init(get_Amplitude())
        Call Alloc_Ham(get_ham_name())
        leap_frog_bulk = .false.
        Call ham%Ham_set()
        ! Test  if  user  has  specified  correct  array  size  for time dependent Hamiltonians
        N_op = Size(OP_V,1)
        do n = 1, N_op
            do i = 1, N_Fl
                if (Op_V(n,i)%get_g_t_alloc()) then
                    if (size(Op_V(n,i)%g_t,1) /= Ltrot) then
                        write(error_unit,*) "Array size of time-dependent coupling Op_V%g_t has to be Ltrot!"
                        CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
                    endif
                endif
            enddo
        enddo
        do n = 1, size(Op_T,1)
            do i = 1, N_Fl
                if (Op_T(n,i)%get_g_t_alloc()) then
                    if (size(Op_T(n,i)%g_t,1) /= Ltrot) then
                        write(error_unit,*) "Array size of time-dependent coupling Op_T%g_t has to be Ltrot!"
                        CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
                    endif
                endif
            enddo
        enddo
         

        ! Test if user has initialized Calc_FL array
        If ( .not. allocated(Calc_Fl)) then
          allocate(Calc_Fl(N_FL))
          Calc_Fl=.True.
        endif
        ! Count number of flavors to be calculated
        N_FL_eff=0
        Do I=1,N_Fl
          if (Calc_Fl(I)) N_FL_eff=N_FL_eff+1
        Enddo
        reconstruction_needed=.false.
        If (N_FL_eff /= N_FL) reconstruction_needed=.true.
        !initialize the flavor map
        allocate(Calc_Fl_map(N_FL_eff),Phase_array(N_FL))
        N_FL_eff=0
        Do I=1,N_Fl
          if (Calc_Fl(I)) then
             N_FL_eff=N_FL_eff+1
             Calc_Fl_map(N_FL_eff)=I
          endif
        Enddo

        if(Projector) then
           if (.not. allocated(WF_R) .or. .not. allocated(WF_L)) then
              write(error_unit,*) "Projector is selected but there are no trial wave functions!"
              CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
           endif
           do nf_eff=1,N_fl_eff
              nf=Calc_Fl_map(nf_eff)
              if (.not. allocated(WF_R(nf)%P) .or. .not. allocated(WF_L(nf)%P)) then
                 write(error_unit,*) "Projector is selected but there are no trial wave functions!"
                 CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
              endif
           enddo
        endif
        
        Call set_default_values_measuring_interval(Thtrot, Ltrot, Projector)

        If ( .not. get_Global_tau_moves() )  then
           ! This  corresponds to the default updating scheme
           call set_Nt_sequential_start(1)
           call set_Nt_sequential_end(Size(OP_V,1))
           call set_N_Global_tau(0)
        else
           !  Gives the possibility to set parameters in the Hamiltonian file
           tmp_Nt_sequential_start = get_Nt_sequential_start()
           tmp_Nt_sequential_end   = get_Nt_sequential_end()
           tmp_N_Global_tau        = get_N_Global_tau()
           Call ham%Overide_global_tau_sampling_parameters(tmp_Nt_sequential_start,tmp_Nt_sequential_end,tmp_N_Global_tau)
           call set_Nt_sequential_start(tmp_Nt_sequential_start)
           call set_Nt_sequential_end(tmp_Nt_sequential_end)
           call set_N_Global_tau(tmp_N_Global_tau)
        endif
        
        call nsigma%make(N_op, Ltrot)
        Do n = 1,N_op
           nsigma%t(n)              = OP_V(n,1)%type
           nsigma%flip_protocol(n)  = OP_V(n,1)%flip_protocol
        Enddo
           
        File_seeds="seeds"
        Call Set_Random_number_Generator(File_seeds,Seed_in)
        !Write(6,*) Seed_in
               
        Call ham%Hamiltonian_set_nsigma(Initial_field)
        if (allocated(Initial_field)) then
           Call nsigma%in(Group_Comm,Initial_field)
           deallocate(Initial_field)
        else
           Call nsigma%in(Group_Comm)
        endif
        Call Hop_mod_init

        IF (ABS(get_CPU_MAX()) > Zero ) call set_NBin(10000000)
        If (get_N_Global_tau() > 0) then
           Call Wrapgr_alloc
        endif
        
#if defined(HDF5)
#if defined(TEMPERING)
        write(file_dat,'(A,I0,A)') "Temp_",igroup,"/data.h5"
#else
        file_dat = "data.h5"
#endif
#if defined(MPI)
        if ( Irank_g == 0 ) then
#endif
          CALL h5open_f(ierr)
          inquire (file=file_dat, exist=file_exists)
          IF (.not. file_exists) THEN
            ! Create HDF5 file
            CALL h5fcreate_f(file_dat, H5F_ACC_TRUNC_F, file_id, ierr)
            call h5ltset_attribute_string_f(file_id, '/', 'program_name', 'ALF', ierr)
            call h5fclose_f(file_id, ierr)
          endif
          call ham%write_parameters_hdf5(file_dat)

#if defined(MPI)
        endif
#endif
#endif


        Call control_init(Group_Comm)
        Call ham%Alloc_obs(get_Ltau())

        If ( mod(Ltrot,get_Nwrap()) == 0  ) then
           Nstm = Ltrot/get_Nwrap()
        else
           nstm = Ltrot/get_Nwrap() + 1
        endif
        allocate ( Stab_nt(0:Nstm) )
        Stab_nt(0) = 0
        do n = 1,Nstm -1
           Stab_nt(n) = get_Nwrap()*n
        enddo

        Stab_nt(Nstm) = Ltrot

      !   Sequential = .true.
        !TODO: check if sequential is done if some fields are discrete (Warning or error termination?)
        if ( get_Langevin() .or.  get_HMC()  ) then
           if ( get_Langevin() ) then
#if defined(MPI)
                if ( Irank_g == 0 ) then
#endif
                    Call check_langevin_schemes_and_variables()
#if defined(TEMPERING)
                    if ( get_N_exchange_steps() > 0 ) then
                        write(output_unit,*) "Langevin mode does not allow tempering updates."
                        write(output_unit,*) "Overwriting N_exchange_steps to 0."
                    endif
#endif
#if defined(MPI)
                endif
#endif
              call set_sequential(.False.)
              call set_HMC(.False.)
              call set_Global_moves(.False.)
              call set_Global_tau_moves(.False.)
#if defined(TEMPERING)
              call set_N_exchange_steps(0)
#endif
           endif
           Call Langevin_HMC%make(get_Langevin(), get_HMC() , get_Delta_t_Langevin_HMC(), get_Max_Force(), get_Leapfrog_Steps())
        else
           Call Langevin_HMC%set_Update_scheme(get_Langevin(), get_HMC() )
        endif
        call lock_QMC_runtime_settings()
        Call check_update_schemes_compatibility()

#if defined(TEMPERING)
        write(file_info,'(A,I0,A)') "Temp_",igroup,"/info"
#else
        file_info = "info"
#endif

#if defined(MPI)
        if ( Irank_g == 0 ) then
#endif
           Open (Unit = 50,file=file_info,status="unknown",position="append")
           Write(50,*) 'Sweeps                              : ', get_NSweep()
           If ( abs(get_CPU_MAX()) < ZERO ) then
              Write(50,*) 'Bins                                : ', get_NBin()
              Write(50,*) 'No CPU-time limitation '
           else
              Write(50,'(" Prog will stop after hours:",2x,F8.4)') get_CPU_MAX()
           endif
           Write(50,*) 'Measure Int.                        : ', get_LOBS_ST(), get_LOBS_EN()
           Write(50,*) 'Stabilization,Wrap                  : ', get_Nwrap()
           Write(50,*) 'Nstm                                : ', NSTM
           Write(50,*) 'Ltau                                : ', get_Ltau()
           Write(50,*) '# of interacting Ops per time slice : ', Size(OP_V,1)
           If ( get_Propose_S0() ) &
                &  Write(50,*) 'Propose Ising moves according to  bare Ising action'
           If ( get_Global_moves() ) Then
              Write(50,*) 'Global moves are enabled   '
              Write(50,*) '# of global moves / sweep :', get_N_Global()
           Endif
           if ( get_sequential() ) then
               If ( get_Global_tau_moves() ) Then
                  Write(50,*) 'Nt_sequential_start: ', get_Nt_sequential_start()
                  Write(50,*) 'Nt_sequential_end  : ', get_Nt_sequential_end()
                  Write(50,*) 'N_Global_tau       : ', get_N_Global_tau()
               else
                  Write(50,*) 'Default sequential updating '
               endif
            endif
           if ( get_Langevin() ) then
              Write(50,*) 'Langevin del_t: ', get_Delta_t_Langevin_HMC()
              Write(50,*) 'Max Force     : ', get_Max_Force()
           endif
           if ( get_HMC() ) then
              Write(50,*) 'HMC del_t     : ', get_Delta_t_Langevin_HMC()
              Write(50,*) 'Leapfrog_Steps: ', get_Leapfrog_steps()
              Write(50,*) 'HMC_Sweeps:     ', get_N_HMC_sweeps()
           endif

           !Write out info  for  amplitude and flip_protocol
           Toggle  = .false.
           Do n = 1,N_op
              if (nsigma%t(n) == 3 .or. nsigma%t(n) == 4)  Toggle = .true.
           Enddo
           if ( Toggle ) then
              Write(50,*) 'Amplitude  for  t=3,4  vertices is  set to: ', get_Amplitude()
           endif
           Toggle  = .false.
           Do n = 1,N_op
              if (nsigma%t(n) == 4)  Toggle = .true.
           Enddo
           If (Toggle)  then 
              Write(50,"('Flip protocal  for  t=4  vertices is  set to')",advance="no")
              do n1  =  1,4 
                 Toggle1 = .false.
                 do n =1,N_op
                    if (nsigma%Flip_protocol(n)  ==  n1  .and. nsigma%t(n) == 4  )  Toggle1 = .true.
                 enddo
                 if  (Toggle1)   Write(50,"(I2,2x)",advance="no")   n1 
              enddo
              Write(50,*)
           endif
           
#if defined(MPI)
           Write(50,*) 'Number of mpi-processes : ', isize_g
           if(use_mpi_shm) Write(50,*) 'Using mpi-shared memory in chunks of ', chunk_size_gb, 'GB.'
#endif
#if defined(GIT)
           Write(50,*) 'This executable represents commit '&
                &      , GIT_COMMIT_HASH , ' of branch ' , GIT_BRANCH , '.'
#endif
#if defined(STAB1)
           Write(50,*) 'STAB1 is defined '
#endif
#if defined(STAB2)
           Write(50,*) 'STAB2 is defined '
#endif
#if defined(STAB3)
           Write(50,*) 'STAB3 is defined '
#endif
#if defined(STABLOG)
           Write(50,*) 'LOG is defined '
#endif
#if defined(QRREF)
           Write(50,*) 'QRREF is defined '
#endif
#if defined(TEMPERING) && !defined(PARALLEL_PARAMS)
           Write(50,*) '# of exchange steps  ', get_N_exchange_steps()
           Write(50,*) 'Tempering frequency  ', get_N_Tempering_frequency()
           Write(50,*) 'Tempering Calc_det   ', get_Tempering_calc_det()
#endif
           close(50)
#if defined(MPI)
        endif
#endif
        
        !Call Test_Hamiltonian
        Allocate ( Test(Ndim,Ndim), GR(NDIM,NDIM,N_FL), GR_Tilde(NDIM,NDIM,N_FL)  )
        ALLOCATE(udvl(N_FL_eff), udvr(N_FL_eff), udvst(NSTM, N_FL_eff))
        do nf_eff = 1, N_FL_eff
           nf=Calc_Fl_map(nf_eff)
           do n = 1, NSTM
              CALL udvst(n,nf_eff)%alloc(ndim)
           ENDDO
           if (Projector) then
              CALL udvl(nf_eff)%init(ndim,'l',WF_L(nf)%P)
              CALL udvr(nf_eff)%init(ndim,'r',WF_R(nf)%P)
              CALL udvst(NSTM, nf_eff)%reset('l',WF_L(nf)%P)
           else
              CALL udvl(nf_eff)%init(ndim,'l')
              CALL udvr(nf_eff)%init(ndim,'r')
              CALL udvst(NSTM, nf_eff)%reset('l')
           endif
        enddo

        DO NST = NSTM-1,1,-1
           NT1 = Stab_nt(NST+1)
           NT  = Stab_nt(NST  )
           !Write(6,*)'Hi', NT1,NT, NST
           CALL WRAPUL(NT1, NT, UDVL)
           Do nf_eff = 1,N_FL_eff
              UDVST(NST, nf_eff) = UDVL(nf_eff)
           ENDDO
        ENDDO
        NT1 = stab_nt(1)
        CALL WRAPUL(NT1, 0, UDVL)



        NVAR = 1
        Phase_array = cmplx(1.d0, 0.d0, kind(0.D0))
        do nf_eff = 1,N_Fl_eff
           nf=Calc_Fl_map(nf_eff)
           CALL CGR(Z, NVAR, GR(:,:,nf), UDVR(nf_eff), UDVL(nf_eff))
           call Op_phase(Z,OP_V,Nsigma,nf)
           Phase_array(nf)=Z
        Enddo
        if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
        Phase=product(Phase_array)
        Phase=Phase**N_SUN
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
                 CALL Exchange_Step(Phase,GR,udvr, udvl,Stab_nt, udvst, get_N_exchange_steps(), get_Tempering_calc_det())
              endif
#endif
              ! Global updates
              If ( get_Global_moves() ) Call Global_Updates(Phase, GR, udvr, udvl, Stab_nt, udvst, get_N_Global() )


              If ( str_to_upper(Langevin_HMC%get_Update_scheme()) == "LANGEVIN" )  then
                 !  Carry out a Langevin update and calculate equal time observables.
                 Call Langevin_HMC%update(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, &
                      &                   get_LOBS_ST(), get_LOBS_EN(), get_Ltau())
                 
                 IF ( get_Ltau() == 1 ) then
                    If (Projector) then 
                       NST = 0 
                       Call Tau_p ( udvl, udvr, udvst, GR, PHASE, NSTM, STAB_NT, NST, get_LOBS_ST(), get_LOBS_EN())
                       call Langevin_HMC%set_L_Forces(.true.)
                    else
                       Call Tau_m( udvst, GR, PHASE, NSTM, get_Nwrap(), STAB_NT, get_LOBS_ST(), get_LOBS_EN() )
                       call Langevin_HMC%set_L_Forces(.true.)
                    endif
                 endif
              endif

              If (  str_to_upper(Langevin_HMC%get_Update_scheme()) == "HMC" )  then
                 if ( get_sequential() ) call Langevin_HMC%set_L_Forces(.False.)
                 Do n=1, get_N_HMC_sweeps()
                     !  Carry out a Langevin update and calculate equal time observables.
                     Call Langevin_HMC%update(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, &
                          &                   get_LOBS_ST(), get_LOBS_EN(), get_Ltau())
                     if (n /= get_N_HMC_sweeps()) then
                        Call Langevin_HMC%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
                             &  get_LOBS_ST(), get_LOBS_EN(), .True. )
                        Call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
                        call Langevin_HMC%set_L_Forces(.true.)
                     endif
                 enddo
                 
                 !Do time-displaced measurements if needed, else set Calc_Obser_eq=.True. for the very first leapfrog ONLY
                 If ( .not. get_sequential() ) then
                    IF ( get_Ltau() == 1 ) then
                       If (Projector) then 
                          NST = 0 
                          Call Tau_p ( udvl, udvr, udvst, GR, PHASE, NSTM, STAB_NT, NST, get_LOBS_ST(), get_LOBS_EN())
                       else
                          Call Tau_m( udvst, GR, PHASE, NSTM, get_Nwrap(), STAB_NT, get_LOBS_ST(), get_LOBS_EN() )
                       endif
                    else
                       Call Langevin_HMC%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
                       &  get_LOBS_ST(), get_LOBS_EN(), .True. )
                       Call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
                    endif
                    call Langevin_HMC%set_L_Forces(.true.)
                 endif
              endif

              If ( get_sequential() )  then 
                 ! Propagation from 1 to Ltrot
                 ! Set the right storage to 1
                 do nf_eff = 1,N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    if (Projector) then
                       CALL udvr(nf_eff)%reset('r',WF_R(nf)%P)
                    else
                       CALL udvr(nf_eff)%reset('r')
                    endif
                 Enddo
                 
                 NST = 1
                 DO NTAU = 0, LTROT-1
                    NTAU1 = NTAU + 1
                    CALL WRAPGRUP(GR,NTAU,PHASE,get_Propose_S0(), get_Nt_sequential_start(), get_Nt_sequential_end(), get_N_Global_tau())
                    
                    If (NTAU1 == Stab_nt(NST) ) then
                       NT1 = Stab_nt(NST-1)
                       CALL WRAPUR(NT1, NTAU1, udvr)
                       Phase_array = cmplx(1.d0, 0.d0, kind(0.D0))
                       Do nf_eff = 1, N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          ! Read from storage left propagation from LTROT to  NTAU1
                          udvl(nf_eff) = udvst(NST, nf_eff)
                          ! Write in storage right prop from 1 to NTAU1
                          udvst(NST, nf_eff) = udvr(nf_eff)
                          NVAR = 1
                          IF (NTAU1 .GT. LTROT/2) NVAR = 2
                          TEST(:,:) = GR(:,:,nf)
                          CALL CGR(Z1, NVAR, GR(:,:,nf), UDVR(nf_eff), UDVL(nf_eff))
                          Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
                          call Op_phase(Z1,OP_V,Nsigma,nf)
                          Phase_array(nf)=Z1
                       ENDDO
                       if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
                       Z=product(Phase_array)
                       Z=Z**N_SUN
                       Call Control_PrecisionP(Z,Phase)
                       Phase = Z
                       NST = NST + 1
                    ENDIF
                    
                    IF (NTAU1.GE. get_LOBS_ST() .AND. NTAU1.LE. get_LOBS_EN() ) THEN
                       !Call  Global_tau_mod_Test(Gr,ntau1)
                       !Stop
                       !write(*,*) "GR before obser sum: ",sum(GR(:,:,1))
                       !write(*,*) "Phase before obser : ",phase
                       Mc_step_weight = 1.d0
                       If (Symm) then
                          Call Hop_mod_Symm(GR_Tilde,GR,ntau1)
                          !reconstruction of NOT calculated block!!!
                          If (reconstruction_needed) Call ham%GR_reconstruction( GR_Tilde )
                          CALL ham%Obser( GR_Tilde, PHASE, Ntau1, Mc_step_weight )
                       else
                          !reconstruction of NOT calculated block!!!
                          If (reconstruction_needed) Call ham%GR_reconstruction( GR )
                          CALL ham%Obser( GR, PHASE, Ntau1, Mc_step_weight  )
                       endif
                    ENDIF
                 ENDDO
                 
                 Do nf_eff = 1,N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    if (Projector) then
                       CALL udvl(nf_eff)%reset('l',WF_L(nf)%P)
                    else
                       CALL udvl(nf_eff)%reset('l')
                    endif
                 ENDDO
                 
                 NST = NSTM-1
                 DO NTAU = LTROT,1,-1
                    NTAU1 = NTAU - 1
                    CALL WRAPGRDO(GR,NTAU, PHASE,get_Propose_S0(),get_Nt_sequential_start(), get_Nt_sequential_end(), get_N_Global_tau())
                    IF (NTAU1.GE. get_LOBS_ST() .AND. NTAU1.LE. get_LOBS_EN() ) THEN
                       !write(*,*) "GR before obser sum: ",sum(GR(:,:,1))
                       !write(*,*) "Phase before obser : ",phase
                       Mc_step_weight = 1.d0
                       If (Symm) then
                          Call Hop_mod_Symm(GR_Tilde,GR,ntau1)
                          !reconstruction of NOT calculated block!!!
                          If (reconstruction_needed) Call ham%GR_reconstruction( GR_Tilde )
                          CALL ham%Obser( GR_Tilde, PHASE, Ntau1, Mc_step_weight )
                       else
                          !reconstruction of NOT calculated block!!!
                          If (reconstruction_needed) Call ham%GR_reconstruction( GR )
                          CALL ham%Obser( GR, PHASE, Ntau1,Mc_step_weight )
                       endif
                    ENDIF
                    IF ( Stab_nt(NST) == NTAU1 .AND. NTAU1.NE.0 ) THEN
                       NT1 = Stab_nt(NST+1)
                       !Write(6,*) 'Wrapul : ', NT1, NTAU1
                       CALL WRAPUL(NT1, NTAU1, udvl)
                       !Write(6,*)  'Write UL, read UR ', NTAU1, NST
                       Phase_array = cmplx(1.d0, 0.d0, kind(0.D0))
                       do nf_eff = 1,N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          ! Read from store the right prop. from 1 to LTROT/NWRAP-1
                          udvr(nf_eff) = udvst(NST, nf_eff)
                          ! WRITE in store the left prop. from LTROT/NWRAP-1 to 1
                          udvst(NST, nf_eff) = udvl(nf_eff)
                          NVAR = 1
                          IF (NTAU1 .GT. LTROT/2) NVAR = 2
                          TEST(:,:) = GR(:,:,nf)
                          CALL CGR(Z1, NVAR, GR(:,:,nf), UDVR(nf_eff), UDVL(nf_eff))
                          Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
                          call Op_phase(Z1,OP_V,Nsigma,nf)
                          Phase_array(nf)=Z1
                       ENDDO
                       if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
                       Z=product(Phase_array)
                       Z=Z**N_SUN
                       Call Control_PrecisionP(Z,Phase)
                       Phase = Z
                       IF( get_Ltau() == 1 .and. Projector .and. Stab_nt(NST)<=THTROT+1 .and. THTROT+1<Stab_nt(NST+1) ) then
                          Call tau_p ( udvl, udvr, udvst, GR, PHASE, NSTM, STAB_NT, NST,  get_LOBS_ST(), get_LOBS_EN() )
                       endif
                       NST = NST -1
                    ENDIF
                 ENDDO
                 
                 !Calculate and compare green functions on time slice 0.
                 NT1 = Stab_nt(0)
                 NT  = Stab_nt(1)
                 CALL WRAPUL(NT, NT1, udvl)
                 
                 do nf_eff = 1,N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    if (Projector) then
                       CALL udvr(nf_eff)%reset('r',WF_R(nf)%P)
                    else
                       CALL udvr(nf_eff)%reset('r')
                    endif
                 ENDDO
                 Phase_array = cmplx(1.d0, 0.d0, kind(0.D0))
                 do nf_eff = 1,N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    TEST(:,:) = GR(:,:,nf)
                    NVAR = 1
                    CALL CGR(Z1, NVAR, GR(:,:,nf), UDVR(nf_eff), UDVL(nf_eff))
                    Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
                    call Op_phase(Z1,OP_V,Nsigma,nf)
                    Phase_array(nf)=Z1
                 ENDDO
                 if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
                 Z=product(Phase_array)
                 Z=Z**N_SUN
                 Call Control_PrecisionP(Z,Phase)
                 Phase = Z
                 NST =  NSTM
                 Do nf_eff = 1,N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    if (Projector) then
                       CALL udvst(NST, nf_eff)%reset('l',WF_L(nf)%P)
                    else
                       CALL udvst(NST, nf_eff)%reset('l')
                    endif
                 enddo
                 
                 IF ( get_Ltau() == 1 .and. .not. Projector ) then
                    Call TAU_M( udvst, GR, PHASE, NSTM, get_Nwrap(), STAB_NT, get_LOBS_ST(), get_LOBS_EN() )
                 endif
              endif

           ENDDO
           Call ham%Pr_obs(get_Ltau())
#if defined(TEMPERING) && !defined(PARALLEL_PARAMS)
           Call Global_Tempering_Pr
#endif

           Call nsigma%out(Group_Comm)

           call system_clock(count_bin_end)
           prog_truncation = .false.
           if ( abs(get_CPU_MAX()) > Zero ) then
              Call make_truncation(prog_truncation,get_CPU_MAX(),count_bin_start,count_bin_end,group_comm)
           endif
           If (prog_truncation) then
              Nbin_eff = nbc
              exit !exit the loop over the bin index, labelled NBC.
           Endif
        Enddo

        ! Deallocate things
        DO nf_eff = 1, N_FL_eff
           CALL udvl(nf_eff)%dealloc
           CALL udvr(nf_eff)%dealloc
           do n = 1, NSTM
              CALL udvst(n,nf_eff)%dealloc
           ENDDO
        ENDDO
        if (Projector) then
           DO nf = 1, N_FL
              CALL WF_clear(WF_R(nf))
              CALL WF_clear(WF_L(nf))
           ENDDO
        endif
        DEALLOCATE(udvl, udvr, udvst)
        DEALLOCATE(GR, TEST, Stab_nt,GR_Tilde)
        if (Projector) DEALLOCATE(WF_R, WF_L)
        If (get_N_Global_tau() > 0) then
           Call Wrapgr_dealloc
        endif
        do nf = 1, N_FL
          do n = 1, size(OP_V,1)
            call Op_clear(Op_V(n,nf),Op_V(n,nf)%N)
          enddo
          do n = 1, size(OP_T,1)
            call Op_clear(Op_T(n,nf),Op_T(n,nf)%N)
          enddo
        enddo

#if defined(MPI)  
        ! Gracefully deallocate all shared MPI memory (thw whole chunks)
        ! irrespective of where they actually have been used
        call deallocate_all_shared_memory
#endif

        Call Control_Print(Group_Comm, Langevin_HMC%get_Update_scheme())

#if defined(MPI)
        If (Irank_g == 0 ) then
#endif
           if ( abs(get_CPU_MAX()) > Zero ) then
#if defined(TEMPERING)
              write(file_info,'(A,I0,A)') "Temp_",igroup,"/info"
#else
              file_info = "info"
#endif
              Open (Unit=50,file=file_info, status="unknown", position="append")
              Write(50,*)' Effective number of bins   : ', Nbin_eff
              Close(50)
           endif
#if defined(MPI)
        endif
#endif
        
        Call Langevin_HMC%clean()
        deallocate(Calc_Fl_map,Phase_array)

         ! Delete the file RUNNING since the simulation finished successfully
#if defined(MPI)
        If (  Irank == 0 ) then
#endif
           open(unit=5, file='RUNNING', status='old')
           close(5, status='delete')
#if defined(MPI)
        endif
#endif

#ifdef MPI
        CALL MPI_FINALIZE(ierr)
#endif

      end Program Main
