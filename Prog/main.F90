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

        Use runtime_error_mod
        Use Operator_mod
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

#ifdef MPI
        Use mpi
#endif
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


        Integer :: Nwrap, NSweep, NBin, NBin_eff,Ltau, NSTM, NT, NT1, NVAR, LOBS_EN, LOBS_ST, NBC, NSW
        Integer :: NTAU, NTAU1
        Real(Kind=Kind(0.d0)) :: CPU_MAX
        Character (len=64) :: file_seeds, file_para, file_dat, file_info, ham_name
        Integer :: Seed_in
        Real (Kind=Kind(0.d0)) , allocatable, dimension(:,:) :: Initial_field

        ! Space for choosing sampling scheme
        Logical :: Propose_S0, Tempering_calc_det
        Logical :: Global_moves, Global_tau_moves
        Integer :: N_Global
        Integer :: Nt_sequential_start, Nt_sequential_end, mpi_per_parameter_set
        Integer :: N_Global_tau
        Logical :: Sequential

#ifdef HDF5
        INTEGER(HID_T) :: file_id
        Logical :: file_exists
#endif
        !  Space for reading in Langevin & HMC  parameters
        Logical                      :: Langevin,  HMC
        Integer                      :: Leapfrog_Steps, N_HMC_sweeps
        Real  (Kind=Kind(0.d0))      :: Delta_t_Langevin_HMC, Max_Force
          
#if defined(TEMPERING)
        Integer :: N_exchange_steps, N_Tempering_frequency
        NAMELIST /VAR_TEMP/  N_exchange_steps, N_Tempering_frequency, mpi_per_parameter_set, Tempering_calc_det
#endif

        NAMELIST /VAR_QMC/   Nwrap, NSweep, NBin, Ltau, LOBS_EN, LOBS_ST, CPU_MAX, &
             &               Propose_S0,Global_moves,  N_Global, Global_tau_moves, &
             &               Nt_sequential_start, Nt_sequential_end, N_Global_tau, &
             &               sequential, Langevin, HMC, Delta_t_Langevin_HMC, &
             &               Max_Force, Leapfrog_steps, N_HMC_sweeps

        NAMELIST /VAR_HAM_NAME/ ham_name

        !  General
        Integer :: Ierr, I,nf, nf_eff, nst, n, N_op
        Complex (Kind=Kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Phase, Z, Z1
        Real    (Kind=Kind(0.d0)) :: ZERO = 10D-8, X, X1
        Real    (Kind=Kind(0.d0)) :: Mc_step_weight

        ! Storage for  stabilization steps
        Integer, dimension(:), allocatable :: Stab_nt 

        ! Space for storage.
        CLASS(UDV_State), Dimension(:,:), ALLOCATABLE :: udvst

        ! For tests
        Real (Kind=Kind(0.d0)) :: Weight, Weight_tot

        ! For the truncation of the program:
        logical                   :: prog_truncation, run_file_exists
        integer (kind=kind(0.d0)) :: count_bin_start, count_bin_end
        
        ! For MPI shared memory
        character(64), parameter :: name="ALF_SHM_CHUNK_SIZE_GB"
        character(64) :: chunk_size_str
        Real    (Kind=Kind(0.d0)) :: chunk_size_gb

#ifdef MPI
        Integer        :: Isize, Irank, Irank_g, Isize_g, color, key, igroup, MPI_COMM_i

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
        mpi_per_parameter_set = 1   ! Default value
        Tempering_calc_det = .true. ! Default value
        OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
        IF (ierr /= 0) THEN
           WRITE(error_unit,*) 'main: unable to open <parameters>',ierr
           CALL Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
        END IF
        READ(5,NML=VAR_TEMP)
        CLOSE(5)
        CALL MPI_BCAST(N_exchange_steps        ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(N_Tempering_frequency   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(mpi_per_parameter_set   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(Tempering_calc_det      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        if ( mod(ISIZE,mpi_per_parameter_set) .ne. 0 ) then
           Write (error_unit,*) "mpi_per_parameter_set is not a multiple of total mpi processes"
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        Call Global_Tempering_setup
#elif !defined(TEMPERING)  && defined(MPI)
        mpi_per_parameter_set = Isize
#elif defined(TEMPERING)  && !defined(MPI)
        Write(error_unit,*) 'Mpi has to be defined for tempering runs'
        CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif
#if defined(PARALLEL_PARAMS) && !defined(TEMPERING)
        Write(error_unit,*) 'TEMPERING has to be defined for PARALLEL_PARAMS'
        CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
#endif

#ifdef MPI
        color = irank/mpi_per_parameter_set
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


#ifdef MPI
#ifdef PARALLEL_PARAMS
        MPI_COMM_i = Group_Comm
        If ( irank_g == 0 ) then
           write(file_para,'(A,I0,A)') "Temp_", igroup, "/parameters"
#else
        MPI_COMM_i = MPI_COMM_WORLD
        If ( Irank == 0 ) then
           file_para = "parameters"
#endif
#else
           file_para = "parameters"
#endif
           ! This is a set of variables that  identical for each simulation.
           Nwrap=0;  NSweep=0; NBin=0; Ltau=0; LOBS_EN = 0;  LOBS_ST = 0;  CPU_MAX = 0.d0
           Propose_S0 = .false. ;  Global_moves = .false. ; N_Global = 0
           Global_tau_moves = .false.; sequential = .true.; Langevin = .false. ; HMC =.false.
           Delta_t_Langevin_HMC = 0.d0;  Max_Force = 0.d0 ; Leapfrog_steps = 0; N_HMC_sweeps = 1
           Nt_sequential_start = 1 ;  Nt_sequential_end  = 0;  N_Global_tau  = 0
           OPEN(UNIT=5,FILE=file_para,STATUS='old',ACTION='read',IOSTAT=ierr)
           IF (ierr /= 0) THEN
              WRITE(error_unit,*) 'main: unable to open <parameters>', file_para, ierr
              CALL Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
           END IF
           READ(5,NML=VAR_QMC)
           REWIND(5)
           READ(5,NML=VAR_HAM_NAME)
           CLOSE(5)
           NBin_eff = NBin
#ifdef MPI
        Endif
        CALL MPI_BCAST(Nwrap                ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(NSweep               ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(NBin                 ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Ltau                 ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(LOBS_EN              ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(LOBS_ST              ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(CPU_MAX              ,1 ,MPI_REAL8    ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Propose_S0           ,1 ,MPI_LOGICAL  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Global_moves         ,1 ,MPI_LOGICAL  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(N_Global             ,1 ,MPI_Integer  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Global_tau_moves     ,1 ,MPI_LOGICAL  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Nt_sequential_start  ,1 ,MPI_Integer  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Nt_sequential_end    ,1 ,MPI_Integer  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(N_Global_tau         ,1 ,MPI_Integer  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(sequential           ,1 ,MPI_LOGICAL  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Langevin             ,1 ,MPI_LOGICAL  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(HMC                  ,1 ,MPI_LOGICAL  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Leapfrog_steps       ,1 ,MPI_Integer  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(N_HMC_sweeps         ,1 ,MPI_Integer  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Max_Force            ,1 ,MPI_REAL8    ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Delta_t_Langevin_HMC ,1 ,MPI_REAL8    ,0,MPI_COMM_i,ierr)

        CALL MPI_BCAST(ham_name             ,64,MPI_CHARACTER,0,MPI_COMM_i,ierr)
#endif
        Call Fields_init()
        Call Alloc_Ham(ham_name)
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
        !  Default values of  measuring interval.
        if (Projector)  then
           if ( LOBS_ST == 0  ) then
              LOBS_ST = Thtrot+1
           else
              If (LOBS_ST < Thtrot+1 ) then
                 Write(error_unit,*) 'Measuring out of dedicating interval, LOBS_ST too small.'
                 CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
              endif
           endif
           if ( LOBS_EN == 0) then
              LOBS_EN = Ltrot-Thtrot
           else
              If (LOBS_EN > Ltrot-Thtrot ) then
                 Write(error_unit,*) 'Measuring out of dedicating interval, LOBS_EN too big.'
                 CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
              endif
           endif
        else
           if ( LOBS_ST == 0  ) then
              LOBS_ST = 1
           endif
           if ( LOBS_EN == 0) then
              LOBS_EN =  Ltrot
           endif
        endif
        If ( .not. Global_tau_moves )  then
           ! This  corresponds to the default updating scheme
           Nt_sequential_start = 1
           Nt_sequential_end   = Size(OP_V,1)
           N_Global_tau        = 0
        else
           !  Gives the possibility to set parameters in the Hamiltonian file
           Call ham%Overide_global_tau_sampling_parameters(Nt_sequential_start,Nt_sequential_end,N_Global_tau)
        endif
        
        call nsigma%make(N_op, Ltrot)
        Do n = 1,N_op
           nsigma%t(n)  = OP_V(n,1)%type
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

        IF (ABS(CPU_MAX) > Zero ) NBIN = 10000000
        If (N_Global_tau > 0) then
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
        Call ham%Alloc_obs(Ltau)

        If ( mod(Ltrot,nwrap) == 0  ) then
           Nstm = Ltrot/nwrap
        else
           nstm = Ltrot/nwrap + 1
        endif
        allocate ( Stab_nt(0:Nstm) )
        Stab_nt(0) = 0
        do n = 1,Nstm -1
           Stab_nt(n) = nwrap*n
        enddo

        Stab_nt(Nstm) = Ltrot

      !   Sequential = .true.
        !TODO: check if sequential is done if some fields are discrete (Warning or error termination?)
        if ( Langevin .or.  HMC  ) then
           if (Langevin) then
#if defined(MPI)
               if ( Irank_g == 0 ) then
#endif
                  if (sequential) then 
                     write(output_unit,*) "Langevin mode does not allow sequential updates."
                     write(output_unit,*) "Overriding Sequential=.True. from parameter files."
                  endif
                  if (HMC) then 
                     write(output_unit,*) "Langevin mode does not allow HMC updates."
                     write(output_unit,*) "Overriding HMC=.True. from parameter files."
                  endif
                  if (Global_moves) then 
                     write(output_unit,*) "Langevin mode does not allow global updates."
                     write(output_unit,*) "Overriding Global_moves=.True. from parameter files."
                  endif
                  if (Global_tau_moves) then 
                     write(output_unit,*) "Langevin mode does not allow global tau updates."
                     write(output_unit,*) "Overriding Global_tau_moves=.True. from parameter files."
                  endif
#if defined(TEMPERING)
                  if ( N_exchange_steps > 0 ) then
                     write(output_unit,*) "Langevin mode does not allow tempering updates."
                     write(output_unit,*) "Overwriting N_exchange_steps to 0."
                  end if
#endif
#if defined(MPI)
               endif
#endif
               Sequential = .False.
               HMC = .False.
               Global_moves = .False.
               Global_tau_moves = .False.
#if defined(TEMPERING)
               N_exchange_steps = 0
#endif
           endif
           Call Langevin_HMC%make(Langevin, HMC , Delta_t_Langevin_HMC, Max_Force, Leapfrog_steps)
        else
           Call Langevin_HMC%set_Update_scheme(Langevin, HMC )
        endif

        if ( .not. Sequential .and. Global_tau_moves) then
           write(output_unit,*) "Warning: Sequential = .False. and Global_tau_moves = .True."
           write(output_unit,*) "in the parameter file. Global tau updates will not occur if"
           write(output_unit,*) "Sequential is set to .False. ."
        endif

        if ( .not. Sequential .and. .not. HMC .and. .not. Langevin .and. .not. Global_moves) then
         write(output_unit,*) "Warning: no updates will occur as Sequential, HMC, Langevin, and"
         write(output_unit,*) "Global_moves are all .False. in the parameter file."
        endif

        if ( Sequential .and. Nt_sequential_end < Nt_sequential_start ) then
         write(output_unit,*) "Warning: Nt_sequential_end is smaller than Nt_sequential_start"
        endif

#if defined(TEMPERING)
        write(file_info,'(A,I0,A)') "Temp_",igroup,"/info"
#else
        file_info = "info"
#endif

#if defined(MPI)
        if ( Irank_g == 0 ) then
#endif
           Open (Unit = 50,file=file_info,status="unknown",position="append")
           Write(50,*) 'Sweeps                              : ', Nsweep
           If ( abs(CPU_MAX) < ZERO ) then
              Write(50,*) 'Bins                                : ', NBin
              Write(50,*) 'No CPU-time limitation '
           else
              Write(50,'(" Prog will stop after hours:",2x,F8.4)') CPU_MAX
           endif
           Write(50,*) 'Measure Int.                        : ', LOBS_ST, LOBS_EN
           Write(50,*) 'Stabilization,Wrap                  : ', Nwrap
           Write(50,*) 'Nstm                                : ', NSTM
           Write(50,*) 'Ltau                                : ', Ltau
           Write(50,*) '# of interacting Ops per time slice : ', Size(OP_V,1)
           If ( Propose_S0 ) &
                &  Write(50,*) 'Propose Ising moves according to  bare Ising action'
           If ( Global_moves ) Then
              Write(50,*) 'Global moves are enabled   '
              Write(50,*) '# of global moves / sweep :', N_Global
           Endif
           if ( sequential ) then
               If ( Global_tau_moves ) Then
                  Write(50,*) 'Nt_sequential_start: ', Nt_sequential_start
                  Write(50,*) 'Nt_sequential_end  : ', Nt_sequential_end
                  Write(50,*) 'N_Global_tau       : ', N_Global_tau
               else
                  Write(50,*) 'Default sequential updating '
               endif
            endif
           if ( Langevin ) then
              Write(50,*) 'Langevin del_t: ', Delta_t_Langevin_HMC
              Write(50,*) 'Max Force     : ', Max_Force
           endif
           if ( HMC ) then
              Write(50,*) 'HMC del_t     : ', Delta_t_Langevin_HMC
              Write(50,*) 'Leapfrog_Steps: ', Leapfrog_Steps
              Write(50,*) 'HMC_Sweeps:     ', N_HMC_sweeps
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
           Write(50,*) '# of exchange steps  ',N_exchange_steps
           Write(50,*) 'Tempering frequency  ',N_Tempering_frequency
           Write(50,*) 'Tempering Calc_det   ',Tempering_calc_det
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

        DO  NBC = 1, NBIN
           ! Here, you have the green functions on time slice 1.
           ! Set bin observables to zero.

           call system_clock(count_bin_start)

           Call ham%Init_obs(Ltau)
#if defined(TEMPERING)
           Call Global_Tempering_init_obs
#endif

           DO NSW = 1, NSWEEP

#if defined(TEMPERING) && !defined(PARALLEL_PARAMS)
              IF (MOD(NSW,N_Tempering_frequency) == 0) then
                 !Write(6,*) "Irank, Call tempering", Irank, NSW, N_exchange_steps
                 CALL Exchange_Step(Phase,GR,udvr, udvl,Stab_nt, udvst, N_exchange_steps, Tempering_calc_det)
              endif
#endif
              ! Global updates
              If (Global_moves) Call Global_Updates(Phase, GR, udvr, udvl, Stab_nt, udvst,N_Global)


              If (  trim(Langevin_HMC%get_Update_scheme()) == "Langevin" )  then
                 !  Carry out a Langevin update and calculate equal time observables.
                 Call Langevin_HMC%update(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, &
                      &                   LOBS_ST, LOBS_EN, LTAU)
                 
                 IF ( LTAU == 1 ) then
                    If (Projector) then 
                       NST = 0 
                       Call Tau_p ( udvl, udvr, udvst, GR, PHASE, NSTM, STAB_NT, NST, LOBS_ST, LOBS_EN)
                       call Langevin_HMC%set_L_Forces(.true.)
                    else
                       Call Tau_m( udvst, GR, PHASE, NSTM, NWRAP, STAB_NT, LOBS_ST, LOBS_EN )
                       call Langevin_HMC%set_L_Forces(.true.)
                    endif
                 endif
              endif

              If (  trim(Langevin_HMC%get_Update_scheme()) == "HMC" )  then
                 if (Sequential) call Langevin_HMC%set_L_Forces(.False.)
                 Do n=1, N_HMC_sweeps
                     !  Carry out a Langevin update and calculate equal time observables.
                     Call Langevin_HMC%update(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, &
                          &                   LOBS_ST, LOBS_EN, LTAU)
                     if (n /= N_HMC_sweeps) then
                        Call Langevin_HMC%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
                             &  LOBS_ST, LOBS_EN, .True. )
                        Call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
                        call Langevin_HMC%set_L_Forces(.true.)
                     endif
                 enddo
                 
                 !Do time-displaced measurements if needed, else set Calc_Obser_eq=.True. for the very first leapfrog ONLY
                 If ( .not. sequential) then
                    IF ( LTAU == 1 ) then
                       If (Projector) then 
                          NST = 0 
                          Call Tau_p ( udvl, udvr, udvst, GR, PHASE, NSTM, STAB_NT, NST, LOBS_ST, LOBS_EN)
                       else
                          Call Tau_m( udvst, GR, PHASE, NSTM, NWRAP, STAB_NT, LOBS_ST, LOBS_EN )
                       endif
                    else
                       Call Langevin_HMC%calc_Forces(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst,&
                       &  LOBS_ST, LOBS_EN, .True. )
                       Call Langevin_HMC_Reset_storage(Phase, GR, udvr, udvl, Stab_nt, udvst)
                    endif
                    call Langevin_HMC%set_L_Forces(.true.)
                 endif
              endif

              If (Sequential)  then 
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
                    CALL WRAPGRUP(GR,NTAU,PHASE,Propose_S0, Nt_sequential_start, Nt_sequential_end, N_Global_tau)
                    
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
                    
                    IF (NTAU1.GE. LOBS_ST .AND. NTAU1.LE. LOBS_EN ) THEN
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
                    CALL WRAPGRDO(GR,NTAU, PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau)
                    IF (NTAU1.GE. LOBS_ST .AND. NTAU1.LE. LOBS_EN ) THEN
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
                       IF( LTAU == 1 .and. Projector .and. Stab_nt(NST)<=THTROT+1 .and. THTROT+1<Stab_nt(NST+1) ) then
                          Call tau_p ( udvl, udvr, udvst, GR, PHASE, NSTM, STAB_NT, NST,  LOBS_ST, LOBS_EN )
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
                 
                 IF ( LTAU == 1 .and. .not. Projector ) then
                    Call TAU_M( udvst, GR, PHASE, NSTM, NWRAP, STAB_NT, LOBS_ST, LOBS_EN )
                 endif
              endif

           ENDDO
           Call ham%Pr_obs(Ltau)
#if defined(TEMPERING)
           Call Global_Tempering_Pr
#endif

           Call nsigma%out(Group_Comm)

           call system_clock(count_bin_end)
           prog_truncation = .false.
           if ( abs(CPU_MAX) > Zero ) then
              Call make_truncation(prog_truncation,CPU_MAX,count_bin_start,count_bin_end,group_comm)
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
        If (N_Global_tau > 0) then
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
           if ( abs(CPU_MAX) > Zero ) then
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
