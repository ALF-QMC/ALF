Program Main

        Use runtime_error_mod
        Use Operator_mod
        Use Lattices_v3
        Use MyMats
        Use Hamiltonian_main
        Use Control
        Use Hop_mod
        Use UDV_State_mod
        Use Fields_mod
        Use WaveFunction_mod
        use iso_fortran_env, only: output_unit, error_unit
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

        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:,:), Allocatable   :: GR
        CLASS(UDV_State), DIMENSION(:,:), ALLOCATABLE :: phi_trial, udvr, phi_0
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:)  , Allocatable   :: Phase_array

        Integer :: N_eqblk, N_blksteps, i_blk, j_step, N_blksteps_eff
        Integer :: Nwrap, itv_pc, itv_Em
        Real(Kind=Kind(0.d0)) :: CPU_MAX
        Character (len=64) :: file_seeds, file_para, file_dat, file_info, ham_name
        Integer :: Seed_in

#ifdef HDF5
        INTEGER(HID_T) :: file_id
        Logical :: file_exists
#endif
        
        NAMELIST /VAR_CPMC/ Nwrap, itv_pc, itv_Em, N_eqblk, N_blksteps, i_blk, j_step, CPU_MAX, N_blksteps

        NAMELIST /VAR_HAM_NAME/ ham_name

        !  General
        Integer :: Ierr, I,nf, nf_eff, nst, n, N_op
        Complex (Kind=Kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Z, Z1
        COMPLEX (Kind=Kind(0.d0)), Dimension(:)  , Allocatable   :: Phase
        Real    (Kind=Kind(0.d0)) :: ZERO = 10D-8, X, X1

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

#ifdef MPI
        color = irank/mpi_per_parameter_set
        key   =  0
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Group_comm, ierr)
        call MPI_Comm_rank(Group_Comm, Irank_g, ierr)
        call MPI_Comm_size(Group_Comm, Isize_g, ierr)
        igroup           = irank/isize_g
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

#ifdef MPI
        MPI_COMM_i = MPI_COMM_WORLD
        If ( Irank == 0 ) then
           file_para = "parameters"
#else
           file_para = "parameters"
#endif

           ! This is a set of variables that  identical for each simulation.
           Nwrap=0; Ltau=0; CPU_MAX = 0.d0;
           OPEN(UNIT=5,FILE=file_para,STATUS='old',ACTION='read',IOSTAT=ierr)
           IF (ierr /= 0) THEN
              WRITE(error_unit,*) 'main: unable to open <parameters>', file_para, ierr
              CALL Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
           END IF
           READ(5,NML=VAR_CPMC)
           REWIND(5)
           READ(5,NML=VAR_HAM_NAME)
           CLOSE(5)
#ifdef MPI
        Endif
        CALL MPI_BCAST(Nwrap                ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(itv_pc               ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(itv_Em               ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(N_eqblk              ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(N_blksteps           ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(CPU_MAX              ,1 ,MPI_REAL8    ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(ham_name             ,64,MPI_CHARACTER,0,MPI_COMM_i,ierr)
#endif

        Call Alloc_Ham(ham_name)
        Call ham%Ham_set()

        N_op = Size(OP_V,1)

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
        allocate(Calc_Fl_map(N_FL_eff),Phase_array(N_FL,N_blk), Phase(N_blk))
        N_FL_eff=0
        Do I=1,N_Fl
          if (Calc_Fl(I)) then
             N_FL_eff=N_FL_eff+1
             Calc_Fl_map(N_FL_eff)=I
          endif
        Enddo
        
        File_seeds="seeds"
        Call Set_Random_number_Generator(File_seeds,Seed_in)

        allocate(nisgma(N_blk))
        do i_blk =1, N_blk
            call nsigma(i_blk)%make(N_op, nwrap)
            Do n = 1,N_op
               nsigma(i_blk)%t(n)  = OP_V(n,1)%type
            Enddo
            Call nsigma(i_blk)%in(Group_Comm)
        enddo
        
        Call Hop_mod_init

        IF (ABS(CPU_MAX) > Zero ) N_blksteps = 10000000

#if defined(HDF5)
        file_dat = "data.h5"
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

        file_info = "info"
#if defined(MPI)
        if ( Irank_g == 0 ) then
#endif
           Open (Unit = 50,file=file_info,status="unknown",position="append")
           If ( abs(CPU_MAX) < ZERO ) then
              Write(50,*) 'N_eqblk                              : ', N_eqblk
              Write(50,*) 'N_blksteps                           : ', N_blksteps
              Write(50,*) 'Nwrap                                : ', Nwrap
              Write(50,*) 'itv_pc                               : ', itv_pc
              Write(50,*) 'itv_Em                               : ', itv_Em
              Write(50,*) 'No CPU-time limitation '
           else
              Write(50,'(" Prog will stop after hours:",2x,F8.4)') CPU_MAX
           endif
           Write(50,*) '# of interacting Ops per time slice : ', Size(OP_V,1)
           
#if defined(MPI)
           Write(50,*) 'Number of mpi-processes : ', isize_g
           if(use_mpi_shm) Write(50,*) 'Using mpi-shared memory in chunks of ', chunk_size_gb, 'GB.'
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

#if defined(MPI)
        endif
#endif

        Call control_init(Group_Comm)
        Call ham%Alloc_obs
        
        Allocate( GR(NDIM,NDIM,N_FL,N_blk) )
        Allocate( phi_trail(N_FL_eff, N_blk) , udvr(N_FL_eff, N_blk), phi_0(N_FL_eff, N_blk))
        
        Phase_array(:,:) = cmplx(1.d0, 0.d0, kind(0.D0))
        
        ! init slater determinant
        do i_blk  = 1, N_blk
            do nf_eff = 1, N_FL_eff
               nf=Calc_Fl_map(nf_eff)
               CALL udvr(nf_eff, i_blk)%init(ndim,'r',WF_R(nf)%P)
               CALL phi_trial(nf_eff, i_blk)%init(ndim,'l',WF_L(nf)%P)
               CALL phi_0(nf_eff, i_blk)%init(ndim,'r',WF_R(nf)%P)
               
               !Carry out U,D,V decomposition.
               CALL udvr(nf_eff, i_blk)%decompose
               CALL phi_trial(nf_eff, i_blk)%decompose
               CALL phi_0(nf_eff, i_blk)%decompose
            enddo

            NVAR = 1
            do nf_eff = 1,N_Fl_eff
               nf=Calc_Fl_map(nf_eff)
               call CGR(Z, NVAR, GR(:,:,nf, i_blk), phi_0(nf_eff, i_blk), phi_trial(nf_eff, i_blk))
               Phase_array(nf, i_blk)=Z
            enddo
            if (reconstruction_needed) call ham%weight_reconstruction(Phase_array(:,i_blk))
            Phase(i_blk)=product(Phase_array(:,i_blk))
            Phase(i_blk)=Phase(i_blk)**N_SUN
        enddo

        Call control_init(Group_Comm)

        !! main loop
        do j_step=1, N_blksteps
            !! also add a function to init Green's function

            call system_clock(count_bin_start)
            call ham%Init_obs(Ltau)
        
            do i_blk =1, N_blk
                
                !! propagate the walkers:
                call stepwlk(Phi, N_wlk, N_sites, w, O, E_blk(i_blk), W_blk(i_blk), H_k, Proj_k_half, flag_mea, Phi_T, N_up, N_par, U, fac_norm, aux_fld);
                if ( mod(j_step,itv_modsvd) .eq. 0 ) call stblz(Phi, N_wlk, O, N_up, N_par); ! re-orthonormalize the walkers
                if ( mod(j_step,itv_pc)     .eq. 0 ) call pop_cntrl(Phi, w, O, N_wlk, N_sites, N_par); ! population control
                
                !! update the exponent of the pre-factor exp(-deltau*(H-E_T))
                if ( mod(j_step,itv_Em)     .eq. 0 ) fac_norm=(real(E_blk(i_blk)/W_blk(i_blk))-0.5*U*N_par)*dtau;
                if ( mod(j_step,itv_Em)     .eq. 0 ) then
                If (reconstruction_needed) Call ham%GR_reconstruction( GR )
                    CALL ham%Obser( GR, PHASE )
                endif
            enddo
            
            call ham%Pr_obs
            call system_clock(count_bin_end)
            prog_truncation = .false.
            if ( abs(CPU_MAX) > Zero ) then
               Call make_truncation(prog_truncation,CPU_MAX,count_bin_start,count_bin_end,group_comm)
            endif
            If (prog_truncation) then
               exit !exit the loop over the step index
            Endif
        
        enddo
        
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
