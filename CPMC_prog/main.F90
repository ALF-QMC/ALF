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
        use gfun_mod
        use set_random
        use stepwlk_mod

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
        CLASS(UDV_State), DIMENSION(:,:), ALLOCATABLE :: phi_trial, phi_0, phi_bp_l, phi_bp_r

        Integer :: N_eqwlk, N_blksteps, i_wlk, j_step, N_blk_eff, i_blk, N_blk, NSTM
        Integer :: Nwrap, itv_pc, Ltau, lmetropolis, ntau_bp, ltrot_bp, nsweep, nwarmup
        Real(Kind=Kind(0.d0)) :: CPU_MAX
        Character (len=64) :: file_seeds, file_para, file_dat, file_info, ham_name
        Integer :: Seed_in

        integer :: mpi_per_parameter_set

#ifdef HDF5
        INTEGER(HID_T) :: file_id
        Logical :: file_exists
#endif
        
        NAMELIST /VAR_CPMC/ Nwrap, ltrot_bp, itv_pc, Ltau, lmetropolis, N_blk, N_blksteps, Nsweep, Nwarmup, CPU_MAX

        NAMELIST /VAR_HAM_NAME/ ham_name

        !  General
        Integer :: Ierr, I,nf, nf_eff, nst, n, N_op, NVAR
        Complex (Kind=Kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Z, Z1, E0_iwlk, Z_2
        Real    (Kind=Kind(0.d0)) :: ZERO = 10D-8, X, X1

        ! Storage for  stabilization steps
        Integer, dimension(:), allocatable :: Stab_nt 
        
        ! Space for storage.
        CLASS(UDV_State), Dimension(:,:,:), ALLOCATABLE :: udvst

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
            write(*,'(a)') ' ===================================================================================='
            write(*,*)
            write(*,'(a)') '          The auxiliary  field quantum monte carlo package                           '
            write(*,*)
            write(*,'(a)') '                 A      FFFFF   QQQ     M   M      CCCC                              '
            write(*,'(a)') '                A A     F      Q   Q   M M M M    C                                  '
            write(*,'(a)') '               AAAAA    FFFFF  Q   Q   M M M M    C                                  '
            write(*,'(a)') '              A     A   F      Q   Q   M M M M    C                                  '
            write(*,'(a)') '             A       A  F       QQQ   M   M   M    CCCC                              '
            write(*,'(a)') '                                  \                                                  '
            write(*,*)
            write(*,*)
            write(*,'(a)') '             written by Zihong Liu ( zihong.liu@tu-dresden.de )                      '
            write(*,*)
            write(*,'(a)') ' ------------------------------------------------------------------------------------'
            write(*,*)

#ifdef MPI
        endif
#endif

#ifdef MPI
        mpi_per_parameter_set = Isize
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
           Nwrap=0; CPU_MAX = 0.d0;
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
        CALL MPI_BCAST(ltrot_bp             ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(itv_pc               ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Ltau                 ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(lmetropolis          ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(N_blk                ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(N_blksteps           ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Nsweep               ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(Nwarmup              ,1 ,MPI_INTEGER  ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(CPU_MAX              ,1 ,MPI_REAL8    ,0,MPI_COMM_i,ierr)
        CALL MPI_BCAST(ham_name             ,64,MPI_CHARACTER,0,MPI_COMM_i,ierr)
#endif
        Call Fields_init()
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
        allocate(Calc_Fl_map(N_FL_eff))
        N_FL_eff=0
        Do I=1,N_Fl
          if (Calc_Fl(I)) then
             N_FL_eff=N_FL_eff+1
             Calc_Fl_map(N_FL_eff)=I
          endif
        Enddo

        File_seeds="seeds"
        Call Set_Random_number_Generator(File_seeds,Seed_in)

        !! To do: modify nsigma structure
        allocate(nsigma_bp(N_wlk))
        allocate(nsigma_qr(N_wlk))
        if ( ltau .eq. 1 ) ltrot_bp = ltrot_bp + ltrot
        do i_wlk =1, N_wlk
            call nsigma_bp(i_wlk)%make(N_op, ltrot_bp)
            call nsigma_qr(i_wlk)%make(N_op, nwrap   )
            Do n = 1,N_op
               nsigma_bp(i_wlk)%t(n)  = OP_V(n,1)%type
               nsigma_qr(i_wlk)%t(n)  = OP_V(n,1)%type
            Enddo
        enddo

        Call Hop_mod_init

        IF (ABS(CPU_MAX) > Zero ) N_blk = 10000000

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
           Write(50,*) 'Nwrap                                : ', Nwrap
           Write(50,*) 'N_blk                                : ', N_blk
           Write(50,*) 'N_blksteps                           : ', N_blksteps
           Write(50,*) 'Nsweep                               : ', Nsweep
           Write(50,*) 'Nwarmup                              : ', Nwarmup
           Write(50,*) 'ltrot_bp                             : ', ltrot_bp
           Write(50,*) 'itv_pc                               : ', itv_pc
           Write(50,*) 'ltau                                 : ', ltau
           Write(50,*) 'lmetropolis                          : ', lmetropolis
           If ( abs(CPU_MAX) < ZERO ) then
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
#if defined(MPI)
        endif
#endif

        Call ham%Alloc_obs(ltau, lmetropolis)

        allocate( gr(ndim,ndim,n_fl,n_grc) )
        allocate( phi_trial(N_FL_eff, N_slat), phi_0   (N_FL_eff, N_wlk))
        allocate( phi_bp_l (N_FL_eff, N_grc ), phi_bp_r(N_FL_eff, N_wlk))

        ! we require ltrot_bp >= nwrap and ltrot_bp <= N_blksteps
        If ( mod(ltrot_bp,nwrap) == 0  ) then
           nstm = ltrot_bp/nwrap
        else
           nstm = ltrot_bp/nwrap + 1
        endif
        Allocate( udvst(NSTM, N_FL_eff, N_grc))
        
        weight_k(:) = 1.d0
        
        ! init slater determinant
        call initial_wlk( phi_trial, phi_0, phi_bp_l, phi_bp_r, udvst, STAB_nt, GR, nwrap )
        call store_phi  ( phi_0, phi_bp_r )

        Call control_init(Group_Comm)

        !! main loop
        ntau_bp = 1 ! ntau_bp is to record the field for back propagation
        NST     = 1
 
        do i_blk=1, N_blk

            call system_clock(count_bin_start)
            
            !! initial obs
            call ham%Init_obs(ltau)
        
            do j_step=1, N_blksteps
                !! population control
                if ( mod(j_step, itv_pc) .eq. 0 ) then
                    call population_control(phi_0, phi_bp_r)
                endif

                !! propagate the walkers:
                call stepwlk_move(Phi_trial, Phi_0, GR, ntau_bp );
                !! QR decomposition for stablization
                if ( ntau_bp .eq. Stab_nt(NST) ) then
                    call re_orthonormalize_walkers(Phi_0, 'U')
                    NST = NST + 1
                endif
                
                ! Measurement and update fac_norm
                if ( ntau_bp .eq. ltrot_bp ) then
                    ntau_bp = 0
                    NST     = 1

                    call backpropagation( GR, phi_bp_l, phi_bp_r, udvst, stab_nt, ltau, lmetropolis, nsweep, nwarmup )
                    
                    !! store phi_0 for the next measurement
                    call store_phi( phi_0, phi_bp_r )

                    !! Update fac_norm
                    call ham%update_fac_norm(GR, j_step+(i_blk-1)*N_blksteps)

                endif
                
                ntau_bp = ntau_bp + 1

            enddo
            !! output print
            call ham%Pr_obs(ltau)

            call seed_vec_out

            call wavefunction_out_hdf5( phi_0 )

            call system_clock(count_bin_end)
            prog_truncation = .false.
            if ( abs(CPU_MAX) > Zero ) then
               Call make_truncation(prog_truncation,CPU_MAX,count_bin_start,count_bin_end,group_comm)
            endif
            If (prog_truncation) then
               exit !exit the loop over the step index
            Endif

        enddo

        deallocate(Calc_Fl_map)

#if defined(MPI)  
        ! Gracefully deallocate all shared MPI memory (thw whole chunks)
        ! irrespective of where they actually have been used
        call deallocate_all_shared_memory
#endif

        Call Control_Print(Group_Comm)

#ifdef MPI
        CALL MPI_FINALIZE(ierr)
#endif
        stop

      end Program Main
