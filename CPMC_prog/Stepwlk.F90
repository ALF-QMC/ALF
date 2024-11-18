    Module stepwlk_mod
        Use Hamiltonian_main 
        Use Operator_mod
        Use Control
        Use Hop_mod
        Use UDV_State_mod
        use gfun_mod
        use upgrade_mod
        
        Contains

        SUBROUTINE initial_wlk( phi_trial, phi_0, phi_bp_l, phi_bp_r, udvst, stab_nt, GR, nwrap )
#ifdef MPI
          Use mpi
#endif
          Implicit none
     
          class(udv_state), Dimension(:,:)  , ALLOCATABLE, INTENT(INOUT) :: phi_trial, phi_0, phi_bp_l, phi_bp_r
          class(udv_state), Dimension(:,:,:), ALLOCATABLE, INTENT(INOUT) :: udvst
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:,:), Allocatable, INTENT(INOUT) :: GR
          INTEGER, dimension(:)    ,allocatable,  INTENT(INOUT) :: stab_nt
          INTEGER, INTENT(IN) :: nwrap

          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, i_grc, NSTM, NST, ltrot_bp, ns
          Integer :: i_st, i_ed, ncslat
          Complex (Kind=Kind(0.d0)) :: overlap_old, overlap_new, Z, Z1,Z2, tot_ene, ZP
          Complex (Kind=Kind(0.d0)) :: tot_c_weight, el_tmp
          complex (Kind=Kind(0.d0)) :: det_Vec(N_FL)
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio, X1, wtmp
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8, tot_re_weight, dz2
          Character (LEN=64) :: FILE_TG, FILE_seeds, file_inst, file_antiinst
          Logical ::   LCONF, LCONF_H5, lconf_inst, lconf_antiinst

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          nstm = Size(udvst, 1)
          ltrot_bp = size(nsigma_bp(1)%f,2)
          tot_ene  = cmplx(0.d0,0.d0,kind(0.d0))
          tot_re_weight = 0.d0

          allocate ( Stab_nt(0:nstm) )
          Stab_nt(0) = 0
          do n = 1,Nstm -1
             Stab_nt(n) = nwrap*n
          enddo
          Stab_nt(Nstm) = ltrot_bp

          do i_wlk = 1, N_wlk
             do nf_eff = 1, N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                CALL phi_0   (nf_eff, i_wlk)%init(ndim,'r',WF_R(nf,1)%P)
                CALL phi_bp_r(nf_eff, i_wlk)%init(ndim,'r',WF_R(nf,1)%P)
             enddo
          enddo
          
          do ns = 1, N_slat
             do nf_eff = 1, N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                call phi_trial(nf_eff, ns)%init(ndim,'l',WF_L(nf,ns)%P)
             enddo
             do i_wlk = 1, N_wlk
                i_grc = ns+(i_wlk-1)*N_slat
                do nf_eff = 1, N_FL_eff
                   nf=Calc_Fl_map(nf_eff)
                   do n = 1, nstm
                      CALL udvst(n,nf_eff, i_grc)%alloc(ndim)
                   ENDDO
                   CALL phi_bp_l(nf_eff, i_grc)%init(ndim,'l',WF_L(nf,ns)%P)
                enddo
             enddo
          enddo
            
          file_inst = 'trial_inst.h5'
          file_antiinst = 'trial_antiinst.h5'
          inquire(file=file_inst, exist=lconf_inst)
          inquire(file=file_antiinst, exist=lconf_antiinst)
          if ( lconf_inst .and. lconf_antiinst ) then
              if ( irank_g .eq. 0 ) write (*,*) "read inst-anti inst config as trial wave function"
              call trial_in_hdf5( phi_0, phi_trial, file_inst, file_antiinst )
          endif

          file_tg = "phiin_0.h5"
          INQUIRE (FILE=file_tg, EXIST=LCONF_H5)
          IF (LCONF_H5) THEN
              if ( irank_g .eq. 0 ) write (*,*) "read input walkers"
              call wavefunction_in_hdf5( phi_0, file_tg )
          endif

          !! initial overlap and green's function
          do i_wlk = 1, N_wlk

             if (weight_k(i_wlk) .le. 0.d0 ) weight_k(i_wlk) = 0.d0

             do ns = 1, N_slat
                 i_grc = ns+(i_wlk-1)*N_slat

                 do nf_eff = 1,N_Fl_eff
                    nf=Calc_Fl_map(nf_eff)
                    call cgrp(Z, GR(:,:,nf,i_grc), phi_0(nf_eff,i_wlk), phi_trial(nf_eff,ns))
                    det_vec(nf_eff) = Z
                 enddo
                 det_Vec(:) = det_Vec(:) * N_SUN
                 if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
                 IF ( .not. LCONF_H5) overlap(i_grc) = sum(det_vec)
             enddo
          enddo

          !! rescale overlap
          call rescale_overlap(overlap)

          !! initial energy
          call ham%update_fac_norm(GR, 0)

          file_seeds = "seedvec_in"
          INQUIRE (FILE=file_seeds, EXIST=LCONF)
          IF (LCONF) THEN
              if ( irank_g .eq. 0 ) write (*,*) "read input seeds"
              call seed_vec_in(file_seeds)
          endif

        END SUBROUTINE initial_wlk

        SUBROUTINE stepwlk_move( phi_trial, phi_0, GR, ntau_bp ) 
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(IN)    :: phi_trial
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:,:), Allocatable, INTENT(INOUT) :: GR
          Integer, INTENT(IN) :: ntau_bp

          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, N_op
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8

          N_op = Size(OP_V,1)

          do i_wlk = 1, N_wlk
             ! update weight by fac_norm
             if ( weight_k(i_wlk) .gt. Zero ) weight_k(i_wlk)=weight_k(i_wlk)*exp(fac_norm)
          enddo

          call half_K_propagation( phi_trial, phi_0, GR )

          !! propagate with interaction
          call int_propagation( phi_0, GR, ntau_bp )
             
          !! Kinetic part exp(-/Delta/tau T/2)
          call half_K_propagation( phi_trial, phi_0, GR )

          !! rescale overlap after each step
          call rescale_overlap(overlap)

        END SUBROUTINE stepwlk_move

        subroutine int_propagation( phi_0, GR, ntau_bp )
          
          Implicit none
     
          class(udv_state), intent(inout), allocatable :: phi_0    (:,:)
          complex (Kind=Kind(0.d0)), intent(inout), allocatable :: GR(:,:,:,:)
          integer, intent(in) :: ntau_bp

          ! local
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, i_st, i_ed
          Integer :: n_op, ns, i_grc
          Complex (Kind=Kind(0.d0)) :: Z, z_alpha
          COMPLEX (Kind=Kind(0.d0)) :: det_Vec(N_FL)
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8, spin
          
          N_op = Size(OP_V,1)
          
          do i_wlk = 1, N_wlk
          
          if ( weight_k(i_wlk) .gt. Zero ) then

             ! upgrade Green's function
             z_alpha = 0.d0
             do n = 1, N_op
                
                N_type = 2
                do ns = 1, N_slat
                   i_grc = ns+(i_wlk-1)*N_slat
                   do nf_eff = 1,N_FL_eff
                      nf=Calc_Fl_map(nf_eff)
                      Call Op_Wrapdo( GR(:,:,nf,i_grc), Op_V(n,nf), 1.d0, Ndim, N_Type, 1)
                   enddo
                enddo

                call Upgrade(GR, n, spin, i_wlk )
                nsigma_bp(i_wlk)%f(n,ntau_bp) = spin

                N_type = 2
                do ns = 1, N_slat
                   i_grc = ns+(i_wlk-1)*N_slat
                   do nf_eff = 1,N_FL_eff
                      nf=Calc_Fl_map(nf_eff)
                      Call Op_Wrapup( Gr(:,:,nf,i_grc), Op_V(n,nf), 1.d0, Ndim, N_Type, 1 )
                   enddo
                enddo

                ! propagate slater determinant
                Do nf_eff = 1,N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    Call Op_mmultR(phi_0(nf_eff,i_wlk)%U,Op_V(n,nf),spin,'n',1)
                    ! store alpha factor in z_alpha
                    z_alpha = z_alpha + &
                        & op_v(n_op,nf)%g*op_v(n_op,nf)%alpha*nsigma_bp(i_wlk)%Phi(n,ntau_bp)
                enddo
             
             enddo

             ! rescale U matrix with z_alpha factor
             Do nf_eff = 1,N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 phi_0(nf_eff,i_wlk)%U(:,:) = exp(z_alpha)*phi_0(nf_eff,i_wlk)%U(:,:) 
             enddo

          endif

          end do
        
        end subroutine int_propagation
        
        subroutine half_K_propagation( phi_trial, phi_0, GR )
          
          implicit none
     
          class(udv_state), intent(in   ), allocatable :: phi_trial(:,:)
          class(udv_state), intent(inout), allocatable :: phi_0    (:,:)
          complex (Kind=Kind(0.d0)), intent(inout), allocatable :: GR(:,:,:,:)
          
          ! local
          integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, i_st, i_ed, ns, i_grc
          complex (Kind=Kind(0.d0)) :: Overlap_old, Overlap_new, Z, sum_o_new, sum_o_old
          complex (Kind=Kind(0.d0)) :: det_Vec(N_FL), log_o_new(n_slat), log_o_old(n_slat)
          real    (Kind=Kind(0.d0)) :: overlap_ratio, re_overlap, re_o_max
          real    (Kind=Kind(0.d0)) :: zero = 1.0E-8

          do i_wlk = 1, N_wlk
          
          if ( weight_k(i_wlk) .gt. zero ) then

             re_o_max = 0.d0
             
             do ns = 1, N_slat
                i_grc = ns+(i_wlk-1)*N_slat
                ! Update Green's function
                do nf_eff = 1,N_Fl_eff
                   nf=Calc_Fl_map(nf_eff)
                   call CGRP(Z, GR(:,:,nf,i_grc), phi_0(nf_eff,i_wlk), phi_trial(nf_eff,ns))
                   det_vec(nf) = Z
                enddo
                det_Vec(:) = det_Vec(:) * N_SUN
                if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
                log_o_old(ns) = sum(det_Vec)
                if ( dble(log_o_old(ns)) .gt. re_o_max ) re_o_max = dble(log_o_old(ns))
             enddo

             !! multi exp(-\Delta\tau T/2)
             Do nf_eff = 1,N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                Call Hop_mod_mmthr_1D2(phi_0(nf_eff,i_wlk)%U,nf,1)
             enddo

             do ns = 1, N_slat
                i_grc = ns+(i_wlk-1)*N_slat
                ! Update Green's function
                do nf_eff = 1,N_Fl_eff
                   nf=Calc_Fl_map(nf_eff)
                   call CGRP(Z, GR(:,:,nf,i_grc), phi_0(nf_eff,i_wlk), phi_trial(nf_eff,ns))
                   det_vec(nf) = Z
                enddo
                det_Vec(:) = det_Vec(:) * N_SUN
                if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
                log_o_new(ns) = sum(det_Vec)
                if ( dble(log_o_new(ns)) .gt. re_o_max ) re_o_max = dble(log_o_new(ns))
             enddo

             sum_o_old = 0.d0
             sum_o_new = 0.d0
             do ns = 1, N_slat
                sum_o_old = sum_o_old + exp(log_o_old(ns)-cmplx(re_o_max,0.d0,kind(0.d0)))
                sum_o_new = sum_o_new + exp(log_o_new(ns)-cmplx(re_o_max,0.d0,kind(0.d0)))
             enddo

             overlap_ratio = sum_o_new/sum_o_old
             re_overlap    = dble( overlap_ratio )
             if ( re_overlap .gt. zero ) then
                 !! upgrade overlap
                 do ns = 1, N_slat
                    i_grc = ns+(i_wlk-1)*N_slat
                    overlap(i_grc) = overlap(i_grc) + &
                        & ( log_o_new(ns) - log_o_old(ns) )
                 enddo
                 weight_k(i_wlk) = weight_k(i_wlk)*re_overlap
             else
                 weight_k(i_wlk) = 0.d0
             endif

          endif

          enddo

        end subroutine half_K_propagation

        subroutine re_orthonormalize_walkers(Phi_0, cop)
          
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0
          Character, Intent(IN)    :: cop
          
          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, N_size, I, i_wlk
          Integer :: ndistance, i_wlk_eff, i_st, i_ed, n_wlk_eff, ns, nrs, i_slat
          real    (Kind=Kind(0.d0)) :: Overlap_ratio, Zero = 1.0E-8, pi=acos(-1.d0)
          real    (Kind=Kind(0.d0)) :: re_o_am, re_o_ph
          Complex (Kind=Kind(0.d0)) :: det_D(N_FL)

          n_wlk_eff = size(phi_0,2)
          ndistance = n_wlk_eff/n_wlk

          do i_wlk_eff = 1, n_wlk_eff

             i_wlk = (i_wlk_eff-1)/ndistance+1     ! index for walker
             ns    = i_wlk_eff-(i_wlk-1)*ndistance ! index for slater det

             if ( weight_k(i_wlk) .gt. Zero ) then

                 N_size = phi_0(1,1)%n_part

                 !Carry out U,D,V decomposition.
                 Do nf_eff = 1, N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    CALL Phi_0(nf_eff,i_wlk_eff)%decompose
                 enddo
                 
                 Det_D = cmplx(0.d0, 0.d0, kind(0.d0))

                 Do nf_eff = 1, N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    DO I = 1,N_size
                       det_D(nf)=det_D(nf)+log(Phi_0(nf_eff,i_wlk_eff)%D(I))
                    ENDDO
                 enddo
                 
                 if (reconstruction_needed) call ham%weight_reconstruction(Det_D)

                 Do nf_eff = 1, N_FL_eff
                    Phi_0(nf_eff,i_wlk_eff)%D(:) = cmplx(1.d0, 0.d0, kind(0.d0))
                 enddo
                 
                 !!!update the overlap when normal propagation
                 !if (cop == 'U') then
                 !    i_slat = (i_wlk_eff-1)*N_slat
                 !    do nrs = 1, N_slat
                 !       i_slat = i_slat + 1
                 !       overlap(i_slat)=overlap(i_slat)-sum(Det_D)
                 !       re_o_am = dble( overlap(i_slat) )
                 !       re_o_ph = mod(aimag( overlap(i_slat) ), 2*pi)
                 !       overlap(i_slat) = dcmplx(re_o_am, re_o_ph)
                 !    enddo
                 !endif

              endif

          enddo

        END SUBROUTINE re_orthonormalize_walkers
        
        SUBROUTINE population_control( phi_0, phi_bp_r ) 
#ifdef MPI
          Use mpi
#endif
          Use Random_wrap
        
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0, phi_bp_r
          
          !Local 
          integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, it_wlk, n_exc,pop_exc(N_wlk_mpi,4)
          integer :: j, it, i_t, i_st, i_ed, nu_wlk, i_src, i_wlk, j_src, j_wlk, n1, n2, nrg, nfrg, ilabel, ncslat
          integer :: i_llim, i_rlim, j_llim, j_rlim
          real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8, d_scal, sum_w, w_count, w_tmp(N_wlk_mpi), weight_mpi(N_wlk_mpi)
          real    (Kind=Kind(0.d0)) :: exp_o_abs(n_slat), exp_o_phase(n_slat), dz2
          complex (Kind=Kind(0.d0)) :: overlap_tmp(N_grc)
          complex (Kind=Kind(0.d0)) :: Z1,Z2,Z3, Z_s_array(N_slat), Z_r_array(N_slat),zp
          Type (Fields)   , dimension(:)  , allocatable :: nsigma_store
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE :: phi_0_m, phi_bp_m

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          !! temporary store
          !! store slater determinant
          pop_exc=0
          allocate(phi_0_m (N_FL_eff,N_wlk))
          allocate(phi_bp_m(N_FL_eff,N_wlk))
          do i_wlk = 1, N_wlk
            do nf_eff = 1, N_FL_eff
               nf=Calc_Fl_map(nf_eff)
               CALL phi_0_m (nf_eff, i_wlk)%alloc(ndim,phi_0(1,1)%n_part)
               CALL phi_bp_m(nf_eff, i_wlk)%alloc(ndim,phi_0(1,1)%n_part)
            enddo
          enddo

          do nf_eff = 1, N_FL_eff
             do i_wlk = 1, N_wlk
                phi_0_m (nf_eff,i_wlk)=phi_0   (nf_eff,i_wlk)
                phi_bp_m(nf_eff,i_wlk)=phi_bp_r(nf_eff,i_wlk)
             enddo
          enddo

          !! store fields
          n1 = size(nsigma_bp(1)%f,1)
          n2 = size(nsigma_bp(1)%f,2)
          allocate(nsigma_store(N_wlk))
          do i_wlk =1, N_wlk
              call nsigma_store(i_wlk)%make(n1, n2)
              nsigma_store(i_wlk)%f = nsigma_bp(i_wlk)%f
              nsigma_store(i_wlk)%t = nsigma_bp(i_wlk)%t
          enddo

          !! store overlap
          overlap_tmp(:)=overlap(:)

          ! population control
          weight_mpi(:)=0.d0
          i_st=irank_g*N_wlk+1
          i_ed=(irank_g+1)*N_wlk

          weight_mpi(i_st:i_ed)=weight_k(:)
          
          CALL MPI_REDUCE(weight_mpi,w_tmp,N_wlk_mpi,MPI_REAL8,MPI_SUM, 0,Group_comm,IERR)
          if (irank_g == 0) weight_mpi=w_tmp
          CALL MPI_BCAST (weight_mpi, N_wlk_mpi, MPI_REAL8, 0,MPI_COMM_WORLD,ierr)
          
          d_scal=dble(N_wlk_mpi)/sum(weight_mpi)
          if (irank_g == 0) sum_w=-ranf_wrap()
          CALL MPI_BCAST (sum_w, 1, MPI_REAL8, 0,MPI_COMM_WORLD,ierr)
          nu_wlk=0

          n_exc=0
          do it_wlk=1, N_wlk_mpi
              i_src=(it_wlk-1)/N_wlk
              i_wlk=it_wlk-N_wlk*i_src

              sum_w=sum_w+weight_mpi(it_wlk)*d_scal
              n=ceiling(sum_w);
              do j=(nu_wlk+1), n
                  j_src=(j-1)/N_wlk
                  j_wlk=j-N_wlk*j_src
                  
                  if ( j_src .ne. i_src ) then
                      n_exc=n_exc+1
                      pop_exc(n_exc,1)=i_src
                      pop_exc(n_exc,2)=i_wlk
                      pop_exc(n_exc,3)=j_src
                      pop_exc(n_exc,4)=j_wlk
                  else
                      if ( irank_g .eq. i_src ) then
                        do nf_eff = 1, N_FL_eff
                            phi_0_m (nf_eff,j_wlk)=phi_0   (nf_eff,i_wlk)
                            phi_bp_m(nf_eff,j_wlk)=phi_bp_r(nf_eff,i_wlk)
                        enddo
                        nsigma_store(j_wlk)%f=nsigma_bp(i_wlk)%f

                        i_llim = 1+(i_wlk-1)*N_slat; i_rlim = i_wlk*N_slat
                        j_llim = 1+(j_wlk-1)*N_slat; j_rlim = j_wlk*N_slat
                        overlap_tmp(j_llim:j_rlim)=overlap(i_llim:i_rlim) 
                      endif
                  endif
              enddo
              nu_wlk=n
          enddo

          nrg = 2+N_fl_eff*6

          do it=1, n_exc
             i_src = pop_exc(it,1); i_wlk = pop_exc(it,2)
             j_src = pop_exc(it,3); j_wlk = pop_exc(it,4)
             if ( irank_g .eq. i_src ) then
                i_llim = 1+(i_wlk-1)*N_slat; i_rlim = i_wlk*N_slat
                Z_s_array(:) = overlap(i_llim:i_rlim)
                
                ilabel = (it-1)*nrg
                call mpi_send(Z_s_array,N_slat,MPI_COMPLEX16,j_src,ilabel,Group_comm,IERR)
                
                do nf_eff = 1, N_FL_eff
                    ilabel = (it-1)*nrg+(nf_eff-1)*6+1
                    call phi_0   (nf_eff,i_wlk)%MPI_send_general(j_src,ilabel, ierr)
                    ilabel = (it-1)*nrg+(nf_eff-1)*6+4
                    call phi_bp_r(nf_eff,i_wlk)%MPI_send_general(j_src,ilabel, ierr)
                enddo
                
                ilabel = (it-1)*nrg+(n_fl_eff)*6+1
                call mpi_send(nsigma_bp(i_wlk)%f,n1*n2,MPI_REAL8,j_src,ilabel,Group_comm,IERR)
             endif
             if ( irank_g .eq. j_src ) then
                ilabel = (it-1)*nrg
                call mpi_recv(Z_r_array,N_slat,MPI_COMPLEX16,i_src,ilabel,Group_comm,STATUS,IERR)
                
                j_llim = 1+(j_wlk-1)*N_slat; j_rlim = j_wlk*N_slat
                overlap_tmp(j_llim:j_rlim) = Z_r_array(:) 
                
                do nf_eff = 1, N_FL_eff
                    ilabel = (it-1)*nrg+(nf_eff-1)*6+1
                    call phi_0_m (nf_eff,j_wlk)%MPI_recv_general(i_src, ilabel, status, ierr)
                    ilabel = (it-1)*nrg+(nf_eff-1)*6+4
                    call phi_bp_m(nf_eff,j_wlk)%MPI_recv_general(i_src, ilabel, status, ierr)
                enddo
                
                ilabel = (it-1)*nrg+(n_fl_eff)*6+1
                call mpi_recv(nsigma_store(j_wlk)%f,n1*n2,MPI_REAL8,i_src,ilabel,Group_comm,STATUS,IERR)
             endif
          enddo

          overlap(:)=overlap_tmp(:)
          ! reset weight
          weight_k(:)=1.d0
          do nf_eff = 1, N_FL_eff
             do i_wlk = 1, N_wlk
                phi_0   (nf_eff,i_wlk)=phi_0_m (nf_eff,i_wlk)
                phi_bp_r(nf_eff,i_wlk)=phi_bp_m(nf_eff,i_wlk)
             enddo
          enddo

          !!! rescale overlap
          !call rescale_overlap(overlap)
          
          do i_wlk =1, N_wlk
              nsigma_bp(i_wlk)%f = nsigma_store(i_wlk)%f
              nsigma_bp(i_wlk)%t = nsigma_store(i_wlk)%t
          enddo

          ! deallocate tmp udv class
          do i_wlk = 1, N_wlk
            do nf_eff = 1, N_FL_eff
               nf=Calc_Fl_map(nf_eff)
               CALL phi_0_m (nf_eff, i_wlk)%dealloc
               CALL phi_bp_m(nf_eff, i_wlk)%dealloc
            enddo
          enddo

          do i_wlk =1, N_wlk
              call nsigma_store(i_wlk)%clear
          enddo

          deallocate(phi_0_m, phi_bp_m)
          deallocate(nsigma_store)

        END SUBROUTINE population_control
        
        SUBROUTINE store_phi( phi_0, phi_bp_r )
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0, phi_bp_r

          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk
          
          do nf_eff = 1, N_FL_eff
             do i_wlk = 1, N_wlk
                phi_bp_r(nf_eff,i_wlk)=phi_0(nf_eff,i_wlk)
                !! phi_bp_r is not used in propagation
                call phi_bp_r(nf_eff,i_wlk)%decompose
             enddo
          enddo

        END SUBROUTINE store_phi

        SUBROUTINE backpropagation( GR_mix, phi_bp_l, phi_bp_r, udvst, Stab_nt, ltau, lmetropolis, nsweep, nwarmup )
#ifdef MPI
          Use mpi
#endif
          Implicit none
     
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:,:), Allocatable, INTENT(INOUT) :: GR_mix
          CLASS(UDV_State), Dimension(:,:)  , ALLOCATABLE, INTENT(INOUT) :: phi_bp_l, phi_bp_r
          CLASS(UDV_State), Dimension(:,:,:), ALLOCATABLE, INTENT(INOUT) :: udvst
          INTEGER, dimension(:), allocatable,  INTENT(IN) :: Stab_nt
          Integer, INTENT(IN) :: ltau, lmetropolis, nsweep, nwarmup

          !Local 
          Complex (Kind=Kind(0.d0)) :: GR_bp(NDIM,NDIM,N_FL,N_grc)
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, ltrot_bp, N_op, nstm, nst, ntau
          Integer :: i_grc, i_st, i_ed, act_mea, ns
          Complex (Kind=Kind(0.d0)) :: z, z_weight, z_sum_overlap, exp_overlap(N_slat)
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
          Real (Kind=Kind(0.d0))    :: Zero = 1.0E-8

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          N_op     = Size(OP_V,1)
          ltrot_bp = size(nsigma_bp(1)%f, 2)
          nstm     = Size(udvst, 1)
        
          !! initialization
          do i_grc = 1, N_grc
             i_wlk = (i_grc-1)/N_slat+1     ! index for walker
             ns    = i_grc-(i_wlk-1)*N_slat ! index for slater det
             do nf_eff = 1, N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                CALL phi_bp_l   (nf_eff, i_grc)%reset('l',wf_l(nf,ns)%P)
                CALL udvst(nstm, nf_eff, i_grc)%reset('l',wf_l(nf,ns)%P)
             enddo
          enddo

          !! backpropagation
          nst= nstm-1
          Do nt = ltrot_bp, 1, -1
             ntau = nt-1
             
             Do i_wlk = 1, N_wlk

                if ( weight_k(i_wlk) .gt. Zero  ) then
                
                Do ns = 1, N_slat
                   i_grc = ns+(i_wlk-1)*N_slat
                   Do nf_eff = 1,N_FL_eff
                      
                      nf=Calc_Fl_map(nf_eff)
                      Call Hop_mod_mmthlc_1D2(phi_bp_l(nf_eff,i_grc)%U,nf,1)

                      Do n = N_op,1,-1
                         Call Op_mmultR(phi_bp_l(nf_eff,i_grc)%U,Op_V(n,nf),nsigma_bp(i_wlk)%f(n,nt),'c',1)
                      enddo
                   
                      Call Hop_mod_mmthlc_1D2(phi_bp_l(nf_eff,i_grc)%U,nf,1)

                   enddo
                enddo

                endif

             Enddo
             
             if ( ntau .eq. stab_nt(nst) .and. ntau .ne. 0 ) then
                 call re_orthonormalize_walkers(phi_bp_l, 'N')
                 Do i_grc = 1, N_grc
                 Do nf_eff = 1,N_FL_eff
                    udvst(nst, nf_eff, i_grc) = phi_bp_l(nf_eff, i_grc)
                 ENDDO
                 ENDDO
                 nst = nst - 1
             endif
             
          Enddo
          
          !! svd at tau = 0
          call re_orthonormalize_walkers(phi_bp_l, 'N')

          !! compute the total weight
          call ham%sum_weight(z_weight)

          !! equal time measurement
          act_mea = 0 + irank
          do i_wlk = 1, N_wlk
             i_st = 1+(i_wlk-1)*N_slat; i_ed = i_wlk*N_slat
             exp_overlap(:) = exp(overlap(i_st:i_ed))
             z_sum_overlap = sum(exp_overlap(:))
             Do ns = 1, N_slat
                i_grc = ns+(i_wlk-1)*N_slat
                do nf_eff = 1,N_Fl_eff
                   nf=Calc_Fl_map(nf_eff)
                   call CGRP(Z, GR_bp(:,:,nf,i_grc), phi_bp_r(nf_eff,i_wlk), phi_bp_l(nf_eff,i_grc))
                enddo
                If (reconstruction_needed) Call ham%GR_reconstruction( GR_bp (:,:,:,i_grc) )
                If (reconstruction_needed) Call ham%GR_reconstruction( GR_mix(:,:,:,i_grc) )
                CALL ham%Obser( GR_bp(:,:,:,i_grc), GR_mix(:,:,:,i_grc), i_wlk, i_grc, z_weight, z_sum_overlap, act_mea )
                act_mea = act_mea + 1
             enddo
          enddo

          if ( lmetropolis .eq. 1 ) then
              ! metripolis sampling
              call metropolis(phi_bp_r, phi_bp_l, udvst, stab_nt, nsweep, nwarmup, ltau)

              if ( ltau .eq. 1 ) then
                 act_mea = 0 + irank
                 do i_wlk = 1, N_wlk
                    call ham%bp_obsert(i_wlk, z_weight, act_mea)
                    act_mea = act_mea + 1
                 enddo
              endif
          else
              !! time dependence measurement if not call metropolis
              if ( ltau .eq. 1 ) call bp_measure_tau(phi_bp_l, phi_bp_r, udvst, stab_nt )
          endif

        end subroutine backpropagation

        subroutine metropolis(phi_bp_r, phi_l_m, udvst, stab_nt, nsweep, nwarmup, ltau)
          
          Implicit none
     
          class(udv_state), dimension(:,:)  , allocatable, intent(in)    :: phi_bp_r
          class(udv_state), dimension(:,:)  , allocatable, intent(inout) :: phi_l_m
          class(udv_state), dimension(:,:,:), allocatable, intent(inout) :: udvst
          integer, dimension(:), allocatable, intent(in) :: stab_nt
          integer, intent(in) :: nsweep, nwarmup, ltau

          !Local 
          integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, N_op, ntau, I, nst, nstm
          integer :: nsw, ltrot_bp, nmea, ltrot_eff, i_slat, ns, i_grc
          Complex (Kind=Kind(0.d0)) :: z, z_weight, detz, z1, z2, zp, ztmp, z_avg, z_sgn_avg, ener_tmp
          real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, overlap_ratio, hs_field
          real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8
          complex (Kind=Kind(0.d0)) :: gr(ndim,ndim,n_fl)
          complex (Kind=Kind(0.d0)) :: det_vec(n_fl), zph1, zph2
          class(udv_state), dimension(:,:), allocatable :: phi_r_m
          complex(Kind=Kind(0.d0)), dimension(:,:,:,:), allocatable :: gtt
          complex(Kind=Kind(0.d0)), dimension(:), allocatable :: overlap_mc

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)
          call mpi_comm_size(mpi_comm_world,isize,ierr)
          call mpi_comm_rank(mpi_comm_world,irank,ierr)
          call mpi_comm_rank(group_comm, irank_g, ierr)
          call mpi_comm_size(group_comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          !! initialization
          n_op     = size(op_v , 1)
          nstm     = size(udvst, 1)
          ltrot_bp = size(nsigma_bp(1)%f, 2)
          ltrot_eff = ltrot

          allocate(overlap_mc(n_grc))
          overlap_mc(:) = overlap(:)

          !! init observables
          if ( ltau .eq. 1 ) call ham%init_obs_mc
          
          !! allocate tmp wavefunction
          allocate(phi_r_m(N_FL_eff,N_wlk))
          allocate(gtt(ndim,ndim,n_fl,n_grc))

          !! compute the total weight
          call ham%sum_weight(z_weight)

          do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             do i_wlk = 1, N_wlk
                call phi_r_m(nf_eff, i_wlk)%init(ndim,'r',phi_bp_r(nf_eff,i_wlk)%U)
             enddo
             do i_grc = 1, N_grc
                i_wlk  = (i_grc-1)/N_slat+1      ! index for walker
                i_slat = i_grc-(i_wlk-1)*N_slat  ! index for slater det
                !! init Green's function
                call cgrp(z, gtt(:,:,nf,i_grc), phi_r_m(nf_eff,i_wlk), phi_l_m(nf_eff,i_grc))
             enddo
          enddo

          do nsw = 1, nsweep

          nst=1
          do ntau = 1, ltrot_bp
 
             do i_wlk = 1, N_wlk

             if ( weight_k(i_wlk) .gt. zero ) then

                 !! Propagate Green's function
                 do ns = 1, N_slat
                    i_grc = ns+(i_wlk-1)*N_slat
                    do nf_eff = 1,N_FL_eff
                       nf=Calc_Fl_map(nf_eff)
                       call hop_mod_mmthr_1d2   (gtt(:,:,nf,i_grc),nf,1)
                       call hop_mod_mmthl_m1_1d2(gtt(:,:,nf,i_grc),nf,1)
                    enddo
                 enddo

                 do n = 1, n_op
                    
                    n_type = 1
                    do ns = 1, N_slat
                       i_grc = ns+(i_wlk-1)*N_slat
                       do nf_eff = 1, N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          hs_field = nsigma_bp(i_wlk)%f(n,ntau) 
                          call op_wrapup(gtt(:,:,nf,i_grc),op_v(n,nf),hs_field,ndim,n_type,ntau)
                       enddo
                    enddo
                     
                    ! metropolis update
                    hs_new = nsigma_bp(i_wlk)%flip(n,ntau)
                    if ( ntau .le. ltrot_eff ) then
                        call upgrade_mc(gtt, n, ntau, hs_new, i_wlk, overlap_mc )
                    endif
       
                    n_type = 2
                    do ns = 1, N_slat
                       i_grc = ns+(i_wlk-1)*N_slat
                       do nf_eff = 1,N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          hs_field = nsigma_bp(i_wlk)%f(n,ntau) 
                          call op_wrapup(gtt(:,:,nf,i_grc),op_v(n,nf),hs_field,ndim,n_type,ntau)
                       enddo
                    enddo

                 enddo

                 do ns = 1, N_slat
                    i_grc = ns+(i_wlk-1)*N_slat
                    do nf_eff = 1,N_FL_eff
                       nf=Calc_Fl_map(nf_eff)
                       call hop_mod_mmthr_1d2   (gtt(:,:,nf,i_grc),nf,1)
                       call hop_mod_mmthl_m1_1d2(gtt(:,:,nf,i_grc),nf,1)
                    enddo
                 enddo
                 
                 !! Propagate wave function
                 do nf_eff = 1,N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    call hop_mod_mmthr_1d2(phi_r_m(nf_eff,i_wlk)%U,nf,1)
                    do n = 1, N_op
                       call op_mmultr(phi_r_m(nf_eff,i_wlk)%U,op_v(n,nf),nsigma_bp(i_wlk)%f(n,ntau),'n',1)
                    enddo
                    call hop_mod_mmthr_1d2(phi_r_m(nf_eff,i_wlk)%U,nf,1)
                 enddo

             endif

             enddo

             !! call svd
             if (  ntau .eq. stab_nt(nst) )  then
                 call re_orthonormalize_walkers(phi_r_m, 'N')

                 do i_wlk = 1, N_wlk

                 if ( weight_k(i_wlk) .gt. zero ) then

                     do ns = 1, N_slat
                        i_grc = ns+(i_wlk-1)*N_slat
                        do nf_eff = 1, N_FL_eff
                           nf=Calc_Fl_map(nf_eff)
                           phi_l_m(nf_eff, i_grc) = udvst(nst, nf_eff, i_grc)
                           !! compute Green's function
                           call cgrp(z, gr(:,:,nf), phi_r_m(nf_eff,i_wlk), phi_l_m(nf_eff,i_grc))
                           det_vec(nf_eff) = z
                           call control_precisiong(gr(:,:,nf),gtt(:,:,nf,i_grc),Ndim)
                           gtt(:,:,nf,i_grc) = gr(:,:,nf)
                        enddo
                        det_Vec(:) = det_Vec(:) * N_SUN
                        zph1 = exp(dcmplx(0.d0, aimag(sum(det_vec))))
                        zph2 = exp(dcmplx(0.d0, aimag(overlap_mc(i_grc))))
                        if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
                        call control_Precisionp(zph1, zph2)
                        overlap_mc(i_grc) = sum(det_vec)
                     enddo

                     !! store on the first slat det storage
                     i_grc = 1+(i_wlk-1)*N_slat
                     do nf_eff = 1, N_FL_eff
                        nf=Calc_Fl_map(nf_eff)
                        udvst(nst, nf_eff, i_grc) = phi_r_m(nf_eff, i_wlk)
                     enddo

                 endif

                 enddo
                 call rescale_overlap(overlap_mc)
                 nst = nst + 1
             endif

          enddo

          do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             do i_grc = 1, N_grc
                i_wlk  = (i_grc-1)/N_slat+1      ! index for walker
                i_slat = i_grc-(i_wlk-1)*N_slat  ! index for slater det
                call phi_l_m(nf_eff, i_grc)%reset('l',wf_l(nf,i_slat)%P)
             enddo
          enddo

          nst = nstm-1
          do ntau = ltrot_bp, 1,-1
             ntau1 = ntau - 1
          
             do i_wlk = 1, N_wlk

             if ( weight_k(i_wlk) .gt. Zero ) then
             
                 !! Propagate Green's function
                 Do ns = 1, N_slat
                    i_grc = ns+(i_wlk-1)*N_slat
                    do nf_eff = 1,N_FL_eff
                       nf=Calc_Fl_map(nf_eff)
                       Call hop_mod_mmthl_1d2   (gtt(:,:,nf,i_grc),nf,1)
                       Call hop_mod_mmthr_m1_1d2(gtt(:,:,nf,i_grc),nf,1)
                    enddo
                 enddo

                 do n = n_op, 1, -1

                    n_type = 2
                    do ns = 1, N_slat
                       i_grc = ns+(i_wlk-1)*N_slat
                       do nf_eff = 1, N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          hs_field = nsigma_bp(i_wlk)%f(n,ntau) 
                          call op_wrapdo(gtt(:,:,nf,i_grc),op_v(n,nf),hs_field,ndim,N_Type,ntau)
                       enddo
                    enddo
                    
                    ! metropolis update
                    hs_new = nsigma_bp(i_wlk)%flip(n,ntau) 
                    if ( ntau .le. ltrot_eff ) then
                        call upgrade_mc(gtt, n, ntau, hs_new, i_wlk, overlap_mc )
                    endif
       
                    n_type = 1
                    hs_field = nsigma_bp(i_wlk)%f(n,ntau)  
                    do ns = 1, N_slat
                       i_grc = ns+(i_wlk-1)*N_slat
                       do nf_eff = 1,N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          call op_wrapdo(gtt(:,:,nf,i_grc),op_v(n,nf),hs_field,ndim,n_Type,ntau)
                       enddo
                    enddo

                 enddo

                 do ns = 1, N_slat
                    i_grc = ns+(i_wlk-1)*N_slat
                    do nf_eff = 1,N_FL_eff
                       nf=Calc_Fl_map(nf_eff)
                       call hop_mod_mmthl_1d2   (gtt(:,:,nf,i_grc),nf,1)
                       call hop_mod_mmthr_m1_1d2(gtt(:,:,nf,i_grc),nf,1)
                    enddo
                 enddo
                 
                 !! Propagate wave function
                 do ns = 1, N_slat
                    i_grc = ns+(i_wlk-1)*N_slat
                    do nf_eff = 1,N_FL_eff
                       nf=Calc_Fl_map(nf_eff)
                       call hop_mod_mmthlc_1d2(phi_l_m(nf_eff,i_grc)%U,nf,1)
                       do n = N_op,1,-1
                          Call op_mmultr(phi_l_m(nf_eff,i_grc)%U,op_v(n,nf),nsigma_bp(i_wlk)%f(n,ntau),'c',1)
                       enddo
                       call hop_mod_mmthlc_1d2(phi_l_m(nf_eff,i_grc)%U,nf,1)
                    enddo
                 enddo

             endif

             enddo

             !! call svd
             if (  ntau1 .eq. stab_nt(nst) .and. ntau1 .ne. 0 )  then
                 call re_orthonormalize_walkers(phi_l_m, 'N')
          
                 do i_wlk = 1, N_wlk
                    if ( weight_k(i_wlk) .gt. Zero ) then
                    
                    !! read phi_r from the 1st slat det storage
                    i_grc = 1+(i_wlk-1)*N_slat
                    do nf_eff = 1, N_FL_eff
                       nf=Calc_Fl_map(nf_eff)
                       phi_r_m(nf_eff, i_wlk) = udvst(nst, nf_eff, i_grc)
                    enddo
                    
                    do ns = 1, N_slat
                       i_grc = ns+(i_wlk-1)*N_slat
                       do nf_eff = 1, N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          !! compute Green's function
                          call cgrp(z, gr(:,:,nf), phi_r_m(nf_eff,i_wlk), phi_l_m(nf_eff,i_grc))
                          det_vec(nf_eff) = z
                          call control_precisiong(gr(:,:,nf),gtt(:,:,nf,i_grc),Ndim)
                          gtt(:,:,nf,i_grc) = gr(:,:,nf)
                          udvst(nst, nf_eff, i_grc) = phi_l_m(nf_eff, i_grc)
                       enddo
                       det_Vec(:) = det_Vec(:) * N_SUN
                       zph1 = exp(dcmplx(0.d0, aimag(sum(det_vec))))
                       zph2 = exp(dcmplx(0.d0, aimag(overlap_mc(i_grc))))
                       if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
                       call control_Precisionp(zph1, zph2)
                       overlap_mc(i_grc) = sum(det_vec)
                    enddo

                    endif
                 enddo
                 call rescale_overlap(overlap_mc)
                 nst = nst - 1
             endif

          enddo

          call re_orthonormalize_walkers(phi_l_m, 'N')

          do i_wlk = 1, N_wlk

             if ( weight_k(i_wlk) .gt. zero ) then
             
             do nf_eff = 1, N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                call phi_r_m(nf_eff, i_wlk)%reset('r',phi_bp_r(nf_eff,i_wlk)%U)
             enddo

             do ns = 1, N_slat
                i_grc = ns+(i_wlk-1)*N_slat
                do nf_eff = 1, N_FL_eff
                   nf=Calc_Fl_map(nf_eff)
                   !! compute Green's function
                   call cgrp(z, gr(:,:,nf), phi_r_m(nf_eff,i_wlk), phi_l_m(nf_eff,i_grc))
                   det_vec(nf_eff) = z
                   call control_precisiong(gr(:,:,nf),gtt(:,:,nf,i_grc),Ndim)
                   gtt(:,:,nf,i_grc) = gr(:,:,nf)
                enddo
                det_Vec(:) = det_Vec(:) * N_SUN
                zph1 = exp(dcmplx(0.d0, aimag(sum(det_vec))))
                zph2 = exp(dcmplx(0.d0, aimag(overlap_mc(i_grc))))
                if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
                call control_Precisionp(zph1, zph2)
                overlap_mc(i_grc) = sum(det_vec)
             enddo
             
             endif
          enddo
          call rescale_overlap(overlap_mc)
          
          do ns = 1, N_slat
             do i_wlk = 1, N_wlk
                i_grc = ns+(i_wlk-1)*N_slat
                do nf_eff = 1, N_FL_eff
                   nf=Calc_Fl_map(nf_eff)
                   call udvst(nstm, nf_eff, i_grc)%reset('l',wf_l(nf,ns)%P)
                enddo
             enddo
          enddo

          if ( (ltau .eq. 1) .and. (nsw .gt. nwarmup) ) call mc_measure_dyn( udvst, gtt, phi_bp_r, stab_nt, overlap_mc )

          enddo
          
          ! deallocate tmp udv class
          do i_wlk = 1, N_wlk
             do nf_eff = 1, N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                call phi_r_m(nf_eff, i_wlk)%dealloc
             enddo
          enddo

          deallocate(phi_r_m)
          deallocate(gtt)
          deallocate(overlap_mc)

        end subroutine metropolis

        subroutine mc_measure_dyn( udvst, gr_in, phi_r_in, stab_nt, overlap_in )
#ifdef MPI
          Use mpi
#endif
          Implicit none
     
          complex(Kind=Kind(0.d0)), dimension(:,:,:,:), allocatable, intent(in) :: gr_in
          complex(Kind=Kind(0.d0)), dimension(:), allocatable, intent(in) :: overlap_in
          class(udv_state), dimension(:,:)  , allocatable, intent(in) :: phi_r_in
          class(udv_state), dimension(:,:,:), allocatable, intent(in) :: udvst
          integer, dimension(:), allocatable, intent(in) :: stab_nt

          !Local 
          Complex (Kind=Kind(0.d0)) :: temp(ndim,ndim), gr(ndim,ndim,n_fl), grc(ndim,ndim,n_fl)
          Complex (Kind=Kind(0.d0)) :: gt0(ndim,ndim,n_fl,n_grc), g00(ndim,ndim,n_fl,n_grc)
          Complex (Kind=Kind(0.d0)) :: gtt(ndim,ndim,n_fl,n_grc), g0t(ndim,ndim,n_fl,n_grc)
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, N_op, ntau, I, nst, nstm
          Integer :: ns, i_grc, act_mea, i_st, i_ed, n_part
          Complex (Kind=Kind(0.d0)) :: Z, Z_weight, DETZ, z_sum_overlap, exp_overlap(N_slat)
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8
          
          class(udv_state), dimension(:,:), allocatable :: phi_r_mea, phi_l_mea

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)
          call mpi_comm_size(MPI_COMM_WORLD,ISIZE,IERR)
          call mpi_comm_rank(MPI_COMM_WORLD,IRANK,IERR)
          call mpi_comm_rank(Group_Comm, irank_g, ierr)
          call mpi_comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          !! initialization
          n_op = size(op_v ,1)
          nstm = size(udvst,1)

          !! allocate tmp wavefunction
          allocate(phi_r_mea(N_FL_eff,N_wlk))
          allocate(phi_l_mea(N_FL_eff,N_grc))

          n_part = phi_r_in(1,1)%n_part

          do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             do i_wlk = 1, N_wlk
                call phi_r_mea(nf_eff, i_wlk)%alloc(ndim,n_part)
                phi_r_mea(nf_eff,i_wlk) = phi_r_in(nf_eff,i_wlk)
                do ns = 1, n_slat
                   i_grc = ns+(i_wlk-1)*N_slat
                   call phi_l_mea(nf_eff, i_grc)%alloc(ndim,n_part)
                enddo
             enddo
          enddo

          gtt = gr_in
          g00 = gtt
          gt0 = gtt
          g0t = gtt
          do i_grc = 1, N_grc
          do nf_eff=1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             do I=1,Ndim
                g0t(I,I,nf,i_grc)=g0t(I,I,nf,i_grc)-1.d0
             enddo
          enddo
          enddo
          
          ntau = 0
          !! call reconstruction of non-calculated flavor blocks
          do i_grc = 1, n_grc
             If (reconstruction_needed) then
                 call ham%gr_reconstruction ( g00(:,:,:,i_grc) )
                 call ham%gr_reconstruction ( gtt(:,:,:,i_grc) )
                 call ham%grt_reconstruction( gt0(:,:,:,i_grc), g0t(:,:,:,i_grc) )
             endif
          enddo
          call ham%obsert_mc(ntau, gt0, g0t, g00, gtt, overlap_in)

          NST=1
          do ntau = 1, ltrot

             do i_wlk = 1, N_wlk
                !! Propagate wave function
                if ( weight_k(i_wlk) .gt. zero ) then
                
                do ns = 1, N_slat
                   i_grc = ns+(i_wlk-1)*N_slat

                   !! Propagate Green's function
                   call propr  (gt0(:,:,:,i_grc),ntau,i_wlk)
                   call proprm1(g0t(:,:,:,i_grc),ntau,i_wlk)
                   call proprm1(gtt(:,:,:,i_grc),ntau,i_wlk)
                   call propr  (gtt(:,:,:,i_grc),ntau,i_wlk)
                
                enddo
                   
                Do nf_eff = 1,N_FL_eff
                   nf=Calc_Fl_map(nf_eff)
                   call Hop_mod_mmthr_1D2(phi_r_mea(nf_eff,i_wlk)%U,nf,1)
                   Do n = 1, N_op
                      call Op_mmultR(phi_r_mea(nf_eff,i_wlk)%U,Op_V(n,nf), nsigma_bp(i_wlk)%f(n,ntau),'n',1)
                   enddo
                   call Hop_mod_mmthr_1D2(phi_r_mea(nf_eff,i_wlk)%U,nf,1)
                enddo

                endif

             enddo

             !! call reconstruction of non-calculated flavor blocks
             do i_grc = 1, n_grc
                If (reconstruction_needed) then
                    call ham%gr_reconstruction ( g00(:,:,:,i_grc) )
                    call ham%gr_reconstruction ( gtt(:,:,:,i_grc) )
                    call ham%grt_reconstruction( gt0(:,:,:,i_grc), g0t(:,:,:,i_grc) )
                endif
             enddo
             call ham%obsert_mc(ntau, gt0, g0t, g00, gtt, overlap_in)
             
             !! call svd
             if (  ntau .eq. stab_nt(nst) )  then
                 call re_orthonormalize_walkers(phi_r_mea, 'N')
          
                 do i_wlk = 1, N_wlk
                
                    if ( weight_k(i_wlk) .gt. zero ) then
                 
                    do ns = 1, N_slat
                       i_grc = ns+(i_wlk-1)*N_slat
                       do nf_eff = 1, N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          phi_l_mea(nf_eff, i_grc) = udvst(nst, nf_eff, i_grc) 
                          call cgrp(detz, gr(:,:,nf), phi_r_mea(nf_eff,i_wlk), phi_l_mea(nf_eff,i_grc))
                          call control_precision_tau(gtt(:,:,nf,i_grc), gr(:,:,nf), Ndim)
                       enddo
                       gtt(:,:,:,i_grc) = gr
                    
                       grc = -gr
                       do nf_eff=1,N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          do I=1,Ndim
                             grc(I,I,nf)=grc(I,I,nf)+1.d0
                          enddo
                       enddo
                  
                       do nf_eff=1,N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          call mmult(temp,gr (:,:,nf),gt0(:,:,nf,i_grc))
                          gt0(:,:,nf,i_grc) = temp
                          call mmult(temp,g0t(:,:,nf,i_grc),grc(:,:,nf))
                          g0t(:,:,nf,i_grc) = temp
                       enddo
                    
                    enddo

                    endif

                 enddo
                 nst = nst + 1
             endif

          enddo
          
          ! deallocate tmp udv class
          do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             do i_wlk = 1, N_wlk
                call phi_r_mea(nf_eff, i_wlk)%dealloc
                do ns = 1, n_slat
                   i_grc = ns+(i_wlk-1)*N_slat
                   call phi_l_mea(nf_eff, i_grc)%dealloc
                enddo
             enddo
          enddo
          
          deallocate(phi_r_mea, phi_l_mea)

        end subroutine mc_measure_dyn

        subroutine bp_measure_tau(phi_bp_l, phi_bp_r, udvst, stab_nt )
#ifdef MPI
          Use mpi
#endif
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_bp_l, phi_bp_r
          CLASS(UDV_State), Dimension(:,:,:), ALLOCATABLE, INTENT(IN) :: udvst
          INTEGER, dimension(:)    ,allocatable,  INTENT(IN) :: Stab_nt

          !Local 
          Complex (Kind=Kind(0.d0)) :: Temp(NDIM,NDIM), GR(NDIM,NDIM,N_FL), GRC(NDIM,NDIM,N_FL)
          Complex (Kind=Kind(0.d0)) :: GT0(NDIM,NDIM,N_FL,N_grc), G00(NDIM,NDIM,N_FL,N_grc)
          Complex (Kind=Kind(0.d0)) :: GTT(NDIM,NDIM,N_FL,N_grc), G0T(NDIM,NDIM,N_FL,N_grc)
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, N_op, ntau, I, nst, nstm
          Integer :: ns, i_grc, act_mea, i_st, i_ed
          Complex (Kind=Kind(0.d0)) :: Z, Z_weight, DETZ, z_sum_overlap, exp_overlap(N_slat)
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          !! initialization
          N_op     = Size(OP_V,1)
          nstm     = Size(udvst, 1)
          
          !! compute the total weight
          call ham%sum_weight(z_weight)
          
          do i_grc = 1, N_grc
             i_wlk = (i_grc-1)/N_slat+1     ! index for walker
             ns    = i_grc-(i_wlk-1)*N_slat ! index for slater det
             do nf_eff = 1, N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                CALL CGRP(DetZ, GR(:,:,nf), phi_bp_r(nf_eff,i_wlk), phi_bp_l(nf_eff,i_grc))
             enddo
             GTT(:,:,:,i_grc) = GR
          enddo

          G00 = GTT
          GT0 = GTT
          G0T = GTT
          do i_grc = 1, N_grc
          do nf_eff=1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             do I=1,Ndim
                G0T(I,I,nf,i_grc)=G0T(I,I,nf,i_grc)-1.d0
             enddo
          enddo
          enddo
          
          ntau = 0
          act_mea = 0 + irank
          do i_wlk = 1, N_wlk
             i_st = 1+(i_wlk-1)*N_slat; i_ed = i_wlk*N_slat
             exp_overlap(:) = exp(overlap(i_st:i_ed))
             z_sum_overlap = sum(exp_overlap(:))
             Do ns = 1, N_slat
                i_grc = ns+(i_wlk-1)*N_slat
                !! call reconstruction of non-calculated flavor blocks
                If (reconstruction_needed) then
                    Call ham%GR_reconstruction ( G00(:,:,:,i_grc) )
                    Call ham%GR_reconstruction ( GTT(:,:,:,i_grc) )
                    Call ham%GRT_reconstruction( GT0(:,:,:,i_grc), G0T(:,:,:,i_grc) )
                endif
                CALL ham%obserT( ntau,GT0(:,:,:,i_grc),G0T(:,:,:,i_grc),G00(:,:,:,i_grc), & 
                    & GTT(:,:,:,i_grc), i_wlk, i_grc, z_weight, z_sum_overlap, act_mea )
                act_mea = act_mea + 1
             enddo
          enddo

          NST=1
          do ntau = 1, ltrot

             do i_wlk = 1, N_wlk
                !! Propagate wave function
                if ( weight_k(i_wlk) .gt. zero ) then
                
                Do ns = 1, N_slat
                   i_grc = ns+(i_wlk-1)*N_slat

                   !! Propagate Green's function
                   CALL PROPR  (GT0(:,:,:,i_grc),ntau,i_wlk)
                   CALL PROPRM1(G0T(:,:,:,i_grc),ntau,i_wlk)
                   CALL PROPRM1(GTT(:,:,:,i_grc),ntau,i_wlk)
                   CALL PROPR  (GTT(:,:,:,i_grc),ntau,i_wlk)
                
                enddo
                   
                Do nf_eff = 1,N_FL_eff
                   
                   nf=Calc_Fl_map(nf_eff)
                   call Hop_mod_mmthr_1D2    (phi_bp_r(nf_eff,i_wlk)%U,nf,1)

                   Do n = 1, N_op
                      call Op_mmultR(phi_bp_r(nf_eff,i_wlk)%U,Op_V(n,nf), nsigma_bp(i_wlk)%f(n,ntau),'n',1)
                   enddo
                
                   call Hop_mod_mmthr_1D2    (phi_bp_r(nf_eff,i_wlk)%U,nf,1)

                enddo

                endif

                !! call reconstruction of non-calculated flavor blocks
                i_st = 1+(i_wlk-1)*N_slat; i_ed = i_wlk*N_slat
                exp_overlap(:) = exp(overlap(i_st:i_ed))
                z_sum_overlap = sum(exp_overlap(:))
                do ns = 1, N_slat
                   i_grc = ns+(i_wlk-1)*N_slat
                   If (reconstruction_needed) then
                       Call ham%GR_reconstruction ( G00(:,:,:,i_grc) )
                       Call ham%GR_reconstruction ( GTT(:,:,:,i_grc) )
                       Call ham%GRT_reconstruction( GT0(:,:,:,i_grc), G0T(:,:,:,i_grc) )
                   endif
                   CALL ham%obserT(ntau,GT0(:,:,:,i_grc),G0T(:,:,:,i_grc),G00(:,:,:,i_grc), & 
                       & GTT(:,:,:,i_grc), i_wlk, i_grc, z_weight, z_sum_overlap, act_mea)
                   act_mea = act_mea + 1
                enddo
             enddo
             
             !! call svd
             if (  ntau .eq. stab_nt(nst) )  then
                 call re_orthonormalize_walkers(phi_bp_r, 'N')
          
                 do i_wlk = 1, N_wlk
                
                    if ( weight_k(i_wlk) .gt. zero ) then
                 
                    do ns = 1, N_slat
                       i_grc = ns+(i_wlk-1)*N_slat
                       do nf_eff = 1, N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          phi_bp_l(nf_eff, i_grc) = udvst(nst, nf_eff, i_grc) 
                          CALL CGRP(DetZ, GR(:,:,nf), phi_bp_r(nf_eff,i_wlk), phi_bp_l(nf_eff,i_grc))
                          Call Control_Precision_tau(GTT(:,:,nf,i_grc), GR(:,:,nf), Ndim)
                       enddo
                       GTT(:,:,:,i_grc) = GR
                    
                       GRC = -GR
                       do nf_eff=1,N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          do I=1,Ndim
                             GRC(I,I,nf)=GRC(I,I,nf)+1.d0
                          enddo
                       enddo
                  
                       do nf_eff=1,N_FL_eff
                          nf=Calc_Fl_map(nf_eff)
                          CALL MMULT(TEMP,GR (:,:,nf),GT0(:,:,nf,i_grc))
                          GT0(:,:,nf,i_grc) = TEMP
                          CALL MMULT(TEMP,G0T(:,:,nf,i_grc),GRC(:,:,nf))
                          G0T(:,:,nf,i_grc) = TEMP
                       enddo
                    
                    enddo

                    endif

                 enddo
                 nst = nst + 1
             endif

          enddo

        end subroutine bp_measure_tau
        
        subroutine rescale_overlap(overlap_in)
          
          Implicit none
          
          complex(Kind=Kind(0.d0)), dimension(:), allocatable, intent(inout) :: overlap_in
        
          integer :: nf, nf_eff, n, m, nt, i_wlk, i_grc, ns
          integer :: i_st, i_ed, ncslat
          real    (Kind=Kind(0.d0)) :: log_o_abs(n_slat), log_o_phase(n_slat), dz2
          real    (kind=kind(0.d0)) :: pi = acos(-1.d0), dre_o, zero = 1.0E-8
          complex (Kind=Kind(0.d0)) :: z1, zp

          do i_wlk = 1, N_wlk
                
             if ( weight_k(i_wlk) .gt. zero ) then

             i_st = 1+(i_wlk-1)*N_slat
             i_ed = i_wlk*N_slat

             ncslat = 0
             do i_grc = i_st, i_ed
                ncslat = ncslat + 1
                log_o_abs  (ncslat) = dble( overlap_in(i_grc) )
                log_o_phase(ncslat) = mod( aimag( overlap_in(i_grc) ), 2.d0*pi )
             enddo

             dz2 = maxval(log_o_abs(:))

             ncslat = 0
             do i_grc = i_st, i_ed
                ncslat = ncslat + 1
                dre_o =  log_o_abs(ncslat) - dz2
                overlap_in(i_grc) = cmplx(dre_o,log_o_phase(ncslat), kind(0.d0))
             enddo

             endif
          enddo

        end subroutine rescale_overlap

        SUBROUTINE PROPR(AIN,NT,i_wlk) 

          ! Ain =       B(NT-1, NT1) 
          ! Aout= Ain = B(NT  , NT1)

          Implicit none
          Complex (Kind=Kind(0.D0)), intent(INOUT) :: Ain(Ndim,Ndim,N_FL)
          Integer, INTENT(IN) :: NT, i_wlk

          !Locals
          Integer :: nf,nf_eff,n 

          Do nf_eff = 1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             Call Hop_mod_mmthr_1D2(Ain(:,:,nf),nf,nt)
             Do n = 1,Size(Op_V,1)
                Call Op_mmultR(Ain(:,:,nf),Op_V(n,nf),nsigma_bp(i_wlk)%f(n,nt),'n',nt)
             ENDDO
             Call Hop_mod_mmthr_1D2(Ain(:,:,nf),nf,nt)
          Enddo

        end SUBROUTINE PROPR


        SUBROUTINE PROPRM1(AIN,NT,i_wlk)

          !Ain = B^{-1}(NT-1, NT1) 
          !Aout= B^{-1}(NT  , NT1)

          Implicit none

          !Arguments 
          Complex (Kind=Kind(0.D0)), intent(Inout) ::  AIN(Ndim, Ndim, N_FL) 
          Integer :: NT, i_wlk

          ! Locals 
          Integer :: nf,nf_eff, n 

          do nf_eff = 1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             Call Hop_mod_mmthl_m1_1D2(Ain(:,:,nf),nf,nt)
             Do n =1,Size(Op_V,1)
                Call Op_mmultL(Ain(:,:,nf),Op_V(n,nf),-nsigma_bp(i_wlk)%f(n,nt),'n',nt)
             Enddo
             Call Hop_mod_mmthl_m1_1D2(Ain(:,:,nf),nf,nt)
          enddo

        END SUBROUTINE PROPRM1
 
        SUBROUTINE seed_vec_in(file_tg)
#ifdef MPI
          Use mpi
#endif
          Implicit none
          
          CHARACTER (LEN=64), INTENT(IN)  :: FILE_TG
          
          ! LOCAL
          INTEGER             :: K, ii, intseed(2), seed_tmp(2)
          CHARACTER (LEN=64)  :: filename
          INTEGER,ALLOCATABLE :: SEED_VEC(:), seed_output(:,:), s_tmp(:)

#if defined(MPI)
          INTEGER        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          filename=file_tg
          
          CALL GET_SEED_LEN(K)
          ALLOCATE(SEED_VEC(K))
          ALLOCATE(s_tmp(K))

          if ( irank_g .eq. 0 ) then
              ALLOCATE(seed_output(K,isize_g))
              OPEN (UNIT = 10, FILE=filename, STATUS='OLD', ACTION='Read')
              do ii = 1, isize_g
                  read(10,*) seed_output(:,ii)
              enddo
              close(10)
          ENDIF

          if ( irank_g .ne. 0 ) then
              call mpi_recv(seed_vec, K, mpi_integer, 0, 1, MPI_COMM_WORLD,STATUS,IERR)
              call RANSET(SEED_VEC)
          else
             do ii = 1, isize_g-1
                s_tmp(:)=seed_output(:,ii+1)
                call mpi_send(s_tmp, K, mpi_integer, ii, 1, MPI_COMM_WORLD,IERR)
             ENDDO
          ENDIF
            
          if ( irank_g .eq. 0 ) then
               SEED_VEC=seed_output(:,1)
               call RANSET(SEED_VEC)
               deallocate(seed_output)
          endif

          deallocate(seed_vec)
          deallocate(s_tmp)

        end SUBROUTINE seed_vec_in

        SUBROUTINE seed_vec_out
#ifdef MPI
          Use mpi
#endif
          Implicit none
          
          ! LOCAL
          INTEGER             :: K, ii, intseed(2), seed_tmp(2)
          CHARACTER (LEN=64)  :: FILE_TG, filename
          INTEGER,ALLOCATABLE :: SEED_VEC(:), seed_output(:,:), s_tmp(:)

#if defined(MPI)
          INTEGER        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif
          filename = "seedvec_out"
          
          CALL GET_SEED_LEN(K)
          ALLOCATE(SEED_VEC(K))
          ALLOCATE(s_tmp(K))
          CALL RANGET(SEED_VEC)

          if ( irank_g .eq. 0 ) then
              ALLOCATE(seed_output(K,isize_g))
              seed_output(:,1)=SEED_VEC(:)
          ENDIF

          if ( irank_g .ne. 0 ) then
              call mpi_send(seed_vec, K, mpi_integer, 0, 1, MPI_COMM_WORLD,IERR)
          else
             do ii = 1, isize_g-1
                call mpi_recv(s_tmp, K, mpi_integer, ii, 1, MPI_COMM_WORLD,STATUS,IERR)
                seed_output(:,ii+1) = s_tmp(:)
             ENDDO
          ENDIF
            
          if ( irank_g .eq. 0 ) then
               OPEN (UNIT = 10, FILE=filename, STATUS='UNKNOWN', ACTION='WRITE')
               do ii = 1, isize_g
                   WRITE(10,*) seed_output(:,ii)
               enddo
               close(10)
               deallocate(seed_output)
          endif

          deallocate(seed_vec)
          deallocate(s_tmp)

        end SUBROUTINE seed_vec_out

#if defined HDF5
        SUBROUTINE wavefunction_out_hdf5( phi_0 )
#ifdef MPI
          Use mpi
#endif
#if defined HDF5
          Use hdf5
          use h5lt
#endif
          Implicit none
     
          class(udv_state), Dimension(:,:), ALLOCATABLE, INTENT(IN) :: phi_0

          ! LOCAL
          CHARACTER (LEN=64) :: FILE_TG, filename
          Complex (Kind=Kind(0.d0)), pointer :: phi0_out(:,:,:,:)
          Complex (Kind=Kind(0.d0)), pointer :: overlap_out(:)
          Real    (Kind=Kind(0.d0)), pointer :: weight_out(:)
          Complex (Kind=Kind(0.d0)), allocatable :: otphi_tmp(:), p0_tmp(:,:,:,:), p1_tmp(:,:,:,:)
          Real    (Kind=Kind(0.d0)), allocatable :: wt_tmp(:)

          INTEGER             :: K, hdferr, rank, nf, nw, n_part, i0, i1, i2, i_st, i_ed, Ndt, ii
          INTEGER             :: i_st2, i_ed2
          INTEGER(HSIZE_T), allocatable :: dims(:), dimsc(:)
          Logical             :: file_exists
          INTEGER(HID_T)      :: file_id, crp_list, space_id, dset_id, dataspace
          Character (len=64)  :: dset_name
          TYPE(C_PTR)         :: dat_ptr

#if defined(MPI)
          INTEGER        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif
          filename = "phiout_0"

          write(filename,'(A,A)') trim(filename), ".h5"
            
          n_part=phi_0(1,1)%n_part
          
          if (irank_g .eq. 0 ) then
              allocate(phi0_out(ndim,n_part,n_fl_eff,n_wlk_mpi))
              allocate(weight_out(N_wlk_mpi), overlap_out(N_grc_mpi))
          endif

          allocate(p0_tmp(ndim,n_part,n_fl_eff,n_wlk))
          allocate(p1_tmp(ndim,n_part,n_fl_eff,n_wlk))
          allocate(wt_tmp(N_wlk))
          allocate(otphi_tmp(N_grc))

          do nf = 1, N_FL_eff
          do nw = 1, N_wlk
              p0_tmp(:,:,nf,nw)=phi_0(nf,nw)%U(:,:)
          enddo
          enddo

          if ( irank_g .ne. 0 ) then
              call mpi_send(overlap    ,N_grc,mpi_complex16, 0, 0, MPI_COMM_WORLD,IERR)
              call mpi_send(weight_k   ,N_wlk,    mpi_real8, 0, 1, MPI_COMM_WORLD,IERR)
              Ndt=N_FL_eff*N_wlk*ndim*n_part
              call mpi_send(p0_tmp     ,  Ndt,mpi_complex16, 0, 2, MPI_COMM_WORLD,IERR)
          else
             do ii = 1, isize_g-1
                i_st=ii*N_wlk+1; i_st2=ii*N_grc+1
                i_ed=(ii+1)*N_wlk; i_ed2=(ii+1)*N_grc

                call mpi_recv(otphi_tmp,N_grc,mpi_complex16, ii, 0, MPI_COMM_WORLD,STATUS,IERR)
                overlap_out(i_st2:i_ed2) = otphi_tmp(:)
                call mpi_recv(wt_tmp,N_wlk,    mpi_real8, ii, 1, MPI_COMM_WORLD,STATUS,IERR)
                weight_out(i_st:i_ed) = wt_tmp(:)
                Ndt=N_FL_eff*N_wlk*ndim*n_part
                call mpi_recv(p1_tmp,  Ndt,mpi_complex16, ii, 2, MPI_COMM_WORLD,STATUS,IERR)
                phi0_out(:,:,:,i_st:i_ed)=p1_tmp
             ENDDO
          ENDIF

          if ( irank_g .eq. 0 ) then
              i_st=1; i_st2=1
              i_ed=N_wlk; i_ed2=N_grc
              overlap_out(i_st2:i_ed2) = overlap(:)
              weight_out (i_st:i_ed) = weight_k(:)
              phi0_out(:,:,:,i_st:i_ed)=p0_tmp
          endif

          if ( irank .eq. 0 ) then

          inquire (file=filename, exist=file_exists)
          IF (.not. file_exists) THEN
              CALL h5open_f(ierr)
              CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdferr)
              
              !Create and write dataset for field overlap
              dset_name = "overlap"
              rank = 2
              allocate( dims(2), dimsc(2) )
              dims  = [2, N_grc_mpi]
              dimsc = dims
              CALL h5screate_simple_f(rank, dims, space_id, hdferr)
              CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
              CALL h5pset_chunk_f(crp_list, rank, dimsc, hdferr)
!!#if defined HDF5_ZLIB
!!                ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
!!              CALL h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
!!#endif
              !Create a dataset using cparms creation properties.
              CALL h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, space_id, &
                              dset_id, hdferr, crp_list )
              dat_ptr = C_LOC(overlap_out(1))
              CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !Close objects
              deallocate( dims, dimsc )
              CALL h5sclose_f(space_id, hdferr)
              CALL h5pclose_f(crp_list, hdferr)
              CALL h5dclose_f( dset_id, hdferr)
              
              !Create and write dataset for real part of weight
              dset_name = "weight_re"
              rank = 2
              allocate( dims(2), dimsc(2) )
              dims  = [1, N_wlk_mpi]
              dimsc = dims
              CALL h5screate_simple_f(rank, dims, space_id, hdferr)
              CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
              CALL h5pset_chunk_f(crp_list, rank, dimsc, hdferr)
!!#if defined HDF5_ZLIB
!!                ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
!!              CALL h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
!!#endif
              !Create a dataset using cparms creation properties.
              CALL h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, space_id, &
                              dset_id, hdferr, crp_list )
              dat_ptr = C_LOC(weight_out(1))
              CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !Close objects
              deallocate( dims, dimsc )
              CALL h5sclose_f(space_id, hdferr)
              CALL h5pclose_f(crp_list, hdferr)
              CALL h5dclose_f( dset_id, hdferr)

              !Create and write dataset for wave function
              dset_name = "wavefunction"
              rank = 5
              allocate( dims(5), dimsc(5) )
              dims  = [2,ndim,n_part,N_FL_eff,N_wlk_mpi]
              dimsc = dims
              CALL h5screate_simple_f(rank, dims, space_id, hdferr)
              CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, hdferr)
              CALL h5pset_chunk_f(crp_list, rank, dimsc, hdferr)
#if defined HDF5_ZLIB
                ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
              CALL h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
#endif
              !Create a dataset using cparms creation properties.
              CALL h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, space_id, &
                              dset_id, hdferr, crp_list )
              dat_ptr = C_LOC(phi0_out(1,1,1,1))
              CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !Close objects
              deallocate( dims, dimsc )
              CALL h5sclose_f(space_id,  hdferr)
              CALL h5pclose_f(crp_list,  hdferr)
              CALL h5dclose_f( dset_id,  hdferr)

              !close file
              CALL h5fclose_f(file_id, hdferr)

         else
              !open file
              CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, hdferr)

              !open and write field overlap
              dset_name = "overlap"
              !Open the  dataset.
              CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
              dat_ptr = C_LOC(overlap_out(1))
              !Write data
              CALL H5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !close objects
              CALL h5dclose_f(dset_id,   hdferr)

              !open and write real weight
              dset_name = "weight_re"
              !Open the  dataset.
              CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
              dat_ptr = C_LOC(weight_out(1))
              !Write data
              CALL H5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !close objects
              CALL h5dclose_f(dset_id,   hdferr)
              
              !open and write real weight
              dset_name = "wavefunction"
              !Open the  dataset.
              CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
              dat_ptr = C_LOC(phi0_out(1,1,1,1))
              !Write data
              CALL H5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !close objects
              CALL h5dclose_f(dset_id,   hdferr)
                
              !close file
              CALL h5fclose_f(file_id, hdferr)
         endif

         endif !irank 0

         if (irank_g .eq. 0 ) then
             deallocate(phi0_out, weight_out, overlap_out)
         endif
         deallocate(p0_tmp, p1_tmp, wt_tmp, otphi_tmp)

        END SUBROUTINE wavefunction_out_hdf5

        SUBROUTINE wavefunction_in_hdf5( phi_0, file_tg )
#ifdef MPI
          Use mpi
#endif
#if defined HDF5
          Use hdf5
          use h5lt
#endif
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0
          CHARACTER (LEN=64), INTENT(IN)  :: FILE_TG

          ! LOCAL
          CHARACTER (LEN=64) :: filename
          Complex (Kind=Kind(0.d0)), pointer :: phi0_out(:,:,:,:)
          Complex (Kind=Kind(0.d0)), pointer :: overlap_out(:)
          Real    (Kind=Kind(0.d0)), pointer :: weight_out(:)
          Complex (Kind=Kind(0.d0)), allocatable :: otphi_tmp(:), p0_tmp(:,:,:,:), p1_tmp(:,:,:,:)
          Real    (Kind=Kind(0.d0)), allocatable :: wt_tmp(:)

          INTEGER             :: K, hdferr, rank, nf, nw, n_part, i0, i1, i2, i_st, i_ed, Ndt, ii
          INTEGER             :: nwalk_in, ngrc_in, i_st2, i_ed2
          integer             :: nf_eff
          INTEGER(HSIZE_T), allocatable :: dims(:), dimsc(:), maxdims(:)
          Logical             :: file_exists
          INTEGER(HID_T)      :: file_id, crp_list, space_id, dset_id, dataspace
          Character (len=64)  :: dset_name
          TYPE(C_PTR)         :: dat_ptr

#if defined(MPI)
          INTEGER        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif
          filename = file_tg

          n_part=phi_0(1,1)%n_part
          
          allocate(p0_tmp(ndim,n_part,n_fl_eff,n_wlk))
          allocate(p1_tmp(ndim,n_part,n_fl_eff,n_wlk))
          allocate(wt_tmp(N_wlk))
          allocate(otphi_tmp(N_grc))
          
          if ( irank .eq. 0 ) then

              !open file
              CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, hdferr)

              !open and read overlap
              dset_name = "overlap"
              !Open the  dataset.
              CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
              !Get dataset's dataspace handle.
              CALL h5dget_space_f(dset_id, dataspace, hdferr)
              !Get dataspace's rank.
              CALL h5sget_simple_extent_ndims_f(dataspace, rank, ierr)
              allocate( dims(rank), maxdims(rank) )
              !Get dataspace's dimensions.
              CALL h5sget_simple_extent_dims_f(dataspace, dims, maxdims, ierr)
              ngrc_in=dims(rank)
              !! allocate !!
              allocate(overlap_out(ngrc_in))
              !!-----------!!
              dat_ptr = C_LOC(overlap_out(1))
              CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !close dataspace
              CALL h5sclose_f(dataspace, hdferr)
              !close objects
              CALL h5dclose_f(dset_id,   hdferr)
              deallocate(dims, maxdims)

              !open and read weight
              dset_name = "weight_re"
              !Open the  dataset.
              CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
              !Get dataset's dataspace handle.
              CALL h5dget_space_f(dset_id, dataspace, hdferr)
              !Get dataspace's rank.
              CALL h5sget_simple_extent_ndims_f(dataspace, rank, ierr)
              allocate( dims(rank), maxdims(rank) )
              !Get dataspace's dimensions.
              CALL h5sget_simple_extent_dims_f(dataspace, dims, maxdims, ierr)
              nwalk_in=dims(rank)
              !! allocate !!
              allocate(phi0_out(ndim,n_part,n_fl_eff,nwalk_in))
              allocate(weight_out(nwalk_in))
              !!-----------!!
              dat_ptr = C_LOC(weight_out(1))
              CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !close dataspace
              CALL h5sclose_f(dataspace, hdferr)
              !close objects
              CALL h5dclose_f(dset_id,   hdferr)
              deallocate(dims, maxdims)
              
              !open and read real weight
              dset_name = "wavefunction"
              !Open the  dataset.
              CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
              dat_ptr = C_LOC(phi0_out(1,1,1,1))
              !Write data
              CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !close objects
              CALL h5dclose_f(dset_id,   hdferr)
                
              !close file
              CALL h5fclose_f(file_id, hdferr)

         endif !irank 0

         if ( irank_g .eq. 0 ) then
            do ii = 1, isize_g-1
               i_st=ii*N_wlk+1  ; i_st2=ii*N_grc+1
               i_ed=(ii+1)*N_wlk; i_ed2=(ii+1)*N_grc
               
               otphi_tmp(:)=overlap_out(i_st2:i_ed2)
               call mpi_send(otphi_tmp,N_grc,mpi_complex16, ii, 0, MPI_COMM_WORLD,IERR)
               wt_tmp(:)=weight_out(i_st:i_ed)
               call mpi_send(wt_tmp,N_wlk,    mpi_real8, ii, 1, MPI_COMM_WORLD,IERR)
               Ndt=N_FL_eff*N_wlk*ndim*n_part
               p1_tmp=phi0_out(:,:,:,i_st:i_ed)
               call mpi_send(p1_tmp,  Ndt,mpi_complex16, ii, 2, MPI_COMM_WORLD,IERR)
            ENDDO
         else
            call mpi_recv(otphi_tmp,N_grc,mpi_complex16, 0, 0, MPI_COMM_WORLD,STATUS,IERR)
            overlap(:) = otphi_tmp(:)
            call mpi_recv(wt_tmp,N_wlk,    mpi_real8, 0, 1, MPI_COMM_WORLD,STATUS,IERR)
            weight_k(:)    = wt_tmp(:)
            Ndt=N_FL_eff*N_wlk*ndim*n_part
            call mpi_recv(p0_tmp,  Ndt,mpi_complex16, 0, 2, MPI_COMM_WORLD,STATUS,IERR)
            do nf = 1, N_FL_eff
            do nw = 1, N_wlk
                phi_0(nf,nw)%U(:,:)=p0_tmp(:,:,nf,nw)
            enddo
            enddo
         ENDIF

         if ( irank_g .eq. 0 ) then
             i_st=1    ; i_st2=1
             i_ed=N_wlk; i_ed2=N_grc
             p0_tmp=phi0_out(:,:,:,i_st:i_ed)
             
             weight_k(:) = weight_out(i_st:i_ed)
             overlap (:) = overlap_out(i_st2:i_ed2)
             do nf = 1, N_FL_eff
             do nw = 1, N_wlk
                 phi_0(nf,nw)%U(:,:)=p0_tmp(:,:,nf,nw)
             enddo
             enddo
         endif

         if (irank_g .eq. 0 ) then
             deallocate(phi0_out, weight_out, overlap_out)
         endif
         deallocate(p0_tmp, p1_tmp, wt_tmp, otphi_tmp)

        END SUBROUTINE wavefunction_in_hdf5

        subroutine trial_in_hdf5( phi_0_r, phi_0_l, file_inst, file_antiinst )
#ifdef MPI
          Use mpi
#endif
#if defined HDF5
          Use hdf5
          use h5lt
#endif
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0_r, phi_0_l
          CHARACTER (LEN=64), INTENT(IN)  :: file_inst, file_antiinst 

          ! LOCAL
          CHARACTER (LEN=64) :: filename
          integer, allocatable :: ipiv(:)
          Complex (kind=kind(0.d0)), pointer :: phi0_in(:,:,:)
          Complex (kind=kind(0.d0)), allocatable :: p1_tmp(:,:,:), p2_tmp(:,:,:), log_zdet(:)
          complex (kind=kind(0.d0)), allocatable, dimension(:,:) :: sMat
          complex (kind=kind(0.d0)) :: alpha, beta, zdet, phase, t_overlap, z_norm, c1, ctmp
          real    (kind=kind(0.d0)) :: d_norm

          INTEGER             :: K, hdferr, rank, nf, nw, n_part, i0, i1, i2, i_st, i_ed, Ndt, ii, nwalk_in
          Integer             :: nf_eff, n_fl_in, ns, i_wlk, info, n
          INTEGER(HSIZE_T), allocatable :: dims(:), dimsc(:), maxdims(:)
          Logical             :: file_exists
          INTEGER(HID_T)      :: file_id, crp_list, space_id, dset_id, dataspace
          Character (len=64)  :: dset_name
          TYPE(C_PTR)         :: dat_ptr

#if defined(MPI)
          INTEGER        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif
          n_part=phi_0_l(1,1)%n_part

          n_fl_in = 1
          
          allocate(p1_tmp(ndim,n_part,n_fl))
          allocate(p2_tmp(ndim,n_part,n_fl))
          
          if ( irank .eq. 0 ) then

              !! inst
              !open file
              CALL h5fopen_f (file_inst, H5F_ACC_RDWR_F, file_id, hdferr)

              !! allocate !!
              allocate(phi0_in(ndim,n_part,n_fl_in))
              !!-----------!!
              !open and read real weight
              dset_name = "wavefunction"
              !Open the  dataset.
              CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
              dat_ptr = C_LOC(phi0_in(1,1,1))
              !Write data
              CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !close objects
              CALL h5dclose_f(dset_id,   hdferr)
                
              !close file
              CALL h5fclose_f(file_id, hdferr)

              !! store input
              do nf_eff = 1, N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 p1_tmp(:,:,nf)=phi0_in(:,:,1)
              enddo

              !! anti-inst
              !open file
              CALL h5fopen_f (file_antiinst, H5F_ACC_RDWR_F, file_id, hdferr)

              !!-----------!!
              !open and read real weight
              dset_name = "wavefunction"
              !Open the  dataset.
              CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
              dat_ptr = C_LOC(phi0_in(1,1,1))
              !Write data
              CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !close objects
              CALL h5dclose_f(dset_id,   hdferr)
                
              !close file
              CALL h5fclose_f(file_id, hdferr)
              
              !! store input
              do nf_eff = 1, N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 p2_tmp(:,:,nf)=phi0_in(:,:,1)
              enddo

         endif !irank 0

         if ( irank_g .eq. 0 ) then
            do ii = 1, isize_g-1
               Ndt=N_FL*ndim*n_part
               call mpi_send(p1_tmp,  Ndt,mpi_complex16, ii, 2, MPI_COMM_WORLD,IERR)
               call mpi_send(p2_tmp,  Ndt,mpi_complex16, ii, 3, MPI_COMM_WORLD,IERR)
            ENDDO
         else
            Ndt=N_FL*ndim*n_part
            call mpi_recv(p1_tmp,  Ndt,mpi_complex16, 0, 2, MPI_COMM_WORLD,STATUS,IERR)
            call mpi_recv(p2_tmp,  Ndt,mpi_complex16, 0, 3, MPI_COMM_WORLD,STATUS,IERR)
         ENDIF

         !! allocate tmp matrix
         allocate(smat(n_part,n_part), ipiv(N_part), log_zdet(N_slat))

         !! test whether overlap has minus sign
         alpha=1.d0
         beta=0.d0
         ctmp = 0.d0
         do nf_eff = 1, N_FL_eff
            nf=Calc_Fl_map(nf_eff)
            call zgemm('C','N',N_part,N_part,Ndim,alpha,p1_tmp(1,1,nf),Ndim,phi_0_r(nf_eff,1)%U(1,1),ndim,beta,smat(1,1),N_part)
            ! ZGETRF computes an LU factorization of a general M-by-N matrix A
            ! using partial pivoting with row interchanges.
            call ZGETRF(N_part, N_part, smat, N_part, ipiv, info)
            ! obtain log of det
            zdet  = 0.d0
            phase = 1.d0
            Do n=1,N_part
               if (ipiv(n).ne.n) then
                  phase = -phase
               endif
               zdet = zdet + log(smat(n,n))
            enddo
            zdet = zdet + log(phase)

            ctmp = ctmp + zdet
         enddo
         z_norm = exp(ctmp)
         d_norm = dble(z_norm)
         if (d_norm .lt. 0.d0) then
            c1 = cmplx(0.d0,1.d0,kind(1.d0))
         else
            c1 = cmplx(1.d0,0.d0,kind(1.d0))
         endif

         !! combine the trial wave function by 
         !! | \Psi_T > = c1*| \Psi_T^1 > + c1^{\dagger}*| \Psi_T^2 >
         do nf_eff = 1, N_FL_eff
            nf=Calc_Fl_map(nf_eff)
            phi_0_l(nf_eff,1)%U(:,:)=p1_tmp(:,:,nf)*c1
            phi_0_l(nf_eff,2)%U(:,:)=p2_tmp(:,:,nf)*conjg(c1)
            WF_L(nf,1)%P(:,:)=p1_tmp(:,:,nf)*c1
            WF_L(nf,2)%P(:,:)=p2_tmp(:,:,nf)*conjg(c1)
         enddo

         !! normalization of overlap <\Psi_T | \phi_k^0>

         alpha=1.d0
         beta=0.d0
         log_zdet(:) = 0.d0
         do ns = 1, n_slat
         do nf_eff = 1, N_FL_eff
            nf=Calc_Fl_map(nf_eff)
            call zgemm('C','N',N_part,N_part,Ndim,alpha,phi_0_l(nf_eff,ns)%U(1,1),Ndim,phi_0_r(nf_eff,1)%U(1,1),ndim,beta,smat(1,1),N_part)
            ! ZGETRF computes an LU factorization of a general M-by-N matrix A
            ! using partial pivoting with row interchanges.
            call ZGETRF(N_part, N_part, smat, N_part, ipiv, info)
            ! obtain log of det
            zdet  = 0.d0
            phase = 1.d0
            Do n=1,N_part
               if (ipiv(n).ne.n) then
                  phase = -phase
               endif
               zdet = zdet + log(smat(n,n))
            enddo
            zdet = zdet + log(phase)

            log_zdet(ns) = log_zdet(ns) + zdet
         enddo
         enddo

         z_norm = exp(log_zdet(1)) + exp(log_zdet(2))
         if (irank_g .eq. 0 ) write(*,*) 'overlap for input slater', z_norm
         !d_norm = dble(z_norm)
         !z_norm = (1.d0/d_norm)**(1.d0/dble(2*N_part))
         z_norm = (1.d0/z_norm)**(1.d0/dble(2*N_part))
         if (irank_g .eq. 0 ) write(*,*) 'renormalized factor for input slater z_norm', z_norm

         do i_wlk  = 1, n_wlk
         do nf_eff = 1, N_FL_eff
            nf=Calc_Fl_map(nf_eff)
            phi_0_r(nf_eff,i_wlk)%U(:,:) = phi_0_r(nf_eff,i_wlk)%U(:,:)*z_norm
         enddo
         enddo

         if (irank_g .eq. 0 ) then
             deallocate(phi0_in)
         endif
         deallocate(p1_tmp)
         deallocate(p2_tmp)
         deallocate(smat, log_zdet, ipiv)

        END SUBROUTINE trial_in_hdf5
#endif
        
    end Module stepwlk_mod
