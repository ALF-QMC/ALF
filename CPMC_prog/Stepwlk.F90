    Module stepwlk_mod
        Use Hamiltonian_main 
        Use Operator_mod
        Use Control
        Use Hop_mod
        Use UDV_State_mod
        use gfun_mod
        use upgrade_mod
        
        Contains

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

          N_op     = Size(OP_V,1)

          do i_wlk = 1, N_wlk

             !! Kinetic part exp(-/Delta/tau T/2)
             if ( weight_k(i_wlk) .gt. Zero ) then

                 ! update weight by fac_norm
                 weight_k(i_wlk)=weight_k(i_wlk)*exp(fac_norm)
                 
                 call half_K_propagation( phi_trial(:,i_wlk), phi_0(:,i_wlk), GR(:,:,:,i_wlk), i_wlk )

             endif

             !! propagate with interaction
             do n = 1, N_op

             if ( weight_k(i_wlk) .gt. Zero ) then

                ! upgrade Green's function
                N_type = 2
                do nf_eff = 1,N_FL_eff
                   nf=Calc_Fl_map(nf_eff)
                   Call Op_Wrapdo( GR(:,:,nf,i_wlk), Op_V(n,nf), 1.d0, Ndim, N_Type,1)
                enddo

                call Upgrade(GR(:,:,:,i_wlk),n,spin, i_wlk )
                nsigma_bp(i_wlk)%f(n,ntau_bp) = spin

                N_type = 2
                do nf_eff = 1,N_FL_eff
                   nf=Calc_Fl_map(nf_eff)
                   Call Op_Wrapup( Gr(:,:,nf,i_wlk), Op_V(n,nf), 1.d0, Ndim, N_Type, 1 )
                enddo

                ! propagate slater determinant
                Do nf_eff = 1,N_FL_eff
                    nf=Calc_Fl_map(nf_eff)
                    Call Op_mmultR(phi_0(nf_eff,i_wlk)%U,Op_V(n,nf),spin,'n',1)
                enddo

             endif

             enddo
             
             !! Kinetic part exp(-/Delta/tau T/2)
             if ( weight_k(i_wlk) .gt. Zero ) then
                 call half_K_propagation( phi_trial(:,i_wlk), phi_0(:,i_wlk), GR(:,:,:,i_wlk), i_wlk )
             endif

          enddo

        END SUBROUTINE stepwlk_move

        subroutine re_orthonormalize_walkers(Phi_0, cop)
          
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0
          Character, Intent(IN)    :: cop
          
          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, N_size, I, i_wlk
          Real    (Kind=Kind(0.d0)) :: Overlap_ratio, Zero = 1.0E-8
          Complex (Kind=Kind(0.d0)) :: Det_D(N_FL)

          do i_wlk = 1, N_wlk
          
          if ( weight_k(i_wlk) .gt. Zero ) then

              N_size = phi_0(1,1)%n_part

              !Carry out U,D,V decomposition.
              Do nf_eff = 1, N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 CALL Phi_0(nf_eff,i_wlk)%decompose
              enddo
              
              Det_D = cmplx(1.d0, 0.d0, kind(0.d0))

              Do nf_eff = 1, N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 DO I = 1,N_size
                    Det_D(nf)=Det_D(nf)*Phi_0(nf_eff,i_wlk)%D(I)
                 ENDDO
              enddo
              
              if (reconstruction_needed) call ham%weight_reconstruction(Det_D)

              Do nf_eff = 1, N_FL_eff
                 Phi_0(nf_eff,i_wlk)%D(:) = cmplx(1.d0, 0.d0, kind(0.d0))
              enddo
              
              ! update the overlap when normal propagation
              if (cop == 'U') Overlap(i_wlk)=Overlap(i_wlk)/product(Det_D)

          endif

          enddo

        END SUBROUTINE re_orthonormalize_walkers
        
        SUBROUTINE initial_wlk( phi_trial, phi_0, phi_bp_l, phi_bp_r, udvst, STAB_Nt, GR, nwrap )
#ifdef MPI
        Use mpi
#endif
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:)  , ALLOCATABLE, INTENT(INOUT) :: phi_trial, phi_0, phi_bp_l, phi_bp_r
          CLASS(UDV_State), Dimension(:,:,:), ALLOCATABLE, INTENT(INOUT) :: udvst
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:,:), Allocatable, INTENT(INOUT) :: GR
          INTEGER, dimension(:)    ,allocatable,  INTENT(INOUT) :: Stab_nt
          INTEGER, INTENT(IN) :: nwrap

          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, i_grc, NSTM, NST, ltrot_bp, ns
          Complex (Kind=Kind(0.d0)) :: Overlap_old, Overlap_new, Z, Z1,Z2, tot_ene, ZP
          Complex (Kind=Kind(0.d0)) :: tot_c_weight, el_tmp
          complex (Kind=Kind(0.d0)) :: det_Vec(N_FL)
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio, X1, wtmp
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8, tot_re_weight
          Character (LEN=64) :: FILE_TG, FILE_seeds
          Logical ::   LCONF, LCONF_H5

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          NSTM = Size(udvst, 1)
          ltrot_bp = size(nsigma_bp(1)%f,2)
          tot_ene    = cmplx(0.d0,0.d0,kind(0.d0))
          tot_re_weight = 0.d0

          allocate ( Stab_nt(0:Nstm) )
          Stab_nt(0) = 0
          do n = 1,Nstm -1
             Stab_nt(n) = nwrap*n
          enddo
          Stab_nt(Nstm) = ltrot_bp

          do i_wlk = 1, N_wlk
             do nf_eff = 1, N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                CALL phi_0   (nf_eff, i_wlk)%init(ndim,'r',WF_R(1,nf)%P)
                CALL phi_bp_r(nf_eff, i_wlk)%init(ndim,'r',WF_R(1,nf)%P)
             enddo
          enddo
          
          do i_grc = 1, N_grc
             i_wlk = (i_grc-1)/N_slat+1      ! index for walker
             ns    = i_grc-(i_wlk-1)*N_slat  ! index for slater det
             do nf_eff = 1, N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                do n = 1, nstm
                   CALL udvst(n,nf_eff, i_grc)%alloc(ndim)
                ENDDO
                CALL phi_trial(nf_eff, i_grc)%init(ndim,'l',WF_L(ns,nf)%P)
                CALL phi_bp_l (nf_eff, i_grc)%init(ndim,'l',WF_L(ns,nf)%P)
             enddo
          enddo
          
          !file_tg = "trial_0.h5"
          !INQUIRE (FILE=file_tg, EXIST=LCONF_H5)
          !IF (LCONF_H5) THEN
          !    if ( irank_g .eq. 0 ) write (*,*) "read input trial"
          !    call trial_in_hdf5( phi_0, phi_trial, file_tg )
          !endif

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
                    call CGRP(Z, GR(:,:,nf,i_grc), phi_0(nf_eff,i_wlk), phi_trial(nf_eff,i_grc))
                    det_vec(nf_eff) = Z
                 enddo
                 det_Vec(:) = det_Vec(:) * N_SUN
                 if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
                 overlap(i_grc) = sum(det_vec)
             enddo
          enddo

          !! initial energy
          tot_ene = 0.d0
          tot_c_weight = 0.d0
          do i_wlk = 1, N_wlk
             
             if (weight_k(i_wlk) .ge. 0.d0 ) then

             Z = 0.d0
             do ns = 1, N_slat
                 i_grc = ns+(i_wlk-1)*N_slat
                 Z = Z + exp(overlap(i_grc))
             enddo
             
             do ns = 1, N_slat
                 i_grc = ns+(i_wlk-1)*N_slat
                 el_tmp  = dble(ham%E0_local(GR(:,:,:,i_grc)))
                 tot_ene = tot_ene + el_tmp*weight_k(i_wlk)*exp(overlap(i_grc))/Z
             enddo
             
             tot_c_weight  = tot_c_weight  + weight_k(i_wlk)

             endif

          enddo
          tot_re_weight = dble(tot_c_weight)
          
          CALL MPI_REDUCE(tot_ene      ,Z1,1,MPI_COMPLEX16,MPI_SUM, 0,Group_comm,IERR)
          CALL MPI_REDUCE(tot_re_weight,X1,1,MPI_REAL8    ,MPI_SUM, 0,Group_comm,IERR)
          
          if (Irank_g == 0 ) then
              Z1 = Z1/X1
              fac_norm= dble(Z1)
          endif
          CALL MPI_BCAST(fac_norm, 1, MPI_REAL8, 0,MPI_COMM_WORLD,ierr)

          file_seeds = "seedvec_in"
          INQUIRE (FILE=file_seeds, EXIST=LCONF)
          IF (LCONF) THEN
              if ( irank_g .eq. 0 ) write (*,*) "read input seeds"
              call seed_vec_in(file_seeds)
          endif

        END SUBROUTINE initial_wlk

        SUBROUTINE population_control( phi_0, phi_bp_r ) 
#ifdef MPI
          Use mpi
#endif
          Use Random_wrap
        
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0, phi_bp_r
          
          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, it_wlk, n_exc,pop_exc(N_wlk_mpi,4)
          Integer :: j, it, i_t, i_st, i_ed, nu_wlk, i_src, i_wlk, j_src, j_wlk, n1, n2, nrg, nfrg, ilabel
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8, d_scal, sum_w, w_count, w_tmp(N_wlk_mpi), weight_mpi(N_wlk_mpi)
          Complex (Kind=Kind(0.d0)) :: overlap_tmp(N_wlk)
          Complex (Kind=Kind(0.d0)) :: Z1,Z2,Z3, Z_s_array(1), Z_r_array(1)
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
               CALL phi_0_m (nf_eff, i_wlk)%init(ndim,'r',WF_R(nf)%P)
               CALL phi_bp_m(nf_eff, i_wlk)%init(ndim,'r',WF_R(nf)%P)
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
          overlap_tmp=overlap

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
                        overlap_tmp(j_wlk)=overlap(i_wlk) 
                        nsigma_store(j_wlk)%f=nsigma_bp(i_wlk)%f
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
                Z_s_array(1) = overlap(i_wlk)
                
                ilabel = (it-1)*nrg
                call mpi_send(Z_s_array,1,MPI_COMPLEX16,j_src,ilabel,Group_comm,IERR)
                
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
                call mpi_recv(Z_r_array,1,MPI_COMPLEX16,i_src,ilabel,Group_comm,STATUS,IERR)
                overlap_tmp(j_wlk) = Z_r_array(1) 
                
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

          overlap=overlap_tmp
          ! reset weight
          weight_k(:)=1.d0
          do nf_eff = 1, N_FL_eff
             do i_wlk = 1, N_wlk
                phi_0   (nf_eff,i_wlk)=phi_0_m (nf_eff,i_wlk)
                phi_bp_r(nf_eff,i_wlk)=phi_bp_m(nf_eff,i_wlk)
             enddo
          enddo
          
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
                call phi_bp_r(nf_eff,i_wlk)%decompose
             enddo
          enddo

        END SUBROUTINE store_phi

        SUBROUTINE backpropagation( GR_mix, phi_bp_l, phi_bp_r, udvst, Stab_nt, ltau )
          
          Implicit none
     
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:,:), Allocatable, INTENT(IN) :: GR_mix
          CLASS(UDV_State), Dimension(:,:)  , ALLOCATABLE, INTENT(INOUT) :: phi_bp_l, phi_bp_r
          CLASS(UDV_State), Dimension(:,:,:), ALLOCATABLE, INTENT(INOUT) :: udvst
          INTEGER, dimension(:)    ,allocatable,  INTENT(IN) :: Stab_nt
          Integer, INTENT(IN) :: ltau

          !Local 
          Complex (Kind=Kind(0.d0)) :: GR_bp(NDIM,NDIM,N_FL,N_wlk)
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, ltrot_bp, N_op, nstm, nst, ntau
          Complex (Kind=Kind(0.d0)) :: Z, Z_weight
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
          Real (Kind=Kind(0.d0))    :: Zero = 1.0E-8

          N_op     = Size(OP_V,1)
          ltrot_bp = size(nsigma_bp(1)%f, 2)
          nstm     = Size(udvst, 1)
        
          !! initialization
          do i_wlk = 1, N_wlk
            do nf_eff = 1, N_FL_eff
               nf=Calc_Fl_map(nf_eff)
               CALL phi_bp_l(nf_eff, i_wlk)%reset('l',WF_L(nf)%P)
               CALL udvst(nstm, nf_eff, i_wlk)%reset('l',WF_L(nf)%P)
            enddo
          enddo

          !! backpropagation
          nst= nstm-1
          Do nt = ltrot_bp, 1, -1
             ntau = nt-1
             
             Do i_wlk = 1, N_wlk

                if ( weight_k(i_wlk) .gt. Zero  ) then
                
                Do nf_eff = 1,N_FL_eff
                   
                   nf=Calc_Fl_map(nf_eff)
                   Call Hop_mod_mmthlc_1D2(phi_bp_l(nf_eff,i_wlk)%U,nf,1)

                   Do n = N_op,1,-1
                      Call Op_mmultR(phi_bp_l(nf_eff,i_wlk)%U,Op_V(n,nf),nsigma_bp(i_wlk)%f(n,nt),'c',1)
                   enddo
                
                   Call Hop_mod_mmthlc_1D2(phi_bp_l(nf_eff,i_wlk)%U,nf,1)

                enddo

                endif

             Enddo
             
             if ( ntau .eq. stab_nt(nst) .and. ntau .ne. 0 ) then
                 call re_orthonormalize_walkers(phi_bp_l, 'N')
                 Do i_wlk = 1, N_wlk
                 Do nf_eff = 1,N_FL_eff
                    udvst(nst, nf_eff, i_wlk) = phi_bp_l(nf_eff, i_wlk)
                 ENDDO
                 ENDDO
                 nst = nst - 1
             endif
             
          Enddo
          
          !! svd at tau = 0
          call re_orthonormalize_walkers(phi_bp_l, 'N')

          !! compute the total weight
          Z_weight = ham%sum_weight

          !! equal time measurement

          do i_wlk = 1, N_wlk
             do nf_eff = 1,N_Fl_eff
                nf=Calc_Fl_map(nf_eff)
                call CGRP(Z, GR_bp(:,:,nf,i_wlk), phi_bp_r(nf_eff,i_wlk), phi_bp_l(nf_eff,i_wlk))
             enddo
             If (reconstruction_needed) Call ham%GR_reconstruction( GR_bp(:,:,:,i_wlk) )
             CALL ham%Obser( GR_bp(:,:,:,i_wlk), GR_mix(:,:,:,i_wlk), i_wlk, Z_weight )
          enddo

          !! time dependence measurement
          if ( ltau .eq. 1 ) then
             call bp_measure_tau(phi_bp_l, phi_bp_r, udvst, stab_nt )
          endif

        END SUBROUTINE backpropagation

        SUBROUTINE bp_measure_tau(phi_bp_l, phi_bp_r, udvst, stab_nt )
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_bp_l, phi_bp_r
          CLASS(UDV_State), Dimension(:,:,:), ALLOCATABLE, INTENT(INOUT) :: udvst
          INTEGER, dimension(:)    ,allocatable,  INTENT(IN) :: Stab_nt

          !Local 
          Complex (Kind=Kind(0.d0)) :: Temp(NDIM,NDIM), GR(NDIM,NDIM,N_FL), GRC(NDIM,NDIM,N_FL)
          Complex (Kind=Kind(0.d0)) :: GT0(NDIM,NDIM,N_FL,N_wlk), G00(NDIM,NDIM,N_FL,N_wlk)
          Complex (Kind=Kind(0.d0)) :: GTT(NDIM,NDIM,N_FL,N_wlk), G0T(NDIM,NDIM,N_FL,N_wlk)
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, N_op, ntau, I, nst, nstm
          Complex (Kind=Kind(0.d0)) :: Z, Z_weight, DETZ
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
          Real (Kind=Kind(0.d0))    :: Zero = 1.0E-8

          !! initialization
          N_op     = Size(OP_V,1)
          nstm     = Size(udvst, 1)
          
          !! compute the total weight
          Z_weight = ham%sum_weight
          
          do i_wlk = 1, N_wlk
             do nf_eff = 1, N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                CALL CGRP(DetZ, GR(:,:,nf), phi_bp_r(nf_eff,i_wlk), phi_bp_l(nf_eff,i_wlk))
             enddo
             GTT(:,:,:,i_wlk) = GR
          enddo

          G00 = GTT
          GT0 = GTT
          G0T = GTT
          do i_wlk = 1, N_wlk
          do nf_eff=1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             do I=1,Ndim
                G0T(I,I,nf,i_wlk)=G0T(I,I,nf,i_wlk)-1.d0
             enddo
          enddo
          enddo
          
          ntau = 0
          do i_wlk = 1, N_wlk
            !! call reconstruction of non-calculated flavor blocks
            If (reconstruction_needed) then
                Call ham%GR_reconstruction ( G00(:,:,:,i_wlk) )
                Call ham%GR_reconstruction ( GTT(:,:,:,i_wlk) )
                Call ham%GRT_reconstruction( GT0(:,:,:,i_wlk), G0T(:,:,:,i_wlk) )
            endif
            CALL ham%obserT(ntau,GT0(:,:,:,i_wlk),G0T(:,:,:,i_wlk),G00(:,:,:,i_wlk), & 
                & GTT(:,:,:,i_wlk), i_wlk, Z_weight)
          enddo

          NST=1
          do ntau = 1, ltrot

             do i_wlk = 1, N_wlk
                !! Propagate wave function
                if ( weight_k(i_wlk) .gt. Zero ) then

                !! Propagate Green's function
                CALL PROPR  (GT0(:,:,:,i_wlk),ntau,i_wlk)
                CALL PROPRM1(G0T(:,:,:,i_wlk),ntau,i_wlk)
                CALL PROPRM1(GTT(:,:,:,i_wlk),ntau,i_wlk)
                CALL PROPR  (GTT(:,:,:,i_wlk),ntau,i_wlk)
                
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
                If (reconstruction_needed) then
                    Call ham%GR_reconstruction ( G00(:,:,:,i_wlk) )
                    Call ham%GR_reconstruction ( GTT(:,:,:,i_wlk) )
                    Call ham%GRT_reconstruction( GT0(:,:,:,i_wlk), G0T(:,:,:,i_wlk) )
                endif
                CALL ham%obserT(ntau,GT0(:,:,:,i_wlk),G0T(:,:,:,i_wlk),G00(:,:,:,i_wlk), & 
                    & GTT(:,:,:,i_wlk), i_wlk, Z_weight)
             enddo
             
             !! call svd
             if (  ntau .eq. stab_nt(nst) )  then
                 call re_orthonormalize_walkers(phi_bp_r, 'N')
          
                 do i_wlk = 1, N_wlk
                 
                    do nf_eff = 1, N_FL_eff
                       nf=Calc_Fl_map(nf_eff)
                       phi_bp_l(nf_eff, i_wlk) = udvst(nst, nf_eff, i_wlk) 
                       CALL CGRP(DetZ, GR(:,:,nf), phi_bp_r(nf_eff,i_wlk), phi_bp_l(nf_eff,i_wlk))
                    enddo
                    GTT(:,:,:,i_wlk) = GR
                    
                    GRC = -GR
                    do nf_eff=1,N_FL_eff
                       nf=Calc_Fl_map(nf_eff)
                       do I=1,Ndim
                          GRC(I,I,nf)=GRC(I,I,nf)+1.d0
                       enddo
                    enddo
                  
                    do nf_eff=1,N_FL_eff
                       nf=Calc_Fl_map(nf_eff)
                       CALL MMULT(TEMP,GR (:,:,nf),GT0(:,:,nf,i_wlk))
                       GT0(:,:,nf,i_wlk) = TEMP
                       CALL MMULT(TEMP,G0T(:,:,nf,i_wlk),GRC(:,:,nf))
                       G0T(:,:,nf,i_wlk) = TEMP
                    enddo

                 enddo
                 nst = nst + 1
             endif

          enddo

        END SUBROUTINE bp_measure_tau
        
        SUBROUTINE half_K_propagation( phi_trial, phi_0, GR, i_wlk )
          
          Implicit none
     
          CLASS(UDV_State), INTENT(IN   ) :: phi_trial(N_FL)
          CLASS(UDV_State), INTENT(INOUT) :: phi_0    (N_FL)
          COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: GR(Ndim,Ndim,N_FL)
          Integer, INTENT(IN) :: i_wlk
          
          ! local
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR
          Complex (Kind=Kind(0.d0)) :: Overlap_old, Overlap_new, log_O_new, log_O_old, Z
          COMPLEX (Kind=Kind(0.d0)) :: det_Vec(N_FL)
          Real    (Kind=Kind(0.d0)) :: overlap_ratio, re_overlap
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8
          
          log_O_old = overlap

          Do nf_eff = 1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             Call Hop_mod_mmthr_1D2(phi_0(nf_eff)%U,nf,1)
          enddo
          
          ! Update Green's function
          do nf_eff = 1,N_Fl_eff
             nf=Calc_Fl_map(nf_eff)
             call CGRP(Z, GR(:,:,nf), phi_0(nf_eff), phi_trial(nf_eff))
             det_vec(nf_eff) = Z
          enddo

          det_Vec(:) = det_Vec(:) * N_SUN
          if (reconstruction_needed) call ham%weight_reconstruction(det_Vec)
          log_O_new = sum(det_Vec) 

          overlap_ratio = exp(log_O_new-log_O_old)
          re_overlap    = dble( overlap_ratio )
          if ( re_overlap .gt. zero ) then
              overlap (i_wlk) = overlap (i_wlk) + (log_o_new - log_o_old)
              weight_k(i_wlk) = weight_k(i_wlk)*re_overlap
          else
              weight_k(i_wlk) = 0.d0
          endif

        END SUBROUTINE half_K_propagation

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
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(IN) :: phi_0

          ! LOCAL
          CHARACTER (LEN=64) :: FILE_TG, filename
          Complex (Kind=Kind(0.d0)), pointer :: phi0_out(:,:,:,:)
          Complex (Kind=Kind(0.d0)), pointer :: overlap_out(:)
          Real    (Kind=Kind(0.d0)), pointer :: weight_out(:)
          Complex (Kind=Kind(0.d0)), allocatable :: otphi_tmp(:), p0_tmp(:,:,:,:), p1_tmp(:,:,:,:)
          Real    (Kind=Kind(0.d0)), allocatable :: wt_tmp(:)

          INTEGER             :: K, hdferr, rank, nf, nw, n_part, i0, i1, i2, i_st, i_ed, Ndt, ii
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
              allocate(weight_out(N_wlk_mpi), overlap_out(N_wlk_mpi))
          endif

          allocate(p0_tmp(ndim,n_part,n_fl_eff,n_wlk))
          allocate(p1_tmp(ndim,n_part,n_fl_eff,n_wlk))
          allocate(wt_tmp(N_wlk))
          allocate(otphi_tmp(N_wlk))

          do nf = 1, N_FL_eff
          do nw = 1, N_wlk
              p0_tmp(:,:,nf,nw)=phi_0(nf,nw)%U(:,:)
          enddo
          enddo
          
          if ( irank_g .ne. 0 ) then
              call mpi_send(overlap    ,N_wlk,mpi_complex16, 0, 0, MPI_COMM_WORLD,IERR)
              call mpi_send(weight_k   ,N_wlk,    mpi_real8, 0, 1, MPI_COMM_WORLD,IERR)
              Ndt=N_FL_eff*N_wlk*ndim*n_part
              call mpi_send(p0_tmp     ,  Ndt,mpi_complex16, 0, 2, MPI_COMM_WORLD,IERR)
          else
             do ii = 1, isize_g-1
                i_st=ii*N_wlk+1
                i_ed=(ii+1)*N_wlk

                call mpi_recv(otphi_tmp,N_wlk,mpi_complex16, ii, 0, MPI_COMM_WORLD,STATUS,IERR)
                overlap_out(i_st:i_ed) = otphi_tmp(:)
                call mpi_recv(wt_tmp,N_wlk,    mpi_real8, ii, 1, MPI_COMM_WORLD,STATUS,IERR)
                weight_out(i_st:i_ed) = wt_tmp(:)
                Ndt=N_FL_eff*N_wlk*ndim*n_part
                call mpi_recv(p1_tmp,  Ndt,mpi_complex16, ii, 2, MPI_COMM_WORLD,STATUS,IERR)
                phi0_out(:,:,:,i_st:i_ed)=p1_tmp
             ENDDO
          ENDIF

          if ( irank_g .eq. 0 ) then
              i_st=1
              i_ed=N_wlk
              weight_out (i_st:i_ed) = weight_k(:)
              overlap_out(i_st:i_ed) = overlap(:)
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
              dims  = [2, N_wlk_mpi]
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
#if defined HDF5_ZLIB
                ! Set ZLIB / DEFLATE Compression using compression level HDF5_ZLIB
              CALL h5pset_deflate_f(crp_list, HDF5_ZLIB, hdferr)
#endif
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

          INTEGER             :: K, hdferr, rank, nf, nw, n_part, i0, i1, i2, i_st, i_ed, Ndt, ii, nwalk_in
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
          allocate(otphi_tmp(N_wlk))
          
          if ( irank .eq. 0 ) then

              !open file
              CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, hdferr)

              !open and read field overlap
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
              nwalk_in=dims(rank)
              !! allocate !!
              allocate(phi0_out(ndim,n_part,n_fl_eff,nwalk_in))
              allocate(weight_out(nwalk_in), overlap_out(nwalk_in))
              !!-----------!!
              dat_ptr = C_LOC(overlap_out(1))
              CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !close dataspace
              CALL h5sclose_f(dataspace, hdferr)
              !close objects
              CALL h5dclose_f(dset_id,   hdferr)
              deallocate(dims, maxdims)


              !open and read real weight
              dset_name = "weight_re"
              !Open the  dataset.
              CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
              dat_ptr = C_LOC(weight_out(1))
              CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !close objects
              CALL h5dclose_f(dset_id,   hdferr)
              
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
               i_st=ii*N_wlk+1
               i_ed=(ii+1)*N_wlk
               
               otphi_tmp(:)=overlap_out(i_st:i_ed)
               call mpi_send(otphi_tmp,N_wlk,mpi_complex16, ii, 0, MPI_COMM_WORLD,IERR)
               wt_tmp(:)=weight_out(i_st:i_ed)
               call mpi_send(wt_tmp,N_wlk,    mpi_real8, ii, 1, MPI_COMM_WORLD,IERR)
               Ndt=N_FL_eff*N_wlk*ndim*n_part
               p1_tmp=phi0_out(:,:,:,i_st:i_ed)
               call mpi_send(p1_tmp,  Ndt,mpi_complex16, ii, 2, MPI_COMM_WORLD,IERR)
            ENDDO
         else
            call mpi_recv(otphi_tmp,N_wlk,mpi_complex16, 0, 0, MPI_COMM_WORLD,STATUS,IERR)
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
             i_st=1
             i_ed=N_wlk
             p0_tmp=phi0_out(:,:,:,i_st:i_ed)
             
             weight_k(:) = weight_out(i_st:i_ed)
             overlap (:) = overlap_out(i_st:i_ed)
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

        SUBROUTINE trial_in_hdf5( phi_0_r, phi_0_l, file_tg )
#ifdef MPI
          Use mpi
#endif
#if defined HDF5
          Use hdf5
          use h5lt
#endif
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0_r, phi_0_l
          CHARACTER (LEN=64), INTENT(IN)  :: FILE_TG

          ! LOCAL
          CHARACTER (LEN=64) :: filename
          Complex (Kind=Kind(0.d0)), pointer :: phi0_out(:,:,:)
          Complex (Kind=Kind(0.d0)), allocatable :: p1_tmp(:,:,:)

          INTEGER             :: K, hdferr, rank, nf, nw, n_part, i0, i1, i2, i_st, i_ed, Ndt, ii, nwalk_in
          Integer             :: nf_eff
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

          n_part=phi_0_l(1,1)%n_part
          
          allocate(p1_tmp(ndim,n_part,n_fl_eff))
          
          if ( irank .eq. 0 ) then

              !open file
              CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, hdferr)

              !! allocate !!
              allocate(phi0_out(ndim,n_part,n_fl))
              !!-----------!!
              !open and read real weight
              dset_name = "wavefunction"
              !Open the  dataset.
              CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
              dat_ptr = C_LOC(phi0_out(1,1,1))
              !Write data
              CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dat_ptr, hdferr)
              !close objects
              CALL h5dclose_f(dset_id,   hdferr)
                
              !close file
              CALL h5fclose_f(file_id, hdferr)

              do nf_eff = 1, N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 p1_tmp(:,:,nf_eff)=phi0_out(:,:,nf)
              enddo

         endif !irank 0

         if ( irank_g .eq. 0 ) then
            do ii = 1, isize_g-1
               Ndt=N_FL_eff*ndim*n_part
               call mpi_send(p1_tmp,  Ndt,mpi_complex16, ii, 2, MPI_COMM_WORLD,IERR)
            ENDDO
         else
            Ndt=N_FL_eff*ndim*n_part
            call mpi_recv(p1_tmp,  Ndt,mpi_complex16, 0, 2, MPI_COMM_WORLD,STATUS,IERR)
         ENDIF

         do nf_eff = 1, N_FL_eff
            nf=Calc_Fl_map(nf_eff)
            do nw = 1, N_wlk
                phi_0_l(nf_eff,nw)%U(:,:)=p1_tmp(:,:,nf_eff)
                phi_0_r(nf_eff,nw)%U(:,:)=p1_tmp(:,:,nf_eff)
            enddo
            WF_L(nf)%P(:,:)=p1_tmp(:,:,nf_eff)
            WF_R(nf)%P(:,:)=p1_tmp(:,:,nf_eff)
         enddo

         if (irank_g .eq. 0 ) then
             deallocate(phi0_out)
         endif
         deallocate(p1_tmp)

        END SUBROUTINE trial_in_hdf5
#endif
        
    end Module stepwlk_mod
