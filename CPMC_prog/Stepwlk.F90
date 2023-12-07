    Module stepwlk_mod
        Use Hamiltonian_main 
        Use Operator_mod
        Use Control
        Use Hop_mod
        Use UDV_State_mod
        use cgr1_mod
        use upgrade_mod
        
        Contains

        SUBROUTINE stepwlk_move( phi_trial, phi_0, GR, phase, phase_alpha, ntau_bp ) 
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(IN)    :: phi_trial
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:,:), Allocatable, INTENT(INOUT) :: GR
          COMPLEX (Kind=Kind(0.d0)), Dimension(:)      , Allocatable, INTENT(INOUT) :: phase, phase_alpha
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
                 weight_k(i_wlk)=weight_k(i_wlk)*exp(fac_norm);
                 
                 call half_K_propagation( phi_trial(:,i_wlk), phi_0(:,i_wlk), GR(:,:,:,i_wlk), phase(i_wlk), i_wlk )

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

                Call Upgrade(GR(:,:,:,i_wlk),n,PHASE(i_wlk), PHASE_alpha(i_wlk), spin, i_wlk )
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
                 call half_K_propagation( phi_trial(:,i_wlk), phi_0(:,i_wlk), GR(:,:,:,i_wlk), phase(i_wlk), i_wlk )
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

              N_size = phi_0(1,1)%ndim

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
              if (cop == 'R') Overlap(i_wlk)=Overlap(i_wlk)/product(Det_D)

          endif

          enddo

        END SUBROUTINE re_orthonormalize_walkers
        
        SUBROUTINE initial_wlk( phi_trial, phi_0, phi_bp_l, phi_bp_r, udvst, STAB_Nt, GR, phase, nwrap )
#ifdef MPI
        Use mpi
#endif
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_trial, phi_0, phi_bp_l, phi_bp_r
          CLASS(UDV_State), Dimension(:,:,:), ALLOCATABLE, INTENT(INOUT) :: udvst
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:,:), Allocatable, INTENT(INOUT) :: GR
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable, INTENT(INOUT) :: phase
          INTEGER, dimension(:)    ,allocatable,  INTENT(INOUT) :: Stab_nt
          INTEGER, INTENT(IN) :: nwrap

          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, NSTM, NST, ltrot_bp
          Complex (Kind=Kind(0.d0)) :: Overlap_old, Overlap_new, Z, Z1, tot_ene
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio, X1
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8, tot_re_weight
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable :: Phase_array

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
          allocate(Phase_array(N_FL))
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
               do n = 1, NSTM
                  CALL udvst(n,nf_eff, i_wlk)%alloc(ndim)
               ENDDO
               CALL phi_trial(nf_eff, i_wlk)%init(ndim,'l',WF_L(nf)%P)
               CALL phi_0(nf_eff, i_wlk)%init(ndim,'r',WF_R(nf)%P)
               CALL phi_bp_l(nf_eff, i_wlk)%init(ndim,'l',WF_L(nf)%P)
               CALL phi_bp_r(nf_eff, i_wlk)%init(ndim,'r',WF_R(nf)%P)
            enddo

            NVAR = 1
            do nf_eff = 1,N_Fl_eff
               nf=Calc_Fl_map(nf_eff)
               call CGR(Z, NVAR, GR(:,:,nf,i_wlk), phi_0(nf_eff,i_wlk), phi_trial(nf_eff,i_wlk))
               Phase_array(nf)=Z
            enddo
            if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
            Phase(i_wlk)=product(Phase_array)
            Phase(i_wlk)=Phase(i_wlk)**N_SUN
            
            tot_ene       = tot_ene       + ham%E0_local(GR(:,:,:,i_wlk))*weight_k(i_wlk)
            tot_re_weight = tot_re_weight + weight_k(i_wlk)

          enddo
          
          CALL MPI_REDUCE(tot_ene      ,Z1,1,MPI_COMPLEX16,MPI_SUM, 0,Group_comm,IERR)
          CALL MPI_REDUCE(tot_re_weight,X1,1,MPI_REAL8    ,MPI_SUM, 0,Group_comm,IERR)
          
          if (Irank_g == 0 ) then
              fac_norm= real(Z1, kind(0.d0))/X1
          endif
          CALL MPI_BCAST(fac_norm, 1, MPI_REAL8, 0,MPI_COMM_WORLD,ierr)

        END SUBROUTINE initial_wlk

        SUBROUTINE population_control( phi_0, phi_bp_r, phase_alpha, phase ) 
#ifdef MPI
          Use mpi
#endif
          Use Random_wrap
        
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0, phi_bp_r
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable, INTENT(INOUT) :: phase_alpha
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable, INTENT(INOUT) :: phase
          
          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, it_wlk, n_exc,pop_exc(N_wlk_mpi,4)
          Integer :: j, it, i_t, i_st, i_ed, nu_wlk, i_src, i_wlk, j_src, j_wlk, n1, n2
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8, d_scal, sum_w, w_count, w_tmp(N_wlk_mpi), weight_mpi(N_wlk_mpi)
          Complex (Kind=Kind(0.d0)) :: Overlap_tmp(N_wlk), phase_alpha_tmp(N_wlk), phase_tmp(N_wlk)
          Complex (Kind=Kind(0.d0)) :: Z1,Z2,Z3
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
                call phi_0_m (nf_eff,i_wlk)%assign(phi_0   (nf_eff,i_wlk))
                call phi_bp_m(nf_eff,i_wlk)%assign(phi_bp_r(nf_eff,i_wlk))
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

          !! store phase and overlap
          Overlap_tmp=Overlap
          phase_alpha_tmp=phase_alpha
          phase_tmp=phase

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
                            call phi_0_m (nf_eff,j_wlk)%assign(phi_0   (nf_eff,i_wlk))
                            call phi_bp_m(nf_eff,j_wlk)%assign(phi_bp_r(nf_eff,i_wlk))
                        enddo
                        Overlap_tmp    (j_wlk)=Overlap    (i_wlk) 
                        phase_tmp      (j_wlk)=phase      (i_wlk) 
                        phase_alpha_tmp(j_wlk)=phase_alpha(i_wlk) 
                        nsigma_store(j_wlk)%f=nsigma_bp(i_wlk)%f
                      endif
                  endif
              enddo
              nu_wlk=n
          enddo

          do it=1, n_exc
             i_src = pop_exc(it,1); i_wlk = pop_exc(it,2)
             j_src = pop_exc(it,3); j_wlk = pop_exc(it,4)
             if ( irank_g .eq. i_src ) then
                call mpi_send(overlap    (i_wlk),1,MPI_COMPLEX16,j_src,0,Group_comm,IERR)
                call mpi_send(phase      (i_wlk),1,MPI_COMPLEX16,j_src,1,Group_comm,IERR)
                call mpi_send(phase_alpha(i_wlk),1,MPI_COMPLEX16,j_src,2,Group_comm,IERR)
                do nf_eff = 1, N_FL_eff
                    call phi_0   (nf_eff,i_wlk)%MPI_send_general(j_src, 3, ierr)
                    call phi_bp_r(nf_eff,i_wlk)%MPI_send_general(j_src, 6, ierr)
                enddo
                call mpi_send(nsigma_bp(i_wlk)%f,n1*n2,MPI_REAL8,j_src,9,Group_comm,IERR)
             endif
             if ( irank_g .eq. j_src ) then
                call mpi_recv(overlap_tmp    (j_wlk),1,MPI_COMPLEX16,i_src,0,Group_comm,STATUS,IERR)
                call mpi_recv(phase_tmp      (j_wlk),1,MPI_COMPLEX16,i_src,1,Group_comm,STATUS,IERR)
                call mpi_recv(phase_alpha_tmp(j_wlk),1,MPI_COMPLEX16,i_src,2,Group_comm,STATUS,IERR)
                do nf_eff = 1, N_FL_eff
                    call phi_0_m (nf_eff,j_wlk)%MPI_recv_general(i_src, 3, status, ierr)
                    call phi_bp_m(nf_eff,j_wlk)%MPI_recv_general(i_src, 6, status, ierr)
                enddo
                call mpi_recv(nsigma_store(j_wlk)%f,n1*n2,MPI_REAL8,i_src,9,Group_comm,STATUS,IERR)
             endif
          enddo

          Overlap=Overlap_tmp
          phase_alpha=phase_alpha_tmp
          phase=phase_tmp
          weight_k(:)=1.d0
          do nf_eff = 1, N_FL_eff
             do i_wlk = 1, N_wlk
                call phi_0   (nf_eff,i_wlk)%assign(phi_0_m (nf_eff,i_wlk))
                call phi_bp_r(nf_eff,i_wlk)%assign(phi_bp_m(nf_eff,i_wlk))
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
                call phi_bp_r(nf_eff,i_wlk)%assign(phi_0(nf_eff,i_wlk))
                call Phi_bp_r(nf_eff,i_wlk)%decompose
             enddo
          enddo

        END SUBROUTINE store_phi

        SUBROUTINE backpropagation( phi_bp_l, phi_bp_r, udvst, phase, phase_alpha, Stab_nt, ltau )
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:)  , ALLOCATABLE, INTENT(INOUT) :: phi_bp_l, phi_bp_r
          CLASS(UDV_State), Dimension(:,:,:), ALLOCATABLE, INTENT(INOUT) :: udvst
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable, INTENT(IN) :: phase_alpha
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable, INTENT(IN) :: phase
          INTEGER, dimension(:)    ,allocatable,  INTENT(IN) :: Stab_nt
          Integer, INTENT(IN) :: ltau

          !Local 
          Complex (Kind=Kind(0.d0)) :: GR_bp(NDIM,NDIM,N_FL,N_wlk)
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, ltrot_bp, N_op, nstm, nst
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
               CALL udvst(NSTM, nf_eff, i_wlk)%reset('l',WF_L(nf)%P)
            enddo
          enddo

          !! backpropagation
          nst= nstm-1
          Do nt = ltrot_bp, 1, -1
             
             Do i_wlk = 1, N_wlk

                if ( weight_k(i_wlk) .gt. Zero ) then
                
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
             
             if ( nt .eq. stab_nt(nst) ) then
                 call re_orthonormalize_walkers(phi_bp_l, 'L')
                 Do i_wlk = 1, N_wlk
                 Do nf_eff = 1,N_FL_eff
                    udvst(nst, nf_eff, i_wlk) = phi_bp_l(nf_eff, i_wlk)
                 ENDDO
                 ENDDO
                 nst = nst - 1
             endif
             
          Enddo
          
          !! svd at tau = 0
          call re_orthonormalize_walkers(phi_bp_l, 'L')

          !! compute the total weight
          Z_weight = ham%sum_weight(PHASE, PHASE_alpha)

          !! equal time measurement
          do i_wlk = 1, N_wlk
             do nf_eff = 1,N_Fl_eff
                nf=Calc_Fl_map(nf_eff)
                NVAR = 1
                call CGR(Z, NVAR, GR_bp(:,:,nf,i_wlk), phi_bp_r(nf_eff,i_wlk), phi_bp_l(nf_eff,i_wlk))
             enddo
             If (reconstruction_needed) Call ham%GR_reconstruction( GR_bp(:,:,:,i_wlk) )
             CALL ham%Obser( GR_bp(:,:,:,i_wlk), PHASE(i_wlk), PHASE_alpha(i_wlk), i_wlk, Z_weight )
          enddo

          !! time dependence measurement
          if ( ltau .eq. 1 ) then
             call bp_measure_tau(phi_bp_l, phi_bp_r, udvst, phase, phase_alpha, stab_nt )
          endif

        END SUBROUTINE backpropagation

        SUBROUTINE bp_measure_tau(phi_bp_l, phi_bp_r, udvst, phase, phase_alpha, stab_nt )
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_bp_l, phi_bp_r
          CLASS(UDV_State), Dimension(:,:,:), ALLOCATABLE, INTENT(INOUT) :: udvst
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable, INTENT(IN) :: phase_alpha
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable, INTENT(IN) :: phase
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
          Z_weight = ham%sum_weight(PHASE, PHASE_alpha)
          
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
                & GTT(:,:,:,i_wlk),PHASE(i_wlk), PHASE_ALPHA(i_wlk), i_wlk, Z_weight)
          enddo

          NST=1
          do ntau = 1, ltrot
             
             !! call svd
             if (  ntau .eq. stab_nt(nst) )  then
                 call re_orthonormalize_walkers(phi_bp_r, 'L')
          
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

             do i_wlk = 1, N_wlk
                !! Propagate Green's function
                CALL PROPR  (GT0(:,:,:,i_wlk),ntau,i_wlk)
                CALL PROPRM1(G0T(:,:,:,i_wlk),ntau,i_wlk)
                CALL PROPRM1(GTT(:,:,:,i_wlk),ntau,i_wlk)
                CALL PROPR  (GTT(:,:,:,i_wlk),ntau,i_wlk)
                
                !! Propagate wave function
                if ( weight_k(i_wlk) .gt. Zero ) then
                
                Do nf_eff = 1,N_FL_eff
                   
                   nf=Calc_Fl_map(nf_eff)
                   Call Hop_mod_mmthr_1D2    (phi_bp_r(nf_eff,i_wlk)%U,nf,1)

                   Do n = 1, N_op
                      Call Op_mmultR(phi_bp_r(nf_eff,i_wlk)%U,Op_V(n,nf), nsigma_bp(i_wlk)%f(n,ntau),'n',1)
                   enddo
                
                   Call Hop_mod_mmthr_1D2    (phi_bp_r(nf_eff,i_wlk)%U,nf,1)

                enddo

                endif

                !! call reconstruction of non-calculated flavor blocks
                If (reconstruction_needed) then
                    Call ham%GR_reconstruction ( G00(:,:,:,i_wlk) )
                    Call ham%GR_reconstruction ( GTT(:,:,:,i_wlk) )
                    Call ham%GRT_reconstruction( GT0(:,:,:,i_wlk), G0T(:,:,:,i_wlk) )
                endif
                CALL ham%obserT(ntau,GT0(:,:,:,i_wlk),G0T(:,:,:,i_wlk),G00(:,:,:,i_wlk), & 
                    & GTT(:,:,:,i_wlk),PHASE(i_wlk), PHASE_ALPHA(i_wlk), i_wlk, Z_weight)
             enddo

          enddo

        END SUBROUTINE bp_measure_tau
        
        SUBROUTINE half_K_propagation( phi_trial, phi_0, GR, phase, i_wlk )
          
          Implicit none
     
          CLASS(UDV_State), INTENT(IN   ) :: phi_trial(N_FL)
          CLASS(UDV_State), INTENT(INOUT) :: phi_0    (N_FL)
          COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: GR(Ndim,Ndim,N_FL)
          COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: phase
          Integer, INTENT(IN) :: i_wlk
          
          ! local
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR
          Complex (Kind=Kind(0.d0)) :: Overlap_old, Overlap_new, log_O_new, log_O_old, Z
          Complex (Kind=Kind(0.d0)) :: phase_new, phase_old
          COMPLEX (Kind=Kind(0.d0)) :: Phase_array(N_FL), Det_Vec(N_FL)
          Real    (Kind=Kind(0.d0)) :: Overlap_ratio
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8
          
          call Compute_overlap(Phase_array, Det_Vec, phi_0, phi_trial)
          Det_Vec(:) = Det_Vec(:) * N_SUN
          if (reconstruction_needed) call ham%weight_reconstruction(Det_Vec    )
          if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
          log_O_old = sum(Det_Vec) 
          Phase_old = product(Phase_array)**dble(N_SUN)

          Do nf_eff = 1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             Call Hop_mod_mmthr_1D2(phi_0(nf_eff)%U,nf,1)
          enddo
          
          ! Update Green's function
          NVAR = 1
          do nf_eff = 1,N_Fl_eff
             nf=Calc_Fl_map(nf_eff)
             call CGR(Z, NVAR, GR(:,:,nf), phi_0(nf_eff), phi_trial(nf_eff))
             Phase_array(nf)=Z
          enddo
          if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)

          call Compute_overlap(Phase_array, Det_Vec, phi_0, phi_trial)
          Det_Vec(:) = Det_Vec(:) * N_SUN
          if (reconstruction_needed) call ham%weight_reconstruction(Det_Vec    )
          if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
          log_O_new = sum(Det_Vec) 
          Phase_new = product(Phase_array)**dble(N_SUN)
          
          Phase=Phase*(Phase_new/Phase_old)

          Overlap_ratio = real(exp(log_O_new-log_O_old), kind(0.d0))

          if ( Overlap_ratio .gt. Zero ) then
              Overlap (i_wlk) = exp(log_O_new)
              weight_k(i_wlk) = weight_k(i_wlk)*Overlap_ratio
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
 

        SUBROUTINE seed_vec_out( Group_Comm )
#ifdef MPI
          Use mpi
#endif
          Implicit none
     
          Integer, INTENT(IN) :: Group_Comm
          
          ! LOCAL
          INTEGER             :: K, ii, intseed(2), seed_tmp(2)
          CHARACTER (LEN=64)  :: FILE_TG, filename
          INTEGER,ALLOCATABLE :: SEED_VEC(:), seed_output(:,:), s_tmp(:,:)

#if defined(MPI)
          INTEGER        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif
          filename = "seedvec_out"
          
          CALL GET_SEED_LEN(K)
          ALLOCATE(SEED_VEC(K))
          ALLOCATE(seed_output(K,isize_g))
          ALLOCATE(s_tmp(K,isize_g))
          seed_output(:,:)=0
          CALL RANGET(SEED_VEC)

          seed_output(:,irank_g+1)=SEED_VEC(:)
          CALL MPI_REDUCE(seed_output,s_tmp,K*isize_g,MPI_Integer,MPI_SUM, 0,Group_comm,IERR)
            
          if ( irank_g .eq. 0 ) then
                OPEN (UNIT = 10, FILE=filename, STATUS='UNKNOWN', ACTION='WRITE')
                do ii = 1, isize_g
                    WRITE(10,*) s_tmp(:,ii)
                enddo
                close(10)
          endif

          deallocate(seed_vec)
          deallocate(seed_output)
          deallocate(s_tmp)

        end SUBROUTINE seed_vec_out

#if defined HDF5
        SUBROUTINE wavefunction_out_hdf5( phi_0, phase_alpha, Group_Comm )
#ifdef MPI
          Use mpi
#endif
#if defined HDF5
          Use hdf5
          use h5lt
#endif
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(IN) :: phi_0
          Complex (Kind=Kind(0.d0)), Dimension(:), Allocatable, INTENT(IN) :: phase_alpha
          Integer, INTENT(IN) :: Group_Comm

          ! LOCAL
          CHARACTER (LEN=64) :: FILE_TG, filename
          Complex (Kind=Kind(0.d0)), pointer :: phi0_out(:,:,:,:)
          Complex (Kind=Kind(0.d0)), pointer :: phasef_out(:)
          Real    (Kind=Kind(0.d0)), pointer :: weight_out(:)
          Complex (Kind=Kind(0.d0)), allocatable :: pf_tmp(:), p0_tmp(:,:,:,:), p1_tmp(:,:,:,:)
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
              allocate(weight_out(N_wlk_mpi), phasef_out(N_wlk_mpi))
          endif

          allocate(p0_tmp(ndim,n_part,n_fl_eff,n_wlk))
          allocate(p1_tmp(ndim,n_part,n_fl_eff,n_wlk))
          allocate(wt_tmp(N_wlk))
          allocate(pf_tmp(N_wlk))

          do nf = 1, N_FL_eff
          do nw = 1, N_wlk
              p0_tmp(:,:,nf,nw)=phi_0(nf,nw)%U(:,:)
          enddo
          enddo
          
          if ( irank_g .ne. 0 ) then
              call mpi_send(phase_alpha,N_wlk,mpi_complex16, 0, 0, MPI_COMM_WORLD,IERR)
              call mpi_send(weight_k   ,N_wlk,    mpi_real8, 0, 1, MPI_COMM_WORLD,IERR)
              Ndt=N_FL_eff*N_wlk*ndim*n_part
              call mpi_send(p0_tmp     ,  Ndt,mpi_complex16, 0, 2, MPI_COMM_WORLD,IERR)
          else
             do ii = 1, isize_g-1
                i_st=ii*N_wlk+1
                i_ed=(ii+1)*N_wlk

                call mpi_recv(pf_tmp,N_wlk,mpi_complex16, ii, 0, MPI_COMM_WORLD,STATUS,IERR)
                phasef_out(i_st:i_ed) = pf_tmp(:)
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
              weight_out(i_st:i_ed) = weight_k(:)
              phasef_out(i_st:i_ed) = phase_alpha(:)
              phi0_out(:,:,:,i_st:i_ed)=p0_tmp
          endif

          if ( irank .eq. 0 ) then

          inquire (file=filename, exist=file_exists)
          IF (.not. file_exists) THEN
              CALL h5open_f(ierr)
              CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdferr)
              
              !Create and write dataset for field phase
              dset_name = "phasef"
              rank = 2
              allocate( dims(2), dimsc(2) )
              dims  = [2, N_wlk_mpi]
              dimsc = dims
              dimsc(rank) = 1
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
              dat_ptr = C_LOC(phasef_out(1))
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
              dimsc(rank) = 1
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
              dims  = [2,ndim,n_part,N_FL,N_wlk_mpi]
              dimsc = dims
              dimsc(rank) = 1
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

              !open and write field phase
              dset_name = "phasef"
              !Open the  dataset.
              CALL h5dopen_f(file_id, dset_name, dset_id, hdferr)
              dat_ptr = C_LOC(phasef_out(1))
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
             deallocate(phi0_out, weight_out, phasef_out)
         endif
         deallocate(p0_tmp, p1_tmp, wt_tmp, pf_tmp)

        END SUBROUTINE wavefunction_out_hdf5
#endif
        
    end Module stepwlk_mod
