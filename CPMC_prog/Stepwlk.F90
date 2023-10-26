    Module stepwlk_mod
        Use Hamiltonian_main 
        Use Operator_mod
        Use Control
        Use Hop_mod
        Use UDV_State_mod
        use cgr1_mod
        use upgrade_mod
        
        Contains

        SUBROUTINE stepwlk_move( phi_trial, phi_0, GR, phase, phase_alpha, N_op, ntau_qr, ntau_bp ) 
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(IN)    :: phi_trial
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:,:), Allocatable, INTENT(INOUT) :: GR
          COMPLEX (Kind=Kind(0.d0)), Dimension(:)      , Allocatable, INTENT(INOUT) :: phase, phase_alpha
          Integer, INTENT(IN) :: N_op, ntau_qr, ntau_bp

          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8

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
                nsigma_qr(i_wlk)%f(n,ntau_qr) = spin
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

              Phi_0(nf_eff,i_wlk)%D(:) = cmplx(1.d0, 0.d0, kind(0.d0))

              ! update the overlap when normal propagation
              if (cop == 'R') Overlap(i_wlk)=Overlap(i_wlk)/product(Det_D)

          endif

          enddo

        END SUBROUTINE re_orthonormalize_walkers
        
        SUBROUTINE initial_wlk( phi_trial, phi_0, GR, phase )
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_trial, phi_0
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:,:), Allocatable, INTENT(INOUT) :: GR
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable, INTENT(INOUT) :: phase

          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk
          Complex (Kind=Kind(0.d0)) :: Overlap_old, Overlap_new, Z
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
          Real (Kind=Kind(0.d0))    :: Zero = 1.0E-8
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable :: Phase_array

          allocate(Phase_array(N_FL))
          tot_ene    = cmplx(0.d0,0.d0,kind(0.d0))
          tot_weight = 0.d0

          do i_wlk = 1, N_wlk
          
            do nf_eff = 1, N_FL_eff
               nf=Calc_Fl_map(nf_eff)
               CALL phi_trial(nf_eff, i_wlk)%init(ndim,'l',WF_L(nf)%P)
               CALL phi_0(nf_eff, i_wlk)%init(ndim,'r',WF_R(nf)%P)
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
            
            tot_ene    = tot_ene    + ham%E0_local(GR(:,:,:,i_wlk))*weight_k(i_wlk)
            tot_weight = tot_weight + weight_k(i_wlk)

          enddo
          fac_norm= real(tot_ene, kind(0.d0))/tot_weight

        END SUBROUTINE initial_wlk

        SUBROUTINE population_control( phi_0, phase_alpha, phase ) 
#ifdef MPI
        Use mpi
#endif
        Use Random_wrap
        
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable, INTENT(INOUT) :: phase_alpha
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable, INTENT(INOUT) :: phase
          
          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, it_wlk 
          Integer :: j, it, i_t, i_st, i_ed, nu_wlk, i_src, i_wlk, j_src, j_wlk
          Real    (Kind=Kind(0.d0)), Dimension(:), Allocatable :: weight_mpi
          Real    (Kind=Kind(0.d0)) :: Zero = 1.0E-8, d_scal, sum_w, w_count, w_tmp(N_wlk)
          Complex (Kind=Kind(0.d0)) :: Overlap_tmp(N_wlk), phase_alpha_tmp(N_wlk), phase_tmp(N_wlk)
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE :: phi_0_m

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          ! temporary store
          allocate(phi_0_m(N_FL_eff,N_wlk))
          do i_wlk = 1, N_wlk
            do nf_eff = 1, N_FL_eff
               nf=Calc_Fl_map(nf_eff)
               CALL phi_0_m(nf_eff, i_wlk)%init(ndim,'r',WF_R(nf)%P)
            enddo
          enddo

          do nf_eff = 1, N_FL_eff
             do i_wlk = 1, N_wlk
                call phi_0_m(nf_eff,i_wlk)%assign(phi_0(nf_eff,i_wlk))
             enddo
          enddo
          Overlap_tmp=Overlap
          phase_alpha_tmp=phase_alpha
          phase_tmp=phase

          ! population control
          allocate(weight_mpi(N_wlk_mpi))
          if (irank_g == 0) then
              weight_mpi(1:N_wlk)=weight_k(:)
          endif
          if (irank_g == 0) then
              do it = 1, isize_g-1
                 i_st=it*N_wlk+1
                 i_ed=(it+1)*N_wlk

                 call mpi_send(weight_k,N_wlk,MPI_REAL8,0 ,it+1024,Group_comm,IERR)
                 call mpi_recv(w_tmp   ,N_wlk,MPI_REAL8,it,it+1024,Group_comm,IERR)
                 weight_mpi(i_st:i_ed)=w_tmp(:)
              enddo
          endif
          CALL MPI_BCAST(weight_mpi, N_wlk_mpi, MPI_REAL8, 0,MPI_COMM_WORLD,ierr)
          
          d_scal=dble(N_wlk_mpi)/sum(weight_mpi)
          sum_w=-ranf_wrap()
          nu_wlk=0
          do it_wlk=1, N_wlk_mpi
              i_src=(it_wlk-1)/N_wlk
              i_wlk=it_wlk-N_wlk*i_src

              sum_w=sum_w+weight_mpi(it_wlk)*d_scal
              n=ceiling(sum_w);
              do j=(nu_wlk+1), n
                  j_src=(j-1)/N_wlk
                  j_wlk=j-N_wlk*j_src
                  if ( j_src .ne. i_src ) then
                      if ( irank_g .eq. j_src ) then
                        call mpi_send(overlap    (i_wlk),1,MPI_COMPLEX16,j_src,j_src+1024,Group_comm,IERR)
                        call mpi_recv(overlap_tmp(j_wlk),1,MPI_COMPLEX16,i_src,j_src+1024,Group_comm,IERR)

                        do nf_eff = 1, N_FL_eff
                            call phi_0_m(nf_eff,j_wlk)%MPI_sendrecv_general(phi_0(nf_eff,i_wlk),j_src, &
                                & j_src, i_src, j_src, status, ierr)
                        enddo

                        call mpi_send(phase_alpha    (i_wlk),1,MPI_COMPLEX16,j_src,j_src+1024,Group_comm,IERR)
                        call mpi_recv(phase_alpha_tmp(j_wlk),1,MPI_COMPLEX16,i_src,j_src+1024,Group_comm,IERR)
                        
                        call mpi_send(phase    (i_wlk),1,MPI_COMPLEX16,j_src,j_src+1024,Group_comm,IERR)
                        call mpi_recv(phase_tmp(j_wlk),1,MPI_COMPLEX16,i_src,j_src+1024,Group_comm,IERR)
                      endif
                  else
                      if ( irank_g .eq. i_src ) then
                        do nf_eff = 1, N_FL_eff
                            call phi_0_m(nf_eff,j_wlk)%assign(phi_0(nf_eff,i_wlk))
                        enddo
                        Overlap_tmp(j_wlk)=Overlap(i_wlk) 
                        phase_alpha_tmp(j_wlk)=phase_alpha(i_wlk) 
                        phase_tmp(j_wlk)=phase(i_wlk) 
                      endif
                  endif
              enddo
              nu_wlk=n
          enddo

          Overlap=Overlap_tmp
          phase_alpha=phase_alpha_tmp
          phase=phase_tmp
          weight_k(:)=1
          do nf_eff = 1, N_FL_eff
             do i_wlk = 1, N_wlk
                call phi_0(nf_eff,i_wlk)%assign(phi_0_m(nf_eff,i_wlk))
             enddo
          enddo

          ! deallocate tmp udv class
          do i_wlk = 1, N_wlk
            do nf_eff = 1, N_FL_eff
               nf=Calc_Fl_map(nf_eff)
               CALL phi_0_m(nf_eff, i_wlk)%dealloc
            enddo
          enddo

          deallocate(weight_mpi)
          deallocate(phi_0_m)

        END SUBROUTINE population_control
        
        SUBROUTINE store_phi( phi_0, phi_bp_r )
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_0, phi_bp_r

          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk
          
          do nf_eff = 1, N_FL_eff
             do i_wlk = 1, N_wlk
                call phi_bp_r(nf_eff,i_wlk)%assign(phi_0(nf_eff,i_wlk))
             enddo
          enddo

        END SUBROUTINE store_phi

        SUBROUTINE backpropagation( phi_bp_l, nwrap )
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:,:), ALLOCATABLE, INTENT(INOUT) :: phi_bp_l
          Integer, INTENT(IN) :: nwrap

          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, i_wlk, ltrot_bp, N_op
          Complex (Kind=Kind(0.d0)) :: Z
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
          Real (Kind=Kind(0.d0))    :: Zero = 1.0E-8

          N_op     = Size(OP_V,1)
          ltrot_bp = size(nsigma_bp(1)%f, 2)
            
          !! initialization
          do i_wlk = 1, N_wlk
            do nf_eff = 1, N_FL_eff
               nf=Calc_Fl_map(nf_eff)
               CALL phi_bp_l(nf_eff, i_wlk)%init(ndim,'l',WF_L(nf)%P)
            enddo
          enddo

          Do nt = ltrot_bp, 1, -1
            
             Do i_wlk = 1, N_wlk
             Do nf_eff = 1,N_FL_eff
                
                nf=Calc_Fl_map(nf_eff)
                Call Hop_mod_mmthl_1D2(phi_bp_l(nf_eff,i_wlk)%U,nf,1)

                Do n = N_op,1,-1
                   Call Op_mmultL(phi_bp_l(nf_eff,i_wlk)%U,Op_V(n,nf),nsigma_bp(i_wlk)%f(n,nt),'n',1)
                enddo
             
                Call Hop_mod_mmthl_1D2(phi_bp_l(nf_eff,i_wlk)%U,nf,1)

             enddo
             Enddo
             
             if ( ( mod(nt, Nwrap) .eq. 0 ) .or. ( nt .eq. 1 ) ) then
                 call re_orthonormalize_walkers(phi_bp_l, 'L')
             endif

          Enddo

        END SUBROUTINE backpropagation
        
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
        
    end Module stepwlk_mod
