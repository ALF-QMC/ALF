    Module stepwlk_mod
        Use Hamiltonian_main 
        Use Operator_mod
        Use Control
        Use Hop_mod
        Use UDV_State_mod
        use cgr1_mod
        
        Contains

        SUBROUTINE stepwlk( phi_trial, phi_0, GR, phase, phase_alpha, i_wlk, N_op, ntau_qr, ntau_bp ) 
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:), ALLOCATABLE, INTENT(IN)    :: phi_trial
          CLASS(UDV_State), Dimension(:), ALLOCATABLE, INTENT(INOUT) :: phi_0
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable, INTENT(INOUT) :: GR
          COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: phase, phase_alpha
          Integer, INTENT(IN) :: N_op, i_wlk, ntau_qr, ntau_bp

          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR
          Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot, Overlap_old, Overlap_new
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
          Real (Kind=Kind(0.d0))    :: Zero = 1.0E-8
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable :: Phase_array
          REAL    (Kind=Kind(0.d0)), Dimension(:), Allocatable :: Det_Vec

          allocate(Phase_array(N_FL), Det_Vec(N_FL))

          if ( weight_k(i_wlk) .gt. Zero ) then
                
              ! update weight by fac_norm
              weight_k(i_wlk)=weight_k(i_wlk)*exp(fac_norm);
            
              call Compute_overlap(Phase_array, Det_Vec, phi_0, phi_trail)
              Det_Vec(:) = Det_Vec(:) * N_SUN
              if (reconstruction_needed) call ham%weight_reconstruction(Det_Vec    )
              if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
              Overlap_old = exp(sum(Det_Vec)) 

              ! Kinetic part
              Do nf_eff = 1,N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 Call Hop_mod_mmthr(phi_0(nf_eff)%U,nf,1)
              enddo
              
              ! Update Green's function
              NVAR = 1
              do nf_eff = 1,N_Fl_eff
                 nf=Calc_Fl_map(nf_eff)
                 call CGR(Z, NVAR, GR(:,:,nf), phi_0(nf_eff), phi_trail(nf_eff))
                 Phase_array(nf)=Z
              enddo
              if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
              Phase=product(Phase_array)
              Phase=Phase**N_SUN

              call Compute_overlap(Phase_array, Det_Vec, phi_0, phi_trail)
              Det_Vec(:) = Det_Vec(:) * N_SUN
              if (reconstruction_needed) call ham%weight_reconstruction(Det_Vec    )
              if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
              Overlap_new = exp(sum(Det_Vec)) 
              
              Phase=product(Phase_array)
              Phase=Phase**N_SUN

              Overlap_ratio = real(Overlap_new/Overlap_old, kind(0.d0))

              if ( Overlap_ratio .gt. Zero ) then
                  Overlap (i_wlk) = Overlap_new
                  weight_k(i_wlk) = weight_k(i_wlk)*Overlap_ratio
              else
                  weight_k(i_wlk) = 0.d0
              endif

          endif

          ! propagate with interaction
          do n = 1, N_op

          if ( weight_k(i_wlk) .gt. Zero ) then
             
             ! upgrade Green's function
             N_type = 2
             do nf_eff = 1,N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                Call Op_Wrapdo( GR(:,:,nf), Op_V(n,nf), 1.d0, Ndim, N_Type,1)
             enddo

             Call Upgrade(GR,n,PHASE, PHASE_alpha, spin, i_wlk )
             nsigma_qr(i_wlk)%f(n,ntau_qr) = spin
             nsigma_bp(i_wlk)%f(n,ntau_bp) = spin

             N_type = 1
             do nf_eff = 1,N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), spin, Ndim, N_Type, 1 )
             enddo

             ! propagate slater determinant
             Do nf_eff = 1,N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 Call Op_mmultR(phi_0(nf_eff)%U,Op_V(n,nf),spin,'n',nt)
             enddo

          endif

          enddo

        END SUBROUTINE stepwlk


        subroutine re_orthonormalize_walkers(Phi_0, i_wlk)
          
          CLASS(UDV_State), Dimension(:), ALLOCATABLE, INTENT(INOUT) :: phi_0
          
          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR, N_size, I
          Real    (Kind=Kind(0.d0)) :: Overlap_ratio, Zero = 1.0E-8
          Complex (Kind=Kind(0.d0)) :: Det_D(N_FL)

          if ( weight_k(i_wlk) .gt. Zero ) then

              N_size = udvl%ndim

              !Carry out U,D,V decomposition.
              Do nf_eff = 1, N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 CALL Phi_0(nf_eff)%decompose
              enddo
              
              Det_D = cmplx(1.d0, 0.d0, kind(0.d0))

              Do nf_eff = 1, N_FL_eff
                 nf=Calc_Fl_map(nf_eff)
                 DO I = 1,N_size
                    Det_D(nf)=Det_D(nf)*Phi_0(nf_eff)%D(I)
                 ENDDO
              enddo
              
              if (reconstruction_needed) call ham%weight_reconstruction(Det_D)

              Overlap(i_wlk)=Overlap(i_wlk)/real(product(Det_D),0.d0,kind(0.d0))
              
              Phi_0(nf_eff)%D(:) = cmplx(1.d0, 0.d0, kind(0.d0))

          endif

        END SUBROUTINE re_orthonormalize_walkers
        
        SUBROUTINE initial_wlk( phi_trial, phi_0, udvr, GR, phase, phase_alpha, i_wlk ) 
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:), ALLOCATABLE, INTENT(INOUT) :: phi_trial, phi_0, udvr
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable, INTENT(INOUT) :: GR
          COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: phase, phase_alpha
          Integer, INTENT(IN) :: i_wlk

          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR
          Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot, Overlap_old, Overlap_new
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new, Overlap_ratio
          Real (Kind=Kind(0.d0))    :: Zero = 1.0E-8
          COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable :: Phase_array
          REAL    (Kind=Kind(0.d0)), Dimension(:), Allocatable :: Det_Vec

          allocate(Phase_array(N_FL), Det_Vec(N_FL))

          do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             CALL udvr(nf_eff)%init(ndim,'r',WF_R(nf)%P)
             CALL phi_trial(nf_eff)%init(ndim,'l',WF_L(nf)%P)
             CALL phi_0(nf_eff)%init(ndim,'r',WF_R(nf)%P)
             
          enddo

          NVAR = 1
          do nf_eff = 1,N_Fl_eff
             nf=Calc_Fl_map(nf_eff)
             call CGR(Z, NVAR, GR(:,:,nf), phi_0(nf_eff), phi_trial(nf_eff))
             Phase_array(nf)=Z
          enddo
          if (reconstruction_needed) call ham%weight_reconstruction(Phase_array)
          Phase=product(Phase_array)
          Phase=Phase**N_SUN

        END SUBROUTINE initial_wlk
        
    end Module stepwlk_mod
