    Module stepwlk_mod
        Use Hamiltonian_main 
        Use Operator_mod
        Use Control
        Use Hop_mod
        Use UDV_State_mod
        use cgr1_mod
        
        Contains

        SUBROUTINE stepwlk( phi_trial, phi_0, GR, i_blk, N_op, ntau ) 
          
          Implicit none
     
          CLASS(UDV_State), Dimension(:), ALLOCATABLE, INTENT(IN)    :: phi_trial
          CLASS(UDV_State), Dimension(:), ALLOCATABLE, INTENT(INOUT) :: phi_0
          COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable, INTENT(INOUT) :: GR
          Integer, INTENT(IN) :: N_op, i_blk
          Integer, INTENT(INOUT) :: NTAU

          !Local 
          Integer :: nf, nf_eff, N_Type, NTAU1, n, m, nt, NVAR
          Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot
          Real    (Kind=Kind(0.d0)) :: S0_ratio, spin, HS_new
          Character (Len=64)        :: Mode
          Logical                   :: Acc, toggle1

          ! Kinetic part
          Do nf_eff = 1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             Call Hop_mod_mmthr(phi_0(nf_eff)%U,nf,1)
          enddo
          
          ! Update Green's function
          NVAR = 1
          do nf_eff = 1,N_Fl_eff
             nf=Calc_Fl_map(nf_eff)
             call CGR(Z, NVAR, GR(:,:,nf, i_blk), phi_0(nf_eff, i_blk), phi_trial(nf_eff, i_blk))
             Phase_array(nf, i_blk)=Z
          enddo
          if (reconstruction_needed) call ham%weight_reconstruction(Phase_array(:,i_blk))
          Phase(i_blk)=product(Phase_array(:,i_blk))
          Phase(i_blk)=Phase(i_blk)**N_SUN

          ! propagate with interaction by upgrade Green's function
          do n = 1, N_op

             N_type = 2
             do nf_eff = 1,N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                Call Op_Wrapdo( GR(:,:,nf), Op_V(n,nf), 1.d0, Ndim, N_Type,ntau)
             enddo

             Call Upgrade(GR,n,ntau,PHASE, spin, S0_ratio )
             nsigma(i_blk)%f(n,ntau) = spin

             N_type = 2
             do nf_eff = 1,N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), 1.d0, Ndim, N_Type, ntau )
             enddo

          enddo

          Ntau = Ntau + 1

        END SUBROUTINE stepwlk
        
    end Module stepwlk_mod
