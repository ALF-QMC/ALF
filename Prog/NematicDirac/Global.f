        ! Functions for Global moves.  These move are not implemented in this example.
        Subroutine Global_move(T0_Proposal_ratio,nsigma_old,size_clust)
          
          !>  The input is the field nsigma declared in this module. This routine generates a 
          !>  global update with  and returns the propability  
          !>  T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)  
          !>   
          Implicit none
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
          Integer, dimension(:,:),  allocatable, intent(in)  :: nsigma_old
          Integer :: size_cluster
          
          call Wolff_cluster_start(size_cluster, T0_Proposal_ratio, nsigma_old)
!           call Geo_cluster_start(size_cluster)
          
          Write(6,*) "cluster finished. Size,Size/N_Spins:", size_cluster, dble(size_cluster)/dble(N_ising*Ltrot)
          size_clust = dble(size_cluster)/dble(N_ising*Ltrot)
!           T0_Proposal_ratio = 1
        End Subroutine Global_move
        
        Subroutine Wolff_cluster_start(size_cluster, T0_Proposal_ratio, nsigma_old)
          implicit none
          Integer, intent(out) :: size_cluster
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio
          Integer, dimension(:,:),  allocatable, intent(in)  :: nsigma_old
          Integer :: n0, nt0, sigma_clust
          Integer, dimension(:), allocatable :: list_R, list_t
!           Real (Kind=Kind(0.d0)) :: Delta_S0
          
          Allocate( list_R(N_ising*Ltrot), list_t(N_ising*Ltrot) )
          
          size_cluster = 0
          n0  = ceiling( dble(N_ising) * RANF_WRAP() )
          nt0 = ceiling( dble(Ltrot)  * RANF_WRAP() )
          sigma_clust = nsigma(n0,nt0)
          
          size_cluster = size_cluster + 1
          list_R(size_cluster) = n0
          list_t(size_cluster) = nt0
          call Wolff_cluster_add(n0, nt0, sigma_clust, size_cluster, list_R, list_t)
          
          T0_Proposal_ratio = Wolff_T0_Proposal_ratio(size_cluster, sigma_clust, nsigma_old, list_R, list_t)
!           Delta_S0 = Delta_S0_global(nsigma_old)
!           Write(6,*) T0_Proposal_ratio*Delta_S0
          
        End Subroutine Wolff_cluster_start
        
        Function Wolff_T0_Proposal_ratio(size_cluster, sigma_clust, nsigma_old, list_R, list_t)
          implicit none
          Real (Kind=Kind(0.d0)) :: Wolff_T0_Proposal_ratio
          Integer, intent(in) :: size_cluster, sigma_clust
          Integer, dimension(:,:), allocatable, intent(in) :: nsigma_old
          Integer, dimension(:), allocatable, intent(in) :: list_R, list_t
          Integer :: x, n, nt, i, n1, nt1, n_bond_space, n_bond_tau
          n_bond_space = 0
          n_bond_tau = 0
          
          do x = 1, size_cluster
            n  = list_R(x)
            nt = list_t(x)
            
            do i = 1,4
              n1 = Ising_nnlist(n,i)
              If ( nsigma(n1,nt) == nsigma_old(n1,nt) ) then
                If ( nsigma(n1,nt) == sigma_clust ) then
                  n_bond_space = n_bond_space+1
                else
                  n_bond_space = n_bond_space-1
                endif
              Endif
            enddo
            
            nt1 = nt +1 
            if (nt1 > Ltrot) nt1 = 1
            If ( nsigma(n,nt1) == nsigma_old(n,nt1) ) then
              if ( nsigma(n,nt1) == sigma_clust ) then
                n_bond_tau = n_bond_tau+1
              else
                n_bond_tau = n_bond_tau-1
              endif
            Endif
          
            nt1 = nt - 1 
            if (nt1 < 1  ) nt1 = Ltrot
            If ( nsigma(n,nt1) == nsigma_old(n,nt1) ) then
              if ( nsigma(n,nt1) == sigma_clust ) then
                n_bond_tau = n_bond_tau+1
              else
                n_bond_tau = n_bond_tau-1
              endif
            Endif
          enddo
          
          Wolff_T0_Proposal_ratio = (1-addProb_space)**(-n_bond_space) * (1-addProb_tau)**(-n_bond_tau)
        
        End Function Wolff_T0_Proposal_ratio
        
        Recursive Subroutine Wolff_cluster_add(n, nt, sigma_clust, size_cluster, list_R, list_t)
          Implicit none
          Integer, intent(in) :: n, nt, sigma_clust
          Integer :: i, n1, nt1
          Integer, intent(inout) :: size_cluster
          Integer, dimension(:), allocatable, intent(inout) :: list_R, list_t
          
          nsigma(n,nt) = -nsigma(n,nt)
          
          do i = 1,4
            n1 = Ising_nnlist(n,i)
            If ( ( nsigma(n1,nt) == sigma_clust ) .and. ( addProb_space > RANF_WRAP() ) ) then
              size_cluster = size_cluster + 1
              list_R(size_cluster) = n1
              list_t(size_cluster) = nt
              call Wolff_cluster_add(n1, nt, sigma_clust, size_cluster, list_R, list_t)
            Endif
          enddo
          
          nt1 = nt +1 
          if (nt1 > Ltrot) nt1 = 1
          If ( ( nsigma(n,nt1) == sigma_clust ) .and. ( addProb_tau > RANF_WRAP() ) ) then
            size_cluster = size_cluster + 1
            list_R(size_cluster) = n
            list_t(size_cluster) = nt1
            call Wolff_cluster_add(n, nt1, sigma_clust, size_cluster, list_R, list_t)
          Endif
          
          nt1 = nt - 1 
          if (nt1 < 1  ) nt1 = Ltrot
          If ( ( nsigma(n,nt1) == sigma_clust ) .and. ( addProb_tau > RANF_WRAP() ) ) then
            size_cluster = size_cluster + 1
            list_R(size_cluster) = n
            list_t(size_cluster) = nt1
            call Wolff_cluster_add(n, nt1, sigma_clust, size_cluster, list_R, list_t)
          Endif
        End Subroutine Wolff_cluster_add
        
        
        Subroutine Geo_cluster_start(size_cluster)
          Implicit none
          Integer, intent(out) :: size_cluster
          Integer :: i1, i1t, i2, i2t, sigma_i1, sigma_i2
          Integer :: j1, j1t, j2, j2t, i
          Logical, allocatable, dimension(:,:) :: Geo_cluster
        
          size_cluster = 0
          
          ! Start by randomly selecting two space-time-points (i1, i1t) and (i2, i2t)
          ! The center of these two points becomes the symmetry center of the move
          i1  = ceiling( dble(Latt%N) * RANF_WRAP() )
          i1t = ceiling( dble(Ltrot)  * RANF_WRAP() )
          i2  = ceiling( dble(Latt%N) * RANF_WRAP() )
          i2t = ceiling( dble(Ltrot)  * RANF_WRAP() )
          
          sigma_i1 = nsigma(i1, i1t)
          sigma_i2 = nsigma(i2, i2t)
          If(sigma_i1 == sigma_i2) return
          
          Allocate (Geo_cluster(Latt%N,Ltrot))
!           Allocate (Geo_cluster_test(Latt%N,Ltrot))
          Geo_cluster(:,:) = .false.
!           Geo_cluster_test(:,:) = .false.
          
          !Defining symmetry
          R_init(1) = Latt%List(i1,1) + Latt%List(i2,1)
          R_init(2) = Latt%List(i1,2) + Latt%List(i2,2)
          if ( R_init(1) < -(L1-1)/2  ) R_init(1) = R_init(1) + L1
          if ( R_init(1) >   L1/2     ) R_init(1) = R_init(1) - L1
          if ( R_init(2) < -(L2-1)/2  ) R_init(2) = R_init(2) + L2
          if ( R_init(2) >   L2/2     ) R_init(2) = R_init(2) - L2
          Write(6,*) "Symmetry-defining vector", R_init(1), R_init(2)

          nsigma(i1, i1t) = sigma_i2
          nsigma(i2, i2t) = sigma_i1
          Geo_cluster(i1, i1t) = .true.
          Geo_cluster(i2, i2t) = .true.
          size_cluster = 2
          
          do i = 1,4
            j1 = Ising_nnlist(i1,i)
            j2 = Find_geo_partner_space(j1)
            call Geo_cluster_tryadd(sigma_i1, j1, i1t, j2, i2t, .false., size_cluster, Geo_cluster)
          enddo
          
          j1t = i1t + 1
          if (j1t > Ltrot) j1t = 1
          j2t = i2t - 1
          if (j2t < 1) j2t = Ltrot
          call Geo_cluster_tryadd(sigma_i1, i1, j1t, i2, j2t, .true., size_cluster, Geo_cluster)

          j1t = i1t - 1
          if (j1t < 1) j1t = Ltrot
          j2t = i2t + 1
          if (j2t > Ltrot) j2t =1
          call Geo_cluster_tryadd(sigma_i1, i1, j1t, i2, j2t, .true., size_cluster, Geo_cluster)
          
          deallocate(Geo_cluster)
        End Subroutine Geo_cluster_start
        
        Integer Function Find_geo_partner_space(j1)
          Implicit none
          Integer, intent(in)  :: j1
          Integer :: R_j2(2)
          Integer :: test(2)
          
          R_j2(1) = R_init(1) - Latt%List(j1,1)
          R_j2(2) = R_init(2) - Latt%List(j1,2)
          if ( R_j2(1) < -(L1-1)/2  ) R_j2(1) = R_j2(1) + L1
          if ( R_j2(1) >   L1/2     ) R_j2(1) = R_j2(1) - L1
          if ( R_j2(2) < -(L2-1)/2  ) R_j2(2) = R_j2(2) + L2
          if ( R_j2(2) >   L2/2     ) R_j2(2) = R_j2(2) - L2
          Find_geo_partner_space = Latt%invlist(R_j2(1),R_j2(2))
          
          test(1) = R_init(1) - (Latt%list(j1,1) + Latt%list(Find_geo_partner_space,1))
          If ( mod(test(1),L1) .ne. 0 ) then
            Write(6,*) "ERROR: this is not zero:", test(1)
          endif
          test(2) = R_init(2) - (Latt%list(j1,2) + Latt%list(Find_geo_partner_space,2))
          If (  mod(test(2),L2) .ne. 0 ) then
            Write(6,*) "ERROR: this is not zero:", test(2)
          endif
        End Function Find_geo_partner_space
        
        
        Recursive Subroutine Geo_cluster_tryadd(sigma_i1, j1, j1t, j2, j2t, in_tau, size_cluster, Geo_cluster)
          Implicit none
          Integer, intent(in) :: sigma_i1, j1, j1t, j2, j2t
          Logical, intent(in) :: in_tau
          Integer, intent(inout) :: size_cluster
          Logical, allocatable, dimension(:,:), intent(inout) :: Geo_cluster
          Integer :: sigma_j1, sigma_j2
          Integer :: i, k1, k1t, k2, k2t
          
          If( Geo_cluster(j1, j1t) ) return
          
          sigma_j1 = nsigma(j1, j1t)
          sigma_j2 = nsigma(j2, j2t)
          If(sigma_j1 == sigma_j2) return

          If( in_tau ) then
            If ( RANF_WRAP() > Geo_addProb_tau   ) return
          else
            If ( RANF_WRAP() > Geo_addProb_space ) return
          Endif
          size_cluster = size_cluster + 2
          
          nsigma(j1, j1t) = sigma_j2
          nsigma(j2, j2t) = sigma_j1
          Geo_cluster(j1, j1t) = .true.
          Geo_cluster(j2, j2t) = .true.
          
          Do i = 1,4
            k1 = Ising_nnlist(j1,i)
            k2 = Find_geo_partner_space(k1)
            call Geo_cluster_tryadd(sigma_j1, k1, j1t, k2, j2t, .false., size_cluster, Geo_cluster)
          enddo

          k1t = j1t + 1
          if (k1t > Ltrot) k1t = 1
          k2t = j2t - 1
          if (k2t < 1) k2t = Ltrot
          call Geo_cluster_tryadd(sigma_j1, j1, k1t, j2, k2t, .true., size_cluster, Geo_cluster)

          k1t = j1t - 1
          if (k1t < 1) k1t = Ltrot
          k2t = j2t + 1
          if (k2t > Ltrot) k2t =1
          call Geo_cluster_tryadd(sigma_j1, j1, k1t, j2, k2t, .true., size_cluster, Geo_cluster)  
          
        End Subroutine Geo_cluster_tryadd
!========================================================================
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none 
          
          !> Arguments
          Integer, dimension(:,:), allocatable, intent(IN) :: Nsigma_old
          !> Local
          Integer :: I, nt, nt1, I1, I2, nc_J, nc_h_p, nc_h_m, N
          
          Delta_S0_global = 1.D0
          nc_J = 0
          nc_h_p = 0
          nc_h_m = 0
          
          If ( Model == "NematicDirac") then
            N = Latt%N
          elseif ( Model == "NematicDirac3") then
            N = Latt%N*2
          elseif ( Model == "NematicDirac2") then
            N = Latt%N
          elseif ( Model == "yyhe") then
            N = Latt%N
          else
            Write(6,*) "Error in Delta_S0_global: Model not yet implemented!"
            Stop
          endif
          
          Do I = 1,N
            If ( I > Latt%N ) then
              I1 = Latt%nnlist(I - Latt%N,1,0) + Latt%N
              I2 = Latt%nnlist(I - Latt%N,0,1) + Latt%N
            else
              I1 = Latt%nnlist(I,1,0)
              I2 = Latt%nnlist(I,0,1)
            endif
            Do nt = 1,Ltrot
               nt1 = nt + 1
               if (nt == Ltrot) nt1 = 1
               if (nsigma(I,nt) == nsigma(I,nt1) ) then 
                 nc_h_p = nc_h_p + 1
               else
                 nc_h_m = nc_h_m + 1
               endif
               if (nsigma_old(I,nt) == nsigma_old(I,nt1) ) then 
                 nc_h_p = nc_h_p - 1
               else
                 nc_h_m = nc_h_m - 1
               endif
               
               nc_J = nc_J + nsigma(I,nt)*nsigma(I1,nt) &
                    &      + nsigma(I,nt)*nsigma(I2,nt) &
                    &      - nsigma_old(I,nt)*nsigma_old(I1,nt) &
                    &      - nsigma_old(I,nt)*nsigma_old(I2,nt)
            enddo
          enddo
          
          Delta_S0_global = ( sinh(Dtau*Ham_h)**nc_h_m ) * (cosh(Dtau*Ham_h)**nc_h_p) &
                  &         * exp( Dtau * Ham_J*real(nc_J,kind(0.d0)))
          
        end Function Delta_S0_global
        
!---------------------------------------------------------------------
        Subroutine Global_move_tau(T0_Proposal_ratio, S0_ratio, &
             &                     Flip_list, Flip_length,Flip_value,ntau)

!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief 
!> On input: 
!> GR(tau,m) as defined in  Global_tau_mod_PlaceGR and the direction of updating scheme
!> direction=u --> You are visiting the time slices from tau = 1  to tau =Ltrot
!> direction=d --> You are visiting the time slices from tau = Ltrot to tau = 1
!> 
!> On input the field configuration is in the array nsigma.
!> On output: 
!> Flip_list   ::  A list of spins that are to be fliped. Refers to the entires  in OP_V
!> Flip_values ::  The values of the fliped spins
!> Flip_length ::  The number of flips. The first Flip_length entries of Flip_list and Flip_values are relevant
!> S0_ratio          = e^( S_0(sigma_new) ) / e^( S_0(sigma) )
!> T0_Proposal_ratio = T0( sigma_new -> sigma ) /  T0( sigma -> sigma_new)  
!> T0_proposal       = T0 ( sigma -> sigma_new )
!--------------------------------------------------------------------
          
          Implicit none 
          Real (Kind= kind(0.d0)), INTENT(INOUT) :: T0_Proposal_ratio,  S0_ratio
          Integer,    allocatable, INTENT(INOUT) :: Flip_list(:), Flip_value(:)
          Integer, INTENT(INOUT) :: Flip_length
          Integer, INTENT(IN)    :: ntau

        end Subroutine Global_move_tau
        
!---------------------------------------------------------------------
        Subroutine  Hamiltonian_set_random_nsigma
          
          ! The user can set the initial configuration
          
          Implicit none
          
          Integer :: I, nt
          
          Do nt = 1,Ltrot
             Do I = 1,Size(OP_V,1)
                nsigma(I,nt)  = 1
                if ( ranf_wrap()  > 0.5D0 ) nsigma(I,nt)  = -1
             enddo
          enddo
          
        end Subroutine Hamiltonian_set_random_nsigma
