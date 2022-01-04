  ! Functions for Global moves.  These move are not implemented in this example.
  Subroutine Global_move(T0_Proposal_ratio,nsigma_old,size_clust)

    !>  The input is the field nsigma declared in this module. This routine generates a
    !>  global update with  and returns the propability
    !>  T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)
    !>
    Implicit none
    Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
!     Integer, dimension(:,:),  allocatable, intent(in)  :: nsigma_old
    type(fields), intent(in) :: nsigma_old
    Integer :: size_cluster

    Integer :: i, nt

    n_global = n_global + 1

    select case( Global_type )
    case ( 'Wolff' )
      call Wolff_cluster_start(size_cluster, T0_Proposal_ratio, nsigma_old)
    case ( 'Geo' )
      call Geo_cluster_start(size_cluster)
    case ( 'switch' )
      size_clust = latt%N*2
      T0_Proposal_ratio = 1
      Do I = 1, latt%N
        Do nt = 1,Ltrot
          nsigma%f(I,nt) = nsigma_old%f(I+latt%N,nt)
          nsigma%f(I+latt%N,nt) = nsigma_old%f(I,nt)
        enddo
      enddo
      if ( mod(n_global,2) == 0 ) nsigma%f(:,:) = -nsigma_old%f(:,:)
    case ( 'flip' )
      T0_Proposal_ratio = 1
      nsigma%f(:,:) = -nsigma_old%f(:,:)
    case default
      write(*,*) "Error in Global_move: Unknown Global_type:", Global_type
      stop
    end select

    write(*,*) Delta_S0_global(nsigma_old) * T0_Proposal_ratio

    !Write(6,*) "cluster finished. Size,Size/N_Spins:", size_cluster, dble(size_cluster)/dble(N_ising*Ltrot)
    size_clust = dble(size_cluster)/dble(N_ising*Ltrot)
    !T0_Proposal_ratio = 1
  End Subroutine Global_move

  Subroutine Wolff_cluster_start(size_cluster, T0_Proposal_ratio, nsigma_old)
    implicit none
    Integer, intent(out) :: size_cluster
    Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio
!     Integer, dimension(:,:),  allocatable, intent(in)  :: nsigma_old
    type(fields), intent(in) :: nsigma_old

    Integer :: n0, nt0, sigma_clust
    Integer, dimension(:), allocatable :: list_R, list_t
    !Real (Kind=Kind(0.d0)) :: Delta_S0

    Allocate( list_R(N_ising*Ltrot), list_t(N_ising*Ltrot) )

    size_cluster = 0
    n0  = ceiling( dble(N_ising) * RANF_WRAP() )
    nt0 = ceiling( dble(Ltrot)   * RANF_WRAP() )
    sigma_clust = nsigma%i(n0,nt0)

    size_cluster = size_cluster + 1
    list_R(size_cluster) = n0
    list_t(size_cluster) = nt0
    call Wolff_cluster_add(n0, nt0, sigma_clust, size_cluster, list_R, list_t)

    T0_Proposal_ratio = Wolff_T0_Proposal_ratio(size_cluster, sigma_clust, nsigma_old, list_R, list_t)
    !Delta_S0 = Delta_S0_global(nsigma_old)
    !Write(6,*) T0_Proposal_ratio*Delta_S0

  End Subroutine Wolff_cluster_start

  Function Wolff_T0_Proposal_ratio(size_cluster, sigma_clust, nsigma_old, list_R, list_t)
    implicit none
    Real (Kind=Kind(0.d0)) :: Wolff_T0_Proposal_ratio
    Integer, intent(in) :: size_cluster, sigma_clust
!     Integer, dimension(:,:), allocatable, intent(in) :: nsigma_old
    type(fields), intent(in) :: nsigma_old
    Integer, dimension(:), allocatable, intent(in) :: list_R, list_t
    Integer :: x, n, nt, i, n1, nt1, n_bond_space, n_bond_tau
    n_bond_space = 0
    n_bond_tau = 0

    do x = 1, size_cluster
      n  = list_R(x)
      nt = list_t(x)

      do i = 1,4
        n1 = Ising_nnlist(n,i)
        If ( nsigma%i(n1,nt) == nsigma_old%i(n1,nt) ) then
          If ( nsigma%i(n1,nt) == sigma_clust ) then
            n_bond_space = n_bond_space+1
          else
            n_bond_space = n_bond_space-1
          endif
        Endif
      enddo

      nt1 = nt +1
      if (nt1 > Ltrot) nt1 = 1
      If ( nsigma%i(n,nt1) == nsigma_old%i(n,nt1) ) then
        if ( nsigma%i(n,nt1) == sigma_clust ) then
          n_bond_tau = n_bond_tau+1
        else
          n_bond_tau = n_bond_tau-1
        endif
      Endif

      nt1 = nt - 1
      if (nt1 < 1  ) nt1 = Ltrot
      If ( nsigma%i(n,nt1) == nsigma_old%i(n,nt1) ) then
        if ( nsigma%i(n,nt1) == sigma_clust ) then
          n_bond_tau = n_bond_tau+1
        else
          n_bond_tau = n_bond_tau-1
        endif
      Endif
    enddo

    Wolff_T0_Proposal_ratio = (1-Wolff_addProb_space)**(-n_bond_space) * (1-Wolff_addProb_tau)**(-n_bond_tau)

  End Function Wolff_T0_Proposal_ratio

  Recursive Subroutine Wolff_cluster_add(n, nt, sigma_clust, size_cluster, list_R, list_t)
    Implicit none
    Integer, intent(in) :: n, nt, sigma_clust
    Integer :: i, n1, nt1
    Integer, intent(inout) :: size_cluster
    Integer, dimension(:), allocatable, intent(inout) :: list_R, list_t

    nsigma%f(n,nt) = -nsigma%f(n,nt)

    do i = 1,4
      n1 = Ising_nnlist(n,i)
      If ( ( nsigma%i(n1,nt) == sigma_clust ) .and. ( Wolff_addProb_space > RANF_WRAP() ) ) then
        size_cluster = size_cluster + 1
        list_R(size_cluster) = n1
        list_t(size_cluster) = nt
        call Wolff_cluster_add(n1, nt, sigma_clust, size_cluster, list_R, list_t)
      Endif
    enddo

    nt1 = nt +1
    if (nt1 > Ltrot) nt1 = 1
    If ( ( nsigma%i(n,nt1) == sigma_clust ) .and. ( Wolff_addProb_tau > RANF_WRAP() ) ) then
      size_cluster = size_cluster + 1
      list_R(size_cluster) = n
      list_t(size_cluster) = nt1
      call Wolff_cluster_add(n, nt1, sigma_clust, size_cluster, list_R, list_t)
    Endif

    nt1 = nt - 1
    if (nt1 < 1  ) nt1 = Ltrot
    If ( ( nsigma%i(n,nt1) == sigma_clust ) .and. ( Wolff_addProb_tau > RANF_WRAP() ) ) then
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
    i1  = ceiling( dble(latt%N) * RANF_WRAP() )
    i1t = ceiling( dble(Ltrot)  * RANF_WRAP() )
    i2  = ceiling( dble(latt%N) * RANF_WRAP() )
    i2t = ceiling( dble(Ltrot)  * RANF_WRAP() )

    sigma_i1 = nsigma%i(i1, i1t)
    sigma_i2 = nsigma%i(i2, i2t)
    If(sigma_i1 == sigma_i2) return

    Allocate (Geo_cluster(latt%N,Ltrot))
    !Allocate (Geo_cluster_test(latt%N,Ltrot))
    Geo_cluster(:,:) = .false.
    !Geo_cluster_test(:,:) = .false.

    !Defining symmetry
    R_init(1) = latt%List(i1,1) + latt%List(i2,1)
    R_init(2) = latt%List(i1,2) + latt%List(i2,2)
    if ( R_init(1) < -(L1-1)/2  ) R_init(1) = R_init(1) + L1
    if ( R_init(1) >   L1/2     ) R_init(1) = R_init(1) - L1
    if ( R_init(2) < -(L2-1)/2  ) R_init(2) = R_init(2) + L2
    if ( R_init(2) >   L2/2     ) R_init(2) = R_init(2) - L2
    Write(6,*) "Symmetry-defining vector", R_init(1), R_init(2)

    nsigma%f(i1, i1t) = sigma_i2
    nsigma%f(i2, i2t) = sigma_i1
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

    R_j2(1) = R_init(1) - latt%List(j1,1)
    R_j2(2) = R_init(2) - latt%List(j1,2)
    if ( R_j2(1) < -(L1-1)/2  ) R_j2(1) = R_j2(1) + L1
    if ( R_j2(1) >   L1/2     ) R_j2(1) = R_j2(1) - L1
    if ( R_j2(2) < -(L2-1)/2  ) R_j2(2) = R_j2(2) + L2
    if ( R_j2(2) >   L2/2     ) R_j2(2) = R_j2(2) - L2
    Find_geo_partner_space = latt%invlist(R_j2(1),R_j2(2))

    test(1) = R_init(1) - (latt%list(j1,1) + latt%list(Find_geo_partner_space,1))
    If ( mod(test(1),L1) .ne. 0 ) then
      Write(6,*) "ERROR: this is not zero:", test(1)
    endif
    test(2) = R_init(2) - (latt%list(j1,2) + latt%list(Find_geo_partner_space,2))
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

    sigma_j1 = nsigma%i(j1, j1t)
    sigma_j2 = nsigma%i(j2, j2t)
    If(sigma_j1 == sigma_j2) return

    If( in_tau ) then
      If ( RANF_WRAP() > Geo_addProb_tau   ) return
    else
      If ( RANF_WRAP() > Geo_addProb_space ) return
    Endif
    size_cluster = size_cluster + 2

    nsigma%f(j1, j1t) = sigma_j2
    nsigma%f(j2, j2t) = sigma_j1
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
  Real (Kind=kind(0.d0)) Function Delta_S0_global(nsigma_old)

    !>  This function computes the ratio:  e^{-S0(nsigma%f)}/e^{-S0(nsigma_old)}
    Implicit none

    !> Arguments
    type(fields), intent(in) :: nsigma_old
!     Integer, dimension(:,:), allocatable, intent(IN) :: Nsigma_old
    !> Local
    Integer :: I, nt, nt1, I1, I2, nc_J, nc_h_p, nc_h_m, N

    Delta_S0_global = 1.D0
    nc_J = 0
    nc_h_p = 0
    nc_h_m = 0

    select case( Model_vers )
    case( 0,1 )
      N = latt%N
    case( 3 )
      N = latt%N*2
    case( 2 )
      N = latt%N
    case default
      Write(6,*) "Error in Delta_S0_global: Model not yet implemented!"
      Stop
    end select

    Do I = 1,N
      If ( I > latt%N ) then
        I1 = latt%nnlist(I - latt%N,1,0) + latt%N
        I2 = latt%nnlist(I - latt%N,0,1) + latt%N
      else
        I1 = latt%nnlist(I,1,0)
        I2 = latt%nnlist(I,0,1)
      endif
      Do nt = 1,Ltrot
          nt1 = nt + 1
          if (nt == Ltrot) nt1 = 1
          if (nsigma%i(I,nt) == nsigma%i(I,nt1) ) then
            nc_h_p = nc_h_p + 1
          else
            nc_h_m = nc_h_m + 1
          endif
          if (nsigma_old%i(I,nt) == nsigma_old%i(I,nt1) ) then
            nc_h_p = nc_h_p - 1
          else
            nc_h_m = nc_h_m - 1
          endif

          nc_J = nc_J + nsigma%i(I,nt)*nsigma%i(I1,nt) &
              &      + nsigma%i(I,nt)*nsigma%i(I2,nt) &
              &      - nsigma_old%i(I,nt)*nsigma_old%i(I1,nt) &
              &      - nsigma_old%i(I,nt)*nsigma_old%i(I2,nt)
      enddo
    enddo

    Delta_S0_global = ( sinh(Dtau*Ham_h)**nc_h_m ) * (cosh(Dtau*Ham_h)**nc_h_p) &
            &         * exp( Dtau * Ham_J*real(nc_J,kind(0.d0)))

  end Function Delta_S0_global
