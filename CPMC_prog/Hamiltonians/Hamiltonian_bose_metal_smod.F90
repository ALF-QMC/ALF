    submodule(Hamiltonian_main) ham_bose_metal_smod

       use Operator_mod
       use WaveFunction_mod
       use Lattices_v3
       use MyMats
       use Random_Wrap
       use Files_mod
       use Matrix
       use Observables
       use Fields_mod
       use Predefined_Hoppings

       implicit none

       type, extends(ham_base) :: ham_bose_metal
       contains
          ! Set Hamiltonian-specific procedures
          procedure, nopass :: Ham_Set
          procedure, nopass :: alloc_obs
          procedure, nopass :: e0_local
          procedure, nopass :: set_xloc
          procedure, nopass :: sum_weight
          procedure, nopass :: update_fac_norm
#ifdef HDF5
          procedure, nopass :: write_parameters_hdf5
#endif
       end type ham_bose_metal

       !#PARAMETERS START# VAR_lattice
       character(len=64) :: Model = "bose_metal"  ! Value not relevant
       character(len=64) :: Lattice_type = "square_anisotropic"
       integer            :: L1 = 6   ! Length in direction a_1
       integer            :: L2 = 6   ! Length in direction a_2
       !#PARAMETERS END#

       !#PARAMETERS START# VAR_Model_Generic
       !Integer              :: N_SUN        = 1        ! Number of colors
       !Integer              :: N_FL         = 2        ! Number of flavors
       !Integer              :: N_hfb        = 1        ! Number of slater on trial wave function
       !Integer              :: N_wlk        = 1        ! Number of walker
       !Integer              :: ltrot        = 10       ! length of imaginary time for dynamical measure
       real(Kind=kind(0.d0)) :: Phi_X = 0.d0     ! Twist along the L_1 direction, in units of the flux quanta
       real(Kind=kind(0.d0)) :: Phi_Y = 0.d0     ! Twist along the L_2 direction, in units of the flux quanta
       logical               :: Bulk = .true.   ! Twist as a vector potential (.T.), or at the boundary (.F.)
       integer               :: N_Phi = 0        ! Total number of flux quanta traversing the lattice
       real(Kind=kind(0.d0)) :: Dtau = 0.1d0    ! Thereby Ltrot=Beta/dtau
       real(Kind=kind(0.d0)) :: Beta = 5.d0     ! Inverse temperature
       logical               :: Checkerboard = .true.   ! Whether checkerboard decomposition is used
       !logical              :: Symm         = .true.   ! Whether symmetrization takes place
       !#PARAMETERS END#

       !#PARAMETERS START# VAR_bose_metal
       real(Kind=kind(0.d0)) :: ham_t     = 1.d0     ! Hopping parameter
       real(Kind=kind(0.d0)) :: ham_alpha = 1.d0     ! Hopping parameter
       real(Kind=kind(0.d0)) :: ham_chem  = 0.d0     ! Chemical potential
       real(Kind=kind(0.d0)) :: ham_U     = 4.d0     ! attractive Hubbard interaction
       integer               :: N_dope    = 0
       !#PARAMETERS END#

       type(Lattice), target :: Latt
       type(Unit_cell), target :: Latt_unit
       type(Hopping_Matrix_type), allocatable :: Hopping_Matrix(:)
       integer, allocatable :: List(:, :), Invlist(:, :)  ! For orbital structure of Unit cell
       Integer, allocatable :: rot_del_list(:,:,:), rot_list(:,:), del_list(:,:,:)
       complex (Kind=Kind(0.d0)), allocatable :: ff_s(:,:,:)

    contains

       module subroutine Ham_Alloc_bose_metal
          allocate (ham_bose_metal::ham)
       end subroutine Ham_Alloc_bose_metal

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_bose_metal_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
       subroutine Ham_Set

#if defined (MPI) || defined(TEMPERING)
          use mpi
#endif
          implicit none

          integer                :: ierr, nf, unit_info
          character(len=64)     :: file_info

          ! L1, L2, Lattice_type, List(:,:), Invlist(:,:) -->  Lattice information
          ! Ham_T, Chem, Phi_X, XB_B, Checkerboard, Symm   -->  Hopping
          ! Interaction                              -->  Model
          ! Simulation type                          -->  Finite  T or Projection  Symmetrize Trotter.

#ifdef MPI
          integer        :: Isize, Irank, irank_g, isize_g, igroup
          integer        :: STATUS(MPI_STATUS_SIZE)
          call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
          call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup = irank/isize_g
#endif

          call read_parameters()

          allocate (weight_k(N_wlk))
          N_grc = N_hfb*N_wlk
          allocate (overlap(N_grc))
          N_wlk_mpi = N_wlk*isize_g
          N_grc_mpi = N_grc*isize_g

          ! Setup the Bravais lattice
          call Ham_Latt

          ! Setup the hopping / single-particle part
          call Ham_Hop

          ! Setup the interaction.
          call Ham_V

          ! Setup the trival wave function, in case of a projector approach
          call Ham_Trial

#ifdef MPI
          if (Irank_g == 0) then
#endif
             File_info = "info"
#if defined(TEMPERING)
             write (File_info, '(A,I0,A)') "Temp_", igroup, "/info"
#endif
             open (newunit=unit_info, file=file_info, status="unknown", position="append")
             write (unit_info, *) '====================================='
             write (unit_info, *) 'Model is      : ', Model
             write (unit_info, *) 'Lattice is    : ', Lattice_type
             write (unit_info, *) '# unit cells  : ', Latt%N
             write (unit_info, *) '# of orbitals : ', Latt_unit%Norb
             write (unit_info, *) 'Flux_1        : ', Phi_X
             write (unit_info, *) 'Flux_2        : ', Phi_Y
             if (Bulk) then
                write (unit_info, *) 'Twist as phase factor in bulk'
             else
                write (unit_info, *) 'Twist as boundary condition'
             end if
             write (unit_info, *) 'Checkerboard  : ', Checkerboard
             write (unit_info, *) 'Symm. decomp  : ', Symm
             write (unit_info, *) 'dtau,Ltrot_eff: ', dtau, Ltrot
             write (unit_info, *) 'N_SUN         : ', N_SUN
             write (unit_info, *) 'N_FL          : ', N_FL
             write (unit_info, *) 'N_wlk         : ', N_wlk
             write (unit_info, *) 'N_wlk_mpi     : ', N_wlk_mpi
             write (unit_info, *) 'N_hfb         : ', N_hfb
             write (unit_info, *) 'N_grc         : ', N_grc
             write (unit_info, *) 'N_grc_mpi     : ', N_grc_mpi
             write (unit_info, *) 't             : ', ham_t
             write (unit_info, *) 'alpha         : ', ham_alpha
             write (unit_info, *) 'Ham_U         : ', Ham_U
             write (unit_info, *) 'Ham_chem      : ', Ham_chem
             write (unit_info, *) 'N_dope        : ', N_dope
             do nf = 1, N_FL
                write (unit_info, *) 'Flavor of:  ', nf
                write (unit_info, *) 'Degen of initial SD walker: ', wf_r(nf)%degen
                write (unit_info, *) 'number of particle for flavor:  ', size(wf_r(nf)%P,2)
             end do
             close (unit_info)
#ifdef MPI
          end if
#endif
       end subroutine Ham_Set

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the  Lattice
!--------------------------------------------------------------------
       subroutine Ham_Latt

          use Predefined_Lattices

          implicit none
          
          integer :: i, i1, no, a0, a1, a2, b0, b1, b2, k, k1, k2, ntype, j1
          integer :: ix, iy, rix, riy, rdx, rdy, ndx, ndy, nsi, nx, ny
          real    (kind=kind(0.d0)) :: pi = acos(-1.d0)
          
          ! Use predefined stuctures or set your own lattice.
          call Predefined_Latt(Lattice_type, L1, L2, Ndim, List, Invlist, Latt, Latt_Unit)
            
          allocate(rot_del_list(latt%n,latt_unit%norb,4))
          allocate(del_list    (latt%n,latt_unit%norb,4))
          allocate(rot_list    (latt%n,latt_unit%norb  ))

          del_list=0
          rot_list=0
          rot_del_list=0

          !! rotation mapping
          do i1 = 1, latt%n
             ix  = latt%list(i1,1); iy = latt%list(i1,2)
             rix = iy; riy = -ix
             rdx = rix - ix; rdy = riy - iy
             ndx = abs(rdx); ndy = abs(rdy)
             nsi = i1
             a1 = sign(1,rdx); a2 = sign(1,rdy)
             do nx = 1, ndx;
                nsi = latt%nnlist(nsi,a1,0)
             enddo
             do ny = 1, ndy;
                nsi = latt%nnlist(nsi,0,a2)
             enddo
             rot_list(i1,1) = nsi
          enddo
          
          !! rotation + delta mapping
          do ntype = 1, 4
             do j1 = 1, latt%n
                select case(ntype)
                case(1)
                    i1 = latt%nnlist(j1, 1, 0)
                case(2)
                    i1 = latt%nnlist(j1,-1, 0)
                case(3)
                    i1 = latt%nnlist(j1, 0, 1)
                case(4)
                    i1 = latt%nnlist(j1, 0,-1)
                end select
                ix  = latt%list(i1,1); iy = latt%list(i1,2)
                rix = iy; riy = -ix
                rdx = rix - ix; rdy = riy - iy
                ndx = abs(rdx); ndy = abs(rdy)
                nsi = i1
                a1 = sign(1,rdx); a2 = sign(1,rdy)
                do nx = 1, ndx;
                   nsi = latt%nnlist(nsi,a1,0)
                enddo
                do ny = 1, ndy;
                   nsi = latt%nnlist(nsi,0,a2)
                enddo
                rot_del_list(j1,1,ntype) = nsi
             enddo
          enddo
          
          !! delta mapping
          do ntype = 1, 4
             do j1 = 1, latt%n
                select case(ntype)
                case(1)
                    i1 = latt%nnlist(j1, 1, 0)
                case(2)
                    i1 = latt%nnlist(j1,-1, 0)
                case(3)
                    i1 = latt%nnlist(j1, 0, 1)
                case(4)
                    i1 = latt%nnlist(j1, 0,-1)
                end select
                del_list(j1,1,ntype) = i1
             enddo
          enddo
            
          !! d wave 1, px wave 2, py wave 3
          allocate(ff_s(4,1,3))
          do ntype = 1, 3
              select case(ntype)
              case (1)
                  ff_s(1,1,ntype) =  1.d0
                  ff_s(2,1,ntype) =  1.d0
                  ff_s(3,1,ntype) = -1.d0
                  ff_s(4,1,ntype) = -1.d0
              case (2)
                  ff_s(1,1,ntype) =  1.d0
                  ff_s(2,1,ntype) = -1.d0
                  ff_s(3,1,ntype) =  0.d0
                  ff_s(4,1,ntype) =  0.d0
              case (3)
                  ff_s(1,1,ntype) =  0.d0
                  ff_s(2,1,ntype) =  0.d0
                  ff_s(3,1,ntype) =  1.d0
                  ff_s(4,1,ntype) = -1.d0
              end select
          enddo

       end subroutine Ham_Latt
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the Hopping
!--------------------------------------------------------------------
       subroutine Ham_Hop

          implicit none

          real(Kind=kind(0.d0)) ::  Ham_Lambda = 0.d0

          real(Kind=kind(0.d0)), allocatable :: Ham_tx_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
               &                                Ham_ty_vec(:), Ham_Lambda_vec(:)
          integer, allocatable ::   N_Phi_vec(:)

          ! Use predefined stuctures or set your own hopping
          integer :: n, nth

          allocate (ham_tx_vec(N_FL), ham_ty_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
               &    N_Phi_vec(N_FL), Ham_Lambda_vec(N_FL))

          ! Here we consider no N_FL  dependence of the hopping parameters.
          Ham_tx_vec = ham_t
          Ham_ty_vec = ham_t
          Ham_Chem_vec = ham_chem
          Phi_X_vec = phi_X
          Phi_Y_vec = phi_Y
          N_Phi_vec = n_phi

          ham_tx_vec(1) = ham_t;
          ham_tx_vec(2) = ham_t*ham_alpha;
          ham_ty_vec(1) = ham_t*ham_alpha;
          ham_ty_vec(2) = ham_t;

          select case (Lattice_type)
          case ("square_anisotropic")
            call set_hopping_parameters_square_anisotropic(Hopping_Matrix, ham_tx_vec, ham_ty_vec, Ham_Chem_vec, &
                   & Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
          end select

          call Predefined_Hoppings_set_OPT(Hopping_Matrix, List, Invlist, Latt, Latt_unit, Dtau, Checkerboard, Symm, OP_T)

          deallocate (ham_tx_vec, ham_ty_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, N_Phi_vec)

       end subroutine Ham_Hop
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the trial wave function
!--------------------------------------------------------------------
       subroutine Ham_Trial()

          use Predefined_Trial

          implicit none

          integer :: N_part, nf
          ! Use predefined stuctures or set your own Trial  wave function
          N_part = ndim/2-n_dope
          call Predefined_TrialWaveFunction(Lattice_type, Ndim, List, Invlist, Latt, Latt_unit, &
               &                            N_part, ham_alpha, N_FL, N_hfb, WF_L, WF_R)

       end subroutine Ham_Trial

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the interaction
!--------------------------------------------------------------------
       subroutine Ham_V

          use Predefined_Int
          implicit none

          integer :: nf, I, I1, I2, nc, J, no, N_ops
          real(Kind=kind(0.d0)) :: X, zero = 1.d-10

          N_ops = 0
          if (abs(ham_u) > zero) N_ops = N_ops + Latt%N*Latt_unit%Norb

          allocate (op_v(n_ops, n_fl))
          nc = 0
          do i1 = 1, latt%N
             do no = 1, latt_unit%norb
                i = invlist(i1, no)
                if (abs(ham_u) > zero) then
                   nc = nc + 1
                   do nf = 1, n_fl
                       call op_make( op_v(nc,nf), 1 )
                       op_v(nc,nf)%p(1) = I
                       
                       op_v(nc,nf)%o(1,1) = cmplx(1.d0, 0.d0, kind(0.D0)) 
                       op_V(nc,nf)%g      = sqrt(cmplx(dtau*ham_u/2.d0, 0.D0, kind(0.D0)))
                       op_v(nc,nf)%alpha  = cmplx(-0.5d0, 0.d0, kind(0.D0))
                       
                       op_v(nc,nf)%type   = 3
                       call op_set(op_v(nc,nf))
                   enddo
                end if
             end do
          end do

          allocate(x_local(size(op_v,1),n_wlk))
          ! initial field shift
          x_local(:,:) = cmplx(0.d0,0.d0,kind(0.d0))

       end subroutine Ham_V

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specifiy the equal time and time displaced observables
!> @details
!--------------------------------------------------------------------
       subroutine alloc_obs(ltau)

          implicit none
          !>  Ltau=1 if time displaced correlations are considered.
          integer, intent(In) :: ltau
          integer    ::  i, N, Nt
          character(len=64) ::  Filename
          character(len=:), allocatable ::  Channel

          ! Scalar observables
          allocate (obs_scal(4))
          do I = 1, size(obs_scal, 1)
             select case (I)
             case (1)
                N = 1; Filename = "Kin"
             case (2)
                N = 1; Filename = "Pot"
             case (3)
                N = 1; Filename = "Part"
             case (4)
                N = 1; Filename = "Ener"
             case default
                write (6, *) ' Error in Alloc_obs '
             end select
             call obser_Vec_make(Obs_scal(I), N, Filename)
          end do

          ! Equal time correlators
          allocate (Obs_eq(4))
          do I = 1, size(Obs_eq, 1)
             select case (I)
             case (1)
                Filename = "Green"
             case (2)
                Filename = "SpinZ"
             case (3)
                Filename = "Den"
             case (4)
                Filename = "swave"
             case default
                write (6, *) ' Error in Alloc_obs '
             end select
             Nt = 1
             Channel = '--'
             call obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
          end do

       end subroutine alloc_obs
       
       !!========================================================================!!
       !!     compute local energy of a given walker
       !!========================================================================!!
       complex(Kind=kind(0.d0)) function e0_local(gr, kappa, kappa_bar)
          implicit none

          complex(Kind=kind(0.d0)), intent(in) :: gr       (ndim, ndim, n_fl)
          complex(Kind=kind(0.d0)), intent(in) :: kappa    (ndim, ndim, n_fl)
          complex(Kind=kind(0.d0)), intent(in) :: kappa_bar(ndim, ndim, n_fl)

          !Local
          complex(Kind=kind(0.d0)) :: grc(ndim, ndim, n_fl), ZK
          complex(Kind=kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP, ZS, ZZ, ZXY
          integer :: I, J, imj, nf, dec, I1, J1, no_I, no_J, n
          real(Kind=kind(0.d0)) :: X

          !! grc_{i,j} = <c^{\dagger}_i c_j>
          do nf = 1, N_FL
             do I = 1, Ndim
                do J = 1, Ndim
                   grc(I, J, nf) = -gr(J, I, nf)
                end do
                grc(I, I, nf) = 1.d0 + grc(I, I, nf)
             end do
          end do

          zkin = cmplx(0.d0, 0.d0, kind(0.d0))
          call Predefined_Hoppings_Compute_Kin(Hopping_Matrix, List, Invlist, Latt, Latt_unit, gr, zkin)
          zkin = zkin*dble(n_sun)
          
          zpot = cmplx(0.d0, 0.d0, kind(0.d0))
          do I = 1, Latt%N
             do no_I = 1, Latt_unit%Norb
                I1 = Invlist(I, no_I)
                !zpot = zpot + 2.d0*gr(i1, i1, 1)*gr(i1, i1, 2) - &
                !    & 2.d0*kappa_bar(i1,i1,1)*kappa(i1,i1,1) & 
                !    & - gr(i1, i1, 1) - gr(i1, i1, 2) + 1.d0
                zpot = zpot + 2.d0*gr(i1, i1, 1)*gr(i1, i1, 2) - 2.d0*kappa_bar(i1,i1,1)*kappa(i1,i1,1)
             end do
          end do
          zpot = zpot*(-ham_u/2.d0)

          e0_local = zpot + zkin

       end function e0_local
       
       !!========================================================================!!
       !!     compute local shift of auxillary field 
       !!========================================================================!!
       complex(Kind=kind(0.d0)) function xbar_loc_compute(gr,nc)
          implicit none

          complex(Kind=kind(0.d0)), intent(in) :: gr(ndim, ndim, n_fl)
          integer, intent(in) :: nc
         
          complex(Kind=kind(0.d0)) :: ztmp
          integer :: nf, i, j, i1, j1, i_grc, n, I0, J0, n_op
         
          ztmp = cmplx(0.d0,0.d0,kind(0.d0))
          do nf = 1,N_FL
             do I = 1,op_v(nc,nf)%N
             do J = 1,op_v(nc,nf)%N
                I1 = op_v(nc,nf)%P(I)
                J1 = op_v(nc,nf)%P(J)

                ztmp = ztmp - op_v(nc,nf)%g*gr(I1,J1,nf)*op_v(nc,nf)%o(I,J)
             enddo
             enddo
             ztmp = ztmp - op_v(nc,nf)%g*op_v(nc,nf)%alpha
          enddo
             
          xbar_loc_compute = ztmp
       
       end function xbar_loc_compute
       
       !!========================================================================!!
       !!     set up the shift of auxillary field for importance sampling
       !!========================================================================!!
       subroutine set_xloc(gr,kappa,kappa_bar)
         Implicit none
          
         complex(Kind=kind(0.d0)), intent(in) :: gr       (ndim, ndim, n_fl, n_grc)
         complex(Kind=kind(0.d0)), intent(in) :: kappa    (ndim, ndim, n_fl, n_grc)
         complex(Kind=kind(0.d0)), intent(in) :: kappa_bar(ndim, ndim, n_fl, n_grc)
         
         ! local
         complex (kind=kind(0.d0)) :: ZK, Zn, ztmp, x_tmp(n_hfb)
         complex (kind=kind(0.d0)) :: ztmp1
         integer :: i_wlk, nf, i, j, i1, j1, i_grc, n, I0, J0, n_op, nc, ns
         real(Kind=kind(0.d0)) :: sign_w, pi = acos(-1.d0), zero = 1.0e-12, re_o_max

         n_op = size(op_v,1)
      
         do i_wlk = 1, n_wlk

            sign_w = cos(aimag(weight_k(i_wlk)))
            if ( sign_w .gt. zero ) then

            do nc = 1, n_op
               
               !! the max real part of the overlap
               re_o_max = 0.d0
               do ns = 1, n_hfb
                  i_grc = ns + (i_wlk - 1)*n_hfb
                  if ( dble(overlap(i_grc)) .gt. re_o_max ) re_o_max = dble(overlap(i_grc))
               enddo

               ztmp = 0.d0
               do ns = 1, N_hfb
                  i_grc = ns + (i_wlk - 1)*N_hfb
                  x_tmp(n_hfb) = & 
                      & xbar_loc_compute(gr(:,:,:,i_grc),nc)*exp(overlap(i_grc)-re_o_max)
                  ztmp = ztmp + exp(overlap(i_grc)-re_o_max)
               end do
               
               x_local(nc,i_wlk) = sum(x_tmp(:))/ztmp
            enddo

            endif
         
         enddo
         
       end subroutine set_xloc

       !!===================================================================!!
       !!     compute the sum of the weight and rescale weight
       !!===================================================================!!
       subroutine sum_weight(z_sum_weight)
          implicit none

          complex(Kind=kind(0.d0)), intent(out) :: z_sum_weight

          !local
          integer :: i_wlk, ii, i_st, i_ed
          real   (Kind=kind(0.d0)) :: X1, tot_re_weight, max_re_w, ang_w, re_lw
          complex(Kind=kind(0.d0)) :: Z1, Z2, wtmp
          
          complex(Kind=kind(0.d0)), allocatable :: weight_mpi(:), w_arr(:)

#ifdef MPI
          integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          integer        :: STATUS(MPI_STATUS_SIZE)
          call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
          call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup = irank/isize_g
#endif
          
          if ( irank_g .eq. 0 ) then
              allocate(weight_mpi(N_wlk_mpi))
              allocate(w_arr     (N_wlk    ))
              i_st = 1
              i_ed = N_wlk
              weight_mpi(i_st:i_ed) = weight_k(:)
          endif
          
          if ( irank_g .eq. 0 ) then
              do ii = 1, isize_g - 1
                 i_st = ii*N_wlk + 1
                 i_ed = (ii + 1)*N_wlk
                 call mpi_recv(w_arr, n_wlk, MPI_COMPLEX16, ii, 1, group_comm, status, ierr)
                 weight_mpi(i_st:i_ed) = w_arr
              enddo
          else
              call mpi_send(weight_k, n_wlk, MPI_COMPLEX16, 0, 1, group_comm, ierr)
          endif
          
          if ( irank_g .eq. 0 ) then
              max_re_w = maxval(dble(weight_mpi))
              deallocate(weight_mpi, w_arr)
          endif
          call MPI_BCAST(max_re_w, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

          Z1 = cmplx(0.d0, 0.d0, kind(0.d0))
          Z2 = cmplx(0.d0, 0.d0, kind(0.d0))
          do i_wlk = 1, N_wlk
             ang_w = aimag(weight_k(i_wlk)) 
             !! Real part of the weight
             if ( cos(ang_w) .gt. 0.d0 ) then 
                 weight_k(i_wlk) = weight_k(i_wlk) - max_re_w
                 re_lw = dble (weight_k(i_wlk))
                 z1 = z1 + exp(re_lw)*cos(ang_w)
             endif
          end do
          call MPI_REDUCE(Z1, Z2, 1, MPI_COMPLEX16, MPI_SUM, 0, Group_comm, IERR)
          call MPI_BCAST (Z2, 1, MPI_COMPLEX16, 0, MPI_COMM_WORLD, ierr)

          z_sum_weight = Z2

       end subroutine sum_weight

       !!===================================================================!!
       !!     update the normalization factor
       !!===================================================================!!
       subroutine update_fac_norm(gr, kappa, kappa_bar, ntw)
          implicit none

          complex(Kind=kind(0.d0)), intent(in) :: gr       (ndim, ndim, n_fl, n_grc)
          complex(Kind=kind(0.d0)), intent(in) :: kappa    (ndim, ndim, n_fl, n_grc)
          complex(Kind=kind(0.d0)), intent(in) :: kappa_bar(ndim, ndim, n_fl, n_grc)
          integer, intent(in) :: ntw

          !local
          integer :: i_wlk, ii, i_st, i_ed, ns, i_grc
          real   (Kind=kind(0.d0)) :: X1, re_w_tmp, max_re_w, ang_w, re_lw, re_o_max
          complex(Kind=kind(0.d0)) :: Z1, Z2, wtmp, el_tmp, Z, tot_ene, zr1, zr2
          character(LEN=64)  :: filename
          complex(Kind=kind(0.d0)), allocatable :: weight_mpi(:), w_arr(:)

#ifdef MPI
          integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
          integer        :: STATUS(MPI_STATUS_SIZE)
          call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
          call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup = irank/isize_g
#endif
          filename = "E_local.dat"

          if ( irank_g .eq. 0 ) then
              allocate(weight_mpi(N_wlk_mpi))
              allocate(w_arr     (N_wlk    ))
              i_st = 1
              i_ed = N_wlk
              weight_mpi(i_st:i_ed) = weight_k(:)
          endif

          if ( irank_g .eq. 0 ) then
              do ii = 1, isize_g - 1
                 i_st = ii*N_wlk + 1
                 i_ed = (ii + 1)*N_wlk
                 call mpi_recv(w_arr, n_wlk, MPI_COMPLEX16, ii, 1, group_comm, status, ierr)
                 weight_mpi(i_st:i_ed) = w_arr
              enddo
          else
              call mpi_send(weight_k, n_wlk, MPI_COMPLEX16, 0, 1, group_comm, ierr)
          endif
          
          if ( irank_g .eq. 0 ) then
              max_re_w = maxval(dble(weight_mpi))
              deallocate(weight_mpi, w_arr)
          endif
          call MPI_BCAST(max_re_w, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      
          Z1 = cmplx(0.d0, 0.d0, kind(0.d0))
          Z2 = cmplx(0.d0, 0.d0, kind(0.d0))
          do i_wlk = 1, N_wlk
             ang_w = aimag(weight_k(i_wlk)) 
             !! Real part of the weight
             if ( cos(ang_w) .gt. 0.d0 ) then
                 weight_k(i_wlk) = weight_k(i_wlk) - max_re_w
                 re_lw = dble (weight_k(i_wlk))
                 
                 !! the max real part of the overlap
                 re_o_max = 0.d0
                 do ns = 1, n_hfb
                    i_grc = ns + (i_wlk - 1)*n_hfb
                    if ( dble(overlap(i_grc)) .gt. re_o_max ) re_o_max = dble(overlap(i_grc))
                 enddo

                 z = 0.d0
                 do ns = 1, N_hfb
                    i_grc = ns + (i_wlk - 1)*N_hfb
                    z = z + exp(overlap(i_grc)-re_o_max)
                 end do

                 !! real part of mix estimated energy
                 tot_ene = cmplx(0.d0, 0.d0, kind(0.d0))
                 do ns = 1, N_hfb
                    i_grc = ns + (i_wlk - 1)*N_hfb
                    el_tmp = dble(ham%E0_local(gr(:,:,:,i_grc), & 
                        & kappa(:,:,:,i_grc), kappa_bar(:,:,:,i_grc)))
                    tot_ene = tot_ene + el_tmp*exp(overlap(i_grc)-re_o_max)/Z
                 end do
                 re_w_tmp = exp(re_lw)*cos(ang_w)
                 z1 = z1 + re_w_tmp
                 z2 = z2 + re_w_tmp*tot_ene
             endif
          end do
          call mpi_reduce(z1, zr1, 1, MPI_COMPLEX16, MPI_SUM, 0, Group_comm, ierr)
          call mpi_reduce(z2, zr2, 1, MPI_COMPLEX16, MPI_SUM, 0, Group_comm, ierr)
          
          if ( irank_g == 0 ) then
             zr2 = zr2/zr1
             fac_norm = dble(zr2)*dtau
             open (UNIT=77, FILE=filename, STATUS='UNKNOWN', position="append")
             write (77, *) ntw*dtau, dble(zr2)
             close (77)
          end if
          call MPI_BCAST(fac_norm, 1, MPI_COMPLEX16, 0, MPI_COMM_WORLD, ierr)

       end subroutine update_fac_norm

    end submodule ham_bose_metal_smod
