    submodule (Hamiltonian_main) ham_Hubbard_smod

      Use Operator_mod
      Use WaveFunction_mod
      Use Lattices_v3
      Use MyMats
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables
      Use Fields_mod
      Use Predefined_Hoppings

      Implicit none
      
      type, extends(ham_base) :: ham_Hubbard
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
        procedure, nopass :: S0
        procedure, nopass :: E0_local
        procedure, nopass :: sum_weight
        procedure, nopass :: update_fac_norm 
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Hubbard

      !#PARAMETERS START# VAR_lattice
      Character (len=64) :: Model = 'Hubbard'  ! Value not relevant
      Character (len=64) :: Lattice_type = 'Square'
      Integer            :: L1 = 6   ! Length in direction a_1
      Integer            :: L2 = 6   ! Length in direction a_2
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_Model_Generic
      !Integer              :: N_SUN        = 1        ! Number of colors
      !Integer              :: N_FL         = 2        ! Number of flavors
      !Integer              :: N_slat       = 1        ! Number of slater on trial wave function
      !Integer              :: N_wlk        = 1        ! Number of walker
      !Integer              :: ltrot        = 10       ! length of imaginary time for dynamical measure
      real(Kind=Kind(0.d0)) :: Phi_X        = 0.d0     ! Twist along the L_1 direction, in units of the flux quanta
      real(Kind=Kind(0.d0)) :: Phi_Y        = 0.d0     ! Twist along the L_2 direction, in units of the flux quanta
      logical               :: Bulk         = .true.   ! Twist as a vector potential (.T.), or at the boundary (.F.)
      Integer               :: N_Phi        = 0        ! Total number of flux quanta traversing the lattice
      real(Kind=Kind(0.d0)) :: Dtau         = 0.1d0    ! Thereby Ltrot=Beta/dtau
      logical               :: Checkerboard = .true.   ! Whether checkerboard decomposition is used
      !logical              :: Symm         = .true.   ! Whether symmetrization takes place
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_Hubbard
      real(Kind=Kind(0.d0)) :: ham_T      = 1.d0     ! Hopping parameter
      real(Kind=Kind(0.d0)) :: Ham_chem   = 0.d0     ! Chemical potential
      real(Kind=Kind(0.d0)) :: Ham_U      = 4.d0     ! Hubbard interaction
      real(Kind=Kind(0.d0)) :: ham_T2     = 1.d0     ! For bilayer systems
      real(Kind=Kind(0.d0)) :: Ham_U2     = 4.d0     ! For bilayer systems
      real(Kind=Kind(0.d0)) :: ham_Tperp  = 1.d0     ! For bilayer systems
      Integer               :: N_dope     = 0        ! Number of doping electrons
      !#PARAMETERS END#

      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:)
      Integer, allocatable :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell

    contains
      
      module Subroutine Ham_Alloc_hubbard
        allocate(ham_Hubbard::ham)
      end Subroutine Ham_Alloc_hubbard

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Hubbard_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
      Subroutine Ham_Set

#if defined (MPI)
          Use mpi
#endif
          Implicit none

          integer                :: ierr, nf, unit_info
          Character (len=64)     :: file_info


          ! L1, L2, Lattice_type, List(:,:), Invlist(:,:) -->  Lattice information
          ! Ham_T, Chem, Phi_X, XB_B, Checkerboard, Symm   -->  Hopping
          ! Interaction                              -->  Model
          ! Simulation type                          -->  Finite  T or Projection  Symmetrize Trotter.
                  
#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          ! From dynamically generated file "Hamiltonian_Hubbard_read_write_parameters.F90"
          call read_parameters()

          allocate( weight_k(N_wlk) )
          N_grc = N_slat*N_wlk
          allocate( overlap(N_grc) )
          N_wlk_mpi=N_wlk*isize_g
          N_grc_mpi=N_grc*isize_g

          ! Setup the Bravais lattice
          Call Ham_Latt

          ! Setup the hopping / single-particle part
          Call Ham_Hop

          ! Setup the interaction.
          call Ham_V

          ! Setup the trival wave function, in case of a projector approach
          Call Ham_Trial()

#ifdef MPI
          If (Irank_g == 0) then
#endif
             File_info = "info"
             Open(newunit=unit_info, file=file_info, status="unknown", position="append")
             Write(unit_info,*) '====================================='
             Write(unit_info,*) 'Model is      : ', Model
             Write(unit_info,*) 'Lattice is    : ', Lattice_type
             Write(unit_info,*) '# unit cells  : ', Latt%N 
             Write(unit_info,*) '# of orbitals : ', Latt_unit%Norb
             Write(unit_info,*) 'Flux_1        : ', Phi_X
             Write(unit_info,*) 'Flux_2        : ', Phi_Y
             If (Bulk) then
                Write(unit_info,*) 'Twist as phase factor in bulk'
             Else
                Write(unit_info,*) 'Twist as boundary condition'
             endif
             Write(unit_info,*) 'HS  couples to z-component of spin'
             Write(unit_info,*) 'Checkerboard  : ', Checkerboard
             Write(unit_info,*) 'Symm. decomp  : ', Symm
             Write(unit_info,*) 'dtau: ', dtau
             Write(unit_info,*) 'ltrot         : ', ltrot
             Write(unit_info,*) 'N_SUN         : ', N_SUN
             Write(unit_info,*) 'N_FL          : ', N_FL
             Write(unit_info,*) 'N_wlk         : ', N_wlk
             Write(unit_info,*) 'N_wlk_mpi     : ', N_wlk_mpi
             Write(unit_info,*) 'N_slat        : ', N_slat
             Write(unit_info,*) 'N_grc         : ', N_grc
             Write(unit_info,*) 'N_grc_mpi     : ', N_grc_mpi
             Write(unit_info,*) 'N_dope        : ', N_dope
             Write(unit_info,*) 't             : ', Ham_T
             Write(unit_info,*) 'Ham_U         : ', Ham_U
             Write(unit_info,*) 't2            : ', Ham_T2
             Write(unit_info,*) 'Ham_U2        : ', Ham_U2
             Write(unit_info,*) 'Ham_tperp     : ', Ham_tperp
             Write(unit_info,*) 'Ham_chem      : ', Ham_chem
             Close(unit_info)
#ifdef MPI
          Endif
#endif
        end Subroutine Ham_Set

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the  Lattice
!--------------------------------------------------------------------
        Subroutine Ham_Latt

          Use Predefined_Lattices

          Implicit none
          ! Use predefined stuctures or set your own lattice.
          Call Predefined_Latt(Lattice_type, L1,L2,Ndim, List,Invlist,Latt,Latt_Unit)

        end Subroutine Ham_Latt
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the Hopping
!--------------------------------------------------------------------
        Subroutine Ham_Hop

          Implicit none

          Real (Kind=Kind(0.d0) ) ::  Ham_Lambda = 0.d0

          Real (Kind=Kind(0.d0) ), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
               &                                  Ham_T2_vec(:),  Ham_Lambda_vec(:)
          Integer, allocatable ::   N_Phi_vec(:)

          ! Use predefined stuctures or set your own hopping
          Integer :: n,nth

          Allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
               &                                   N_Phi_vec(N_FL), Ham_Lambda_vec(N_FL) )

          ! Here we consider no N_FL  dependence of the hopping parameters.
          Ham_T_vec      = Ham_T
          Ham_Tperp_vec  = Ham_Tperp
          Ham_Chem_vec   = Ham_Chem
          Phi_X_vec      = Phi_X
          Phi_Y_vec      = Phi_Y
          Ham_T2_vec     = Ham_T2
          Ham_Lambda_vec = Ham_Lambda
          N_Phi_vec      = N_Phi

          Select case (Lattice_type)
          Case ("Square")
             Call  Set_Default_hopping_parameters_square(Hopping_Matrix,Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                  &                                      Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
          Case ("Honeycomb")
             Ham_Lambda = 0.d0
             Call  Set_Default_hopping_parameters_honeycomb(Hopping_Matrix, Ham_T_vec, Ham_Lambda_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                  &                                         Bulk,  N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
          end Select

          Call  Predefined_Hoppings_set_OPT(Hopping_Matrix,List,Invlist,Latt,  Latt_unit,  Dtau, Checkerboard, Symm, OP_T )

          Deallocate (Ham_T_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
               &                                   N_Phi_vec,  Ham_Lambda_vec )

        end Subroutine Ham_Hop
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the trial wave function
!--------------------------------------------------------------------
        Subroutine Ham_Trial()

          Use Predefined_Trial

          Implicit none

          Integer :: N_part, nf
          ! Use predefined stuctures or set your own Trial  wave function
          N_part = Ndim/2-N_dope

          Call Predefined_TrialWaveFunction(Lattice_type ,Ndim,  List,Invlist,Latt, Latt_unit, &
               &                            N_part, N_FL, N_slat, wf_l, wf_r)
          overlap(:) = cmplx(1.d0, 0.d0, kind(0.d0))

        end Subroutine Ham_Trial

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the interaction
!--------------------------------------------------------------------
        Subroutine Ham_V

          Use Predefined_Int
          Implicit none

          Integer :: nf, I, I1, I2,  nc,  J, no,  N_ops
          Real (Kind=Kind(0.d0)) :: X,  Zero = 1.D-10
          Real (Kind=Kind(0.d0)), allocatable :: Ham_U_vec(:)

          Allocate (Ham_U_vec(Latt_unit%Norb))

          N_ops = 0
          Ham_U_vec(:)  = Ham_U
          If (abs(Ham_U ) > Zero ) N_ops = N_ops + Latt%N*Latt_unit%Norb
          
          Allocate(Op_V(N_ops,N_FL))
          Ham_U_vec = Ham_U_vec/real(N_SUN,kind(0.d0))
          nc = 0
          Do I1 = 1,Latt%N
             do no = 1, Latt_unit%Norb
                I = invlist(I1,no)
                if (abs(Ham_U_vec(no)) > Zero ) then
                   nc = nc + 1
                   Call Predefined_Int_U_MZ( OP_V(nc,1), OP_V(nc,2), I,  DTAU, Ham_U_vec(no) )
                endif
             enddo
          enddo

          Deallocate (Ham_U_vec)

        end Subroutine Ham_V


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specifiy the equal time and time displaced observables
!> @details
!--------------------------------------------------------------------
        Subroutine  Alloc_obs(Ltau)

          Implicit none
          !>  Ltau=1 if time displaced correlations are considered.
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Nt
          Character (len=64) ::  Filename
          Character (len=2)  ::  Channel

          ! Scalar observables
          Allocate ( Obs_scal(6) )
          Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
                N = 1;   Filename = "Kin"
             case (2)
                N = 1;   Filename = "Pot"
             case (3)
                N = 1;   Filename = "Part"
             case (4)
                N = 1;   Filename = "Ener"
             case (5)
                N = ndim*ndim*n_fl;   Filename ="grc"
             case (6)
                N = ndim*ndim*n_fl;   Filename ="mixgrc"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo
             
          Allocate ( Obs_eq(5) )
          Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Filename = "Green"
             case (2)
                Filename = "SpinZ"
             case (3)
                Filename = "SpinXY"
             case (4)
                Filename = "SpinT"
             case (5)
                Filename = "Den"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Nt = 1
             Channel = '--'
             Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
          enddo
          
          If (Ltau == 1) then
             ! Time-displaced correlators
             Allocate ( Obs_tau(5) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Channel = 'P' ; Filename = "Green"
                case (2)
                   Channel = 'PH'; Filename = "SpinZ"
                case (3)
                   Channel = 'PH'; Filename = "SpinXY"
                case (4)
                   Channel = 'PH'; Filename = "SpinT"
                case (5)
                   Channel = 'PH'; Filename = "Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = Ltrot+1
                Channel = 'T0'
                Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             enddo
          endif

        End Subroutine Alloc_obs

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes equal time observables
!> @details
!> @param [IN] Gr   Complex(:,:,:)
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!> @param [IN] Ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!-------------------------------------------------------------------
        subroutine Obser(GR,GR_mix,i_wlk,i_grc,sum_w,sum_o,act_mea)

          Use Predefined_Obs

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR_mix(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: sum_w, sum_o
          Integer, Intent(IN) :: i_wlk, i_grc, act_mea

          !Local
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK, invsumw, zone, ztmp, z_ol
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, ZZ, ZXY, ZW, Re_ZW, Z_fac
          Integer :: I,J, imj, nf, dec, I1, J1, no_I, no_J,n, nc
          Real    (Kind=Kind(0.d0)) :: X

          Z_ol  = exp(overlap(i_grc))/sum_o
          ZW    = cmplx(weight_k(i_wlk),0.d0,kind(0.d0))/sum_w
          Z_fac = Z_ol*ZW
          
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                Enddo
                GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >

          ! Compute scalar observables.
          if ( act_mea .eq. 0 ) then
              do I = 1,Size(Obs_scal,1)
                 Obs_scal(I)%N         =  Obs_scal(I)%N + 1
                 Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + 1.d0
              enddo
          endif

          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Call Predefined_Hoppings_Compute_Kin(Hopping_Matrix,List,Invlist, Latt, Latt_unit, GRC, ZKin)
          Zkin = Zkin* dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *Z_fac


          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Latt%N
             do no_I = 1,Latt_unit%Norb
                I1 = Invlist(I,no_I)
                ZPot = ZPot + Grc(i1,i1,1) * Grc(i1,i1,2)* ham_U
             enddo
          Enddo
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot *Z_fac


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf)
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)

          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * Z_fac

          Obs_scal(4)%Obs_vec(1)  =    Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*Z_fac

          nc = 0
          Do nf = 1,N_FL
             Do I = 1,Ndim
             Do J = 1,Ndim
                nc = nc + 1
                ztmp = grc(i,j,nf)
                Obs_scal(5)%Obs_vec(nc) = Obs_scal(5)%Obs_vec(nc) + ztmp*Z_fac
             Enddo
             Enddo
          Enddo

          nc = 0
          Do nf = 1,N_FL
             Do I = 1,Ndim
             Do J = 1,Ndim
                nc = nc + 1
                zone = cmplx(0.d0,0.d0,kind(0.d0))
                if ( I .eq. J ) zone = cmplx(1.d0,0.d0,kind(0.d0))
                Ztmp = zone-GR_mix(J,I,nf)
                Obs_scal(6)%Obs_vec(nc) = Obs_scal(6)%Obs_vec(nc) + ztmp*Z_fac
             Enddo
             Enddo
          Enddo

          ! Compute the standard two-point correlations
          if ( act_mea .eq. 0 ) then
              Do I = 1, Size(Obs_eq,1)
                 obs_eq(I)%N        = obs_eq(I)%N + 1
                 obs_eq(I)%Ave_sign = obs_eq(I)%Ave_sign + 1.d0
              enddo
          endif
          
          ! Standard two-point correlations
          Call Predefined_Obs_eq_Green_measure ( Latt, Latt_unit, List, GR, GRC, N_SUN, act_mea, Z_fac, Obs_eq(1) )
          Call Predefined_Obs_eq_SpinMz_measure( Latt, Latt_unit, List, GR, GRC, N_SUN, act_mea, Z_fac, Obs_eq(2),Obs_eq(3),Obs_eq(4) )
          Call Predefined_Obs_eq_Den_measure   ( Latt, Latt_unit, List, GR, GRC, N_SUN, act_mea, Z_fac, Obs_eq(5) )

        end Subroutine Obser

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes time displaced  observables
!> @details
!> @param [IN] NT, Integer
!> \verbatim
!>  Imaginary time
!> \endverbatim
!> @param [IN] GT0, GTT, G00, GTT,  Complex(:,:,:)
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
!>  G00(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(0  )>
!>  GTT(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(tau)>
!> \endverbatim
!-------------------------------------------------------------------
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, i_wlk,i_grc,sum_w,sum_o,act_mea)

          Use Predefined_Obs

          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: sum_w, sum_o
          Integer, Intent(IN) :: i_wlk, act_mea
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY, ZW, Re_ZW, Z_fac, Z_ol
          Real    (Kind=Kind(0.d0)) :: X
          Integer :: IMJ, I, J, I1, J1, no_I, no_J

          Z_ol  = exp(overlap(i_grc))/sum_o
          ZW    = cmplx(weight_k(i_wlk),0.d0,kind(0.d0))/sum_w
          Z_fac = Z_ol*ZW

          ! Compute the dynamic two-point correlations
          if ( (act_mea .eq. 0) .and. (nt .eq. 0) ) then
              Do I = 1, Size(Obs_tau,1)
                 obs_tau(I)%N        = obs_tau(I)%N + 1
                 obs_tau(I)%Ave_sign = obs_tau(I)%Ave_sign + 1.d0
              enddo
          endif

          ! Standard two-point correlations

          Call Predefined_Obs_tau_Green_measure  ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT, N_SUN, act_mea, Z_fac, Obs_tau(1) )
          Call Predefined_Obs_tau_SpinMz_measure ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT, N_SUN, act_mea, Z_fac, Obs_tau(2),&
               &                                   Obs_tau(3), Obs_tau(4) )
          Call Predefined_Obs_tau_Den_measure    ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT, N_SUN, act_mea, Z_fac, Obs_tau(5) )

        end Subroutine OBSERT

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Single spin flip S0 ratio
!> @details
!> S0=exp(-S0(new))/exp(-S0(old)) where the new configuration correpsonds to the old one up to
!> a spin flip of Operator n on time slice nt
!> @details
!--------------------------------------------------------------------
      Real (Kind=Kind(0.d0)) function S0(n,nt,Hs_new)
        Implicit none
        !> Operator index
        Integer, Intent(IN) :: n
        !> Time slice
        Integer, Intent(IN) :: nt
        !> New local field on time slice nt and operator index n
        Real (Kind=Kind(0.d0)), Intent(In) :: Hs_new

        Integer :: nt1,I

        S0 = 1.d0

      end function S0

      Complex (Kind=Kind(0.d0)) function E0_local(GR)
        Implicit none

        Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          
        !Local
        Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
        Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, ZZ, ZXY
        Integer :: I,J, imj, nf, dec, I1, J1, no_I, no_J,n
        Real    (Kind=Kind(0.d0)) :: X
          
        Do nf = 1,N_FL
           Do I = 1,Ndim
              Do J = 1,Ndim
                 GRC(I, J, nf) = -GR(J, I, nf)
              Enddo
              GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
           Enddo
        Enddo
          
        Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
        Call Predefined_Hoppings_Compute_Kin(Hopping_Matrix,List,Invlist, Latt, Latt_unit, GRC, ZKin)
        Zkin = Zkin* dble(N_SUN)
          
        ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
        Do I = 1,Latt%N
           do no_I = 1,Latt_unit%Norb
              I1 = Invlist(I,no_I)
              ZPot = ZPot + Grc(i1,i1,1) * Grc(i1,i1,2)* ham_U
           enddo
        Enddo

        E0_local = (Zpot + ZKin)*dtau

      end function E0_local

      Complex (Kind=Kind(0.d0)) function sum_weight
        Implicit none
        
        !local
        Integer                   :: i_wlk
        Real    (Kind=Kind(0.d0)) :: X1, tot_re_weight
        Complex (Kind=Kind(0.d0)) :: Z1, Z2, ZP, wtmp

#ifdef MPI
        Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
        Integer        :: STATUS(MPI_STATUS_SIZE)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        call MPI_Comm_rank(Group_Comm, irank_g, ierr)
        call MPI_Comm_size(Group_Comm, isize_g, ierr)
        igroup           = irank/isize_g
#endif
        
        Z1 = cmplx(0.d0, 0.d0, kind(0.d0))
        Z2 = cmplx(0.d0, 0.d0, kind(0.d0))
        do i_wlk  = 1, N_wlk
            ZP      = 1.d0
            wtmp    = cmplx(weight_k(i_wlk), 0.d0, kind(0.d0))
            Z1 = Z1 + wtmp*ZP
        enddo
        CALL MPI_REDUCE(Z1,Z2,1,MPI_COMPLEX16,MPI_SUM,0,Group_comm    ,IERR)
        CALL MPI_BCAST (Z2,   1,MPI_COMPLEX16        ,0,MPI_COMM_WORLD,ierr)
        
        sum_weight = Z2
        
      end function sum_weight
      
      Subroutine update_fac_norm(GR, ntw)
        Implicit none
         
        Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL,N_wlk)
        Integer                  , INTENT(IN) :: ntw
        
        !local
        Integer                   :: i_wlk, ns, i_grc
        Real    (Kind=Kind(0.d0)) :: X1, tot_re_weight
        Complex (Kind=Kind(0.d0)) :: tot_ene, Z1, Z2, tot_c_weight, ZP, wtmp, el_tmp, Z
        CHARACTER (LEN=64)  :: filename 

#ifdef MPI
        Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
        Integer        :: STATUS(MPI_STATUS_SIZE)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        call MPI_Comm_rank(Group_Comm, irank_g, ierr)
        call MPI_Comm_size(Group_Comm, isize_g, ierr)
        igroup           = irank/isize_g
#endif
        filename = "E_local.dat"

        !! initial energy
        tot_ene      = 0.d0
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
            OPEN (UNIT = 77, FILE=filename, STATUS='UNKNOWN', position="append")
            write(77,*) ntw*dtau,  dble(Z1)/dtau
            close(77)
        endif
        CALL MPI_BCAST(fac_norm, 1, MPI_REAL8, 0,MPI_COMM_WORLD,ierr)
        
      end subroutine update_fac_norm
        
    end submodule ham_Hubbard_smod
