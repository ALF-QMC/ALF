!  Copyright (C) 2022 The ALF project
!
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version


!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This File is a template for defining new models. 
!> One can define a new model class by copying this file, replacing alle occurences
!> of Nematic_Dirac_demo by the Hamiltonian name, populating the subroutines below as needed
!> adding the Hamiltonian name to the file Prog/Hamiltonians.list.

!> @details
!> The public variables of this module are the following
!>
!>
!> @param [public] OP_V
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable
!> List of operators of type=1,2 and 3 describing the sequence of interactions on a time slice.
!> The first index runs over this sequence. The second corresponds to the flavor index.  \endverbatim
!>
!> @param [public] OP_T
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable
!> Sequence of  operators  accounting for the  hopping on a  time slice. This can include  various
!> checkerboard decompositions. The first index runs over this sequence. The second corresponds to
!> the flavor index. \endverbatim
!> *  The progagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n}  \f$.  That is
!> first the hopping and then the potential energy.
!>
!>@param [public] WF_L
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Left trial wave function.  \endverbatim
!>
!> @param [public] WF_R
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Right trial wave function.   For both wave functions the index runs over the flavor index. \endverbatim
!>
!> @param [public]  nsigma
!> \verbatim Type(Fields)
!> Contains all auxiliary fields in the variable f(:,:). The first index runs through the operator
!> sequence. The second through the time slices.   \endverbatim
!
!> @param [public]  Ndim
!> \verbatim Integer
!> Total number of orbitals. e.g. # unit cells * # orbitals per unit cell.  \endverbatim
!
!> @param [public]  N_FL
!> \verbatim Integer
!> # of flavors.  Propagation is block diagonal in flavors.  \endverbatim
!
!> @param [public]  N_SUN
!> \verbatim Integer
!> # of colors.  Propagation is color independent.  \endverbatim
!>
!> @param [public] Ltrot
!> \verbatim Integer
!> Available measurment interval in units of Delta Tau. \endverbatim
!>
!> @param [public] Thtrot
!>  \verbatim Integer
!> Effective projection parameter in units of Delta Tau.  (Only relevant if projective option is turned on) \endverbatim
!>
!> @param [public] Projector
!> \verbatim Logical
!> Flag for projector. If true then the total number of time slices will correspond to Ltrot + 2*Thtrot \endverbatim
!>
!> @param [public] Group_Comm
!> \verbatim Integer
!> Defines MPI communicator  \endverbatim
!
!> @param [public] Symm
!> \verbatim Logical  \endverbatim
!> If set to true then the green functions will be symmetrized
!> before being  sent to the Obser, ObserT subroutines.
!> In particular, the transformation,  \f$ \tilde{G} =  e^{-\Delta \tau T /2 } G e^{\Delta \tau T /2 } \f$
!> will be carried out  and \f$ \tilde{G} \f$  will be sent to the Obser and ObserT subroutines.  Note that
!> if you want to use this  feature, then you have to be sure the hopping and interaction terms are decomposed
!> symmetrically. If Symm is true, the propagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=N_T}^{1}e^{T_n/2} \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n/2}  \f$
!>
!>
!> You still have to add some docu for the other private variables in this module.
!>
!--------------------------------------------------------------------

    submodule (Hamiltonian_main) ham_Nematic_Dirac_demo_smod

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
      Use LRC_Mod

      Implicit none
      
      type, extends(ham_base) :: ham_Nematic_Dirac_demo
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
        procedure, nopass :: S0
        ! procedure, nopass :: ##PROCEDURE_NAME##  ! Some other procedure defined in ham_base
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Nematic_Dirac_demo

      !#PARAMETERS START# VAR_Nematic_Dirac
      !Integer :: N_SUN = 2       ! SU(N) symmetry
      real(Kind=Kind(0.d0)) :: Dtau = 0.1d0    ! Imaginary time step size
      real(Kind=Kind(0.d0)) :: beta = 10.d0    ! Reciprocal temperature
      Integer  :: L1 = 4          ! Size of lattice in a1 direction
      Integer  :: L2 = 4          ! Size of lattice in a2 direction
      real(Kind=Kind(0.d0)) :: ham_t = 1.d0    ! Hopping amplitude of fermions
      real(Kind=Kind(0.d0)) :: Ham_h = 3.d0    ! Ising transverse field
      real(Kind=Kind(0.d0)) :: Ham_J = 1.d0    ! Ferromagnetic Ising interaction
      real(Kind=Kind(0.d0)) :: Ham_xi = 1.d0   ! Coupling Ising spins <-> fermions
      real(Kind=Kind(0.d0)) :: Ham_chem = 0.d0 ! Chemical potential
      !#PARAMETERS END#

      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:)
      Integer, allocatable :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell

      !> Second unit cell
      Type (Unit_cell), target :: latt_unit_ising

      !> Constanst for the Ising action
      Integer                :: N_ising
      Real (Kind=Kind(0.d0)) :: DW_Ising_tau(-1:1), DW_Ising_Space(-1:1)
      Integer, allocatable   :: Ising_nnlist(:,:)

      !>    Private variables for observing Z_x_ising
      Real (Kind=Kind(0.d0)) :: eq_x_ising, neq_x_ising

      !> Counting measurements
      Integer :: n_measure = 0


    contains
      
      module Subroutine Ham_Alloc_Nematic_Dirac_demo
        allocate(ham_Nematic_Dirac_demo::ham)
      end Subroutine Ham_Alloc_Nematic_Dirac_demo

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Nematic_Dirac_demo_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
      Subroutine Ham_Set

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif
          Implicit none

          integer                :: ierr, nf, unit_info
          Character (len=64)     :: file_info


#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

          ! From dynamically generated file "Hamiltonian_Nematic_Dirac_demo_read_write_parameters.F90"
          call read_parameters()

          ! Set fixed parameters of Main Hamiltonian
          N_FL = 1
          Ltrot = nint(beta/dtau)
          Projector = .false.
          Thtrot = 0
          Symm = .false.

          ! Setup the Bravais lattice
          call Ham_Latt

          ! Setup the hopping / single-particle part
          call Ham_Hop

          ! Setup the interaction.
          call Ham_V
          call Setup_Ising_action

#ifdef MPI
          If (Irank_g == 0) then
#endif
             File_info = "info"
#if defined(TEMPERING)
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif
             Open(newunit=unit_info, file=file_info, status="unknown", position="append")
             Write(unit_info,*) '====================================='
             Write(unit_info,*) 'Model is            : ', 'Nematic_Dirac'
             Write(unit_info,*) 'L1                  : ', L1
             Write(unit_info,*) 'L2                  : ', L2
             Write(unit_info,*) 'N_SUN               : ', N_SUN
             Write(unit_info,*) 'ham_t               : ', ham_t
             Write(unit_info,*) 'dtau                : ', dtau
             Write(unit_info,*) 'beta                : ', beta
             Write(unit_info,*) 'Ham_h               : ', Ham_h
             Write(unit_info,*) 'Ham_J               : ', Ham_J
             Write(unit_info,*) 'Ham_xi              : ', Ham_xi
             Write(unit_info,*) 'Ham_chem            : ', Ham_chem
             Close(unit_info)
#ifdef MPI
          Endif
#endif
        end Subroutine Ham_Set

        subroutine ham_latt()
            implicit none
  
            Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
            Integer :: I, nc, no
  
            latt_unit%Norb = 2
            latt_unit%N_coord = 4
            allocate(latt_unit%orb_pos_p(2, 2))
            latt_unit%orb_pos_p(1, :) = [0.D0,  0.D0]
            latt_unit%orb_pos_p(2, :) = [0.D0,  1.D0]
  
            latt_unit_ising%Norb = 1
            latt_unit_ising%N_coord = 2
            allocate(latt_unit_ising%orb_pos_p(1, 2))
            latt_unit_ising%orb_pos_p(1, :) = [0.D0,  0.D0]
  
            a1_p(:) = [1.D0,  1.D0]
            a2_p(:) = [1.D0, -1.D0]
            L1_p    = dble(L1)*a1_p
            L2_p    = dble(L2)*a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt)
            If ( L1 == 1 .or. L2 == 1 ) then
              Write(error_unit,*) 'One dimensional systems not implemented'
              error stop
            endif
  
            ! Ndim is needed by the Monte Carlo core!
            Ndim = Latt%N*latt_unit%Norb
  
            ! This is for the orbital structure.
            Allocate (List(Ndim,2), Invlist(Latt%N, latt_unit%Norb))
            nc = 0
            Do I = 1, Latt%N
              Do no = 1, latt_unit%Norb
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc
              Enddo
            Enddo
          end subroutine ham_latt


        subroutine ham_hop()
          implicit none

          Integer :: I, I1, J1, nc1
          complex (Kind=Kind(0.d0)) :: t_temp

          allocate(Op_T(1, N_FL))
          Call Op_make(Op_T(1, 1), Ndim)

          DO I = 1, Latt%N
            I1 = Invlist(I,1)
            Do nc1 = 1, 4
              select case (nc1)
              case (1)
                J1 = invlist(Latt%nnlist(I, 0, 0), 2)
                t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
              case (2)
                J1 = invlist(Latt%nnlist(I, 0, 1), 2)
                t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
              case (3)
                J1 = invlist(Latt%nnlist(I, -1, 1), 2)
                t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
              case (4)
                J1 = invlist(Latt%nnlist(I, -1, 0), 2)
                t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
              end select
              Op_T(1, 1)%O(I1,J1) = t_temp
              Op_T(1, 1)%O(J1,I1) = conjg(t_temp)
            Enddo
          Enddo

          Do I = 1,Ndim
             Op_T(1,1)%P(i) = i
             Op_T(1,1)%O(i,i) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
          Enddo
          Op_T(1,1)%g = -Dtau
          Op_T(1,1)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
          Call Op_set(Op_T(1,1))
        end subroutine ham_hop

        subroutine Ham_V()
          implicit none
          integer :: I, J1, nc1
          complex (Kind=Kind(0.d0)) :: t_temp

          Allocate(Op_V(Latt%N, N_FL))
          do I = 1, Latt%N
            call Op_make(Op_V(I, 1), 5)
            Op_V(I, 1)%P(1) = Invlist(I,1)
            Do nc1 = 1, 4
              select case (nc1)
              case (1)
                J1 = invlist(Latt%nnlist(I, 0, 0), 2)
                t_temp = -cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
              case (2)
                J1 = invlist(Latt%nnlist(I, 0, 1), 2)
                t_temp =  cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
              case (3)
                J1 = invlist(Latt%nnlist(I, -1, 1), 2)
                t_temp =  cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
              case (4)
                J1 = invlist(Latt%nnlist(I, -1, 0), 2)
                t_temp = -cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
              end select
              Op_V(I, 1)%P(nc1+1) = J1
              Op_V(I, 1)%O(1   ,nc1+1) = t_temp
              Op_V(I, 1)%O(nc1+1,1   ) = conjg(t_temp)
            Enddo
            Op_V(I, 1)%g     = cmplx(-dtau*ham_xi,0.d0, kind(0.D0))
            Op_V(I, 1)%alpha = cmplx(0.d0,0.d0, kind(0.D0))
            Op_V(I, 1)%type  = 1  ! Set to Ising type interaction
            Call Op_set( Op_V(I, 1) )
          Enddo

        end subroutine Ham_V
        
        Subroutine Setup_Ising_action()

          ! This subroutine sets up lists and arrays to enable an
          ! an efficient calculation of  S0(n,nt)

          Implicit none
          Integer :: I

          allocate(Ising_nnlist(Latt%N,4))
          N_ising = Latt%N
          do I = 1, Latt%N
            Ising_nnlist(I,1) = Latt%nnlist(I, 1, 0)
            Ising_nnlist(I,2) = Latt%nnlist(I, 0, 1)
            Ising_nnlist(I,3) = Latt%nnlist(I,-1, 0)
            Ising_nnlist(I,4) = Latt%nnlist(I, 0,-1)
          enddo

          ! exp(-S0(new))/exp(-S0(old)) of one Ising bond. Index is sa*sb
          DW_Ising_tau  ( 1) = tanh(Dtau*Ham_h)
          DW_Ising_tau  (-1) = 1.D0/DW_Ising_tau(1)
          DW_Ising_Space( 1) = exp(-2.d0*Dtau*Ham_J)
          DW_Ising_Space(-1) = exp( 2.d0*Dtau*Ham_J)
        End Subroutine Setup_Ising_action

!===================================================================================
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
        complex(Kind=Kind(0.d0)), Intent(In) :: Hs_new

        Integer :: nt1,I
        S0 = 1.d0
        If ( Op_V(n,1)%type == 1 ) then
            do i = 1,4
              S0 = S0*DW_Ising_space(nsigma%i(n,nt)*nsigma%i(Ising_nnlist(n,i),nt))
            enddo
            nt1 = nt +1
            if (nt1 > Ltrot) nt1 = 1
            S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt1))
            nt1 = nt - 1
            if (nt1 < 1  ) nt1 = Ltrot
            S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt1))
            If (S0 < 0.d0) Write(6,*) 'S0 : ', S0
        endif

      end function S0


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes the ratio exp(S0(new))/exp(S0(old))
!>
!> @details
!> This function computes the ratio \verbatim  e^{-S0(nsigma)}/e^{-S0(nsigma_old)} \endverbatim
!> @param [IN] nsigma_old,  Type(Fields)
!> \verbatim
!>  Old configuration. The new configuration is stored in nsigma.
!> \endverbatim
!-------------------------------------------------------------------
      Real (Kind=kind(0.d0)) Function Delta_S0_global(nsigma_old)
      
      !>  This function computes the ratio:  e^{-S0(nsigma%f)}/e^{-S0(nsigma_old)}
      Implicit none
  
      !> Arguments
      type(fields), intent(in) :: nsigma_old
  
      !> Local
      Integer :: I, nt, nt1, I1, I2, nc_J, nc_h_p, nc_h_m
  
      Delta_S0_global = 1.D0
      nc_J = 0    ! Number of parallel minus number of anti-parallel Ising bonds in space
      nc_h_p = 0  ! Number of parallel Ising bonds in time
      nc_h_m = 0  ! Number of anti-parallel Ising bonds in time
  
      Do I = 1, latt%N
        I1 = latt%nnlist(I,1,0)
        I2 = latt%nnlist(I,0,1)
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
          Character (len=:), allocatable ::  Channel


          ! Scalar observables
          Allocate ( Obs_scal(5) )
          Do I = 1,Size(Obs_scal,1)
            select case (I)
            case (1)
              N = 1;   Filename ="Part"
            case (2)
              N = 1;   Filename ="ising_z"
            case (3)
              N = 3;   Filename ="Kin_Pot_E"
            case (4)
              N = 1;   Filename ="ising_x"
              eq_x_ising  = tanh(Dtau*Ham_h)
              neq_x_ising = 1/tanh(Dtau*Ham_h)
            case (5)
              N = 3;   Filename ="m"
            case default
              Write(6,*) ' Error in Alloc_obs scal'
            end select
            call Obser_Vec_make(Obs_scal(I), N, Filename)
          enddo
        
        
          ! Equal time correlators
          Allocate ( Obs_eq(6) )
          Channel = '--'
          Do I = 1,Size(Obs_eq,1)
            select case (I)
            case (1)
              Filename ="IsingX"
              Nt = 1
              call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_ising, Channel, dtau)
            case (2)
              Filename ="IsingZ"
              Nt = 1
              call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_ising, Channel, dtau)
            case (3)
              Filename ="Green"
              Nt = 1
              call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
            case (4)
              Filename ="Den"
              Nt = 1
              call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
            case (5)
              Filename ="IsingZT"
              Nt = Ltrot+1
              call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_ising, Channel, dtau)
            case (6)
              Filename ="IsingXT"
              Nt = Ltrot+1
              call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_ising, Channel, dtau)
            case default
              Write(6,*) ' Error in Alloc_obs eq'
            end select
          enddo
        
          If (Ltau == 1) then
             ! Equal time correlators
             Allocate ( Obs_tau(2) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                  Filename ="Green"
                case (2)
                  Filename ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs tau'
                end select
                Nt = Ltrot+1
                Channel = '--'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             enddo
          endif

        End Subroutine Alloc_obs

  ! Generic funtion for kinetic energy
          function E_kin(GRC)
              Implicit none
              Complex (Kind=Kind(0.d0)), intent(in) :: GRC(Ndim,Ndim,N_FL)
            
              Complex (Kind=Kind(0.d0)) :: E_kin
            
              Integer :: n, nf, I, J
            
              E_kin = cmplx(0.d0, 0.d0, kind(0.D0))
              Do n  = 1,Size(Op_T,1)
                  Do nf = 1,N_FL
                    Do I = 1,Size(Op_T(n,nf)%O,1)
                        Do J = 1,Size(Op_T(n,nf)%O,2)
                          E_kin = E_kin + Op_T(n,nf)%O(i, j)*Grc( Op_T(n,nf)%P(I), Op_T(n,nf)%P(J), nf )
                        ENddo
                    Enddo
                  Enddo
              Enddo
              E_kin = E_kin * dble(N_SUN)
            end function E_kin
        

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
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!> @param [IN] Ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!-------------------------------------------------------------------
        subroutine Obser(GR,Phase,Ntau, Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer,                   INTENT(IN) :: Ntau
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Local
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)) :: ZP, ZS
          Integer :: I, J, nf
          ! Add local variables as needed
          Complex (Kind=Kind(0.d0)) :: Z_z_ising, Z_x_ising, Z_m
          Complex (Kind=Kind(0.d0)) :: Zkin, Zpot, Zrho
          Integer :: nc1, imj, nt1, nt, dnt, Ntau1, n, I1

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

          ZS = ZS*Mc_step_weight
          
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
          Do I = 1,Size(Obs_scal,1)
            Obs_scal(I)%N         =  Obs_scal(I)%N + 1
            Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
         Enddo
  
          Z_z_ising = cmplx(0.d0, 0.d0, kind(0.D0))
          Z_x_ising = cmplx(0.d0, 0.d0, kind(0.D0))
          Z_m       = cmplx(0.d0, 0.d0, kind(0.D0))

          nt1 = Ntau+1; if ( Ntau == Ltrot ) nt1 = 1

          ! Particle number
          Zrho = cmplx(0.d0, 0.d0, kind(0.D0))
          do I = 1, Ndim
            Zrho = Zrho + Grc(i, i, 1)
          enddo
          Zrho = Zrho * dble(N_SUN)
          Obs_scal(1)%Obs_vec(1) = Obs_scal(1)%Obs_vec(1) + Zrho *ZP*ZS
        
          ! Kinetic energy
          Zkin = E_kin(GRC)
        
          ! Potential energy
          Zpot = cmplx(0.d0, 0.d0, kind(0.D0))
          do I = 1, size(Op_V, 1)
            I1 = Op_V(I,1)%P(1)
            do nc1 = 2, size( Op_V(I,1)%O, 1 )
              Zpot = Zpot + nsigma%f(I,Ntau) * GRC(I1, Op_V(I,1)%P(nc1) ,1) * Op_V(I,1)%O(1,nc1)
              Zpot = Zpot + nsigma%f(I,Ntau) * GRC(Op_V(I,1)%P(nc1), I1 ,1) * Op_V(I,1)%O(nc1,1)
            enddo
          enddo
          Zpot = Zpot * dble(N_SUN) * ham_t * ham_xi
          Obs_scal(3)%Obs_vec(1) = Obs_scal(3)%Obs_vec(1) + Zkin *ZP*ZS
          Obs_scal(3)%Obs_vec(2) = Obs_scal(3)%Obs_vec(2) + Zpot *ZP*ZS
          Obs_scal(3)%Obs_vec(3) = Obs_scal(3)%Obs_vec(3) + Zkin+Zpot *ZP*ZS
        
        
          do I = 1,Latt%N
            Z_z_ising = Z_z_ising + nsigma%f(I,Ntau)
            Z_m = Z_m + nsigma%f(1,Ntau)*nsigma%f(I,Ntau)
        
            if ( nsigma%f(I,nt1) == nsigma%f(I,Ntau) ) then
              Z_x_ising = Z_x_ising + eq_x_ising
            else
              Z_x_ising = Z_x_ising + neq_x_ising
            endif
          enddo
          Z_z_ising = Z_z_ising/Latt%N
          Z_x_ising = Z_x_ising/Latt%N
          Z_m = Z_m/Latt%N

          Obs_scal(2)%Obs_vec(1) = Obs_scal(2)%Obs_vec(1) + Z_z_ising *ZP*ZS
          Obs_scal(4)%Obs_vec(1) = Obs_scal(4)%Obs_vec(1) + Z_x_ising *ZP*ZS
          Obs_scal(5)%Obs_vec(1) = Obs_scal(5)%Obs_vec(1) + Z_m *ZP*ZS
          Obs_scal(5)%Obs_vec(2) = Obs_scal(5)%Obs_vec(2) + Z_m**2 *ZP*ZS
          Obs_scal(5)%Obs_vec(3) = Obs_scal(5)%Obs_vec(3) + Z_m**4 *ZP*ZS
      


          ! Compute equal-time correlations
          ! counting up correlation functions
          ! DO I = 1,Size(Obs_eq,1)
          DO I = 1,2
            Obs_eq(I)%N        = Obs_eq(I)%N + 1
            Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO
        
          ! Compute Ising X-X and Z-Z correlation functions
          nt1 = Ntau+1; if ( Ntau == Ltrot ) nt1 = 1
          Do I = 1,Latt%N
            if ( nsigma%f(I,nt1) == nsigma%f(I,Ntau) ) then
              Z_x_ising = eq_x_ising
            else
              Z_x_ising = neq_x_ising
            endif
            Obs_eq(1)%Obs_Latt0(1) = Obs_eq(1)%Obs_Latt0(1) + Z_x_ising * ZP*ZS
            ! Improved estimator: <s_z> = 0 
            ! Obs_eq(2)%Obs_Latt0(1) = Obs_eq(2)%Obs_Latt0(1) + nsigma%f(I,Ntau) * ZP*ZS
            Do J = 1,Latt%N
              imj = Latt%imj(I,J)
              Obs_eq(2)%Obs_Latt(imj,1,1,1) = Obs_eq(2)%Obs_Latt(imj,1,1,1) + nsigma%f(I,Ntau) * nsigma%f(J,Ntau) * ZP*ZS
              if ( nsigma%f(J,nt1) == nsigma%f(J,Ntau) ) then
                Obs_eq(1)%Obs_Latt(imj,1,1,1) = Obs_eq(1)%Obs_Latt(imj,1,1,1) + Z_x_ising * eq_x_ising  * ZP*ZS
              else
                Obs_eq(1)%Obs_Latt(imj,1,1,1) = Obs_eq(1)%Obs_Latt(imj,1,1,1) + Z_x_ising * neq_x_ising * ZP*ZS
              endif
            enddo
          enddo
 
          ! Compute Green-function
          call Predefined_Obs_eq_Green_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(3) )

          ! Compute Density-Density correlations
          call Predefined_Obs_eq_Den_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(4) )
       
         ! Computing time-displaced X-X and Z-Z correlation functions
          n_measure = n_measure + 1
          if ( n_measure == Ltrot ) then
            n_measure = 0
        
            DO I = 5,6
              Obs_eq(I)%N        = Obs_eq(I)%N + 1
              Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
            ENDDO
          
            Z_x_ising = cmplx(0.d0, 0.d0, kind(0.D0))
            Ntau1 = Ntau+1; if ( Ntau == Ltrot ) Ntau1 = 1
          
            do dnt= 0, Ltrot
              nt = Ntau + dnt
              if ( nt > Ltrot ) nt = nt - Ltrot
              nt1 = nt+1; if ( nt == Ltrot ) nt1 = 1
              do I = 1,Latt%N
                if ( nsigma%f(I,nt1) == nsigma%f(I,nt) ) then
                  Z_x_ising = eq_x_ising
                else
                  Z_x_ising = neq_x_ising
                endif
                ! Improved estimator: <s_z> = 0 
                ! Obs_eq(5)%Obs_Latt0(1) = Obs_eq(5)%Obs_Latt0(1) + nsigma%f(I,nt) * ZP*ZS
                Obs_eq(6)%Obs_Latt0(1) = Obs_eq(6)%Obs_Latt0(1) + Z_x_ising    * ZP*ZS
                do J = 1,Latt%N
                  imj = Latt%imj(I,J)
                  Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + nsigma%f(I,Ntau) * nsigma%f(J,nt) * ZP*ZS
                  if ( nt == Ntau .and. I == J ) then
                    Obs_eq(6)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(6)%Obs_Latt(imj,dnt+1,1,1) + 1  * ZP*ZS
                  elseif ( nsigma%f(J,Ntau1) == nsigma%f(J,Ntau) ) then
                    Obs_eq(6)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(6)%Obs_Latt(imj,dnt+1,1,1) + Z_x_ising * eq_x_ising  * ZP*ZS
                  else
                    Obs_eq(6)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(6)%Obs_Latt(imj,dnt+1,1,1) + Z_x_ising * neq_x_ising * ZP*ZS
                  endif
                enddo
              enddo
            enddo
       
          endif

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
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!-------------------------------------------------------------------
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE,  Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: ZP, ZS
          ! Add local variables as needed

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS * Mc_step_weight

          ! Compute observables
          call Predefined_Obs_tau_Green_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(1) )
          call Predefined_Obs_tau_Den_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(2) )

        end Subroutine OBSERT
        
    end submodule ham_Nematic_Dirac_demo_smod
