!  Copyright (C) 2016 - 2020 The ALF project
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
!> This module defines the  Hamiltonian and observables.  Here, we have included a
!> set of predefined Hamiltonians. They include the Hubbard and SU(N) tV models
!> on honeycomb, pi-flux and square lattices.

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


    submodule (Hamiltonian_main) ham_Hubbard_Plain_Vanilla_smod

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
      use runtime_error_mod


      Implicit none
      
      type, extends(ham_base) :: ham_Hubbard_Plain_Vanilla
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
        procedure, nopass :: weight_reconstruction
        procedure, nopass :: GR_reconstruction
        procedure, nopass :: GRT_reconstruction
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Hubbard_Plain_Vanilla
      
      !#PARAMETERS START# VAR_lattice
      Character (len=64) :: Model = ''  ! Value irrelevant
      Character (len=64) :: Lattice_type = 'Square'  ! Possible Values: 'Square'
      Integer            :: L1 = 6   ! Length in direction a_1
      Integer            :: L2 = 6   ! Length in direction a_2
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_Hubbard_Plain_Vanilla
      real(Kind=Kind(0.d0)) :: ham_T    = 1.d0      ! Hopping parameter
      real(Kind=Kind(0.d0)) :: Ham_chem = 0.d0      ! Chemical potential
      real(Kind=Kind(0.d0)) :: Ham_U    = 4.d0      ! Hubbard interaction
      real(Kind=Kind(0.d0)) :: Dtau     = 0.1d0     ! Thereby Ltrot=Beta/dtau
      real(Kind=Kind(0.d0)) :: Beta     = 5.d0      ! Inverse temperature
      !logical              :: Projector= .false.   ! Whether the projective algorithm is used
      logical               :: Adiabatic= .false.   ! If  true,  and projector  true then adiabatic  switching on of  U. 
      real(Kind=Kind(0.d0)) :: Theta    = 10.d0     ! Projection parameter
      !logical              :: Symm     = .true.    ! Whether symmetrization takes place
      Integer               :: N_part   = -1        ! Number of particles in trial wave function. If N_part < 0 -> N_part = L1*L2/2
      !#PARAMETERS END#

      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit
      INTEGER :: nf_calc, nf_reconst

    contains
      
      module Subroutine Ham_Alloc_Hubbard_Plain_Vanilla
        allocate(ham_Hubbard_Plain_Vanilla::ham)
      end Subroutine Ham_Alloc_Hubbard_Plain_Vanilla

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Hubbard_Plain_Vanilla_read_write_parameters.F90"

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

          ! From dynamically generated file "Hamiltonian_Hubbard_Plain_Vanilla_read_write_parameters.F90"
          call read_parameters()
          

          If (L1 == 1) then
             Write(6,*) 'For  one-dimensional lattices set L2=1'
             Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif

          if ( (.not. projector) .and. adiabatic) then
            write(output_unit,*) "Adiabatic mode is only implemented for projective code."
            write(output_unit,*) "Overriding Adiabatic=.True. from parameter files."
          endif

          if (N_part < 0) N_part = L1*L2/2
          Ltrot  = nint(beta/dtau)
          Thtrot = 0
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot = Ltrot+2*Thtrot
          N_SUN        = 1
          N_FL         = 2

          ! Setup the Bravais lattice
          Call  Ham_Latt

          ! Setup the hopping / single-particle part
          Call  Ham_Hop

          ! Setup the interaction.
          call Ham_V

          ! Setup the trival wave function, in case of a projector approach
          if (Projector) Call Ham_Trial()

#ifdef MPI
          If (Irank_g == 0) then
#endif
             File_info = "info"
#if defined(TEMPERING)
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif

             Open(newunit=unit_info, file=file_info, status="unknown", position="append")
             Write(unit_info,*) '====================================='
             Write(unit_info,*) 'Model is      : Hubbard_Plain_Vanilla'
             Write(unit_info,*) 'Lattice is    : ', Lattice_type
             Write(unit_info,*) 'L1            : ', L1
             Write(unit_info,*) 'L2            : ', L2
             Write(unit_info,*) '# of orbitals : ', Ndim
             Write(unit_info,*) 'Symm. decomp  : ', Symm
             if (Projector) then
                Write(unit_info,*) 'Projective version'
                Write(unit_info,*) 'Theta         : ', Theta
                Write(unit_info,*) 'Tau_max       : ', beta
                Write(unit_info,*) '# of particles: ', N_part
                Write(unit_info,*) 'Adiabatic switching on of U: ', Adiabatic
             else
                Write(unit_info,*) 'Finite temperture version'
                Write(unit_info,*) 'Beta          : ', Beta
             endif
             Write(unit_info,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
             Write(unit_info,*) 't             : ', Ham_T
             Write(unit_info,*) 'Ham_U         : ', Ham_U
             Write(unit_info,*) 'Ham_chem      : ', Ham_chem
             If  (  Ham_U >=0.d0  .and.   Ham_chem   == 0.d0 )    then
                Write(unit_info,*) 'Assuming particle hole symmetry' 
             endif
             If  (  Ham_U <  0.d0  )    then
                Write(unit_info,*) 'Assuming time  reversal  symmetry' 
             endif
             if (Projector) then
                Do nf = 1,N_FL
                   Write(unit_info,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                   Write(unit_info,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
                enddo
             endif
             Close(unit_info)
#ifdef MPI
          Endif
#endif

          ! Particle-hole  symmetry   for repulsive  U  only if  chemical potential vanishes
          ! Time  reversal  symmetry  for  attractive U
          If ( (Ham_U >= 0.d0 .and.   Ham_chem  == 0.d0)  .or.  Ham_U < 0.d0  )    then
             allocate(Calc_Fl(N_FL))
             nf_calc=2
             nf_reconst=1
             Calc_Fl(nf_calc)=.True.
             Calc_Fl(nf_reconst)=.False.
          endif
          
          
        End Subroutine Ham_Set

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the  Lattice
!--------------------------------------------------------------------
        Subroutine Ham_Latt


          Implicit none
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          
          If (Lattice_Type /=  "Square")  then
             Write(6,*) 'The plain vanilla Hubbard model is only defined for the square lattice'
             Call Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          Endif
          
          Latt_unit%Norb    = 1
          Latt_unit%N_coord = 2
          allocate(Latt_unit%Orb_pos_p(Latt_unit%Norb,2))
          Latt_unit%Orb_pos_p(1, :) = [0.d0, 0.d0]

          a1_p(1) =  1.0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0  ; a2_p(2) =  1.d0
          L1_p    =  dble(L1)*a1_p
          L2_p    =  dble(L2)*a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
         
          Ndim = Latt%N*Latt_unit%Norb
          
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

          Integer :: nf , I, Ix, Iy
          allocate(Op_T(1,N_FL))
          do nf = 1,N_FL
             Call Op_make(Op_T(1,nf),Ndim)
             Do I = 1,Latt%N
                Ix = Latt%nnlist(I,1,0)
                Op_T(1,nf)%O(I,  Ix) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                Op_T(1,nf)%O(Ix, I ) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                If ( L2 > 1 ) then
                   Iy = Latt%nnlist(I,0,1)
                   Op_T(1,nf)%O(I,  Iy) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                   Op_T(1,nf)%O(Iy, I ) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                endif
                Op_T(1,nf)%O(I,  I ) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                Op_T(1,nf)%P(i) = i
             Enddo
             Op_T(1,nf)%g      = -Dtau
             Op_T(1,nf)%alpha  =  cmplx(0.d0,0.d0, kind(0.D0))
             Call Op_set(Op_T(1,nf))
          enddo


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

          Integer                              :: nf, Ix, Iy, I, n
          Real (Kind=Kind(0.d0)), allocatable  :: H0(:,:),  U0(:,:), E0(:)
          Real (Kind=Kind(0.d0))               :: Pi = acos(-1.d0), Delta = 0.01d0

          Allocate(WF_L(N_FL),WF_R(N_FL))
          do nf=1,N_FL
             Call WF_alloc(WF_L(nf),Ndim,N_part)
             Call WF_alloc(WF_R(nf),Ndim,N_part)
          enddo


          Allocate(H0(Ndim,Ndim),  U0(Ndim, Ndim),  E0(Ndim) )
          H0 = 0.d0; U0 = 0.d0;  E0=0.d0
          Do I = 1,Latt%N
             Ix = Latt%nnlist(I,1,0)
             H0(I,  Ix) = -Ham_T*(1.d0   +   Delta*cos(Pi*real(Latt%list(I,1) + Latt%list(I,2),Kind(0.d0))))
             H0(Ix, I ) = -Ham_T*(1.d0   +   Delta*cos(Pi*real(Latt%list(I,1) + Latt%list(I,2),Kind(0.d0))))
             If (L2  > 1 ) Then
                Iy = Latt%nnlist(I,0,1)
                H0(I,  Iy) = -Ham_T *(1.d0  -   Delta)
                H0(Iy, I ) = -Ham_T *(1.d0  -   Delta)
             Endif
          Enddo
          Call  Diag(H0,U0,E0)
!!$          Do I = 1,Ndim
!!$             Write(6,*) I,E0(I)
!!$          Enddo
          Do nf = 1,N_FL
             do n=1,N_part
                do I=1,Ndim
                   WF_L(nf)%P(I,n)=U0(I,n)
                   WF_R(nf)%P(I,n)=U0(I,n)
                enddo
             enddo
             WF_L(nf)%Degen = E0(N_part+1) - E0(N_part)
             WF_R(nf)%Degen = E0(N_part+1) - E0(N_part)
          enddo

          Deallocate(H0,  U0,  E0 )

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

          Integer :: nf, I, nt
          Real (Kind=Kind(0.d0)) :: X


          Allocate(Op_V(Ndim,N_FL))

          do nf = 1,N_FL
             do i  = 1, Ndim
                Call Op_make(Op_V(i,nf), 1)
             enddo
          enddo

          Do nf = 1,N_FL
             X = 1.d0
             if (nf == 2)  X = -1.d0
             Do i = 1,Ndim
                Op_V(i,nf)%P(1)   = I
                Op_V(i,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
                If  (Adiabatic)   then
                   Allocate(OP_V(i,nf)%g_t(Ltrot))
                   Op_V(i,nf)%g_t  = X*SQRT(CMPLX(DTAU*ham_U/2.d0, 0.D0, kind(0.D0)))
                   do  nt = 1, Thtrot
                      Op_V(i,nf)%g_t(nt)            = X*SQRT(CMPLX(DTAU*dble(nt)/dble(thtrot)*ham_U/2.d0, 0.D0, kind(0.D0)))
                      Op_V(i,nf)%g_t(Ltrot-(nt-1))  = X*SQRT(CMPLX(DTAU*dble(nt)/dble(thtrot)*ham_U/2.d0, 0.D0, kind(0.D0)))
                   enddo
                else
                   Op_V(i,nf)%g      = X*SQRT(CMPLX(DTAU*ham_U/2.d0, 0.D0, kind(0.D0)))
                endif
                Op_V(i,nf)%alpha  = cmplx(-0.5d0, 0.d0, kind(0.D0))
                Op_V(i,nf)%type   = 2
                Call Op_set( Op_V(i,nf) )
             Enddo
          Enddo

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
          Character (len=:), allocatable ::  Channel


          ! Scalar observables
          Allocate ( Obs_scal(4) )
          Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
                N = 1;   Filename ="Kin"
             case (2)
                N = 1;   Filename ="Pot"
             case (3)
                N = 1;   Filename ="Part"
             case (4)
                N = 1;   Filename ="Ener"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo

          ! Equal time correlators
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
                Write(6,*) "Error in Alloc_obs"
             end select
             Nt = 1
             Channel = "--"
             Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
          enddo

          If (Ltau == 1) then
             ! time-displaced correlators
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
                Nt = Ltrot+1-2*Thtrot
                If(Projector) Channel = 'T0'
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
          Integer, INTENT(IN)          :: Ntau
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Local
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, ZZ, ZXY, ZDen
          Integer :: I,J, imj, nf,  Ix, Iy
          Real    (Kind=Kind(0.d0)) :: X

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


          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Zkin = Zkin* dble(N_SUN)
          Do I = 1,Latt%N
             Ix = Latt%nnlist(I,1,0)
             Zkin = Zkin  + GRC(I,Ix,1)  + GRC(Ix,I,1)  &
                  &       + GRC(I,Ix,2)  + GRC(Ix,I,2)
          Enddo
          If (L2 > 1) then
             Do I = 1,Latt%N
                Iy = Latt%nnlist(I,0,1)
                Zkin = Zkin + GRC(I,Iy,2)  + GRC(Iy,I,2)   &
                     &      + GRC(I,Iy,1)  + GRC(Iy,I,1)
             Enddo
          Endif
          Zkin = Zkin*cmplx(-Ham_T,0.d0,Kind(0.d0))
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Ndim
             ZPot = ZPot + Grc(i,i,1) * Grc(i,i, 2)
          Enddo
          Zpot = Zpot*ham_U
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do I = 1,Ndim
             Zrho = Zrho + Grc(i,i,1) +  Grc(i,i,2)
          enddo
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS
          Obs_scal(4)%Obs_vec(1)  =    Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS

          Do I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N         =  Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign  =  Obs_eq(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo

          Do I = 1,Latt%N
             Do J = 1,Latt%N
                imj  = latt%imj(I,J)
                ZXY  = GRC(I,J,1) * GR(I,J,2) +  GRC(I,J,2) * GR(I,J,1)
                ZZ   = GRC(I,J,1) * GR(I,J,1) +  GRC(I,J,2) * GR(I,J,2)    + &
                       (GRC(I,I,2) - GRC(I,I,1))*(GRC(J,J,2) - GRC(J,J,1))

                ZDen = (GRC(I,I,1) + GRC(I,I,2)) * (GRC(J,J,1) + GRC(J,J,2)) + &
                     &  GRC(I,J,1) * GR(I,J,1)   +  GRC(I,J,2) * GR(I,J,2)
                
                Obs_eq(1)%Obs_Latt(imj,1,1,1) =  Obs_eq(1)%Obs_Latt(imj,1,1,1) + (GRC(I,J,1) + GRC(I,J,2))*ZP*ZS
                Obs_eq(2)%Obs_Latt(imj,1,1,1) =  Obs_eq(2)%Obs_Latt(imj,1,1,1) +  ZZ  *ZP*ZS
                Obs_eq(3)%Obs_Latt(imj,1,1,1) =  Obs_eq(3)%Obs_Latt(imj,1,1,1) +  ZXY *ZP*ZS
                Obs_eq(4)%Obs_Latt(imj,1,1,1) =  Obs_eq(4)%Obs_Latt(imj,1,1,1) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                Obs_eq(5)%Obs_Latt(imj,1,1,1) =  Obs_eq(5)%Obs_Latt(imj,1,1,1) +  ZDen * ZP * ZS


             enddo
             Obs_eq(5)%Obs_Latt0(1) = Obs_eq(5)%Obs_Latt0(1) + (GRC(I,I,1) + GRC(I,I,2)) *  ZP*ZS
          enddo




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
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE, Mc_step_weight )

          Use Predefined_Obs

          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY, ZDEN
          Real    (Kind=Kind(0.d0)) :: X
          Integer :: IMJ, I, J, I1, J1, no_I, no_J

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS*Mc_step_weight

          If (NT == 0 ) then
             Do I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N         =  Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign  =  Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             Enddo
          Endif
          Do I = 1,Latt%N
             Do J = 1,Latt%N
                imj  = latt%imj(I,J)

                ZZ   =      (GTT(I,I,1) -  GTT(I,I,2) ) * ( G00(J,J,1)  -  G00(J,J,2) )   &
                     &    -  G0T(J,I,1) * GT0(I,J,1)  -  G0T(J,I,2) * GT0(I,J,2)
                ZXY  =    -  G0T(J,I,1) * GT0(I,J,2)  -  G0T(J,I,2) * GT0(I,J,1)


                ZDen =   (cmplx(2.d0,0.d0,kind(0.d0)) -  GTT(I,I,1) - GTT(I,I,2) ) * &
                     &   (cmplx(2.d0,0.d0,kind(0.d0)) -  G00(J,J,1) - G00(J,J,2) )   &
                     &   -G0T(J,I,1) * GT0(I,J,1)  -  G0T(J,I,2) * GT0(I,J,2)

                Obs_tau(1)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(1)%Obs_Latt(imj,NT+1,1,1) + (GT0(I,J,1) + GT0(I,J,2))*ZP*ZS
                Obs_tau(2)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(2)%Obs_Latt(imj,NT+1,1,1) +  ZZ  *ZP*ZS
                Obs_tau(3)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(3)%Obs_Latt(imj,NT+1,1,1) +  ZXY *ZP*ZS
                Obs_tau(4)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(4)%Obs_Latt(imj,NT+1,1,1) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                Obs_tau(5)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(5)%Obs_Latt(imj,NT+1,1,1) +  ZDen * ZP * ZS


             enddo
             Obs_tau(5)%Obs_Latt0(1) = Obs_tau(5)%Obs_Latt0(1) + &
                  &                   (cmplx(2.d0,0.d0,kind(0.d0)) -  GTT(I,I,1) - GTT(I,I,2))  *  ZP*ZS
          enddo



        end Subroutine OBSERT

!--------------------------------------------------------------------
!> @brief
!> Reconstructs dependent flavors of the configuration's weight.
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!--------------------------------------------------------------------
      subroutine weight_reconstruction(weight)
         implicit none
         complex (Kind=Kind(0.d0)), Intent(inout) :: weight(:)

         weight(nf_reconst) = conjg(Weight(nf_calc))  

      end subroutine weight_reconstruction


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reconstructs dependent flavors of equal time Greens function
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!> @param [INOUT] Gr   Complex(:,:,:)
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!-------------------------------------------------------------------
      subroutine GR_reconstruction(GR)

         Implicit none

         Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: GR(Ndim,Ndim,N_FL)
         Integer :: I,J,imj
         real (kind=kind(0.d0)) :: X, ZZ

         If  (Ham_U  >= 0.d0)  then 
            Do J = 1,Ndim
               Do I = 1,Ndim
                  X=-1.0
                  imj = latt%imj(I,J)
                  if (mod(Latt%list(imj,1)+Latt%list(imj,2),2)==0) X=1.d0
                  !  write(*,*) Latt%list(I,:),Latt%list(J,:),mod(Latt%list(imj,1)+Latt%list(imj,2),2), X
                  ZZ=0.d0
                  if (I==J) ZZ=1.d0
                  GR(I,J,nf_reconst) = ZZ - X*GR(J,I,nf_calc)
               Enddo
            Enddo
         else
            Do J = 1,Ndim
               Do I = 1,Ndim
                  GR(I,J,nf_reconst) = Conjg(GR(I,J,nf_calc))
               Enddo
            Enddo
         Endif
      end Subroutine GR_reconstruction


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reconstructs dependent flavors of time displaced Greens function G0T and GT0
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!> @param [INOUT] GT0, G0T,  Complex(:,:,:)
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
!> \endverbatim
!-------------------------------------------------------------------
      Subroutine GRT_reconstruction(GT0, G0T)
         Implicit none

         Complex (Kind=Kind(0.d0)), INTENT(INOUT) :: GT0(Ndim,Ndim,N_FL), G0T(Ndim,Ndim,N_FL)
         Integer :: I,J,imj
         real (kind=kind(0.d0)) :: X, ZZ

         If (Ham_U >= 0.d0)  then
            Do J = 1,Latt%N
               Do I = 1,Latt%N
                  X=-1.0
                  imj = latt%imj(I,J)
                  if (mod(Latt%list(imj,1)+Latt%list(imj,2),2)==0) X=1.d0
                  G0T(I,J,nf_reconst) = -X*conjg(GT0(J,I,nf_calc))
                  GT0(I,J,nf_reconst) = -X*conjg(G0T(J,I,nf_calc))
               enddo
            enddo
         else
            Do J = 1,Latt%N
               Do I = 1,Latt%N
                  G0T(I,J,nf_reconst) = conjg(G0T(I,J,nf_calc))
                  GT0(I,J,nf_reconst) = conjg(GT0(I,J,nf_calc))
               enddo
            enddo
         endif
      end Subroutine GRT_reconstruction

    end submodule ham_Hubbard_Plain_Vanilla_smod
