!  Copyright (C) 2016 - 2023 The ALF project
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
!> This module defines the  Hamiltonian and observables  for the Z2 lattice gauge model coupled to
!> fermionic and Z2 matter.
!--------------------------------------------------------------------



    submodule (Hamiltonian_main) ham_Z2_Matter_smod

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
      Use Predefined_Obs

      Implicit none
      
      type, extends(ham_base) :: ham_Z2_Matter
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
        procedure, nopass :: Global_move_tau
        procedure, nopass :: Hamiltonian_set_nsigma
        procedure, nopass :: Overide_global_tau_sampling_parameters
        procedure, nopass :: Get_Delta_S0_global
        procedure, nopass :: S0
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Z2_Matter

      !#PARAMETERS START# VAR_lattice
      Character (len=64) :: Model = 'Z2_Matter'  ! Possible Values: 'Z2_Matter'
      Character (len=64) :: Lattice_type = 'Square'  ! Possible Values: 'Square'
      Integer            :: L1 = 6   ! Length in direction a_1
      Integer            :: L2 = 6   ! Length in direction a_2
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_Z2_Matter
      !Integer              :: N_SUN = 2
      real(Kind=Kind(0.d0)) :: ham_T    = 1.d0   ! Hopping for fermions
      real(Kind=Kind(0.d0)) :: ham_TZ2  = 1.d0   ! Hopping for orthogonal fermions
      real(Kind=Kind(0.d0)) :: Ham_chem = 0.d0   ! Chemical potential for fermions
      real(Kind=Kind(0.d0)) :: Ham_U    = 0.d0  ! Hubbard for fermions
      real(Kind=Kind(0.d0)) :: Ham_J    = 1.d0   ! Hopping Z2 matter fields
      real(Kind=Kind(0.d0)) :: Ham_K    = 1.d0   ! Plaquette term for gauge fields
      real(Kind=Kind(0.d0)) :: Ham_h    = 1.d0   ! sigma^x-term for matter
      real(Kind=Kind(0.d0)) :: Ham_g    = 1.d0   ! tau^x-term for gauge
      real(Kind=Kind(0.d0)) :: Dtau     = 0.1d0  ! Thereby Ltrot=Beta/dtau
      real(Kind=Kind(0.d0)) :: Beta     = 10.d0  ! Inverse temperature
      !logical              :: Projector = .false.  ! Whether the projective algorithm is used
      real(Kind=Kind(0.d0)) :: Theta    = 0.d0      ! Projection parameter
      Integer               :: N_part   = -1        ! Number of particles in trial wave function. If N_part < 0 -> N_part = L1*L2/2
      !#PARAMETERS END#
      
      Type (Lattice),        target :: Latt
      Type (Unit_cell),      target :: Latt_unit
      Logical                :: One_dimensional
      Integer, allocatable   :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell
      real (Kind=Kind(0.d0)) :: Zero = 1.D-10

      !>    Storage for the Ising action
      Real (Kind=Kind(0.d0)) :: DW_Ising_tau(-1:1), DW_Ising_Space(-1:1), DW_Ising_Flux(-1:1,-1:1)
      Real (Kind=Kind(0.d0)) :: DW_Matter_tau(-1:1), DW_Ising_Matter(-1:1)
      Integer, allocatable   :: Field_list(:,:,:), Field_list_inv(:,:)

    contains
      
      module Subroutine Ham_Alloc_Z2_Matter
        allocate(ham_Z2_Matter::ham)
      end Subroutine Ham_Alloc_Z2_Matter

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Z2_Matter_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian.  Called by main.
!--------------------------------------------------------------------
      Subroutine Ham_Set

#if defined (MPI) || defined(TEMPERING)
          use mpi
#endif
          Implicit none

          integer :: ierr, nf, unit_info
          Character (len=64) :: file_info
          
#ifdef MPI
          Integer        :: Isize, Irank, igroup, irank_g, isize_g
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
          !if ( irank_g == 0 )   write(6,*) "Mpi Test", igroup, isize_g
#endif

          ! From dynamically generated file "Hamiltonian_Z2_Matter_read_write_parameters.F90"
          call read_parameters()
          
          Call Ham_latt

          if ( str_to_upper(Model) == "Z2_MATTER" ) then
             N_FL = 1
             If ( str_to_upper(Lattice_type)  /= "SQUARE" ) then
                Write(error_unit,*) "Ham_set: Z2_Matter is only implemented for a square lattice"
                CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
             Endif
          else
             Write(error_unit,*) "Ham_set: Model not yet implemented!"
             CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif
          if (N_part < 0) N_part = L1*L2/2
          If (Abs(Ham_T) < Zero ) then
              Ham_J = 0.d0 ! Matter-Ising interction
              Ham_h = 0.d0
          endif
          If (Abs(Ham_TZ2) < Zero ) then
              Ham_J = 0.d0 ! Matter-Ising interction
              Ham_K = 0.d0 ! Flux
              Ham_g = 0.d0
          endif

          Call Ham_hop
          Ltrot = nint(beta/dtau)
          Thtrot = 0
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot = Ltrot+2*Thtrot
          
          If  ( str_to_upper(Model) == "Z2_MATTER" )  Call Setup_Ising_action_and_field_list
          
          call Ham_V
           
          if (Projector) Call Ham_Trial()

#if defined(MPI)
           If (irank_g == 0 ) then
#endif
              File_info = "info"
#if defined(TEMPERING)
              write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif
              Open(newunit=unit_info, file=file_info, status="unknown", position="append")
              Write(unit_info,*) '====================================='
              Write(unit_info,*) 'Model is      : Z2_Matter'
              Write(unit_info,*) 'Lattice is    : ', Lattice_type
              Write(unit_info,*) '# of orbitals : ', Ndim
              if (Projector) then
                 Write(unit_info,*) 'Projective version'
                 Write(unit_info,*) 'Theta         : ', Theta
                 Write(unit_info,*) 'Tau_max       : ', beta
                 Write(unit_info,*) '# of particles: ', N_part
              else
                 Write(unit_info,*) 'Finite temperture version'
                 Write(unit_info,*) 'Beta          : ', Beta
              endif
              Write(unit_info,*) 'dtau,Ltrot    : ', dtau,Ltrot
              Write(unit_info,*) 'N_SUN         : ', N_SUN
              Write(unit_info,*) 'N_FL          : ', N_FL
              If (Abs(Ham_T) < Zero) then
                 Write(unit_info,*) 't_Z2          : ', Ham_TZ2
                 Write(unit_info,*) 'g_Z2          : ', Ham_g
                 Write(unit_info,*) 'K_Gauge       : ', Ham_K
              elseif (Abs(Ham_TZ2) < Zero) then
                 Write(unit_info,*) 't_fermion     : ', Ham_T
                 Write(unit_info,*) 'h_Matter      : ', Ham_h
              else
                 Write(unit_info,*) 't_Z2          : ', Ham_TZ2
                 Write(unit_info,*) 'g_Z2          : ', Ham_g
                 Write(unit_info,*) 'K_Gauge       : ', Ham_K
                 Write(unit_info,*) 't_fermion     : ', Ham_T
                 Write(unit_info,*) 'h_Matter      : ', Ham_h
                 Write(unit_info,*) 'J_Gauge_Z2    : ', Ham_J
              endif
              Write(unit_info,*) 'Ham_chem      : ', Ham_chem
              Write(unit_info,*) 'Ham_U         : ', Ham_U
              if (Projector) then
                 Do nf = 1,N_FL
                    Write(unit_info,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                    Write(unit_info,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
                 enddo
              endif
              close(unit_info)
#if defined(MPI)
           endif
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
          If ( L1 == 1 .or. L2 == 1 ) then
             Write(error_unit,*) 'Ham_Latt: One dimensional systems are not included '
             CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif
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


          Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:)

          Real (Kind=Kind(0.d0) ), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:),&
               &                                  Phi_Y_vec(:),  Ham_T2_vec(:),  Ham_Lambda_vec(:)
          Integer, allocatable ::   N_Phi_vec(:)

          Logical ::  Bulk = .False.,  Checkerboard = .False.

          Allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), &
               &    Phi_X_vec(N_FL), Phi_Y_vec(N_FL), N_Phi_vec(N_FL), Ham_Lambda_vec(N_FL) )

          ! Here we consider no N_FL  dependence of the hopping parameters.
          Ham_T_vec      = 0.d0
          Ham_Tperp_vec  = 0.d0
          Ham_Chem_vec   = Ham_Chem
          Phi_X_vec      = 0.d0
          Phi_Y_vec      = 0.d0
          Ham_T2_vec     = 0.d0
          Ham_Lambda_vec = 0.d0
          N_Phi_vec      = 0

          Call  Set_Default_hopping_parameters_square(Hopping_Matrix,Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
               &                                      Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
          Call  Predefined_Hoppings_set_OPT(Hopping_Matrix,List,Invlist,Latt,  Latt_unit,  Dtau, Checkerboard, Symm, OP_T )

          Deallocate (Ham_T_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
               &                                   N_Phi_vec,  Ham_Lambda_vec )


        end Subroutine Ham_Hop
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

          Integer :: nf, I, I1, I2,  nc, nc1,  J, N_Field_type, N_ops
          Real (Kind=Kind(0.d0)) :: X

          N_ops = size(Field_list_inv,1)
          Allocate(Op_V(N_ops, N_FL))


          !Field_list_inv(nc,1) = I
          !Field_list_inv(nc,2) = n_orientation
          !Field_list_inv(nc,3) = N_Field_type
          Do nc = 1, N_ops
             N_Field_type = Field_list_inv(nc,3)
             select case (N_Field_type)
             case (3 ) ! Hubbard
                I = Field_list_inv(nc,1)
                do nf = 1,N_FL
                   Call Predefined_Int_U_SUN( OP_V(nc,nf), I, N_SUN, DTAU, Ham_U  )
                enddo
             case (1 ) ! Z2_Gauge
                I = Field_list_inv(nc,1)
                select case ( Field_list_inv(nc,2) )
                case (1)
                   I1 = Latt%nnlist(I,1,0)
                case (2)
                   I1 = Latt%nnlist(I,0,1)
                end select
                do nf = 1,N_FL
                   Call Predefined_Int_Ising_SUN( OP_V(nc,nf), I, I1, DTAU, -Ham_TZ2  )
                enddo
             case (2 ) ! Bond_Matter
                I = Field_list_inv(nc,1)
                select case ( Field_list_inv(nc,2) )
                case (1)
                   I1 = Latt%nnlist(I,1,0)
                case (2)
                   I1 = Latt%nnlist(I,0,1)
                end select
                do nf = 1,N_FL
                   Call Predefined_Int_Ising_SUN( OP_V(nc,nf), I, I1, DTAU, -Ham_T  )
                enddo
             case (4 ) ! Site Matter
                I = Field_list_inv(nc,1)
                do nf = 1,N_FL
                   Call OP_Make( Op_V(nc,nf),1)
                   Op_V(nc,nf)%P(1)   = I1
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0,0.d0, kind(0.D0))
                   Op_V(nc,nf)%g      = cmplx(0.d0,0.d0, kind(0.D0))
                   Op_V(nc,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 1
                   Call Op_set( Op_V(nc,1) )
                enddo
             end select
          end Do

        end Subroutine Ham_V

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
             H0(I,  Ix) = -(1.d0   +   Delta*cos(Pi*real(Latt%list(I,1) + Latt%list(I,2),Kind(0.d0))))
             H0(Ix, I ) = -(1.d0   +   Delta*cos(Pi*real(Latt%list(I,1) + Latt%list(I,2),Kind(0.d0))))
             If (L2  > 1 ) Then
                Iy = Latt%nnlist(I,0,1)
                H0(I,  Iy) = -(1.d0  -   Delta)
                H0(Iy, I ) = -(1.d0  -   Delta)
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
!> Single spin flip S0 ratio
!> @details
!> S0=exp(-S0(new))/exp(-S0(old)) where the new configuration correpsonds to the old one up to
!> a spin flip of Operator n on time slice nt
!> @details
!--------------------------------------------------------------------
        Real (Kind=Kind(0.d0)) function S0(n,nt,Hs_new)
          Implicit none
          Integer, Intent(IN) :: n,nt
          complex (Kind=Kind(0.d0)), Intent(In) :: Hs_new

          !Local
          Integer :: nt1,I, F1,F2,I1,I2,I3,  n_orientation, n_m

          !> Ratio for local spin-flip  of gauge field only.

          S0 = 1.d0
          If ( Abs(Ham_TZ2) > Zero ) then

             !Field_list_inv(nc,1) = I1
             !Field_list_inv(nc,2) = n_orientation
             !Field_list_inv(nc,3) = N_Field_type

             If (Field_list_inv(n,3) == 1 ) then
                
                If (Abs(Ham_T) > Zero ) then
                   I              = Field_list_inv(n,1)
                   n_orientation  = Field_list_inv(n,2)
                   n_m            = Field_list(I,n_orientation,2)
                   S0 = S0* DW_Ising_Matter(nsigma%i(n,nt)*nsigma%i(n_m,nt) )  ! Coupling to matter field.
                endif

                If (Projector) then
                   if   (nt == Ltrot)  then
                      S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt-1))
                   elseif ( nt == 1 ) then
                      S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt+1))
                   else
                      S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt+1))*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt-1))
                   endif
                else
                   nt1 = nt +1
                   if (nt1 > Ltrot) nt1 = 1
                   S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt1))
                   nt1 = nt - 1
                   if (nt1 < 1  ) nt1 = Ltrot
                   S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt1))
                endif
                ! Magnetic flux term
                I1 = Field_list_inv(n,1)
                if ( Field_list_inv(n,2) == 1 ) then
                   !     I2
                   !     I1 I3
                   I2 = Latt%nnlist(I1,0,1 )
                   I3 = Latt%nnlist(I1,1,0 )
                   F1 = nsigma%i(n,nt)*nsigma%i(Field_list(I1,2,1),nt)* nsigma%i(Field_list(I2,1,1),nt)*nsigma%i(Field_list(I3,2,1),nt)
                   !     I1
                   !     I2 I3
                   I2 = Latt%nnlist(I1,0,-1)
                   I3 = Latt%nnlist(I1,1,-1)
                   F2 = nsigma%i(n,nt)*nsigma%i(Field_list(I2,1,1),nt)* nsigma%i(Field_list(I2,2,1),nt)*nsigma%i(Field_list(I3,2,1),nt)
                else
                   !    I3
                   !    I2  I1
                   I2 = Latt%nnlist(I1,-1,0 )
                   I3 = Latt%nnlist(I1,-1,1 )
                   F1 = nsigma%i(n,nt)*nsigma%i(Field_list(I2,1,1),nt)* nsigma%i(Field_list(I2,2,1),nt)*nsigma%i(Field_list(I3,1,1),nt)
                   !    I2
                   !    I1  I3
                   I2 = Latt%nnlist(I1,0,1)
                   I3 = Latt%nnlist(I1,1,0)
                   F2 = nsigma%i(n,nt)*nsigma%i(Field_list(I1,1,1),nt)* nsigma%i(Field_list(I2,1,1),nt)*nsigma%i(Field_list(I3,2,1),nt)
                endif
                S0 = S0*DW_Ising_Flux(F1,F2)
             else
                S0 = 1.d0
             endif

          endif

        end function S0

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
!> The move will be carried out with prbablity  T0 ( sigma -> sigma_new ).   If T0 ( sigma -> sigma_new ) > Ranf
!> then T0_Proposal_ratio  will be initialized. Otherwise the latter quantity is set to zero.
!--------------------------------------------------------------------
        Subroutine Global_move_tau(T0_Proposal_ratio, S0_ratio, &
             &                     Flip_list, Flip_length,Flip_value,ntau)


          Implicit none
          Real (Kind= kind(0.d0)),INTENT(OUT) :: T0_Proposal_ratio, S0_ratio
          Integer                ,INTENT(OUT) :: Flip_list(:)
          Complex (Kind= Kind(0.d0)),INTENT(out) :: Flip_value(:)
          Integer, INTENT(OUT) :: Flip_length
          Integer, INTENT(IN)    :: ntau


          !Local
          Integer                   ::  ns , nc, n_op, n_op1, ntau_p1, ntau_m1, I, n
          Integer, allocatable      ::  Isigma1(:),Isigma2(:),Isigma3(:)
          Real  (Kind = Kind(0.d0)) ::  S0_Matter, T0_Proposal

          ! Write(6,*) 'In GLob_move', m,direction,ntau, size(Flip_list,1), Size(Flip_value,1), Flip_list(1)
          ! Ising from n_op = 1,Latt_unit%N_coord*Ndim
          ! Hubbard from n_op = Latt_unit%N_coord*Ndim +1, Size(OP_V,1) = Latt_unit%N_coord*Ndim +  Ndim
          ! Write(6,*) 'Global_move_tau ' , S0(Flip_list(1),ntau)

          Allocate (Isigma1(Latt%N), Isigma2(Latt%N), Isigma3(Latt%N) )

          I  =  nranf(Latt%N)
          Flip_length = 4
          S0_Matter = 1.d0
          do n = 1,4
             select case(n)
             case (1)
                n_op  = Field_list(I,1,2)
                if ( Abs(Ham_TZ2) > Zero) then
                   n_op1 = Field_list(I,1,1)
                   S0_Matter = S0_Matter*DW_Ising_Matter( nsigma%i(n_op,ntau) *  nsigma%i(n_op1,ntau) )
                endif
             case (2)
                n_op  = Field_list(I,2,2)
                if ( Abs(Ham_TZ2) > Zero) then
                   n_op1 = Field_list(I,2,1)
                   S0_Matter = S0_Matter*DW_Ising_Matter( nsigma%i(n_op,ntau) *  nsigma%i(n_op1,ntau) )
                endif
             case (3)
                n_op  = Field_list(latt%nnlist(I,-1,0),1,2)
                if ( Abs(Ham_TZ2) > Zero) then
                   n_op1 = Field_list(latt%nnlist(I,-1,0),1,1)
                   S0_Matter = S0_Matter*DW_Ising_Matter( nsigma%i(n_op,ntau) *  nsigma%i(n_op1,ntau) )
                endif
             case (4)
                n_op  = Field_list(latt%nnlist(I,0,-1),2,2)
                if ( Abs(Ham_TZ2) > Zero) then
                   n_op1 = Field_list(latt%nnlist(I,0,-1),2,1)
                   S0_Matter = S0_Matter*DW_Ising_Matter( nsigma%i(n_op,ntau) *  nsigma%i(n_op1,ntau) )
                endif
             case default
                Write(error_unit,*) 'Global_move_tau: Error'
                CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
             end select
             Flip_list(n)  = n_op
             Flip_value(n) = nsigma%flip(n_op,ntau)  
          enddo
          If ( I == Latt%N )   then
             Flip_length   = 5
             n             = 5
             n_op          = Field_list(Latt%N,3,4)
             Flip_list (n) = n_op
             Flip_value(n) = nsigma%flip(n_op,ntau)  
          endif

          If (Projector) then
             if ( ntau == Ltrot ) then
                Call Hamiltonian_set_Z2_matter(Isigma2,ntau   )
                Call Hamiltonian_set_Z2_matter(Isigma3,ntau-1 )
                S0_Matter = S0_Matter* DW_Matter_tau( Isigma2(I)*Isigma3(I) )
             elseif ( ntau == 1 ) then
                Call Hamiltonian_set_Z2_matter(Isigma1,ntau + 1)
                Call Hamiltonian_set_Z2_matter(Isigma2,ntau    )
                S0_Matter = S0_Matter*DW_Matter_tau ( Isigma1(I)*Isigma2(I) )
             else
                Call Hamiltonian_set_Z2_matter(Isigma1,ntau +1 )
                Call Hamiltonian_set_Z2_matter(Isigma2,ntau    )
                Call Hamiltonian_set_Z2_matter(Isigma3,ntau -1 )
                S0_Matter = S0_Matter*DW_Matter_tau ( Isigma1(I)*Isigma2(I) ) * DW_Matter_tau( Isigma2(I)*Isigma3(I) )
             endif
          else
             ntau_p1 = ntau + 1
             if (ntau == Ltrot) ntau_p1 = 1
             ntau_m1 = ntau -1
             if (ntau == 1    ) ntau_m1 = Ltrot
             Call Hamiltonian_set_Z2_matter(Isigma1,ntau_p1)
             Call Hamiltonian_set_Z2_matter(Isigma2,ntau   )
             Call Hamiltonian_set_Z2_matter(Isigma3,ntau_m1)
             !  Check the dynamics and the ergodicity
             S0_Matter = S0_Matter*DW_Matter_tau ( Isigma1(I)*Isigma2(I) ) * DW_Matter_tau( Isigma2(I)*Isigma3(I) )
          endif

          T0_Proposal       =  1.d0 - 1.d0/(1.d0+S0_Matter)
          !  Move acceptance probability.
          If ( T0_Proposal > Ranf_wrap() )  then
             T0_Proposal_ratio =  1.d0 / S0_Matter
             !T0_Proposal       =  1.d0
             !T0_Proposal_ratio =  1.d0
          else
             T0_Proposal_ratio = 0.d0
          endif
          S0_ratio          =  S0_Matter

          Deallocate (Isigma1,Isigma2, Isigma3)

        end Subroutine Global_move_tau

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
        Real (Kind=kind(0.d0)) Function Get_Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none

          !> Arguments
          type (Fields),  Intent(IN)  :: nsigma_old
          !> Local
          Integer :: I,n,n1,n2,n3,n4,nt,nt1, nc_F, nc_J, nc_h_p, nc_h_m, n1_m, n4_m
          Real (Kind=kind(0.d0)) :: exp_delta_S0


          exp_delta_S0 = 1.d0
          If ( Model == "Z2_Matter" ) then
             nc_F = 0
             nc_J = 0
             nc_h_p = 0
             nc_h_m = 0
             Do I = 1,Latt%N
                n1   = Field_list(I,1,1)
                n1_m = Field_list(I,1,2)
                n2   = Field_list(Latt%nnlist(I,1,0),2,1)
                n3   = Field_list(Latt%nnlist(I,0,1),1,1)
                n4   = Field_list(I,2,1)
                n4_m = Field_list(I,2,2)
                do nt = 1,Ltrot
                   nt1 = nt +1
                   if (nt == Ltrot) nt1 = 1
                   if (nsigma%i(n1,nt) == nsigma%i(n1,nt1) ) then
                      nc_h_p = nc_h_p + 1
                   else
                      nc_h_m = nc_h_m + 1
                   endif
                   if (nsigma_old%i(n1,nt) == nsigma_old%i(n1,nt1) ) then
                      nc_h_p = nc_h_p - 1
                   else
                      nc_h_m = nc_h_m - 1
                   endif

                   if (nsigma%i(n4,nt) == nsigma%i(n4,nt1) ) then
                      nc_h_p = nc_h_p + 1
                   else
                      nc_h_m = nc_h_m + 1
                   endif
                   if (nsigma_old%i(n4,nt) == nsigma_old%i(n4,nt1) ) then
                      nc_h_p = nc_h_p - 1
                   else
                      nc_h_m = nc_h_m - 1
                   endif

                   nc_F = nc_F + nsigma%i    (n1,nt)*nsigma%i    (n2,nt)*nsigma%i    (n3,nt)*nsigma%i    (n4,nt)  &
                        &      - nsigma_old%i(n1,nt)*nsigma_old%i(n2,nt)*nsigma_old%i(n3,nt)*nsigma_old%i(n4,nt)

                   nc_J = nc_J + nsigma%i(n1,nt)*nsigma%i(n1_m,nt) + &
                        &        nsigma%i(n4,nt)*nsigma%i(n4_m,nt) - &
                        &        nsigma_old%i(n1,nt)*nsigma_old%i(n1_m,nt) - &
                        &        nsigma_old%i(n4,nt)*nsigma_old%i(n4_m,nt)

                enddo
             enddo
             exp_delta_S0 = ( sinh(Dtau*Ham_h)**nc_h_m ) * (cosh(Dtau*Ham_h)**nc_h_p) * &
                  &            exp( -Dtau*(Ham_K*real(nc_F,kind(0.d0)) + Ham_J*real(nc_J,kind(0.d0))))
             Get_Delta_S0_global = log(exp_delta_S0) ! This can be done better
          endif
        end Function Get_Delta_S0_global

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> This routine sets storage to estimate Ising action as well as the list and types of fields (HS or Ising)
!> so as to know if the field nsimg%i(nc,nt) corresponds to a HS field for the U term, an Ising gauge field,
!> or a Z_2 matter field.
!>  Field_list_inv(nc,1) = I              ! Postion on lattice
!>  Field_list_inv(nc,2) = n_orientation  ! Orientation 1=a_x, 2 = a_y , 3 = no-orientation for on-site interaction.
!>  Field_list_inv(nc,3) = Field_type     ! 1 = gauge Field, 2 = Bond matter field, 3 =  HS for Hubbard, 4 = Matter field at a given site.

!--------------------------------------------------------------------
        Subroutine Setup_Ising_action_and_field_list

          ! This subroutine sets up lists and arrays so as to enable an
          ! an efficient calculation of  S0(n,nt)

          Integer :: nc, nth, n, n1, n2, n3, n4, I, I1, n_orientation, Ix, Iy, N_Field_type,  N_Pos
          Integer :: N_ops
          Real (Kind=Kind(0.d0)) :: X_p(2)


          N_ops = 0
          If (Abs(Ham_U)   > Zero )   N_ops = N_ops + Latt%N                          !  Hubbard
          If (Abs(Ham_TZ2) > Zero )   N_ops = N_ops + Latt%N*Latt_unit%N_coord        !  Z2 gauge fields
          If (Abs(Ham_T  ) > Zero )   N_ops = N_ops + Latt%N*Latt_unit%N_coord + 1    !  Matter fields.

          ! Setup list of bonds for the square lattice.
          Allocate ( Field_list(Latt%N,3,4),  Field_list_inv(N_ops,3) )
          nc = 0
          If (Abs(Ham_U)   > Zero )  then
             DO I = 1,Latt%N
                nc = nc + 1
                N_Pos         = I
                N_orientation = 3
                N_Field_type  = 3
                Field_list_inv(nc,1) = N_Pos
                Field_list_inv(nc,2) = N_Orientation
                Field_list_inv(nc,3) = N_Field_type
                Field_list(N_pos,n_orientation,N_field_type) = nc
             Enddo
          Endif
          If (Abs(Ham_TZ2)   > Zero )  then
             N_Field_type = 1
             DO I = 1,Latt%N
                Ix = Latt%list(I,1)
                Iy = Latt%list(I,2)
                if (mod(Ix + Iy,2) == 0 ) then
                   do n = 1,4
                      nc = nc + 1
                      select case (n)
                      case (1)
                         I1 = I                  ;  n_orientation  = 1
                      case (2)
                         I1 = I                  ;  n_orientation  = 2
                      case (3)
                         I1 = latt%nnlist(I,-1,0);  n_orientation  = 1
                      case (4)
                         I1 = latt%nnlist(I,0,-1);  n_orientation  = 2
                      case default
                         Write(6,*) ' Error in Setup_Ising_action '
                      end select
                      Field_list(I1,n_orientation,N_Field_type) = nc
                      Field_list_inv(nc,1) = I1
                      Field_list_inv(nc,2) = n_orientation
                      Field_list_inv(nc,3) = N_Field_type
                      ! The bond is given by  I1, I1 + a_(n_orientation).
                   enddo
                endif
             Enddo
          Endif
          If (Abs(Ham_T)   > Zero )  then
             N_Field_type = 2
             DO I = 1,Latt%N
                Ix = Latt%list(I,1)
                Iy = Latt%list(I,2)
                if (mod(Ix + Iy,2) == 0 ) then
                   do n = 1,4
                      nc = nc + 1
                      select case (n)
                      case (1)
                         I1 = I                  ;  n_orientation  = 1
                      case (2)
                         I1 = I                  ;  n_orientation  = 2
                      case (3)
                         I1 = latt%nnlist(I,-1,0);  n_orientation  = 1
                      case (4)
                         I1 = latt%nnlist(I,0,-1);  n_orientation  = 2
                      case default
                         Write(6,*) ' Error in Setup_Ising_action '
                      end select
                      Field_list(I1,n_orientation,N_Field_type) = nc
                      Field_list_inv(nc,1) = I1
                      Field_list_inv(nc,2) = n_orientation
                      Field_list_inv(nc,3) = N_Field_type
                      ! The bond is given by  I1, I1 + a_(n_orientation).
                   enddo
                endif
             Enddo
             nc = nc + 1
             I = Latt%N
             n_orientation = 3
             N_Field_type  = 4
             Field_list(I,n_orientation,N_Field_type) = nc
             Field_list_inv(nc,1) = I
             Field_list_inv(nc,2) = n_orientation
             Field_list_inv(nc,3) = N_Field_type
          Endif

          !Test
          !Do I = 1,Latt%N
          !   Write(6,*)
          !   Write(6,*) Latt%list(I,1), Latt%list(I,2), I
          !   Write(6,*) Field_list(I,1), Field_list(I,2), Field_list( latt%nnlist(I,-1,0),1), Field_list( latt%nnlist(I,0,-1),2)
          !Enddo

          DW_Ising_tau = 1.d0
          If (Ham_g > Zero ) then
             DW_Ising_tau  ( 1) = tanh(Dtau*Ham_g)
             DW_Ising_tau  (-1) = 1.D0/DW_Ising_tau(1)
          endif
          DO n = -1,1,2
             do n1 = -1,1,2
                DW_Ising_Flux(n,n1) = exp( Dtau*Ham_K*(dble(n) +  dble(n1) ))/ exp(  -Dtau*Ham_K*(dble(n) +  dble(n1) ))
             enddo
          enddo
          DW_Ising_Matter( 1) = exp( 2.d0*Dtau*Ham_J)
          DW_Ising_Matter(-1) = exp(-2.d0*Dtau*Ham_J)
          DW_Matter_tau = 1.d0
          If (Ham_h > Zero ) then
             DW_Matter_tau  ( 1) = tanh(Dtau*Ham_h)
             DW_Matter_tau  (-1) = 1.D0/DW_Matter_tau(1)
          endif

        End Subroutine Setup_Ising_action_and_field_list

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
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Nt
          Character (len=64) ::  Filename
          Character (len=:), allocatable ::  Channel

          ! Scalar observables
          Allocate ( Obs_scal(4) )
          Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
                N = 1;   Filename ="Part"
             case (2)
                N = 2;   Filename ="Flux"
             case (3)
                N = 2;   Filename ="X"
             case (4)
                N = 1;   Filename ="Q"
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
                Filename ="Greenf"
             case (2)
                Filename ="SpinZ"
             case (3)
                Filename ="Den"
             case (4)
                Filename ="Green"
             case (5)
                Filename ="Q"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Nt = 1
             Channel = '--'
             Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
          enddo

          If (Ltau == 1) then
             ! Time-displaced correlators
             Allocate ( Obs_tau(3) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Channel = 'P' ; Filename ="Green"
                case (2)
                   Channel = 'PH'; Filename ="SpinZ"
                case (3)
                   Channel = 'PH'; Filename ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = Ltrot+1
                If(Projector) Channel = 'T0'
                Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             enddo
          endif

        end Subroutine Alloc_obs

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
!--------------------------------------------------------------------
        Subroutine Obser(GR,Phase,Ntau, Mc_step_weight)

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight


          !Local
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin_mat, ZPot_mat, Z, ZP,ZS, Z1, Z2, ZN
          Complex (Kind=Kind(0.d0)) :: ZQ, ZSTAR, ZQT, ZQTT
          Integer :: I,J, imj, nf, dec, I1, I2,I3,I4, J1,J2,J3,J4, no_I, no_J,  iFlux_tot,  &
               &     no, no1, ntau1, ntau2, L_Vison, L_Wilson, n, nx,ny
          Real (Kind=Kind(0.d0)) :: X_ave, X, XI1,XI2,XI3,XI4, X_p(2)
          Integer,  allocatable  :: Isigma(:), Isigmap1(:)
          Integer ::  IB_x, IB_y, Ix, Iy

          Real (Kind=Kind(0.d0)) :: X_star_i, X_star_j,  X_star_ij

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS*Mc_step_weight

          ZN =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))

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

          If ( abs(Ham_T) > Zero ) then
             ntau1 = ntau + 1
             If (ntau == Ltrot )  ntau1 = 1
             Allocate ( Isigma(Latt%N), Isigmap1(Latt%N) )
             Call Hamiltonian_set_Z2_matter(Isigma  ,ntau  )
             Call Hamiltonian_set_Z2_matter(Isigmap1,ntau1 )
             
             iFlux_tot = 0
             Do I = 1, Ndim
                iFlux_tot = iFlux_tot + iFlux(I,Ntau,2)
             Enddo
             Obs_scal(2)%Obs_vec(2)  =   Obs_scal(2)%Obs_vec(2) + cmplx(dble(iFlux_tot),0.d0,kind(0.d0))*ZP*ZS

             X_ave = 0.d0
             Do I = 1,Latt%N
                X_ave = X_ave + tau_x(I,ntau, Isigma, Isigmap1)
             Enddo
             Obs_scal(3)%Obs_vec(2)  =  Obs_scal(3)%Obs_vec(2) + cmplx(X_ave,0.d0,kind(0.d0)) * ZP*ZS
             
          endif

          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Latt%N
                Zrho = Zrho + Grc(i,i,nf)
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zrho * ZP*ZS


          If ( abs(Ham_TZ2) > Zero ) then
             iFlux_tot = 0
             Do I = 1, Ndim
                iFlux_tot = iFlux_tot + iFlux(I,Ntau,1)
             Enddo
             Obs_scal(2)%Obs_vec(1)  =   Obs_scal(2)%Obs_vec(1) + cmplx(dble(iFlux_tot),0.d0,kind(0.d0))*ZP*ZS
             
             X_ave = 0.d0
             Do I = 1,Latt%N
                do no = 1,2
                   X_ave = X_ave + sigma_x(i,no,ntau)
                Enddo
             Enddo
             Obs_scal(3)%Obs_vec(1)  =  Obs_scal(3)%Obs_vec(1) + cmplx(X_ave,0.d0,kind(0.d0)) * ZP*ZS

          Endif

          ! Constraint.
          do I  = 1, Latt%N
             Z1 = cmplx(1.d0,0.d0,kind(0.d0)) - cmplx(2.d0,0.d0,kind(0.d0))*GRC(I,I,1)
             Z1 = Z1**(N_SUN)
             Z1 = Z1 * cmplx( star_sigma_x(I,ntau)*tau_x(I,ntau, Isigma, Isigmap1) ,0.d0,Kind(0.d0))
             Obs_scal(4)%Obs_vec(1)  =  Obs_scal(4)%Obs_vec(1) + Z1*ZP*ZS
          enddo


          ! Green function for electron.
          Obs_eq(1)%N        = Obs_eq(1)%N + 1
          Obs_eq(1)%Ave_sign = Obs_eq(1)%Ave_sign + real(ZS,kind(0.d0))
          If ( abs(Ham_T) > Zero ) then
             Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
             Do I1 = 1,Latt%N
                Do J1 = 1,Latt%N
                   imj = latt%imj(I1,J1)
                   ! Green_fermion
                   Z1 = cmplx(real(Isigma(I1)*Isigma(J1), kind(0.d0)), 0.d0,kind(0.d0))
                   Obs_eq(1)%Obs_Latt(imj,1,1,1) =  Obs_eq(1)%Obs_Latt(imj,1,1,1) + &
                        &               Z * Z1*GRC(I1,J1,1) *  ZP*ZS
                enddo
             enddo
          endif

          ! Compute spin-spin, Green, and den-den correlation functions
          Call Predefined_Obs_eq_SpinSUN_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(2) )
          Call Predefined_Obs_eq_Den_measure    ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(3) )
          Call Predefined_Obs_eq_Green_measure  ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(4) )

          !  Constraint  correlation
          Obs_eq(5)%N        = Obs_eq(5)%N + 1
          Obs_eq(5)%Ave_sign = Obs_eq(5)%Ave_sign + real(ZS,kind(0.d0))
          Do I = 1,Latt%N
             Do J = 1,Latt%N
                imj = latt%imj(I,J)
                if ( i == j ) then
                   Z1 = cmplx(1.d0,0.d0,kind(0.d0))
                else
                   Z1 =   (cmplx(1.d0,0.d0,kind(0.d0)) - cmplx(2.d0,0.d0,kind(0.d0))*GRC(I,I,1)) *  &
                        & (cmplx(1.d0,0.d0,kind(0.d0)) - cmplx(2.d0,0.d0,kind(0.d0))*GRC(J,J,1)) +  &
                        &  cmplx(4.d0,0.d0,kind(0.d0)) * GRC(I,J,1)*GR(I,J,1)
                   Z1 = Z1**(N_SUN)
                   Z1 = Z1 * cmplx(tau_x_c(I,J,ntau,Isigma, Isigmap1) * star_sigma_x_c(i,j,ntau) ,0.d0,kind(0.d0))
                endif
                Obs_eq(5)%Obs_Latt(imj,1,1,1) =  Obs_eq(5)%Obs_Latt(imj,1,1,1) + Z1*ZP*ZS
             Enddo
             Z1 = cmplx(1.d0,0.d0,kind(0.d0)) - cmplx(2.d0,0.d0,kind(0.d0))*GRC(I,I,1)
             Z1 = Z1**(N_SUN)
             Z1 = Z1 * cmplx(tau_x(I,ntau, Isigma, Isigmap1)*star_sigma_x(i,ntau),0.d0,kind(0.d0))
             Obs_eq(5)%Obs_Latt0(1)  = Obs_eq(5)%Obs_Latt0(1)  + Z1*ZP*ZS
          Enddo

          If (Abs(Ham_T) > Zero )  Deallocate ( Isigma, Isigmap1)
 
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
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE,Mc_step_weight)
          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS
          Integer :: IMJ, I, J, I1, J1, no_I, no_J, NT1

          NT1 = NT
          If (NT == 0 ) NT1 = LTROT
          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

          ZS = ZS*Mc_step_weight

!!$          If (NT == 0 ) then
!!$             DO I = 1,Size(Obs_tau,1)
!!$                Obs_tau(I)%N = Obs_tau(I)%N + 1
!!$                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
!!$             ENDDO
!!$          endif

          Call Predefined_Obs_tau_Green_measure  ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(1) )
          Call Predefined_Obs_tau_SpinSUN_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(2) )
          Call Predefined_Obs_tau_Den_measure    ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(3) )
          
        end Subroutine OBSERT

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Returns the flux on a plaquette. I is the left-bottom corner.
!>
!--------------------------------------------------------------------
      Integer Function  iFlux(I,nt,nb_type)

        Implicit none

        Integer, INTENT(IN) :: I,nt, nb_type

        ! Local
        Integer :: n1,n2,n3,n4

        !   I3  I2
        !   I   I1
        n1  = Field_list(I,1,nb_type)
        n2  = Field_list(Latt%nnlist(I,1,0),2,nb_type)
        n3  = Field_list(Latt%nnlist(I,0,1),1,nb_type)
        n4  = Field_list(I,2,nb_type)
        iFlux =   nsigma%i(n1,nt)*nsigma%i(n2,nt)*nsigma%i(n3,nt)*nsigma%i(n4,nt)

      end Function iFlux

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> The user can set the initial field.
!>
!> @details
!> @param[OUT] Initial_field Real(:,:)
!> \verbatim
!>  Upon entry Initial_field is not allocated. If alloacted then it will contain the
!>  the initial field
!> \endverbatim
!--------------------------------------------------------------------
      Subroutine  Hamiltonian_set_nsigma(Initial_field)

        ! The user can set the initial configuration

        Implicit none

        Complex (Kind=Kind(0.d0)), allocatable, dimension(:,:), Intent(INOUT) :: Initial_field

        ! Local
        Integer :: I,nc, I1, nt, n_orientation, N_ops
        Integer, allocatable::  Isigma(:), Isigma1(:)


        N_ops = size(Field_list_inv,1)

        Allocate  (Initial_field(N_ops, Ltrot) )
        allocate  (Isigma(Latt%N), Isigma1(Latt%N) )

        Initial_field = 0.d0
        If ( Abs(Ham_U) > Zero ) then
           do nt = 1,Ltrot
              do I = 1,Latt%N
                 nc = Field_list(I,3,3)
                 Initial_field(nc,nt) = cmplx(1.D0, 0.d0,Kind(0.d0))
                 if ( ranf_wrap()  > 0.5D0 ) Initial_field(nc,nt)  = cmplx(-1.D0,0.d0,Kind(0.d0))
              enddo
           enddo
        endif
        If ( Abs(Ham_TZ2) > Zero ) then
           !  Start with a pi-flux state.
           Do nt = 1,Ltrot
              Do I = 1, Latt%N
                 if (mod( Latt%list(i,1) + latt%list(i,2), 2 ) == 0 ) then
                    Initial_field(Field_list(I,1,1),nt) =  cmplx( 1.d0, 0.d0, Kind(0.d0))
                    Initial_field(Field_list(I,2,1),nt) =  cmplx(-1.d0, 0.d0, Kind(0.d0))
                 else
                    Initial_field(Field_list(I,1,1),nt) =  cmplx(1.d0, 0.d0, Kind(0.d0))
                    Initial_field(Field_list(I,2,1),nt) =  cmplx(1.d0, 0.d0, Kind(0.d0))
                 endif
              Enddo
           Enddo
        endif
        If ( Abs(Ham_T) > Zero ) then
           Do nt = 1,Ltrot
              Do I = 1,Latt%N
                 Isigma(I) = 1
                 if ( ranf_wrap()  > 0.5D0 ) Isigma(I)  = -1
              enddo
              Do I = 1,Latt%N
                 Do n_orientation = 1,2
                    nc = Field_list(I,n_orientation,2)
                    if (  n_orientation == 1 )  I1 = latt%nnlist(I,1,0)
                    if (  n_orientation == 2 )  I1 = latt%nnlist(I,0,1)
                    Initial_field(nc,nt) = cmplx(real(Isigma(I)*Isigma(I1), kind(0.d0)),0.d0,Kind(0.d0))
                 enddo
              Enddo
              Initial_field(Field_list(Latt%N,3,4),nt) = cmplx(real(Isigma(Latt%N), kind(0.d0)),0.d0,Kind(0.d0))
              do nc = 1,size(Initial_field,1)
                 nsigma%f(nc,nt) = Initial_field(nc,nt)
              enddo
              Call Hamiltonian_set_Z2_matter(Isigma1,nt)
              Do nc = 1,Latt%N
                 if ( Isigma(nc) .ne.  Isigma1(nc)  ) then
                    Write(error_unit,*) 'Error in Hamiltonian_set_Z2_matter'
                    CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
                 endif
              enddo
           enddo
        endif

        deallocate (Isigma, Isigma1)



      end Subroutine Hamiltonian_set_nsigma

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Given the the HS fields nsigma  (mu^{z}_{i,j}, tau^z_{i=Latt%N}) the routine computes
!> the site matter fields tau^{z}_i
!>
!> @details
!--------------------------------------------------------------------
      Subroutine  Hamiltonian_set_Z2_matter(Isigma,nt)

        ! On input :  Link variables  nsigma(:,nt)
        ! On output:  The Z2_matter fields Isigma on the time slice.

        Implicit none

        Integer, Intent(IN)                  :: nt
        Integer, allocatable, INTENT(INOUT)  :: Isigma(:)

        !Local
        Integer :: I, I1, nx, ny

        Isigma(Latt%N) = nsigma%i( Field_list(Latt%N,3,4), nt )
        I = Latt%N
        do nx = 1,L1
           do ny = 1,L2
              I1 = latt%nnlist(I,0,1)
              Isigma(I1)  = Isigma(I)*nsigma%i(Field_list(I,2,2),nt)
              !Write(6,*) Latt%list(I,1), Latt%list(I,2), ' -> ', Latt%list(I1,1), Latt%list(I1,2)
              I = I1
           enddo
           I1          = latt%nnlist(I,1,0)
           Isigma(I1)  = Isigma(I)*nsigma%i(Field_list(I,1,2),nt)
           !Write(6,*) Latt%list(I,1), Latt%list(I,2), ' -> ', Latt%list(I1,1), Latt%list(I1,2)
           I = I1
        enddo

      end Subroutine Hamiltonian_set_Z2_matter

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> This routine allows to user to  determine the global_tau sampling parameters at run time
!> It is especially usefull if these parameters are dependent on other parameters.
!>
!> @details
!> \endverbatim
!--------------------------------------------------------------------
      Subroutine Overide_global_tau_sampling_parameters(Nt_sequential_start,Nt_sequential_end,N_Global_tau)

        Implicit none
        Integer, Intent(INOUT) :: Nt_sequential_start,Nt_sequential_end, N_Global_tau


        Nt_sequential_start = 1
        Nt_sequential_end   = 0
        If (abs(Ham_U  ) > Zero ) Nt_sequential_end = Nt_sequential_end + Latt%N
        If (abs(Ham_TZ2) > Zero ) Nt_sequential_end = Nt_sequential_end + Latt%N*Latt_unit%N_coord
        N_Global_tau = 0
        if (abs(Ham_T) > Zero )  N_Global_tau        = Latt%N/4

      end Subroutine Overide_global_tau_sampling_parameters


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> This function returns the expectation value of the operator tau_x_i on time slice nt
!>
!>
!> @details
!> \endverbatim
!--------------------------------------------------------------------
      Real (Kind=Kind(0.d0)) function tau_x(i,nt, Isigma, Isigmap1)
        
        Implicit none

        Integer, Intent(IN) ::  i, nt
        Integer, Intent(IN) :: Isigma(:), Isigmap1(:)
        
        Real ( Kind =Kind(0.d0) ) :: X
        Integer :: I3, I4, nt1

        tau_x = 1.d0
        If  (Abs(Ham_T) > Zero ) then
           X     =   DW_Matter_tau( Isigma(I)*Isigmap1(I) )
           If  (Abs(Ham_TZ2) > Zero )  then
              !      I2
              !  I3  I  I1
              !      I4
              I3 = Latt%nnlist(I,-1, 0)
              I4 = Latt%nnlist(I, 0,-1)
              X  =   X *  DW_Ising_Matter( nsigma%i(Field_list(I ,1,2),nt) * nsigma%i(Field_list(I ,1,1),nt) ) * &
                   &      DW_Ising_Matter( nsigma%i(Field_list(I ,2,2),nt) * nsigma%i(Field_list(I ,2,1),nt) ) * &
                   &      DW_Ising_Matter( nsigma%i(Field_list(I3,1,2),nt) * nsigma%i(Field_list(I3,1,1),nt) ) * &
                   &      DW_Ising_Matter( nsigma%i(Field_list(I4,2,2),nt) * nsigma%i(Field_list(I4,2,1),nt) )
           Endif
           tau_x  = X
        endif
      end function tau_x

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> This function returns the expectation value of the operator <tau_x_i tau_x_j> on time slice nt
!>
!>
!> @details
!> \endverbatim
!--------------------------------------------------------------------
      Real (Kind=Kind(0.d0)) function tau_x_c(i,j,nt, Isigma, Isigmap1)
        
        Implicit none

        Integer, Intent(IN) ::  i,j, nt
        Integer, Intent(IN) :: Isigma(:), Isigmap1(:)
        
        Real ( Kind =Kind(0.d0) ) :: X
        Integer :: I3, I4, J3, J4

        !Write(6,*) 'In tau_x_c'

        if ( i == j ) then
           tau_x_c = 1.d0
        else
           tau_x_c = 1.d0
           X       = 1.d0
           If  (Abs(Ham_T) > Zero ) then
              X     =  DW_Matter_tau( Isigma  (I)*Isigmap1(I)) * DW_Matter_tau( Isigma  (J)*Isigmap1(J))
              If  (Abs(Ham_TZ2) > Zero )  then
                 I3 = Latt%nnlist(I,-1, 0)
                 I4 = Latt%nnlist(I, 0,-1)
                 J3 = Latt%nnlist(J,-1, 0)
                 J4 = Latt%nnlist(J, 0,-1)
                 If (J == Latt%nnlist(I,-1,0) ) then
                    !   I - J  = a_1
                    !
                    !       J2  I2
                    !   J3  J   I  I1
                    !       J4  I4
                    !
                    X  =   X *  DW_Ising_Matter( nsigma%i(Field_list(I ,1,2),nt) * nsigma%i(Field_list(I ,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(I ,2,2),nt) * nsigma%i(Field_list(I ,2,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(I4,2,2),nt) * nsigma%i(Field_list(I4,2,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J ,2,2),nt) * nsigma%i(Field_list(J ,2,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J3,1,2),nt) * nsigma%i(Field_list(J3,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J4,2,2),nt) * nsigma%i(Field_list(J4,2,1),nt) )
                 elseif (J == Latt%nnlist(I,1,0) ) then
                    !   I - J  = - a_1
                    !
                    !       I2  J2
                    !   I3  I   J  J1
                    !       I4  J4
                    !
                    X  =   X *  DW_Ising_Matter( nsigma%i(Field_list(I ,2,2),nt) * nsigma%i(Field_list(I ,2,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(I3,1,2),nt) * nsigma%i(Field_list(I3,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(I4,2,2),nt) * nsigma%i(Field_list(I4,2,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J ,1,2),nt) * nsigma%i(Field_list(J ,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J ,2,2),nt) * nsigma%i(Field_list(J ,2,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J4,2,2),nt) * nsigma%i(Field_list(J4,2,1),nt) )
                 elseif (J == Latt%nnlist(I,0,-1) ) then
                    !   I - J  =  a_2
                    !
                    !           I2
                    !       I3  I  I1
                    !       J3  J  J1
                    !           J4
                    !
                    X  =   X *  DW_Ising_Matter( nsigma%i(Field_list(I ,1,2),nt) * nsigma%i(Field_list(I ,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(I ,2,2),nt) * nsigma%i(Field_list(I ,2,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(I3,1,2),nt) * nsigma%i(Field_list(I3,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J ,1,2),nt) * nsigma%i(Field_list(J ,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J3,1,2),nt) * nsigma%i(Field_list(J3,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J4,2,2),nt) * nsigma%i(Field_list(J4,2,1),nt) )
                 elseif (J ==  Latt%nnlist(I,0,1) ) then
                    !   I - J  =  -a_2
                    !           J2
                    !       J3  J  J1
                    !       I3  I  I1
                    !           I4
                    !
                    X  =   X *  DW_Ising_Matter( nsigma%i(Field_list(I ,1,2),nt) * nsigma%i(Field_list(I ,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(I3,1,2),nt) * nsigma%i(Field_list(I3,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(I4,2,2),nt) * nsigma%i(Field_list(I4,2,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J ,1,2),nt) * nsigma%i(Field_list(J ,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J ,2,2),nt) * nsigma%i(Field_list(J ,2,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J4,2,2),nt) * nsigma%i(Field_list(J4,2,1),nt) )
                 else
                    X  =   X *  DW_Ising_Matter( nsigma%i(Field_list(I ,1,2),nt) * nsigma%i(Field_list(I ,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(I ,2,2),nt) * nsigma%i(Field_list(I ,2,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(I3,1,2),nt) * nsigma%i(Field_list(I3,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(I4,2,2),nt) * nsigma%i(Field_list(I4,2,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J ,1,2),nt) * nsigma%i(Field_list(J ,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J ,2,2),nt) * nsigma%i(Field_list(J ,2,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J3,1,2),nt) * nsigma%i(Field_list(J3,1,1),nt) ) * &
                         &      DW_Ising_Matter( nsigma%i(Field_list(J4,2,2),nt) * nsigma%i(Field_list(J4,2,1),nt) )
                 Endif
              endif
              tau_x_c  = X
           endif
        endif
        !Write(6,*) 'Out tau_x_c'
      end function tau_x_c

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> This function returns the expectation value of the operator sigma_x_(i,i + a_n_orientation)  on time slice nt
!>
!>
!> @details
!> \endverbatim
!--------------------------------------------------------------------
      Real (Kind=Kind(0.d0)) function sigma_x(i,n_orientation,nt)
        
        Implicit none

        Integer, Intent(IN) ::  i, nt, n_orientation

        Real ( Kind =Kind(0.d0) ) :: X
        Integer :: F1, F2, nt1, nt2


        sigma_x = 1.d0
        X     = 1.d0
        nt1   = nt + 1
        If (nt == Ltrot )  nt1 = 1
        nt2   = nt + 2
        If (nt2 >  Ltrot )  nt2 = nt2 - Ltrot
        !Write(6,*) 'In sigma_x', nt1, nt2, I, n_orientation
        If  (Abs(Ham_TZ2) > Zero ) then
           X = X *  DW_Ising_tau( nsigma%i(Field_list(I ,n_orientation,1),nt )*nsigma%i(Field_list(I ,n_orientation,1),nt1) )
           if ( n_orientation == 1 ) then
              F1 = iFlux(i                  , nt,1)
              F2 = iFlux(latt%nnlist(i,0,-1), nt,1)
           else
              F1 = iFlux(i                  , nt,1)
              F2 = iFlux(latt%nnlist(i,-1,0), nt,1)
           endif
           X  = X * DW_Ising_Flux(F1,F2)
           If (Abs(Ham_T) > Zero )  then
              X = X * DW_Ising_Matter( nsigma%i(Field_list(I ,n_orientation,2),nt) * nsigma%i(Field_list(I ,n_orientation,1),nt) )
           endif
           sigma_x  = X
        endif
        !Write(6,*) 'Out sigma_x'
        
      end function sigma_x

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> This function returns the expectation value of the operator <star_sigma_x(i) >
!> where star_sigma_x(j) = sigma^x(i,i+a_x) sigma^x(i,i-a_x) sigma^x(i,i+a_y) sigma^x(i,i-a_y)
!>
!> @details
!> \endverbatim
!--------------------------------------------------------------------
      Real (Kind=Kind(0.d0)) function star_sigma_x(i,nt)
        
        Implicit none

        Integer, Intent(IN) ::  i, nt

        Real ( Kind =Kind(0.d0) ) :: X
        Integer :: I3, I4, ntp1

        !Write(6,*) 'In star_sigma_x'
        star_sigma_x = 1.d0
        If  (Abs(Ham_TZ2) > Zero ) then
           X     = 1.d0
           ntp1   = nt + 1
           If ( nt == Ltrot )  ntp1 = 1
           !         I2
           !      I3  I  I1
           !         I4
           !
           I3 = Latt%nnlist(I,-1, 0)
           I4 = Latt%nnlist(I, 0,-1)
           X =     DW_Ising_tau  ( nsigma%i(Field_list(I ,1,1),nt)*nsigma%i(Field_list(I ,1,1),ntp1)  ) * &
                &  DW_Ising_tau  ( nsigma%i(Field_list(I ,2,1),nt)*nsigma%i(Field_list(I ,2,1),ntp1)  ) * &
                &  DW_Ising_tau  ( nsigma%i(Field_list(I4,2,1),nt)*nsigma%i(Field_list(I4,2,1),ntp1)  ) * &
                &  DW_Ising_tau  ( nsigma%i(Field_list(I3,1,1),nt)*nsigma%i(Field_list(I3,1,1),ntp1)  )
           ! Flux remains invariant.
           If (Abs(Ham_T) > Zero )  then
              X = X * DW_Ising_Matter( nsigma%i(Field_list(I ,1,2),nt) * nsigma%i(Field_list(I ,1,1),nt) ) * &
                      DW_Ising_Matter( nsigma%i(Field_list(I ,2,2),nt) * nsigma%i(Field_list(I ,2,1),nt) ) * &
                      DW_Ising_Matter( nsigma%i(Field_list(I3,1,2),nt) * nsigma%i(Field_list(I3,1,1),nt) ) * &
                      DW_Ising_Matter( nsigma%i(Field_list(I4,2,2),nt) * nsigma%i(Field_list(I4,2,1),nt) )
           endif
           star_sigma_x  = x
        endif
        !Write(6,*) 'Out star_sigma_x'
        
      end function star_sigma_x

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> This function returns the expectation value of the operator <star_sigma_x(i) * star_sigma_x(j) >
!> where star_sigma_x(j) = sigma^x(i,i+a_x) sigma^x(i,i-a_x) sigma^x(i,i+a_y) sigma^x(i,i-a_y)
!>
!> @details
!> \endverbatim
!--------------------------------------------------------------------
      Real (Kind=Kind(0.d0)) function star_sigma_x_c(i,j,nt)
        
        Implicit none

        Integer, Intent(IN) ::  i, j, nt

        Integer :: I3, I4, J3,J4,  nt1
        Real (Kind=Kind(0.d0)) :: X

        nt1 = nt + 1
        if ( nt == Ltrot ) nt1 = 1

        if ( I == J ) then
           star_sigma_x_c = 1.d0
        elseif ( Abs(Ham_TZ2) < Zero ) then
           star_sigma_x_c = 1.d0
        else
           !      I2
           !  I3  I  I1
           !      I4
           !      J2
           !  J3  J  J1
           !      J4
           I3 = Latt%nnlist(I,-1, 0)
           I4 = Latt%nnlist(I, 0,-1)
           J3 = Latt%nnlist(J,-1, 0)
           J4 = Latt%nnlist(J, 0,-1)
           If (J == Latt%nnlist(I,-1,0) ) then
              !
              !       J2  I2
              !   J3  J   I  I1
              !       J4  I4
              !
              X       = DW_Ising_tau( nsigma%i(Field_list(I ,1,1),nt)*nsigma%i(Field_list(I ,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(I ,2,1),nt)*nsigma%i(Field_list(I ,2,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(I4,2,1),nt)*nsigma%i(Field_list(I4,2,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J ,2,1),nt)*nsigma%i(Field_list(J ,2,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J3,1,1),nt)*nsigma%i(Field_list(J3,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J4,2,1),nt)*nsigma%i(Field_list(J4,2,1),nt1) )
           elseif (J == Latt%nnlist(I,1,0) ) then
              !
              !       I2  J2
              !   I3  I   J  J1
              !       I4  J4
              !
              X       = DW_Ising_tau( nsigma%i(Field_list(I ,2,1),nt)*nsigma%i(Field_list(I ,2,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(I3,1,1),nt)*nsigma%i(Field_list(I3,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(I4,2,1),nt)*nsigma%i(Field_list(I4,2,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J ,1,1),nt)*nsigma%i(Field_list(J ,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J ,2,1),nt)*nsigma%i(Field_list(J ,2,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J4,2,1),nt)*nsigma%i(Field_list(J4,2,1),nt1) )
           elseif (J == Latt%nnlist(I,0,-1) ) then
              !
              !           I3
              !       I2  I  I1
              !       J3  J  J1
              !           J4
              !
              X   =     DW_Ising_tau( nsigma%i(Field_list(I ,1,1),nt)*nsigma%i(Field_list(I ,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(I ,2,1),nt)*nsigma%i(Field_list(I ,2,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(I3,1,1),nt)*nsigma%i(Field_list(I3,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J ,1,1),nt)*nsigma%i(Field_list(J ,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J3,1,1),nt)*nsigma%i(Field_list(J3,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J4,2,1),nt)*nsigma%i(Field_list(J4,2,1),nt1) )
           elseif (J ==  Latt%nnlist(I,0,1) ) then
              !
              !           J3
              !       J2  J  J1
              !       I3  I  I1
              !           I4
              !
              X   =     DW_Ising_tau( nsigma%i(Field_list(I ,1,1),nt)*nsigma%i(Field_list(I ,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(I3,1,1),nt)*nsigma%i(Field_list(I3,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(I4,2,1),nt)*nsigma%i(Field_list(I4,2,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J ,1,1),nt)*nsigma%i(Field_list(J ,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J ,2,1),nt)*nsigma%i(Field_list(J ,2,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J3,1,1),nt)*nsigma%i(Field_list(J3,1,1),nt1) )
           else
              X   =     DW_Ising_tau( nsigma%i(Field_list(I ,1,1),nt)*nsigma%i(Field_list(I ,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(I ,2,1),nt)*nsigma%i(Field_list(I ,2,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(I3,1,1),nt)*nsigma%i(Field_list(I3,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(I4,2,1),nt)*nsigma%i(Field_list(I4,2,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J ,1,1),nt)*nsigma%i(Field_list(J ,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J ,2,1),nt)*nsigma%i(Field_list(J ,2,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J3,1,1),nt)*nsigma%i(Field_list(J3,1,1),nt1) )* &
                   &    DW_Ising_tau( nsigma%i(Field_list(J4,2,1),nt)*nsigma%i(Field_list(J4,2,1),nt1) )
           endif
           If  (Abs(Ham_T) > Zero )  then
              If (J == Latt%nnlist(I,-1,0) ) then
                 !   I - J  = a_1
                 !
                 !       J2  I2
                 !   J3  J   I  I1
                 !       J4  I4
                 !
                 X  =   X *  DW_Ising_Matter( nsigma%i(Field_list(I ,1,2),nt) * nsigma%i(Field_list(I ,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(I ,2,2),nt) * nsigma%i(Field_list(I ,2,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(I4,2,2),nt) * nsigma%i(Field_list(I4,2,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J ,2,2),nt) * nsigma%i(Field_list(J ,2,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J3,1,2),nt) * nsigma%i(Field_list(J3,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J4,2,2),nt) * nsigma%i(Field_list(J4,2,1),nt) )
              elseif (J == Latt%nnlist(I,1,0) ) then
                 !   I - J  = - a_1
                 !
                 !       I2  J2
                 !   I3  I   J  J1
                 !       I4  J4
                 !
                 X  =   X *  DW_Ising_Matter( nsigma%i(Field_list(I ,2,2),nt) * nsigma%i(Field_list(I ,2,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(I3,1,2),nt) * nsigma%i(Field_list(I3,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(I4,2,2),nt) * nsigma%i(Field_list(I4,2,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J ,1,2),nt) * nsigma%i(Field_list(J ,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J ,2,2),nt) * nsigma%i(Field_list(J ,2,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J4,2,2),nt) * nsigma%i(Field_list(J4,2,1),nt) )
              elseif (J == Latt%nnlist(I,0,-1) ) then
                 !   I - J  =  a_2
                 !
                 !           I2
                 !       I3  I  I1
                 !       J3  J  J1
                 !           J4
                 !
                 X  =   X *  DW_Ising_Matter( nsigma%i(Field_list(I ,1,2),nt) * nsigma%i(Field_list(I ,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(I ,2,2),nt) * nsigma%i(Field_list(I ,2,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(I3,1,2),nt) * nsigma%i(Field_list(I3,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J ,1,2),nt) * nsigma%i(Field_list(J ,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J3,1,2),nt) * nsigma%i(Field_list(J3,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J4,2,2),nt) * nsigma%i(Field_list(J4,2,1),nt) )
              elseif (J ==  Latt%nnlist(I,0,1) ) then
                 !   I - J  =  -a_2
                 !           J2
                 !       J3  J  J1
                 !       I3  I  I1
                 !           I4
                 !
                 X  =   X *  DW_Ising_Matter( nsigma%i(Field_list(I ,1,2),nt) * nsigma%i(Field_list(I ,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(I3,1,2),nt) * nsigma%i(Field_list(I3,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(I4,2,2),nt) * nsigma%i(Field_list(I4,2,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J ,1,2),nt) * nsigma%i(Field_list(J ,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J ,2,2),nt) * nsigma%i(Field_list(J ,2,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J4,2,2),nt) * nsigma%i(Field_list(J4,2,1),nt) )
              else
                 X  =   X *  DW_Ising_Matter( nsigma%i(Field_list(I ,1,2),nt) * nsigma%i(Field_list(I ,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(I ,2,2),nt) * nsigma%i(Field_list(I ,2,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(I3,1,2),nt) * nsigma%i(Field_list(I3,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(I4,2,2),nt) * nsigma%i(Field_list(I4,2,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J ,1,2),nt) * nsigma%i(Field_list(J ,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J ,2,2),nt) * nsigma%i(Field_list(J ,2,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J3,1,2),nt) * nsigma%i(Field_list(J3,1,1),nt) ) * &
                      &      DW_Ising_Matter( nsigma%i(Field_list(J4,2,2),nt) * nsigma%i(Field_list(J4,2,1),nt) )
              Endif
              
           endif
           star_sigma_x_c = X
        endif
        
      end function star_sigma_x_c


      end submodule ham_Z2_Matter_smod
