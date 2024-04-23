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


    submodule (Hamiltonian_main) ham_Anderson_smod

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
      
      type, extends(ham_base) :: ham_Anderson
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: E_loc
        procedure, nopass :: weight_loc
        procedure, nopass :: ObserT
        procedure, nopass :: weight_reconstruction
        procedure, nopass :: GR_reconstruction
        procedure, nopass :: GRT_reconstruction
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Anderson

      !#PARAMETERS START# VAR_lattice
      Character (len=64) :: Model = ''  ! Value not relevant
      Character (len=64) :: Lattice_type = 'Bilayer_square'  ! Possible values: 'Bilayer_square', 'Bilayer_honeycomb'
      Integer            :: L1 = 6   ! Length in direction a_1
      Integer            :: L2 = 6   ! Length in direction a_2
      !#PARAMETERS END#


      !#PARAMETERS START# VAR_Anderson
      real(Kind=Kind(0.d0)) :: ham_T    = 1.d0      !Hopping parameter
      real(Kind=Kind(0.d0)) :: Ham_chemc = 0.d0     !Chemical potential on c-orbitals
      real(Kind=Kind(0.d0)) :: Ham_chemf = 0.d0     !Chemical potential on  f-orbials
      real(Kind=Kind(0.d0)) :: Ham_Uc   = 0.d0      !Hubbard interaction  on  c-orbitals Uc
      real(Kind=Kind(0.d0)) :: ham_Uf   = 2.d0      !Hubbard interaction  on  f-orbials  Uf
      real(Kind=Kind(0.d0)) :: Ham_Vcf    = 2.d0    !Hybridasation in periodic Anderson model
      real(Kind=Kind(0.d0)) :: Dtau         = 0.1d0    ! Thereby Ltrot=Beta/dtau
      real(Kind=Kind(0.d0)) :: Beta         = 5.d0     ! Inverse temperature
      !logical              :: Symm         = .fasle.   ! Whether symmetrization takes place
      !logical              :: Projector    = .false.  ! Whether the projective algorithm is used
      logical               :: Adiabatic= .false.   ! If  true,  and projector  true then adiabatic  switching on of  U. 
      real(Kind=Kind(0.d0)) :: Theta        = 10.d0    ! Projection parameter
      Integer               :: N_part   = -1        ! Number of particles in trial wave function. If N_part < 0 -> N_part = L1*L2/2 
      !#PARAMETERS END#

      Type (Lattice),       Target  :: Latt
      Type (Unit_cell),     Target  :: Latt_unit
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:)
      Integer, allocatable   :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell
      
      Type (Unit_cell), Target  :: Latt_unit_f  ! Unit cell for f  correlation functions
      Type (Unit_cell), Target  :: Latt_unit_c  ! Unit cell for c  correlation functions
      Integer, allocatable      :: List_c(:,:)
      Integer, allocatable      :: List_f(:,:)
      Integer  ::  nf_calc,  nf_reconst  

    contains
      
      module Subroutine Ham_Alloc_Anderson
        allocate(ham_Anderson::ham)
      end Subroutine Ham_Alloc_Anderson

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Anderson_read_write_parameters.F90"

   
! NAMELIST /VAR_VAFQMC/ Del_BCS,  mc_AFM, mf_AFM, Vcf_And, chi_ff 

!real(Kind=Kind(0.d0)) :: Del_BCS    = 0.d0   !Varitional bcs  mean-field parameters for trial wave-functions ! This could be a vector  Del_f,  Del_c,  Del_fc
!real(Kind=Kind(0.d0)) :: mc_AFM    = 0.d0    !Varitional AFM order on c-orbitals  for trial wave-functions
!real(Kind=Kind(0.d0)) :: mf_AFM    = 0.d0    !Varitional  AFM order on f-orbitals for trial wave-functions
!real(Kind=Kind(0.d0)) :: Vcf_And    = 0.d0   !Varitional  hybridisation parameters for trial wave-functions
!real(Kind=Kind(0.d0)) :: chi_ff    = 0.d0   !Varitional  spinon hopping for trial wave-functions


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

          ! From dynamically generated file "Hamiltonian_Kondo_read_write_parameters.F90"
          call read_parameters()
          
          !If ( .not. ( Lattice_type == "Bilayer_square" .or.  Lattice_type == "Bilayer_honeycomb") ) then
          !   Write(error_unit,*) "The Kondo Hamiltonian is only defined for bilayer lattices"
          !   CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          !endif

          if ( (.not. projector) .and. adiabatic) then
            write(output_unit,*) "Adiabatic mode is only implemented for projective code."
            write(output_unit,*) "Overriding Adiabatic=.True. from parameter files."
          endif


          if (N_part < 0) N_part = L1*L2/2
          Ltrot = nint(beta/dtau)
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot = Ltrot+2*Thtrot

          N_SUN        = 1
          N_FL         = 2

          !IF ( N_FL > 1 ) then
          !   Write(error_unit,*) 'For the Kondo systems, N_FL has  to be equal to unity'
          !   CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          !Endif



          ! Setup the Bravais lattice
          Call  Ham_Latt

          ! Setup the hopping / single-particle part
          Call  Ham_Hop

          ! Setup the interaction.
          call Ham_V

          ! Setup the trival wave function, in case of a projector approach
          !if (Projector) Call Ham_Trial()

#ifdef MPI
          If (Irank_g == 0) then
#endif
             File_info = "info"
#if defined(TEMPERING)
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif
             Open(newunit=unit_info, file=file_info, status="unknown", position="append")
             Write(unit_info,*) '====================================='
             Write(unit_info,*) 'Model is      : Anderson'
             Write(unit_info,*) 'Lattice is    : ', Lattice_type
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
             Write(unit_info,*) 'N_SUN         : ',   N_SUN
             Write(unit_info,*) 'N_FL          : ', N_FL
             Write(unit_info,*) 't             : ', Ham_T
             Write(unit_info,*) 'Ham_Uc        : ', Ham_Uc
             Write(unit_info,*) 'Ham_Uf        : ', Ham_Uf
             Write(unit_info,*) 'Ham_Vcf       : ', Ham_Vcf
             Write(unit_info,*) 'Ham_chemc     : ', Ham_chemc
             Write(unit_info,*) 'Ham_chemf     : ', Ham_chemf

             If  (Ham_Uc >=0.d0  .and.   Ham_chemc   == 0.d0) then
               Write(unit_info,*) 'Assuming particle hole symmetry' 
            endif

            If  (Ham_Uc <  0.d0 ) then
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
         If ( (Ham_Uf >= 0.d0 .and.   Ham_chemf  == 0.d0)  .or.  Ham_Uf < 0.d0  )    then
            allocate(Calc_Fl(N_FL))
            nf_calc=2
            nf_reconst=1
            Calc_Fl(nf_calc)=.True.
            Calc_Fl(nf_reconst)=.False.
         endif

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
          Integer :: n, nc,I, no
          
          ! Use predefined stuctures or set your own lattice.
          
          Call Predefined_Latt(Lattice_type, L1,L2,Ndim, List,Invlist,Latt,Latt_Unit)

          !  Setup lattices for f-and c-sites.
          Select case (Lattice_type)
          Case ("Bilayer_square")
             Latt_Unit_f%Norb       = 1
             Latt_Unit_f%N_coord    = 2
             Allocate (Latt_Unit_f%Orb_pos_p(1,3))
             Latt_Unit_f%Orb_pos_p(1,:) =  0.d0
             Latt_Unit_f%Orb_pos_p(1,3) = -1.d0

             Latt_Unit_c%Norb       = 1
             Latt_Unit_c%N_coord    = 2
             Allocate (Latt_Unit_c%Orb_pos_p(1,3))
             Latt_Unit_c%Orb_pos_p(1,:) =  0.d0
             Latt_Unit_c%Orb_pos_p(1,3) =  0.d0

          Case ("Bilayer_honeycomb")
             Latt_Unit_f%Norb    = 2
             Latt_Unit_f%N_coord = 3
             Allocate (Latt_Unit_f%Orb_pos_p(2,3))
             Latt_unit_f%Orb_pos_p = -1.d0
             do n = 1,2
                Latt_Unit_f%Orb_pos_p(1,n) = 0.d0
                Latt_Unit_f%Orb_pos_p(2,n) = (Latt%a2_p(n) - 0.5D0*Latt%a1_p(n) ) * 2.D0/3.D0
             Enddo

             Latt_Unit_c%Norb    = 2
             Latt_Unit_c%N_coord = 3
             Allocate (Latt_Unit_c%Orb_pos_p(2,3))
             Latt_unit_c%Orb_pos_p = 0.d0
             do n = 1,2
                Latt_Unit_c%Orb_pos_p(1,n) = 0.d0
                Latt_Unit_c%Orb_pos_p(2,n) = (Latt%a2_p(n) - 0.5D0*Latt%a1_p(n) ) * 2.D0/3.D0
             Enddo

          end Select
          
          Allocate (List_f(Ndim,2))  !  For measuring only on  f-lattice
          List_f = 0
          Do I = 1,Latt%N
             Do no = 1,Latt_Unit_f%Norb
                list_f(Invlist(I,no + Latt_Unit%Norb/2),1)  =  I
                list_f(Invlist(I,no + Latt_Unit%Norb/2),2)  =  no
             Enddo
          Enddo
          
          Allocate (List_c(Ndim,2) ) !  For measuring only on  c-lattice
          List_c = 0 
          nc = 0
          Do I = 1,Latt%N
             Do no = 1,Latt_Unit_c%Norb
                List_c(Invlist(I,no ),1) = I
                List_c(Invlist(I,no ),2) = no
             Enddo
          Enddo
          
        end Subroutine Ham_Latt
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the Hopping
!--------------------------------------------------------------------
        Subroutine Ham_Trial()

         Use Predefined_Trial

         Implicit none

         Integer :: N_part, nf
         ! Use predefined stuctures or set your own Trial  wave function
         N_part = Ndim/2
         Call Predefined_TrialWaveFunction(Lattice_type ,Ndim,  List,Invlist,Latt, Latt_unit, &
              &                            N_part, N_FL,  WF_L, WF_R)

       end Subroutine Ham_Trial
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Copy and  set Operators
!> @details
!--------------------------------------------------------------------

       Subroutine Ham_Hop
         Implicit none
         ! Use predefined structures or set your own hopping
         Real (Kind=Kind(0.d0) ) ::  Ham_Lambda = 0.d0
         Integer :: I, Ix, Iy, nt, Nc_ops, Nf_ops
         Integer :: nth, I1, I2, no, n
         Integer :: nc_f, nc_tot, nc, nf
         Logical :: Checkerboard ,Bulk

         Real (Kind=Kind(0.d0) ), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
         &                                  Ham_T2_vec(:),  Ham_Lambda_vec(:)
         Integer, allocatable ::   N_Phi_vec(:)


         Real (Kind=Kind(0.d0)) :: Factor
         Real (Kind=Kind(0.d0)), allocatable :: hi(:), ti(:), lambdai(:)
         Type (Operator), dimension(:,:), allocatable :: Op_T_c, OP_V_cf
         ! We need hopping along conduction electrons + Non-interacting periodic Anderson term
         ! Conduction electron hopping 
         !call time_dependent_vmfp(hi, ti, lambdai)
       


         Allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
         &                                   N_Phi_vec(N_FL), Ham_Lambda_vec(N_FL) )

         Ham_T_vec      = Ham_T
         Ham_Tperp_vec  = 0.d0
         Ham_Chem_vec   = Ham_Chemc
         Phi_X_vec      = 0.d0 
         Phi_Y_vec      = 0.d0 
         Ham_T2_vec     = 0.d0
         Ham_Lambda_vec = Ham_Lambda
         N_Phi_vec      = 0
         Checkerboard = .false. 
         Bulk = .true.

        ! Select case (Lattice_type)
         !Case ("Bilayer_square")
         !   Call  Set_Default_hopping_parameters_Bilayer_square(Hopping_Matrix,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
         !      &                                              Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
         !      &                                              List, Invlist, Latt, Latt_unit )

         !Case ("Bilayer_honeycomb")
          !  Call  Set_Default_hopping_parameters_Bilayer_honeycomb(Hopping_Matrix,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
          !     &                                                 Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
           !    &                                                 List, Invlist, Latt, Latt_unit )

         !end Select

         
        !Call  Predefined_Hoppings_set_OPT(Hopping_Matrix,List,Invlist,Latt,  Latt_unit,  Dtau, Checkerboard, Symm, OP_T_c )

         !Deallocate (Ham_T_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
          !    &                                   N_Phi_vec,  Ham_Lambda_vec )



         !Ndim  = Latt%N * Latt_Unit%Norb
         !Nc_ops  = Ndim 
         Nc_ops  = Latt%N 
         allocate(Op_T_c(1, N_FL))
         do nf = 1, N_FL
             Call Op_make(Op_T_c(1, nf), Nc_ops)

             Do I = 1,  1, Latt%N!Nc_ops
                 Ix = Latt%nnlist(I, 1, 0)
                 Op_T_c(1, nf)%O(I, Ix) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                 Op_T_c(1, nf)%O(Ix, I) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                 Iy = Latt%nnlist(I, 0, 1)
                 Op_T_c(1, nf)%O(I, Iy) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                 Op_T_c(1, nf)%O(Iy, I) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                 
                 Op_T_c(1, nf)%O(I, I) = cmplx(-Ham_chemc, 0.d0, kind(0.D0))
                 Op_T_c(1, nf)%P(I) = I
             Enddo
             If (Adiabatic) then
                 Allocate(Op_T_c(1, nf)%g_t(Ltrot))
                 Op_T_c(1, nf)%g_t = -dtau
                 do nt = 1, Thtrot
                     !Op_T_c(1, nf)%g_t(nt) = - ti(nt)
                     Op_T_c(1, nf)%g_t(nt) = -dtau * dble(nt) / dble(thtrot)
                     Op_T_c(1, nf)%g_t(Ltrot - (nt - 1)) = -dtau * dble(nt) / dble(thtrot)
                 enddo
             else
                 Op_T_c(1, nf)%g = -dtau
             endif
             Op_T_c(1, nf)%alpha = cmplx(0.d0, 0.d0, kind(0.D0))
             Call Op_set(Op_T_c(1, nf))
         enddo
       

         ! Non-interacting periodic Anderson term
         Factor = 1.d0
        ! Nf_ops = Latt%N * Latt_Unit%Norb
         Nf_ops = Latt%N 
         Allocate(Op_V_cf(Nf_ops, N_FL))

         do nf = 1, N_FL
             do nc = 1, Nf_ops
         !        !do nc = Nf_ops+1,  2*Nf_ops  
                 Call Op_make(Op_V_cf(nc, nf), 2)
             enddo
       
             nc = 0!Latt%N * Latt_Unit%Norb

            Do I = 1, Latt%N
            ! Do I = 1,  Nf_ops
                 Do no = 1, Latt_unit%Norb / 2
                     I1 = Invlist(I, no)
                     I2 = Invlist(I, no + Latt_unit%Norb / 2)
                     nc = nc + 1
                     Op_V_cf(nc, nf)%O(1, 2) = cmplx(-Ham_Vcf, 0.d0, kind(0.d0))
                     Op_V_cf(nc, nf)%O(2, 1) = cmplx(-Ham_Vcf, 0.d0, kind(0.d0))
                     Op_V_cf(nc, nf)%O(1, 1) = cmplx(-Ham_chemf, 0.d0, kind(0.d0))
                     Op_V_cf(nc, nf)%P(1) = I1
                     Op_V_cf(nc, nf)%P(2) = I2
                     !Op_V_cf(nc, nf)%g = -Dtau * Factor
       
                     If (Adiabatic) then
                         Allocate(Op_V_cf(nc, nf)%g_t(Ltrot))
                         Op_V_cf(nc, nf)%g_t = -dtau
                         do nt = 1, Thtrot
                             !Op_V_cf(nc, nf)%g_t(nt) = - ti(nt)
                             Op_V_cf(nc, nf)%g_t(nt) = -dtau * dble(nt) / dble(thtrot)
                             Op_V_cf(nc, nf)%g_t(Ltrot - (nt - 1)) = -dtau * dble(nt) / dble(thtrot)
                         enddo
                     else
                         Op_V_cf(nc, nf)%g = -Dtau * Factor
                     endif
                 enddo
             enddo
         enddo

       
  
         nc_tot = size(OP_V_cf, 1) + size(OP_T_c, 1)
         Allocate(OP_T(nc_tot, N_FL))

         print*, 'OP_V_cf, OP_T_c, OP_T , nc_tot=',size(OP_V_cf, 1),  size(OP_T_c, 1),  size(OP_T, 1), size(OP_V_cf, 2),  size(OP_T_c, 2),  size(OP_T, 2),nc_tot



         do nf = 1, N_FL
             nc = 0
             do n = 1, size(OP_V_cf, 1)
                 nc = nc + 1
                 Call Op_make(Op_T(nc, nf), 2)
                 Call Copy_Op(OP_T(nc, nf), OP_V_cf(n, nf))
             enddo
             do n = 1, size(OP_T_c, 1)
                 nc = nc + 1
                 Call Op_make(Op_T(nc, nf), 2)
                 Call Copy_Op(OP_T(nc, nf), OP_T_c(n, nf))
             enddo
         enddo
       
         ! Clear memory
         do nf = 1, N_FL
             Do n = 1, size(OP_T_c, 1)
                 Call Op_clear(OP_T_c(n, nf), 2)
             enddo
             Do n = 1, size(OP_V_cf, 1)
                 Call Op_clear(OP_V_cf(n, nf), 2)
             enddo
         enddo

         deallocate(OP_T_c, OP_V_cf)
       end Subroutine Ham_Hop
       
!> @author
!> ALF Collaboration
!>
!> @brief
!> Copy and  set Operators
!> @details
!--------------------------------------------------------------------
         Subroutine Copy_Op(OP1,OP2)
            Implicit none
            Type (Operator),     Intent(INOUT) ::  OP1 
            Type (Operator),     Intent(IN)    ::  OP2
            OP1%P     =  OP2%P
            OP1%O     =  OP2%O
            OP1%g     =  OP2%g
            OP1%alpha =  OP2%alpha
            Call Op_set(Op1)
          end Subroutine Copy_Op
   !--------------------------------------------------------------------
     
   
   
!> @author
!> ALF Collaboration
!>
!> @brief
!> Copy and  set Operators
!> @details
!--------------------------------------------------------------------
          Subroutine Ham_V
            Implicit none
            Integer :: nf, I, I1, I2, nc, no, N_ops, nt
            Real (Kind=Kind(0.d0)) :: X, Zero=1.D-10
            Real (Kind=Kind(0.d0)), allocatable :: hi(:), ti(:), lambdai(:)
            ! setting dtau as a varitional parameter for VAFQMC 
            ! call time_dependent_vmfp(hi, ti, lambdai)
          
            N_ops = 0
            if (abs(Ham_Uc)  > Zero) N_ops = N_ops + Latt%N*Latt_Unit%Norb/2
            if (abs(Ham_Uf) > Zero)  N_ops = N_ops + Latt%N*Latt_Unit%Norb/2
        
            print*, 'N_ops ,  Latt%N*Latt_Unit%Norb/2,  Latt%N , Latt_Unit%Norb, N_FL =', &
                    N_ops ,  Latt%N*Latt_Unit%Norb/2, Latt%N, Latt_Unit%Norb, N_FL
            


            Ndim  = Latt%N * Latt_Unit%Norb
            N_ops  = Ndim 

            Allocate(Op_V(N_ops,N_FL))
            do nf = 1,N_FL
                do i = 1, N_ops
                    Call Op_make(Op_V(i,nf), 1)
                end do
            end do
        
            nc = 0
            if (abs(Ham_Uc) > Zero) then
                do nf = 1,N_FL
                    X = 1.d0
                    if (nf == 2) X = -1.d0
        
                    Do I = 1,Latt%N
                        do no = 1, Latt_unit%Norb/2
                            I1 = invlist(I,no)
                            nc = nc + 1
                            Op_V(nc,nf)%P(1) = I1
                            Op_V(nc,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
                            If (Adiabatic) then
                                Allocate(OP_V(nc,nf)%g_t(Ltrot))
                                ! call time_dependent_vmfp(hi, ti, lambdai)
                                Op_V(nc,nf)%g_t = X*SQRT(CMPLX(dtau*ham_Uc/2.d0, 0.D0, kind(0.D0)))
                                do nt = 1, Thtrot
                                    ! Op_V(nc,nf)%g_t(nt) = X*SQRT(hi(nt)*Ham_Uf/2.d0, 0.D0, kind(0.D0)))
                                    Op_V(nc,nf)%g_t(nt) = X*SQRT(CMPLX(Dtau*dble(nt)/dble(thtrot)*ham_Uc/2.d0, 0.D0, kind(0.D0)))
                                    Op_V(nc,nf)%g_t(Ltrot-(nt-1)) = X*SQRT(CMPLX(Dtau*dble(nt)/dble(thtrot)*ham_Uc/2.d0, 0.D0, kind(0.D0)))
                                end do
                            else
                                Op_V(nc,nf)%g = X*SQRT(CMPLX(dtau*Ham_Uc/2.d0, 0.D0, kind(0.D0)))
                            endif
                            Op_V(nc,nf)%alpha = cmplx(-0.5d0, 0.d0, kind(0.D0))
                            Op_V(nc,nf)%type = 2
                            Call Op_set( Op_V(nc,nf) )
                        end do
                    End Do
                end do
            End if
        
            if (abs(Ham_Uf) > Zero) then
                do nf = 1,N_FL
                    X = 1.d0
                    if (nf == 2) X = -1.d0
                    Do I = 1,Latt%N
                        do no = Latt_unit%Norb/2 + 1, Latt_unit%Norb
                            I1 = invlist(I,no)
                            nc = nc + 1
                            print*, 'nf, I, no, I1, nc = ',nf, I, no, I1, nc
                            Op_V(nc,nf)%P(1) = I1
                            Op_V(nc,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
                            If (Adiabatic) then
                                Allocate(OP_V(nc,nf)%g_t(Ltrot))
                                ! call time_dependent_vmfp(hi, ti, lambdai)
                                Op_V(nc,nf)%g_t = X*SQRT(CMPLX(dtau*ham_Uf/2.d0, 0.D0, kind(0.D0)))
                                do nt = 1, Thtrot
                                    ! Op_V(nc,nf)%g_t(nt) = X*SQRT(hi(nt)*Ham_Uf/2.d0, 0.D0, kind(0.D0)))
                                    Op_V(nc,nf)%g_t(nt) = X*SQRT(CMPLX(dtau*dble(nt)/dble(thtrot)*Ham_Uf/2.d0, 0.D0, kind(0.D0)))
                                    Op_V(nc,nf)%g_t(Ltrot-(nt-1)) = X*SQRT(CMPLX(dtau*dble(nt)/dble(thtrot)*Ham_Uf/2.d0, 0.D0, kind(0.D0)))
                                end do
                            else
                                Op_V(nc,nf)%g = X*SQRT(CMPLX(dtau*Ham_Uf/2.d0, 0.D0, kind(0.D0)))
                            endif
                            Op_V(nc,nf)%alpha = cmplx(-0.5d0, 0.d0, kind(0.D0))
                            Op_V(nc,nf)%type = 2
                            Call Op_set( Op_V(nc,nf))
                        end do
                    end Do
                end do
            End if
        
        End Subroutine Ham_V
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
          Allocate ( Obs_scal(5) )
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
             case (5)
                N = 1;   Filename ="Constraint"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo

          ! Equal time correlators
          Allocate ( Obs_eq(4) )
          Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Filename = "Green"
             case (2)
                Filename = "SpinZ"
             case (3)
                Filename = "Den"
             case (4)
                Filename = "Dimer"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Nt = 1
             Channel = '--'
             if (I == 4 ) then
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_f, Channel, dtau)
             else
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit  , Channel, dtau)
             endif
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
                   Channel = 'PH'; Filename = "Den"
                case (4)
                   Channel = 'P' ; Filename = "Greenf"
                case (5)
                   Channel = 'PH'; Filename = "Dimer"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = Ltrot+1-2*Thtrot
                If(Projector) Channel = 'T0'
                if (I == 4 .or.  I == 5 ) then
                   Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_f, Channel, dtau)
                elseif ( I == 1 )  then
                   Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_c, Channel, dtau)
                else
                   Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
                endif
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
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, Zhubc, ZCon, ZJ, Z, ZP,ZS, ZZ, ZXY
          Integer :: I,J, no, n, I_c,I_f, nf, J_c, J_f, no_I, no_J, imj
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
          Call Predefined_Hoppings_Compute_Kin(Hopping_Matrix,List,Invlist, Latt, Latt_unit, GRC, ZKin)
          Zkin = Zkin* dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          Z     = cmplx(real(N_SUN,kind(0.d0)),0.d0,Kind(0.d0))
          ZHubc = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Latt%N
             Do no = 1, Latt_unit%Norb/2
                I_c = invlist(I,no)
                ZHubc =  ZHubc +  Z*( GRC(I_c,I_c,1) - 0.5d0)**2 +  GRC(I_c,I_c,1)* GR(I_c,I_c,1)
             Enddo
          Enddo
          Zhubc = Ham_Uc*Zhubc

          ZJ  = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Latt%N
             Do no = 1, Latt_unit%Norb/2
                I_c  = Invlist(I,no) 
                I_f  = Invlist(I,no + Latt_Unit%Norb/2) 
                ZJ = ZJ +  Z*2.d0*GRC(I_c,I_f,1)* GRC(I_f,I_c,1) +  GRC(I_c,I_c,1)* GR(I_f,I_f,1) + &
                     &     GR(I_c,I_c,1)* GRC(I_f,I_f,1)
             Enddo
          Enddo
          ZJ = -Ham_Vcf*ZJ/2.d0 +  Real(Latt%N* Latt_unit%Norb/8,kind(0.d0))*Ham_Vcf


          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + ( Zhubc + ZJ )*ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf)
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS
          Obs_scal(4)%Obs_vec(1)  =    Obs_scal(4)%Obs_vec(1) + (Zkin + Zhubc + ZJ )*ZP*ZS


          ZCon = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Latt%N
             Do no = 1, Latt_unit%Norb/2  ! Latt_unit%Norb/2 +1 , Latt_unit%Norb
                I_f  = Invlist(I,no + Latt_Unit%Norb/2) 
                ZCon =  ZCon +  Z*( GRC(I_f,I_f,1) - 0.5d0)**2 +  GRC(I_f,I_f,1)* GR(I_f,I_f,1)
             Enddo
          Enddo
          Obs_scal(5)%Obs_vec(1)  =    Obs_scal(5)%Obs_vec(1) + ZCon*ZP*ZS


          ! Standard two-point correlations
          Call Predefined_Obs_eq_Green_measure  ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(1) )
          Call Predefined_Obs_eq_SpinSUN_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(2) )
          Call Predefined_Obs_eq_Den_measure    ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(3) )


          ! Dimer correlations
          obs_eq(4)%N        = obs_eq(4)%N + 1
          obs_eq(4)%Ave_sign = obs_eq(4)%Ave_sign + real(ZS,kind(0.d0))
          Do I = 1,Latt%N
             do no_I  = 1, Latt_unit%Norb / 2
                I_c = Invlist(I,no_I)
                I_f = Invlist(I,no_I + Latt_unit%Norb/2 ) 
                Do J = 1,Latt%N
                   Imj = latt%imj(I,J)
                   do no_J  = 1, Latt_unit%Norb / 2
                      J_c = Invlist(J,no_J)
                      J_f = Invlist(J,no_J + Latt_unit%Norb / 2 )
                      Z  = Predefined_Obs_dimer_eq(I_c,I_f,J_c,J_f, GR, GRC, N_SUN, N_FL) 
                      obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) + Z*ZP*ZS
                   enddo
                enddo
                Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  &
                     &  Predefined_Obs_dimer0_eq(I_c,I_f, GR, N_SUN, N_FL) * ZP*ZS
             enddo
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

        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE,Mc_step_weight)
          Use Predefined_Obs
          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          
          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS, ZZ, ZXY
          Real    (Kind=Kind(0.d0)) :: X
          Integer :: IMJ, I_c, I_f, J_c, J_f, I,J, no_I, no_J

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS*Mc_step_weight

          ! Standard two-point correlations

          Call Predefined_Obs_tau_Green_measure  ( Latt, Latt_unit_c, List_c, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(1) )
          Call Predefined_Obs_tau_SpinSUN_measure( Latt, Latt_unit  , List  , NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(2) )
          Call Predefined_Obs_tau_Den_measure    ( Latt, Latt_unit  , List  , NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(3) )

          ! Greenf correlations
          If (NT == 0 ) then
             obs_tau(4)%N        = obs_tau(4)%N + 1
             obs_tau(4)%Ave_sign = obs_tau(4)%Ave_sign + real(ZS,kind(0.d0))
             obs_tau(5)%N        = obs_tau(5)%N + 1
             obs_tau(5)%Ave_sign = obs_tau(5)%Ave_sign + real(ZS,kind(0.d0))
          endif
          Do I = 1,Latt%N
             do no_I  = 1, Latt_unit%Norb / 2
                I_c = Invlist(I,no_I)
                I_f = Invlist(I,no_I + Latt_unit%Norb / 2 ) 
                Do J = 1,Latt%N
                   Imj = latt%imj(I,J)
                   do no_J  = 1, Latt_unit%Norb / 2
                      J_c = Invlist(J,no_J)
                      J_f = Invlist(J,no_J + Latt_unit%Norb / 2 )
                      Z  = Predefined_Obs_Cotunneling(I_c, I_f, J_c, J_f,  GT0,G0T,G00,GTT, N_SUN, N_FL) 
                      obs_tau(4)%Obs_Latt(imj,NT+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,NT+1,no_I,no_J) + Z*ZP*ZS
                      Z  = Predefined_Obs_dimer_tau(I_c, I_f, J_c, J_f, GT0,G0T,G00,GTT, N_SUN, N_FL) 
                      obs_tau(5)%Obs_Latt(imj,NT+1,no_I,no_J) =  Obs_tau(5)%Obs_Latt(imj,NT+1,no_I,no_J) + Z*ZP*ZS
                   enddo
                enddo
                Z = Predefined_Obs_dimer0_eq(I_c,I_f, GTT, N_SUN, N_FL)
                Obs_tau(5)%Obs_Latt0(no_I) =  Obs_tau(5)%Obs_Latt0(no_I) +  Z*ZP*ZS
             enddo
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

         If  (Ham_Uf  >= 0.d0)  then 
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

         If (Ham_Uf >= 0.d0)  then
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



!(1)local energy
!Requirements of varitional Monte Carlo
!e_n({\boldsymbol \sigma}^\prime,{\boldsymbol \sigma}) = 
!\frac{ \big \langle \psi_\text{MF} \big| U^\dagger_n ({\boldsymbol \sigma}^\prime)  H U_n(\boldsymbol \sigma)\big| \psi_\text{MF} 
! \big\rangle}{W_n({\boldsymbol \sigma}^\prime, {\boldsymbol \sigma})}
         Complex (kind=kind(0.d0)) function E_loc(GR)
            implicit none
            Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
            !Local
            Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL)
            Complex (Kind=Kind(0.d0)) :: Zkin, ZPot,  ZHubc, ZKondo, Z, ZP,ZS
            Integer :: I,J, no, n, I_c,I_f, nf, J_c, J_f, no_I, no_J, imj

            Do nf = 1,N_FL
               Do I = 1,Ndim
                  Do J = 1,Ndim
                     GRC(I, J, nf) = -GR(J, I, nf)
                  Enddo
                  GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
               Enddo
            Enddo
            ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >

            !Kinetic energy
            Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
            Call Predefined_Hoppings_Compute_Kin(Hopping_Matrix,List,Invlist, Latt, Latt_unit, GRC, ZKin)
            Zkin = Zkin* dble(N_SUN)

         
            Z     = cmplx(real(N_SUN,kind(0.d0)),0.d0,Kind(0.d0))
            ZHubc = cmplx(0.d0, 0.d0, kind(0.D0))
            Do I = 1,Latt%N
               Do no = 1, Latt_unit%Norb/2
                  I_c = invlist(I,no)
                  ZHubc =  ZHubc +  Z*( GRC(I_c,I_c,1) - 0.5d0)**2 +  GRC(I_c,I_c,1)* GR(I_c,I_c,1)
               Enddo
            Enddo
            Zhubc = Ham_Uc*Zhubc

            ZKondo  = cmplx(0.d0, 0.d0, kind(0.D0))
            Do I = 1,Latt%N
               Do no = 1, Latt_unit%Norb/2
                  I_c  = Invlist(I,no) 
                  I_f  = Invlist(I,no + Latt_Unit%Norb/2) 
                  ZKondo = ZKondo +  Z*2.d0*GRC(I_c,I_f,1)* GRC(I_f,I_c,1) +  GRC(I_c,I_c,1)* GR(I_f,I_f,1) + &
                        &     GR(I_c,I_c,1)* GRC(I_f,I_f,1)
               Enddo
            Enddo
            ZKondo = -Ham_Vcf*ZKondo/2.d0 +  Real(Latt%N* Latt_unit%Norb/8,kind(0.d0))*Ham_Vcf
            ! Total local energy
            !Kinetic energy -> Zkin
            !Potentail energy -> Zhubc + ZJ
            E_loc =  (Zkin + Zhubc + ZKondo)!
         end function E_loc



!(2)Weight in QMC simulations
!W_n( {\boldsymbol \sigma}^\prime, {\boldsymbol \sigma}) 
!= \big\langle \psi_\text{MF} \big| U^\dagger_n({\boldsymbol \sigma^\prime}) U_n({\boldsymbol \sigma}) \big| \psi_\text{MF}  \big\rangle
!Complex (kind=kind(0.d0)) function weight_loc(phase, weight)
      Complex (kind=kind(0.d0)) function weight_loc(weight)
         implicit none
         complex (Kind=Kind(0.d0)), Intent(in) :: weight(:)
         Integer :: i,j
         !weight_loc =   Weight(nf_calc)*abs(Weight(nf_calc))
         weight_loc =  cmplx(0.d0, 0.d0, kind(0.D0))
      end function weight_loc


      subroutine time_dependent_vmfp(hi, ti, lambdai)
         implicit none
         Real (Kind=Kind(0.d0)) :: gamma,gmin,gmax,fun_val
         ! Local variables
         integer :: i, itr, maxitr, ntl
         Real (Kind=Kind(0.d0)) :: cost, eps, error, U
         Real (Kind=Kind(0.d0)), allocatable :: hi(:), ti(:), lambdai(:)
         !real, external:: fun

         !ALF-> representations
         !Ltrot  = nint(beta/dtau)
         !Thtrot = nint(theta/dtau)!(T=0)
         !Ltrot = Ltrot+2*Thtrot
         !ham_U->U 
         !ntl->n 
         !dtau->dtau 
         !beta->beta 
         ntl = Ltrot
         U = ham_Uf
         allocate(hi(ntl), ti(ntl+1), lambdai(ntl))
         ! Calculate threshhold
         cost = beta / dtau / 2.0
         eps = 1.0e-14
         maxitr = 100
         error = 1.0
         if (cost > real( ntl)) then
            gmin = 0.0
            gmax = 1.0
            itr = 0 
            do while (error > eps .and. itr < maxitr)
               gamma = (gmin + gmax) / 2.0
               fun_val= fun(gamma,  ntl)
               if (fun_val > cost) then
                     gmin = gamma
               else
                     gmax = gamma
               end if
               error = abs(fun_val - cost)
               itr = itr + 1
               ! print *, ' itr,  error =', itr, error, gmin, gmax
            end do
         
            ! Calculate hi
            hi(ntl) = dtau
            do i =  ntl - 1, 1, -1
               hi(i) = hi(i + 1) / gamma
            end do
         else
            gamma = 1.0
            hi = beta / real( ntl) / 2.0
         endif

         do i=1, ntl
            !print*,'i,hi(i)=',i,hi(i)
            lambdai(i) = sqrt(U * hi(i))
         enddo

         ti(1) = hi(1) / 2.0
         do i = 1,  ntl-1
            ti(i+1) = (hi(i+1) + hi(i)) / 2.0
         end do

         ti(ntl + 1) = hi(ntl)
      end subroutine time_dependent_vmfp

      ! Function to fix constraint on time-dependent coupling
      Real (Kind=Kind(0.d0)) function fun(p, nn)
         integer :: nn
         Real (Kind=Kind(0.d0)) :: p
         Real (Kind=Kind(0.d0)) :: cost
         cost = p**(nn - 1)
         fun = (1.d0 - p * cost) / (1.d0 - p) / cost
      end function fun

    end submodule ham_Anderson_smod
