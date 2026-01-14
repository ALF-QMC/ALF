!  Copyright (C) 2022-2023 The ALF project
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
!> of ##NAME## by the Hamiltonian name, populating the subroutines below as needed
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
 
    submodule (Hamiltonian_main) ham_Spin_Peierls_smod

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
      Use Predefined_Lattices 

      Implicit none
      
      type, extends(ham_base) :: ham_Spin_Peierls
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
        procedure, nopass :: weight_reconstruction
        procedure, nopass :: GR_reconstruction
        procedure, nopass :: GRT_reconstruction
        procedure, nopass :: S0
        procedure, nopass :: Hamiltonian_set_nsigma
        ! procedure, nopass :: ##PROCEDURE_NAME##  ! Some other procedure defined in ham_base
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Spin_Peierls

      !#PARAMETERS START# VAR_lattice
      Character (len=64) :: Model = 'Spin_Peierls'  ! Value not relevant
      Character (len=64) :: Lattice_type = 'Square'
      Integer            :: L1 = 6   ! Length in direction a_1
      Integer            :: L2 = 6   ! Length in direction a_2
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_Model_Generic
      !Integer              :: N_SUN        = 2        ! Number of colors
      !Integer              :: N_FL         = 1        ! Number of flavors
      real(Kind=Kind(0.d0)) :: Phi_X        = 0.d0     ! Twist along the L_1 direction, in units of the flux quanta
      real(Kind=Kind(0.d0)) :: Phi_Y        = 0.d0     ! Twist along the L_2 direction, in units of the flux quanta
      logical               :: Bulk         = .true.   ! Twist as a vector potential (.T.), or at the boundary (.F.)
      Integer               :: N_Phi        = 0        ! Total number of flux quanta traversing the lattice
      real(Kind=Kind(0.d0)) :: Dtau         = 0.1d0    ! Thereby Ltrot=Beta/dtau
      real(Kind=Kind(0.d0)) :: Beta         = 5.d0     ! Inverse temperature
      logical               :: Checkerboard = .true.   ! Whether checkerboard decomposition is used
      !logical              :: Symm         = .true.   ! Whether symmetrization takes place
      !logical              :: Projector    = .false.  ! Whether the projective algorithm is used
      real(Kind=Kind(0.d0)) :: Theta        = 10.d0    ! Projection parameter
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_Spin_Peierls
      real(Kind=Kind(0.d0)) :: Ham_Jx       = 0.d0
      real(Kind=Kind(0.d0)) :: Ham_Jy       = 0.d0
      real(Kind=Kind(0.d0)) :: Ham_U        = 0.d0
      real(Kind=Kind(0.d0)) :: Ham_h        = 0.d0
      real(Kind=Kind(0.d0)) :: Ham_g_factor = 2.d0
      real(Kind=Kind(0.d0)) :: Ham_Lambda   = 0.d0
      real(Kind=Kind(0.d0)) :: Ham_Omega0   = 0.d0
      !#PARAMETERS END#
      
      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit, Latt_phi_unit
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:)
      Integer, allocatable :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell
      Integer  ::  nf_calc,  nf_reconst  !  For  flavor  symmetry

      Logical  ::   SU2_Symm =.false.
      Real (Kind=Kind(0.d0))  ::  Ham_M,  Ham_k
      Integer, allocatable  ::   listb(:,:),  invlistb(:,:) 
      ! listb(i:1...LQ,1,N_coord) 

    contains
      
      module Subroutine Ham_Alloc_Spin_Peierls
        allocate(ham_Spin_Peierls::ham)
      end Subroutine Ham_Alloc_Spin_Peierls

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Spin_Peierls_read_write_parameters.F90"

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

!         ! From dynamically generated file "Hamiltonian_Spin_Peierls_read_write_parameters.F90"
          call read_parameters()
! 
          Ltrot = nint(beta/dtau)
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot = Ltrot+2*Thtrot
! 
!         Setup the Bravais lattice
          call Ham_Latt
!         Setup the hopping / single-particle part
          call Ham_Hop
!         Setup the interaction.
          call Ham_V
! 
!           ! Setup the trival wave function, in case of a projector approach
!           if (Projector) Call Ham_Trial()

#ifdef MPI
          If (Irank_g == 0) then
#endif
             File_info = "info"
#if defined(TEMPERING)
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif
             Open(newunit=unit_info, file=file_info, status="unknown", position="append")
             
             Write(unit_info,*) '====================================='
             Write(unit_info,*) 'Antiferromagnet'
             Write(unit_info,*) 'L1            : ', L1
             Write(unit_info,*) 'L2            : ', L2
             if (Projector) then
                Write(unit_info,*) 'Projective version'
                Write(unit_info,*) 'Theta         : ', Theta
                Write(unit_info,*) 'Tau_max       : ', beta
             else
                Write(unit_info,*) 'Finite temperture version'
                Write(unit_info,*) 'Beta          : ', Beta
             endif
             Write(unit_info,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
             Write(unit_info,*) 'N_SUN         : ', N_SUN
             Write(unit_info,*) 'N_FL          : ', N_FL
             Write(unit_info,*) 'J_x           : ', Ham_Jx
             Write(unit_info,*) 'J_y           : ', Ham_Jy
             Write(unit_info,*) 'Ham_U         : ', Ham_U
             Write(unit_info,*) 'Ham_h         : ', Ham_h
             Write(unit_info,*) 'Ham_g         : ', Ham_g_factor
             Write(unit_info,*) 'Ham_Lambda    : ', Ham_Lambda
             Write(unit_info,*) 'Ham_omega0    : ', Ham_omega0
             
             Close(unit_info)
#ifdef MPI
          Endif
#endif

          Ham_k  =  1.d0/ (2.d0*Ham_Lambda)  !Lambda  =  g^2/2k = 1/2k
          Ham_M  =  Ham_k/(Ham_Omega0**2)    !Omega_0 =  sqrt{k/M}
          
          If ( str_to_upper(Lattice_type) /=  "SQUARE" ) then 
             write(error_unit,*)   " The Spin Peierls code, only runs for the square lattice"   
             CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__) 
          endif
          
          If  ( Ham_h <=  1.D-8 )  SU2_Symm = .true.
          If (SU2_Symm   .and. N_FL  .ne. 1 .and.  N_SUN .ne. 2 )  then  
             write(error_unit,*)   " SU(2)  symmetry is present    "
             write(error_unit,*)   " N_FL   has  to be  equal  to 1 "
             write(error_unit,*)   " N_SUN  has  to be  equal  to 2 "
             CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif
          If (.not. SU2_Symm   .and. N_FL  .ne. 2 .and.  N_SUN .ne. 1 )  then  
             write(error_unit,*)   " SU(2)  symmetry is not  present    "
             write(error_unit,*)   " N_FL   has  to be  equal  to 2 "
             write(error_unit,*)   " N_SUN  has  to be  equal  to 1 "
             CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif

          If  (.not.SU2_Symm)  then
             write(error_unit,*)   " SU(2)  spin  symmetry is   required  "
             CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif

          ! Use  particle-hole  symmetry between the two flavors
          If  (.not. SU2_symm )    then
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
!> Sets  the lattice
!> @details
!--------------------------------------------------------------------
        Subroutine  Ham_Latt
          
          Implicit none
          Real (Kind=Kind(0.d0)) ::  a1_p(2), a2_p(2), L1_p(2), L2_p(2) 
          Integer :: nc, no, I 

          Call  Predefined_Latt(Lattice_type, L1, L2, Ndim, List, Invlist, Latt, Latt_Unit )

          Latt_phi_unit%Norb = Latt_Unit%N_coord
          Allocate (Latt_phi_unit%Orb_pos_p(2,2))
          Latt_phi_unit%Orb_pos_p = 0.d0
          
          
        end Subroutine Ham_Latt
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the hopping. In this case this is just  the magnetic  field.
!> @details
!--------------------------------------------------------------------
        Subroutine  Ham_Hop
          Implicit none

          Integer :: nf, n 
          Real (Kind=Kind(0.d0)):: X


          if  (SU2_Symm)   then
             allocate(Op_T(1,N_FL))
             Do nf = 1,N_FL
                x = -1.d0
                If  (nf == 1)  x = 1.d0
                Call Op_make(Op_T(1,nf),1)
                Op_T(1,nf)%P(1)   = 1
                Op_T(1,nf)%O(1,1) = cmplx(1.d0,0.d0,kind(0.D0)) 
                Op_T(1,nf)%g      = cmplx(0.d0,0.d0,kind(0.d0))
                Call Op_set( Op_T(1,nf) )
             enddo
          else
             allocate(Op_T(Ndim,N_FL))
             Do nf = 1,N_FL
                x = -1.d0
                If  (nf == 1)  x = 1.d0
                do n = 1,Ndim 
                   Call Op_make(Op_T(n,nf),1)
                   Op_T(n,nf)%P(1)   = n
                   Op_T(n,nf)%O(1,1) = cmplx(1.d0,0.d0, kind(0.D0)) 
                   Op_T(n,nf)%g      = dtau*x*ham_g_factor*ham_h
                   Call Op_set( Op_T(n,nf) )
                enddo
             enddo
          endif
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
          
          Integer :: nf, I, I1, I2,  nc, nc1,  J, N_op,  Ix, Iy
          Integer :: n, nst,nen, n_fam, no  
          Real (Kind=Kind(0.d0)) :: X,  J_Heis, X_p(2) 



          If  ( L2 ==  1 )  then
             N_op =  Latt%N   +    L1*Latt_Unit%N_coord
             N_fam = 2
          else
             N_op =  Latt%N  +  L1*L2*Latt_unit%N_coord
             N_fam = 4
          endif
          Allocate (Listb(Latt%N,Latt_Unit%N_coord), Invlistb(N_op,2))

          Allocate(Op_V(N_op,N_FL))
          do nf = 1,N_FL
             nc = 0
             do i  = 1, Latt%N 
                nc = nc + 1
                Call Op_make(Op_V(nc,nf),1)
             enddo
             do i  = Latt%N  + 1, N_op
                nc = nc + 1
                Call Op_make(Op_V(nc,nf),2)
             enddo
          enddo
          do nf = 1,N_FL
             nc = 0
             do I  = 1, Latt%N 
                nc = nc + 1
                !Call Predefined_Int_U_SUN( OP_V(nc,nf), I, N_SUN, DTAU, Ham_U  )
                Op_V(nc,nf)%P(1) = I
                Op_V(nc,nf)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                Op_V(nc,nf)%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
                Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0))) 
                Op_V(nc,nf)%type   = 2
                Call Op_set( Op_V(nc,nf) )
             Enddo
             do n = 1,2
                Select case (n)
                case(1)
                   nst = 1;  J_heis=Ham_Jx;  no = 1
                case(2)
                   nst = 2;  J_heis=Ham_Jx;  no = 1
                end Select
                Do Ix = nst, L1, 2
                   do Iy  = 1, L2
                      x_p  = dble(Ix)*latt%a1_p   +  dble(Iy)*Latt%a2_p
                      I =  Inv_R(x_p, Latt) 
                      nc = nc + 1
                      I1 = Invlist(            I     ,1)
                      I2 = Invlist(Latt%nnlist(I,1,0),1)
                      !Call Predefined_Int_V_SUN( OP_V(nc,nf), I1, I2, 1, DTAU, J_Heis/4.d0  )
                      Op_V(nc,nf)%P(1) = I1
                      Op_V(nc,nf)%P(2) = I2
                      Op_V(nc,nf)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
                      Op_V(nc,nf)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
                      Op_V(nc,nf)%g     = SQRT(CMPLX(DTAU*J_heis/4.d0, 0.D0, kind(0.D0))) 
                      Op_V(nc,nf)%alpha = cmplx(0.d0, 0.d0, kind(0.D0))
                      Op_V(nc,nf)%type  = 4
                      Op_V(nc,nf)%flip_protocol = 2
                      Call Op_set( Op_V(nc,nf) )
                      Listb   (I,no)  =  nc
                      Invlistb(nc,1)  =  I
                      Invlistb(nc,2)  =  no
                   enddo
                enddo
             enddo
             do n = 3,N_Fam
                Select case (n)
                case(3)
                   nst = 1;  J_heis=Ham_Jy;  no = 2
                case(4)
                   nst = 2;  J_heis=Ham_Jy;  no = 2  
                end Select
                do Ix  = 1, L1
                   Do Iy = nst, L2, 2
                      x_p  = dble(Ix)*latt%a1_p   +  dble(Iy)*Latt%a2_p
                      I =  Inv_R(x_p, Latt) 
                      nc = nc + 1
                      I1 = Invlist(            I     ,1)
                      I2 = Invlist(Latt%nnlist(I,0,1),1)
                      !Call Predefined_Int_V_SUN( OP_V(nc,nf), I1, I2, 1, DTAU, J_Heis/4.d0  )
                      Op_V(nc,nf)%P(1) = I1
                      Op_V(nc,nf)%P(2) = I2
                      Op_V(nc,nf)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
                      Op_V(nc,nf)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
                      Op_V(nc,nf)%g     = SQRT(CMPLX(DTAU*J_heis/4.d0, 0.D0, kind(0.D0))) 
                      Op_V(nc,nf)%alpha = cmplx(0.d0, 0.d0, kind(0.D0))
                      Op_V(nc,nf)%type  = 4
                      Op_V(nc,nf)%flip_protocol = 2
                      Call Op_set( Op_V(nc,nf) )
                      Listb   (I,no)  =  nc
                      Invlistb(nc,1)  =  I
                      Invlistb(nc,2)  =  no
                   enddo
                enddo
             enddo
          enddo
          
          !Write(6,*) nc, n_op
          
        end Subroutine Ham_V
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
          Complex (Kind=Kind(0.d0)), Intent(In) :: Hs_new

          ! Local
          Real (Kind=Kind(0.d0) )  ::   S_old,  S_new,  J_Heis
          Integer                  ::   ntp1, ntm1

          S0 = 1.d0
          if  (nsigma%t(n) == 4 ) then
             J_Heis =  Ham_Jx
             if  ( invlistb(n,2)   ==  2 )  J_Heis  =  Ham_Jy
             ntp1  = nt + 1
             if ( ntp1 > Ltrot )  ntp1  =  1
             ntm1  = nt - 1 
             if ( ntm1 ==  0   )  ntm1  =  Ltrot
             
             S_old  = Ham_M * ( (aimag(nsigma%f(n,ntp1)-nsigma%f(n,nt)))**2  + (aimag(nsigma%f(n,nt)-nsigma%f(n,ntm1)))**2)/(2.d0*Dtau) + &
                  &   Ham_k * Dtau * ( aimag(nsigma%f(n,nt)) + J_Heis/(4.d0*Ham_k) ) **2/2.d0
             
             S_new  = Ham_M * ( (aimag(nsigma%f(n,ntp1)- Hs_new       ))**2  + (aimag(Hs_new        -nsigma%f(n,ntm1)))**2)/(2.d0*Dtau) + &
                  &   Ham_k * Dtau * ( aimag(Hs_new) +  J_Heis/(4.d0*Ham_k) )**2/2.d0

             S0     = exp(-S_new +  S_old )
             
          endif
        end function S0
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


!           ! Scalar observables
           Allocate ( Obs_scal(3) )
           Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
               N = 1;   Filename = "Pot"
             case (2)
               N = 1;   Filename = "Part"
             case (3)
               N = 2;   Filename = "PhiXY"
             case default
               Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
           enddo

           ! Equal time correlators
           Allocate ( Obs_eq(3) )
           Do I = 1,Size(Obs_eq,1)
              select case (I)
              case (1)
                 Filename = "SpinZ"
              case (2)
                 Filename = "Phi"
              case (3)
                 Filename = "Dimer"
              case default
                 Write(6,*) ' Error in Alloc_obs '
              end select
              Nt = 1
              Channel = '--'
              If (I == 1) then
                 Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
              else
                 Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_phi_unit, Channel, dtau)
              endif
           enddo
           If (Ltau == 1) then
              ! Time-displaced correlators
              Allocate ( Obs_tau(3) )
              Do I = 1,Size(Obs_tau,1)
                 select case (I)
                 case (1)
                    Channel = 'PH' ; Filename = "SpinZ"
                 case (2)
                    Channel = 'PH' ; Filename = "Phi"
                 case (3)
                    Channel = 'PH' ; Filename = "Dimer"
                 case default
                    Write(6,*) ' Error in Alloc_obs '
                 end select
                 Nt = Ltrot+1-2*Thtrot
                 If(Projector) Channel = 'T0'
                 If (I == 1) then
                    Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
                 else
                    Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_phi_unit, Channel, dtau)
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
          Integer,                   INTENT(IN) :: Ntau
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Local
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)) :: ZP, ZS, Zrho,  ZPot,Zmag, Z_phi_x,  Z_phi_y
          Integer :: I, I1, J, J1, nf, no_I,  no_J, nc,  nc1, imj
          ! Add local variables as needed

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

          Z_phi_x  = cmplx(0.d0, 0.d0, kind(0.D0))
          Z_phi_y  = cmplx(0.d0, 0.d0, kind(0.D0))
          do I  = 1,  Latt%N
             nc= Listb(I,1)  
             Z_phi_x  =  Z_phi_x +  aimag(nsigma%f(nc,ntau) )
          enddo
          Z_Phi_x =  Z_Phi_x/dble(Latt%N)
          If  ( L2 >  1 )  then 
             do I  = 1,  Latt%N
                nc= Listb(I,2)  
                Z_phi_y  =  Z_phi_y +  aimag(nsigma%f(nc,ntau) )
             enddo
             Z_Phi_y =  Z_Phi_y/dble(Latt%N)
          endif
          Obs_scal(3)%Obs_vec(1)  =  Obs_scal(3)%Obs_vec(1) + Z_phi_x * ZP*ZS
          Obs_scal(3)%Obs_vec(2)  =  Obs_scal(3)%Obs_vec(2) + Z_phi_y * ZP*ZS
          
          Zrho = cmplx(0.d0, 0.d0, kind(0.D0))
          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          Zmag = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1,Ndim
             ZPot = ZPot + Grc(i,i,1) * Grc(i,i,1)
             ZRho = ZRho + Grc(i,i,1) + Grc(i,i,1)
          Enddo
          Zpot = Zpot 
          ZRho = ZRho*real(N_SUN,kind(0.d0))
          Obs_scal(1)%Obs_vec(1)  =  Obs_scal(1)%Obs_vec(1) + Zpot * ZP*ZS
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + ZRho * ZP*ZS
          
          Call Predefined_Obs_eq_SpinSUN_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(1) )
          !Phonon and  dimer correlations
          Do  I  = 2,3
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          enddo
          Do I  = 1, Latt%N
             Do  No_I = 1, Latt_phi_unit%Norb
                nc  = listb(I,no_I)
                If (no_I == 1) then
                   I1 = latt%nnlist(I,1,0)
                else
                   I1 = latt%nnlist(I,0,1)
                endif
                Do  J  =  1, Latt%N
                   DO no_J  =  1, Latt_phi_unit%Norb
                      nc1  = listb(J,no_J)
                      If (no_J == 1) then
                         J1 = latt%nnlist(J,1,0)
                      else
                         J1 = latt%nnlist(J,0,1)
                      endif
                      imj  = latt%imj(I,J)
                      Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
                           &    Aimag(nsigma%f(nc,ntau)) * Aimag(nsigma%f(nc1,ntau)) *ZP*ZS
                      Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
                           &    Predefined_Obs_dimer_eq(I,I1,J,J1, GR, GRC, N_SUN, N_FL)*ZP*ZS
                   Enddo
                Enddo
                Obs_eq(2)%Obs_Latt0(no_I) = Obs_eq(2)%Obs_Latt0(no_I)  + Aimag(nsigma%f(nc,ntau)) *ZP*ZS
                Obs_eq(3)%Obs_Latt0(no_I) = Obs_eq(3)%Obs_Latt0(no_I) +  &
                     &  Predefined_Obs_dimer0_eq(I,I1, GR, N_SUN, N_FL) * ZP*ZS

             Enddo
          Enddo

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
          Integer :: I, I1, J, J1, No_I, No_J, imj,  nc, nc1, nt_st
          Real (Kind=Kind(0.d0)) :: X
          ! Add local variables as needed

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS * Mc_step_weight

          ! Compute observables
          Call Predefined_Obs_tau_SpinSUN_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(1) )

          If (NT == 0 ) then
             Do  I  =  2,3 
                Obs_tau(I)%N        = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + real(ZS,kind(0.d0))
             enddo
          endif
          Do I  = 1, Latt%N
             Do  No_I = 1, Latt_phi_unit%Norb
                nc  = listb(I,no_I)
                If (no_I == 1) then
                   I1 = latt%nnlist(I,1,0)
                else
                   I1 = latt%nnlist(I,0,1)
                endif
                Do  J  =  1, Latt%N
                   DO no_J  =  1, Latt_phi_unit%Norb
                      nc1  = listb(J,no_J)
                      If (no_J == 1) then
                         J1 = latt%nnlist(J,1,0)
                      else
                         J1 = latt%nnlist(J,0,1)
                      endif
                      imj  = latt%imj(I,J)
                      X  = 0.d0
                      do  nt_st =  1,Ltrot
                         X = X  + Aimag(nsigma%f(nc,nt_st)) * Aimag(nsigma%f(nc1,NPBC_beta(NT+nt_st,Ltrot)))
                      enddo
                      X = X/dble(Ltrot)
                      Obs_tau(2)%Obs_Latt(imj,NT+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,NT+1,no_I,no_J) + &
                           &     cmplx(X,0.d0,kind(0.d0))*ZP*ZS           
                      Obs_tau(3)%Obs_Latt(imj,NT+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,NT+1,no_I,no_J) + &
                           &     Predefined_Obs_dimer_tau(I, I1, J, J1, GT0,G0T,G00,GTT, N_SUN, N_FL)*ZP*ZS
                   Enddo
                Enddo
                Obs_tau(2)%Obs_Latt0(no_I) = Obs_tau(2)%Obs_Latt0(no_I)  + Aimag(nsigma%f(nc,NPBC_beta(NT+1,Ltrot) )) *ZP*ZS
                Obs_tau(3)%Obs_Latt0(no_I) = Obs_tau(3)%Obs_Latt0(no_I)  + &
                     &                       Predefined_Obs_dimer0_eq(I,I1, GTT, N_SUN, N_FL)*ZP*ZS
             Enddo
          Enddo
          
        end Subroutine OBSERT
!--------------------------------------------------------------------
!> @brief
!> Periodic  boundary  conditions
!> @details
!--------------------------------------------------------------------
        Integer  function  NPBC_beta(I,M)
          implicit none

          Integer,  Intent(IN) ::  I, M

          NPBC_beta = I 
          If (I  > M)  NPBC_beta =  I - M
          If (I  < 1)  NPBC_beta =  I + M
          
        end function NPBC_beta

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
          real (kind=kind(0.d0)) :: XI,XJ, ZZ
          
          Do J = 1,Ndim
             XJ = 1.d0
             if (List(J,2) == 1 )  XJ = -1.d0
             Do I = 1,Ndim
                XI=1.0
                if (List(I,2) == 1 )  XI = -1.d0
                ZZ=0.d0
                if (I==J) ZZ=1.d0
                GR(I,J,nf_reconst) = ZZ - XI*XJ*conjg(GR(J,I,nf_calc))
             Enddo
          Enddo
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
        real (kind=kind(0.d0)) :: XI,XJ, ZZ
        
        Do J = 1,NDIM
           XJ = 1.d0
           if (List(J,2) == 1 )  XJ = -1.d0
           Do I = 1,NDIM 
              XI=1.0
              if (List(I,2) == 1 )  XI = -1.d0
              G0T(I,J,nf_reconst) = -XI*XJ*conjg(GT0(J,I,nf_calc))
              GT0(I,J,nf_reconst) = -XI*XJ*conjg(G0T(J,I,nf_calc))
           enddo
         enddo
       end Subroutine GRT_reconstruction

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
             Implicit none

             Complex (Kind=Kind(0.d0)), allocatable, dimension(:,:), Intent(INOUT) :: Initial_field

             Integer ::  N_op,  nt,  nc  ,  I, no
             Logical ::  Test_Dimer  =.False. 

             If  (Test_Dimer)  then 
                N_op  = Size(Op_V,1)
                Allocate ( Initial_field(N_op,Ltrot) )
                
                Initial_field  =  cmplx(0.d0,0.d0,kind(0.d0)) 
                do no  = 1,N_op
                   do nt = 1,Ltrot
                      Initial_field(no,nt)  =  cmplx(1.d0,0.d0,kind(0.d0))
                      If  (ranf_wrap() > 0.5d0)  Initial_field(no,nt)   =  cmplx(-1.d0,0.d0,kind(0.d0))
                   enddo
                enddo
                
                Do  I  =  1,L1
                   nc = Listb(I,1)
                   If  (mod(I,2)  == 0 )  then 
                      do  nt  = 1,Ltrot
                         Initial_field(nc,nt)  =   Initial_field(nc,nt)  + cmplx(0.d0, 1.d0,kind(0.d0))
                      enddo
                   else
                      do  nt  = 1,Ltrot
                         Initial_field(nc,nt)  =   Initial_field(nc,nt)  + cmplx(0.d0,-1.d0,kind(0.d0))
                      enddo
                   endif
                enddo
             endif
             
           end Subroutine Hamiltonian_set_nsigma
        
      end submodule ham_Spin_Peierls_smod
