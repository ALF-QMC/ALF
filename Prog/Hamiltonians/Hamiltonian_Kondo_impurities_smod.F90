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

    submodule (Hamiltonian_main) ham_Kondo_impurities_smod

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
      Use Predefined_Int
      Use LRC_Mod

      Implicit none
      
      type, extends(ham_base) :: ham_Kondo_impurities
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
        procedure, nopass :: Ham_Latt
        procedure, nopass :: Ham_Hop
        procedure, nopass :: weight_reconstruction
        procedure, nopass :: GR_reconstruction
        procedure, nopass :: GRT_reconstruction
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Kondo_impurities

      !#PARAMETERS START# VAR_lattice
      Character (len=64) :: Model = ''  ! Value not relevant
      Character (len=64) :: Lattice_type = 'Square'  ! Possible values: 'Bilayer_square', 'Bilayer_honeycomb'
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
      logical               :: Checkerboard = .false.   ! Whether checkerboard decomposition is used
      !logical              :: Symm         = .true.   ! Whether symmetrization takes place
      !logical              :: Projector    = .false.  ! Whether the projective algorithm is used
      real(Kind=Kind(0.d0)) :: Theta        = 10.d0    ! Projection parameter
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_Kondo_impurities
      real(Kind=Kind(0.d0))     :: ham_T        = 1.d0        ! Hopping parameter
      real(Kind=Kind(0.d0))     :: Ham_T_alterm = 0.0         ! alter-magnetic hopping 
      real(Kind=Kind(0.d0))     :: Ham_chem     = 0.d0        ! Chemical potential
      real(Kind=Kind(0.d0))     :: Ham_mag      = 0.d0        ! Magnetic  field
      real(Kind=Kind(0.d0))     :: Ham_U        = 0.d0        ! Hubbard 
      Character (len=64)        :: d_wave_rep   = "DX2Y2"     ! representation for the alter-magnetic term.
      Character (len=64)        :: Ham_Imp_Kind = "ANDERSON"  ! Kondo           /  Anderson 
      Integer                   :: Ham_N_Imp    = 0           ! # of impurities 
      !#PARAMETERS END#
      !VAR_impurities
      real(Kind=Kind(0.d0)), allocatable :: Imp_t(:,:),  Imp_Jz(:,:), Imp_V(:,:,:)
      !

      Type (Lattice),       target :: Latt_c
      Type (Unit_cell),     target :: Latt_unit_c
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix_c(:)
      Integer, allocatable :: List_c(:,:), Invlist_c(:,:)  

      Type (Lattice),       target :: Latt_f
      Type (Unit_cell),     target :: Latt_unit_f
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix_f(:)  !   Use  this  structure  to  set  the bonds. 
      Integer, allocatable :: List_f(:,:), Invlist_f(:,:)  

      Integer, allocatable ::  Invlist_cf_ff(:,:), GreenPsi_List(:,:)
      real(Kind=Kind(0.d0)), allocatable :: V_cf_ff(:)

      Integer, allocatable :: Invlist_ff_z(:,:)
      Real(Kind=Kind(0.d0)), allocatable :: V_ff_z(:)

      Integer ::  x_shift = 0, y_shift = 0

      !  For  flavor  symmetry   in case  of  particle-hole symmtry
      Integer, allocatable   ::  Prefactor(:)
      Logical  :: Particle_hole
      Integer  ::  nf_calc,  nf_reconst  
      
      NAMELIST /VAR_impurities/    Imp_t, Imp_V, Imp_Jz

    contains
      
      module Subroutine Ham_Alloc_Kondo_impurities
        allocate(ham_Kondo_impurities::ham)
      end Subroutine Ham_Alloc_Kondo_impurities

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Kondo_impurities_read_write_parameters.F90"

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

          integer                :: ierr, nf, unit_info,  unit_para, n, m, n_imp_max
          Character (len=64)     :: file_info,  file_para


#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
#endif

!         !From dynamically generated file "Hamiltonian_Kondo_impurities_read_write_parameters.F90"
          call read_parameters()

!         Position of impurities 
          Allocate (Imp_t(Ham_N_imp,Ham_N_imp), Imp_Jz(Ham_N_imp,Ham_N_Imp), Imp_V(Ham_N_imp,L1,L2) )
          Imp_t = 0.d0;  Imp_V  = 0.d0;  Imp_Jz = 0.d0 
          
#ifdef MPI
          X_Shift =  Irank_g
          Y_Shift =  Irank_g/L1
          !Write(6,*)  Irank, X_Shift, Y_Shift
          If (Irank_g == 0 ) then
#endif          
#ifdef TEMPERING
             write(file_para,'(A,I0,A)') "Temp_", igroup, "/parameters"
#else
             file_para = "parameters"
#endif
             OPEN(NEWUNIT=unit_para, FILE=file_para, STATUS='old', ACTION='read', IOSTAT=ierr)
             READ(unit_para, NML=VAR_impurities)
             CLOSE(unit_para)
#ifdef MPI
          Endif
          !Broadcast parameters to all MPI tasks
          N =  Ham_N_imp*Ham_N_imp
          CALL MPI_BCAST(Imp_t   ,  N,MPI_REAL8    ,0,Group_Comm,ierr)
          CALL MPI_BCAST(Imp_Jz  ,  N,MPI_REAL8    ,0,Group_Comm,ierr)
          N = Ham_N_imp*L1*L2
          CALL MPI_BCAST(Imp_V  ,  N,MPI_REAL8    ,0,Group_Comm,ierr)
#endif
          IF ( .not.( (N_SUN ==2  .and.  N_FL ==1)    .or. ( N_SUN ==1  .and.  N_FL ==2) )  ) THEN
             WRITE(error_unit,*) 'Code  runs  for  N_SUN=2,N_FL =1 or  N_SUN=1,N_FL =2'
             CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif

          IF ( abs(Ham_T_alterm) > 1.D-10 .and. N_FL ==1 )  then
             WRITE(error_unit,*) 'For the Altermagnetic  band, N_FL = 2 is required'
             CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif

           If ( (mod(L1,4) /=0  .or.  mod(L2,4) /=0 ) .and. Checkerboard .and. abs(Ham_T_alterm) > 1.D-10  ) then
            Write(error_unit,*) 'The checkerboard decompostion for the altermagnetic term is only defined for lattices with L1,L2 = 4n'
            Write(error_unit,*) 'Consider setting Checkerboard = .false. or using L1,L2=4n'
            CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif

          ierr = 0
          do  n = 1, Ham_N_imp
             do  m = 1, Ham_N_imp
                if ( abs(Imp_Jz(m,n))   >  1.D-10 )      ierr  = 1
             enddo
          enddo
          If  (abs(Ham_mag)  >  1.D-10 )  ierr = 1
          if (ierr == 1 .and.  N_FL ==1 )  then
             WRITE(error_unit,*) 'Code  needs  N_FL =2, and  N_FL = 1 is set'
             CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif
          
         
          Ltrot = nint(beta/dtau)
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot = Ltrot+2*Thtrot
            
!         ! Setup the Bravais lattice and lists
          call Ham_Latt
          
!  
!         ! Setup the hopping / single-particle part
          call Ham_Hop
! 
!           ! Setup the interaction.
          call Ham_V

          

          if  (N_FL == 2)  then
             ! Setup  prefactor  for falvor  symmetry
             Call  Set_Prefactor(Particle_hole)
             If  (Particle_hole)  then
                allocate(Calc_Fl(N_FL))
                nf_calc=2
                nf_reconst=1
                Calc_Fl(nf_calc)=.True.
                Calc_Fl(nf_reconst)=.False.
             endif
          endif
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
             Write(unit_info,*) 'Model is      :  Kondo_impurities'
             Write(unit_info,*) 'Lattice is    : ', Lattice_type
             Write(unit_info,*) '# c-orbs      : ', Latt_c%N*Latt_unit_c%norb 
             Write(unit_info,*) '# f-orbs      : ', Latt_f%N*Latt_unit_f%norb
             Write(unit_info,*) '# orbs        : ', Ndim
             Write(unit_info,*) '# t           : ', Ham_t
             Write(unit_info,*) '# t_am        : ', Ham_T_alterm
             Write(unit_info,*) '# Symm. t_am  : ', d_wave_rep
             Write(unit_info,*) '# U           : ', Ham_U
             Write(unit_info,*) '# Flux        : ', N_Phi
             Write(unit_info,*) 'Checkerboard  : ', Checkerboard
             Write(unit_info,*) 'Symm. decomp  : ', Symm
             If  (N_FL == 2)  Write(unit_info,*) 'Particle_hole : ', Particle_hole
             if (Projector) then
                Write(unit_info,*) 'Projective version'
                Write(unit_info,*) 'Theta         : ', Theta
                Write(unit_info,*) 'Tau_max       : ', beta
             else
                Write(unit_info,*) 'Finite temperture version'
                Write(unit_info,*) 'Beta          : ', Beta
             endif
             Write(unit_info,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
             Write(unit_info,*) 'N_SUN         : ',   N_SUN
             Write(unit_info,*) 'N_FL          : ', N_FL
             if (Projector) then
                Do nf = 1,N_FL
                   Write(unit_info,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                   Write(unit_info,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
                enddo
             endif
             Write(Unit_info,*) "Number of impurities : ",  Ham_N_imp
             If (str_to_upper(Ham_Imp_Kind) =="ANDERSON")   then   
                Write(Unit_info,*) "Anderson Impurities"
             else
                Write(Unit_info,*) "Kondo  Impurities"
             endif
             
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
!> Specifies the lattice
!> @details
!--------------------------------------------------------------------
        Subroutine  Ham_Latt

          Implicit  none
          
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2), ic_p(2), X
          Integer ::  nc, I, no, Ic, i1, i2, N_bonds, N_bonds_z, nb, n, m

          Latt_Unit_c%Norb      = 1
          Allocate (Latt_unit_c%Orb_pos_p(1,2))
          Latt_Unit_c%Orb_pos_p(1,:) = 0.d0
          a1_p(1) =  1.0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0  ; a2_p(2) =  1.d0
          L1_p    =  dble(L1)*a1_p
          L2_p    =  dble(L2)*a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt_c )

          
          Latt_Unit_f%Norb   = Ham_N_Imp
          Allocate (Latt_unit_f%Orb_pos_p(Ham_N_Imp,2))
          Do n = 1, Ham_N_Imp
            Do i1 = 1, L1
               Do i2 = 1, L2
                  If  ( abs(Imp_V(n,i1,i2)) > 1.0D-10 )  then
                     Latt_Unit_f%Orb_pos_p(n,1) =  i1
                     Latt_Unit_f%Orb_pos_p(n,2) =  i2
                  Endif
               enddo
            enddo
          enddo
          ! The  explicit  position of the f-impurities is  not  required  at  this  point.
          a1_p(1) =  1.0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0  ; a2_p(2) =  1.d0
          L1_p    =  a1_p
          L2_p    =  a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt_f )
          
          Ndim  =  Latt_c%N*Latt_Unit_c%Norb  +  Latt_f%N*Latt_Unit_f%Norb
          !  Setup  the lists
          Allocate (List_c(Ndim,2), Invlist_c(Latt_c%N,Latt_Unit_c%Norb))
          Allocate (List_f(Ndim,2), Invlist_f(Latt_f%N,Latt_Unit_f%Norb))
          List_c  = 0;  List_f = 0; Invlist_c = 0;  Invlist_f = 0
          Nc = 0
          Do I = 1,Latt_c%N
             Do no = 1,Latt_Unit_c%Norb
                nc = nc + 1
                List_c   (nc,1) = I
                List_c   (nc,2) = no
                Invlist_c(I,no) = nc
             Enddo 
          Enddo
          Do I = 1,Latt_f%N
             Do no = 1,Latt_Unit_f%Norb
                nc = nc + 1
                List_f   (nc,1) = I
                List_f   (nc,2) = no
                Invlist_f(I,no) = nc 
             Enddo
          Enddo


          ! Bond  list of  coupling between  c-f orbitals  and  ff orbitals
          N_bonds   =  0
          N_bonds_z =  0 
          do  no  =  1, Latt_Unit_f%Norb
             do i1 = 1,L1
                do  i2 = 1,L2
                   if  (abs(Imp_V(no,i1,i2)) > 1.0D-10 )   N_bonds =  N_bonds + 1
                enddo
             enddo
          enddo
          Do n  =  1,  Ham_N_imp
             Do m  =  n+1, Ham_N_imp   
                If  (abs(Imp_t (n,m)) > 0.d0 .or. abs(Imp_t (m,n)) > 0.d0  )  N_bonds   = n_Bonds + 1
                If  (abs(Imp_Jz(n,m)) > 0.d0 .or. abs(Imp_Jz(m,n)) > 0.d0  )  N_bonds_z = n_Bonds_z + 1
             enddo
          enddo
          Allocate (Invlist_cf_ff(N_bonds  ,2), V_cf_ff(N_bonds  ) )
          Allocate (Invlist_ff_z (N_bonds_z,2), V_ff_z (N_bonds_z) )
          nc  = 0 
          do  no  =  1, Latt_Unit_f%Norb
             do i1 = 1,L1
                do  i2 = 1,L2
                   if  (abs(Imp_V(no,i1,i2)) > 1.0D-10 )    then
                      nc  = nc  + 1
                      ic_p = dble(i1+x_shift)*Latt_c%a1_p + dble(i2+y_shift)*Latt_c%a2_p
                      Ic   = Inv_R(ic_p,Latt_c)
                      Invlist_cf_ff(nc,1) = Invlist_f(1,no)     !  f-orbital
                      Invlist_cf_ff(nc,2) = Invlist_c(Ic,1)     !  c-orbital
                      V_cf_ff(nc)       = Imp_V(no,i1,i2)     !  Amplitude 
                   endif
                enddo
             enddo
          enddo
          Do n  =  1,  Ham_N_imp
             Do m  =  n+1, Ham_N_imp   
                If (abs(Imp_t(n,m)) > 0.d0 .or. abs(Imp_t(m,n)) > 0.d0  )   Then
                   nc= nc + 1
                   X = Imp_t(n,m)
                   if  ( abs(Imp_t(m,n)) >  abs(X) )  X  =  Imp_t(m,n)
                   Invlist_cf_ff(nc,1) = Invlist_f(1,n)     !  f-orbital
                   Invlist_cf_ff(nc,2) = Invlist_f(1,m)     !  f-orbital
                   V_cf_ff(nc)       = X                  !  Amplitude 
                endif
             enddo
          enddo
          
          nc = 0
          Do n  =  1,  Ham_N_imp
             Do m  =  n+1, Ham_N_imp   
                If (abs(Imp_Jz(n,m)) > 0.d0 .or. abs(Imp_Jz(m,n)) > 0.d0  )   Then
                   nc= nc + 1
                   X = Imp_Jz(n,m)
                   if  ( abs(Imp_Jz(m,n)) >  abs(X) )  X  =  Imp_Jz(m,n)
                   Invlist_ff_z(nc,1) = Invlist_f(1,n)     !  f-orbital
                   Invlist_ff_z(nc,2) = Invlist_f(1,m)     !  f-orbital
                   V_ff_z(nc)         = X                  !  Amplitude 
                endif
             enddo
          enddo
          
          !Testing
          !Do nb  =1, N_bonds
          !   Write(6,"(I4,2x,I4,2x,I4,2x,F14.7)") nb,Invlist_cf_ff(nb,1), Invlist_cf_ff(nb,2), V_cf_ff(nb)
          !enddo
          !Write(6,*)
          !Do nb  =1, N_bonds_z
          !   Write(6,"(I4,2x,I4,2x,I4,2x,F14.7)") nb,Invlist_ff_z(nb,1), Invlist_ff_z(nb,2), V_ff_z(nb)
          !enddo
          
        end Subroutine Ham_Latt
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the  hopping
!> @details
!--------------------------------------------------------------------
        Subroutine  Ham_Hop

          Implicit  none

          Real (Kind=Kind(0.d0) ), allocatable :: Ham_T_vec(:), Ham_T_vecam(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:)
          Integer, allocatable ::   N_Phi_vec(:)

          Integer :: nc_f,nc_tot, n, m,  nc, nf,  nb, N_bonds
          Real (Kind=Kind(0.d0)) ::  Factor, X 
          Type (Operator),     dimension(:,:), allocatable :: Op_T_c, OP_T_cf_ff
          Logical ::  Symm1
          


          Select case  (str_to_upper(Lattice_type))
          Case("SQUARE")
             N_phi  =  1
          Case ("PI_FLUX")
             N_phi  =  Latt_c%N/2
          case default
             Write(error_unit,*) 'Lattice  has  to be  set to Pi_flux or  Square'
             CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          end Select


          Allocate (Ham_T_vec(N_FL),  Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL), N_Phi_vec(N_FL) )

          ! Here we consider no N_FL  dependence of the hopping parameters.
          Ham_T_vec      = Ham_T
          Ham_Chem_vec   = Ham_Chem
          Phi_X_vec      = Phi_X
          Phi_Y_vec      = Phi_Y
          N_Phi_vec      = N_Phi

          If (abs(Ham_T_alterm) > 1.d-10 ) then
            Allocate (Ham_T_vecam(N_FL))
            Ham_T_vecam(1)    =  Ham_T_alterm
            Ham_T_vecam(2)    = -Ham_T_alterm
            Call  Set_Default_hopping_parameters_square_am(Hopping_Matrix_c,Ham_T_vec,Ham_T_vecam,Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                &                                      Bulk, N_Phi_vec, N_FL, List_c, Invlist_c, Latt_c, Latt_unit_c,d_wave_rep  )
            deallocate (Ham_T_vecam)
          else
            Call  Set_Default_hopping_parameters_square(Hopping_Matrix_c,Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
               &                                      Bulk, N_Phi_vec, N_FL, List_c, Invlist_c, Latt_c, Latt_unit_c )
          endif
          If  ( str_to_upper(Ham_Imp_Kind) == "ANDERSON"  )  then
             Call  Predefined_Hoppings_set_OPT(Hopping_Matrix_c,List_c,Invlist_c,Latt_c,  Latt_unit_c,  Dtau, Checkerboard, Symm, OP_T_c )
             
             Factor =  1.d0 
             If ( Symm )  Factor = 0.5d0
             N_bonds  = Size(V_cf_ff,1) 
             Allocate(Op_T_cf_ff(N_bonds,N_FL))
             do nf = 1,N_FL
                do nc = 1, N_bonds
                   Call Op_make(Op_T_cf_ff(nc,nf),2)
                enddo
                Do nc  =  1, N_bonds
                   Op_T_cf_ff(nc,nf)%P(1) = Invlist_cf_ff(nc,1) 
                   Op_T_cf_ff(nc,nf)%P(2) = Invlist_cf_ff(nc,2) 
                   Op_T_cf_ff(nc,nf)%O(1,2) = cmplx(V_cf_ff(nc),0.d0,kind(0.d0))
                   Op_T_cf_ff(nc,nf)%O(2,1) = cmplx(V_cf_ff(nc),0.d0,kind(0.d0))
                   Op_T_cf_ff(nc,nf)%g = -Dtau * Factor
                   Op_T_cf_ff(nc,nf)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                   !Call Op_set(Op_T(nc,nf))
                enddo
             enddo
             nc_tot =  size(OP_T_cf_ff,1) + size(OP_T_c,1)
             If ( Symm )   nc_tot  =  2*size(OP_T_cf_ff,1) + size(OP_T_c,1)
             Allocate (OP_T(nc_tot,N_FL))
             do nf  = 1,N_FL
                nc  = 0
                do n =  1, size(OP_T_cf_ff,1)
                   nc  = nc  + 1
                   Call Op_make(Op_T(nc,nf),2)
                   Call Copy_Op(OP_T(nc,nf), OP_T_cf_ff(n,nf))
                enddo
                do n =  1, size(OP_T_c,1)
                   nc  = nc  + 1
                   Call Op_make(Op_T(nc,nf),2)
                   Call Copy_Op(OP_T(nc,nf), OP_T_c(n,nf))
                enddo
                If  (symm) then
                   do n =  size(OP_T_cf_ff,1),1,-1
                      nc  = nc  + 1
                      Call Op_make(Op_T(nc,nf),2)
                      Call Copy_Op(OP_T(nc,nf), OP_T_cf_ff(n,nf))
                   enddo
                endif
             enddo
             !  Clear  memory
             do  nf  = 1,N_FL
                Do  n =  1,size(OP_T_c,1)
                   Call  Op_clear(OP_T_c(n,nf), 2)
                enddo
                Do  n =  1,size(OP_T_cf_ff,1)
                   Call  Op_clear(OP_T_cf_ff(n,nf), 2)
                enddo
             enddo
             deallocate (OP_T_c, OP_T_cf_ff)
             !Write(6,*) "Size OP_T,1: ", Size(OP_T,1)
          else
             Call  Predefined_Hoppings_set_OPT(Hopping_Matrix_c,List_c,Invlist_c,Latt_c,  Latt_unit_c,  Dtau, Checkerboard, Symm, OP_T )
          endif
          
          Deallocate (Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, N_Phi_vec )
          
        end Subroutine Ham_Hop

!--------------------------------------------------------------------
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
!> Sets  the interaction.
!> @details
!> The  Hubbard  interaction  commutes  with the   Heisenberg and
!> Kondo  interactions, such  that we  can pull it out.
!> For  the  other  terms  we  use  a  symmettic    decompostion. 
!> e^{-\Delta \tau H_J^{even}/2} e^{-\Delta \tau H_J^{odd}/2}  e^{-Delta \tau  H_JK}
!> e^{-\Delta \tau H_J^{odd}/2} e^{-\Delta \tau H_J^{even}/2} e^{-\Delta \tau H_U}
!--------------------------------------------------------------------
        Subroutine  Ham_V

          Implicit  none

          Integer :: N_op_U ,  N_op_K, N_op_H,  N_op,  n,  nf, nc
          Integer :: I, J 


          If  (str_to_upper(Ham_Imp_Kind) == "ANDERSON" )  then
             N_op_U   =  Latt_unit_f%Norb  !  Hubbard
             N_op     =  N_op_U
             Allocate(Op_V(N_op,N_FL))
             Do  nf   =  1,N_FL
                nc = 0
                ! Hubbard
                Do  n = 1,Latt_unit_f%Norb
                   nc = nc + 1
                   I = Invlist_f(1,n)  !  f-orbital 
                   Call Predefined_Int_U_SUN( OP_V(nc,nf), I, 1, DTAU, Ham_U/2.d0  )
                enddo
             enddo
          else
             ! Set the Kondo 
             N_op  = Latt_unit_f%Norb +   Size(V_cf_ff,1)  +  Size(V_ff_z,1)
             If  (Symm)   N_op  = Latt_unit_f%Norb +   2*Size(V_cf_ff,1) + 2*Size(V_ff_z,1)
             Allocate(Op_V(N_op,N_FL))
             nc = 0
             ! Exchange
             Do n  =  1, Size(V_cf_ff,1)
                nc = nc + 1
                I = Invlist_cf_ff(n,1)
                J = Invlist_cf_ff(n,2)
                Do  nf =  1,N_FL
                   If  (Symm)  then 
                      Call Predefined_Int_V_SUN( OP_V(nc,nf), I, J, 1, DTAU, V_cf_ff(n)/8.d0  )
                   else
                      Call Predefined_Int_V_SUN( OP_V(nc,nf), I, J, 1, DTAU, V_cf_ff(n)/4.d0  )
                   endif
                enddo
             enddo
             do n = 1,Size(V_ff_z,1)
                nc = nc + 1
                I = Invlist_ff_z(n,1)
                J = Invlist_ff_z(n,2)
                If  (Symm)  then 
                   Call Predefined_Int_Jz( OP_V(nc,1), OP_V(nc,2), I, J, DTAU, V_ff_z(n)/2.d0  )
                else
                   Call Predefined_Int_Jz( OP_V(nc,1), OP_V(nc,2), I, J, DTAU, V_ff_z(n)  )
                endif
             enddo
             ! Hubbard
             Do  n = 1,Latt_unit_f%Norb
                nc = nc + 1
                I = Invlist_f(1,n)  !  f-orbital 
                Do  nf =  1,N_FL
                   Call Predefined_Int_U_SUN( OP_V(nc,nf), I, 1, DTAU, Ham_U/2.d0  )
                enddo
             enddo
             If  (Symm) then 
                do n = Size(V_ff_z,1),1, -1
                   nc = nc + 1
                   I = Invlist_ff_z(n,1)
                   J = Invlist_ff_z(n,2)
                   Call Predefined_Int_Jz( OP_V(nc,1), OP_V(nc,2), I, J, DTAU, V_ff_z(n)/2.d0  )
                enddo
                Do n  =  Size(V_cf_ff,1),1, -1
                   nc = nc + 1
                   I = Invlist_cf_ff(n,1)
                   J = Invlist_cf_ff(n,2) 
                   Do  nf =  1,N_FL
                      Call Predefined_Int_V_SUN( OP_V(nc,nf), I, J, 1, DTAU, V_cf_ff(n)/8.d0  )
                   enddo
                enddo
             endif
          endif

          
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
          Integer    ::  i, N, Nt,  N_files, ix, iy, nc, Ic, nf
          Character (len=64) ::  Filename
          Character (len=64), allocatable ::  Filename_lst(:)
          Character (len=2)  ::  Channel
          Real (Kind=Kind(0.d0))  ::  ic_p(2)

          If ( str_to_upper(Ham_Imp_Kind) == "KONDO" )  then
             
             Filename ="GreenPsi_"
             N_Files = 0
             Do n  = 1,Ham_N_imp
                Do  ix = 1,L1
                   Do iy = 1,L2
                      If  (abs(Imp_V(n,ix,iy)) > 0.d0 )  N_Files =  N_Files  + N_FL
                   Enddo
                Enddo
             Enddo
             Allocate  (Filename_lst(N_Files))
             Allocate  (GreenPsi_List(N_Files,3))
             nc = 0
             Do n  = 1,Ham_N_imp
                Do  ix = 1,L1
                   Do iy = 1,L2
                      If  (abs(Imp_V(n,ix,iy)) > 0.d0 )  Then
                         Do nf = 1,N_FL
                           nc = nc + 1
                           If (N_FL == 1)  then 
                              Write(Filename_lst(nc), '(A,I0,A,I0,A,I0)') trim(Filename),n,"_",ix,"_",iy
                           else
                              Select case(nf) 
                              case(1)
                                Write(Filename_lst(nc), '(A,I0,A,I0,A,I0,A)') trim(Filename),n,"_",ix,"_",iy,"_up"
                              case(2)
                                 Write(Filename_lst(nc), '(A,I0,A,I0,A,I0,A)') trim(Filename),n,"_",ix,"_",iy,"_do"
                              end select
                           endif
                           ic_p = dble(ix+ x_shift)*Latt_c%a1_p + dble(iy + y_shift)*Latt_c%a2_p
                           Ic   = Inv_R(ic_p,Latt_c)
                           GreenPsi_List(nc,1)  = Invlist_f(1,n)  ! f-electron
                           GreenPsi_List(nc,2)  = Invlist_c(Ic,1) ! c-electron
                           GreenPsi_List(nc,3)  = nf              ! Spin  nf = 1 is spin up,  nf = 2 is spin do
                         enddo
                      endif
                   Enddo
                Enddo
             Enddo
          Else
             N_Files = Ham_N_Imp
             Allocate  (Filename_lst(N_Files))
             Filename ="Greend_"
             Do nc  = 1,Ham_N_imp
                Write(Filename_lst(nc), '(A,I0)') trim(Filename),nc
             Enddo
          Endif
          !Write(6,*)  'Number of files: ', N_Files
          !Do  nc  = 1,  N_Files
          !   Write(6,*)  Filename_lst(nc)
          !enddo
          
!         Scalar observables
          Allocate ( Obs_scal(1) )
          Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
                N = 1;   Filename = "Pot"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo
          
!         Equal time correlators
          If  (N_SUN == 1)  then
             Allocate ( Obs_eq(3) )
             Do I = 1,Size(Obs_eq,1)
                select case (I)
                case (1)
                   Filename = "SpinZ"
                case (2)
                   Filename = "SpinXY"
                case (3)
                   Filename = "SpinT"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = 1
                Channel = '--'
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt_f, Latt_unit_f, Channel, dtau)
             enddo
          else
             Allocate ( Obs_eq(1) )
             Do I = 1,Size(Obs_eq,1)
                select case (I)
                case (1)
                   Filename = "SpinZ"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = 1
                Channel = '--'
                Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt_f, Latt_unit_f, Channel, dtau)
             enddo
          endif
! 
          If (Ltau == 1) then
             ! Time-displaced correlators
             If  (N_FL == 2)  then
                Allocate ( Obs_tau(3+N_Files) )
             Else
                Allocate ( Obs_tau(1+N_Files) )
             endif
             Do I = 1,N_Files
                Nt = Ltrot+1-2*Thtrot
                Channel = 'P'
                If (Projector) Channel = 'T0'
                Filename = Filename_lst(I)
                Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt_f, Latt_unit_c, Channel, dtau)
             enddo
             do I  =  N_files + 1,  size( Obs_tau, 1 )
                nc =   I - N_Files
                select case (nc)
                case (1)
                   Channel = 'PH' ; Filename = "SpinZ"
                case (2)
                   Channel = 'PH' ; Filename = "SpinXY"
                case (3)
                   Channel = 'PH' ; Filename = "SpinT"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = Ltrot+1-2*Thtrot
                If(Projector) Channel = 'T0'
                Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt_f, Latt_unit_f, Channel, dtau)
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
          Complex (Kind=Kind(0.d0)) :: GT0_tmp(Ndim,Ndim,N_FL),  GTT_tmp(Ndim,Ndim,N_FL), G0T_tmp(Ndim,Ndim,N_FL)

          Complex (Kind=Kind(0.d0)) :: ZP, ZS
          Complex (Kind=Kind(0.d0)) :: Z_pot, Z
          Integer :: I, I_c, I_f,  J, J_c, J_f,  nf, n, imj, no_I, no_J
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
          ! Compute scalar observables.
          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo

          Z_pot  = cmplx(0.d0,0.d0, kind(0.D0))
          do I  =  1, Latt_f%N
             I_f  =  invlist_f(I,1)  !  f-orbital
             If  (N_FL  == 1 )  then
                Z_pot  = Z_pot  +  Grc(I_f,I_f,1)*Grc(I_f,I_f,1)
             else
                Z_pot  = Z_pot  +  Grc(I_f,I_f,1)*Grc(I_f,I_f,2)
             endif
          enddo
          Obs_scal(1)%Obs_vec(1)  =  Obs_scal(1)%Obs_vec(1) + Z_pot *ZP* ZS


          If (N_FL == 2 )  then
             Call Predefined_Obs_eq_SpinMz_measure( Latt_f, Latt_unit_f, List_f,  GR, GRC, N_SUN, ZS, ZP, &
                  & Obs_eq(1), Obs_eq(2), Obs_eq(3) )
          else
             Call Predefined_Obs_eq_SpinSUN_measure( Latt_f, Latt_unit_f, List_f,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(1) )
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
          Complex (Kind=Kind(0.d0)) :: Z
          Integer ::  I, I_c,I_f, J,J_c, J_f,   imj, N, N_Files, spin
          
          
          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS * Mc_step_weight


          If (N_FL == 1 )  N_files  =  size(Obs_tau,1) - 1
          If (N_FL == 2 )  N_files  =  size(Obs_tau,1) - 3

          Do I = 1, N_Files
             If (NT == 0 ) then
                Obs_tau(I)%N         =  Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign  =  Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             Endif
             If ( str_to_upper(Ham_Imp_Kind)  == "ANDERSON" )   Then
                I_f  =  Invlist_f(1,I)
                If ( N_FL == 1 ) Obs_tau(I)%Obs_Latt(1,NT+1,1,1) =  Obs_tau(I)%Obs_Latt(1,NT+1,1,1) + &
                     &     GT0(I_f,I_f,1)   * dble(N_SUN) * ZP*ZS
                If ( N_FL == 2 ) Obs_tau(I)%Obs_Latt(1,NT+1,1,1) =  Obs_tau(I)%Obs_Latt(1,NT+1,1,1) + &
                     &  ( GT0(I_f,I_f,1) + GT0(I_f,I_f,2) ) * ZP*ZS
             else
                I_f = GreenPsi_List(I,1) ! f-electron
                I_c = GreenPsi_List(I,2) ! c-electron
                If (N_FL == 1) then
                  Obs_tau(I)%Obs_Latt(1,NT+1,1,1) =  Obs_tau(I)%Obs_Latt(1,NT+1,1,1) + &
                        &  Predefined_Obs_Cotunneling(I_c, I_f, I_c, I_f,  GT0,G0T,G00,GTT, N_SUN, N_FL     ) * ZP*ZS
                else
                  spin = GreenPsi_List(I,3) 
                  Obs_tau(I)%Obs_Latt(1,NT+1,1,1) =  Obs_tau(I)%Obs_Latt(1,NT+1,1,1) + &
                        &  Predefined_Obs_Cotunneling(I_c, I_f, I_c, I_f,  GT0,G0T,G00,GTT, N_SUN, N_FL,spin) * ZP*ZS 
                endif
             endif
          enddo
          
          If (N_FL == 2 )  then
             Call Predefined_Obs_tau_SpinMz_measure( Latt_f, Latt_unit_f, List_f, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, &
                  & Obs_tau(N_Files+1), Obs_tau(N_files+2), Obs_tau(N_Files+3) )
          else
             Call Predefined_Obs_tau_SpinSUN_measure( Latt_f, Latt_unit_f, List_f, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(N_Files+1) )
          endif
          

          ! Compute observables
!!$          
!!$          Do I  = 1, Latt_f%N
!!$             I_f  =   Invlist_f(I,1) !  f-orbital 
!!$             I_c  =   Invlist_f(I,2) !  c-orbital 
!!$             Do J = 1,  Latt_f%N
!!$                J_f  =   Invlist_f(J,1) !  f-orbital 
!!$                J_c  =   Invlist_f(J,2) !  c-orbital 
!!$                imj = latt_f%imj(I,J)
!!$                Z =  - G0T(J_f,I_f,1) * GT0(I_f,J_f,1) * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
!!$                Obs_tau(1)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(1)%Obs_Latt(imj,NT+1,1,1) + Z*ZP*ZS
!!$                Z =  Predefined_Obs_Cotunneling(I_c, I_f, J_c, J_f,  GT0,G0T,G00,GTT, N_SUN, N_FL)  
!!$                Obs_tau(2)%Obs_Latt(imj,NT+1,1,1) =  Obs_tau(2)%Obs_Latt(imj,NT+1,1,1) + Z*ZP*ZS
!!$             Enddo
!!$          Enddo

        end Subroutine OBSERT


!-------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief
!> Sets  the  prefactor  for  flavor  symmetry
!--------------------------------------------------------------------

        Subroutine  Set_Prefactor(Particle_Hole)

          Implicit none

          Logical, Intent(out) ::  Particle_Hole
          
          Integer i,  ix, iy, nc, nc1, no, no1,  i_f, i_c
          Real (Kind=Kind(0.d0))  ::  ic_p(2), X

          Allocate(Prefactor(Ndim))
          Particle_hole = .true.
          Prefactor  = 0
          Do  i =  1, Latt_c%N
             ix  = Latt_c%list(i,1) 
             iy  = Latt_c%list(i,2)
             nc  = invlist_c(i,1)
             Prefactor(nc)  = -1
             if (mod(ix+iy,2)  == 0 )   Prefactor(nc)  = 1
          enddo

          nc  = 0
          do no  =  1,Latt_unit_f%Norb
             do ix  = 1,L1
                do iy = 1,L2
                   If  ( abs(Imp_V(no,ix,iy)) >  1.D-10 ) then
                      nc = nc + 1
                      i_f = invlist_f(1,no)
                      ic_p = dble(ix + x_shift)*Latt_c%a1_p + dble(iy + y_shift)*Latt_c%a2_p
                      I    = Inv_R(ic_p,Latt_c)
                      i_c =  Invlist_c(I,1)
                      If  ( str_to_upper(Ham_Imp_Kind) == "KONDO" )  then
                         I =  -Prefactor(i_c)*nint(Imp_V(no,ix,iy)/abs(Imp_V(no,ix,iy)))
                      else
                         I =  -Prefactor(i_c)
                      endif
                      If ( Prefactor(i_f) ==  0 )  then
                         Prefactor(i_f)  = I 
                      else
                         if  ( Prefactor(i_f) /= I ) then
                            Particle_hole = .False.
                         endif
                      endif
                   endif
                enddo
             enddo
          enddo
          if  (nc   == 0 )   Prefactor(Invlist_f(1,1)) = 1   !   Set  the overall factor
          do no  =  1,Latt_unit_f%Norb
             nc = Invlist_f(1,no)
             If (Prefactor(nc) /= 0 )  then 
                do no1  =  1,Latt_unit_f%Norb
                   If  ( abs(Imp_t(no1,no)) > 1.D-10 )  then
                      If  ( str_to_upper(Ham_Imp_Kind) == "KONDO" )  then
                         I = - nint(Imp_t(no1,no)/abs(Imp_t(no1,no)) ) * Prefactor(nc)
                      else
                         I = - Prefactor(nc)
                      endif
                      nc1 = Invlist_f(1,no1)
                      if (  Prefactor(nc1) ==  0 )  then
                         Prefactor(nc1) = I
                      elseif ( Prefactor(nc1) /=  I  ) then
                         Particle_hole= .false.
                      endif
                   endif
                   If  ( abs(Imp_t(no,no1)) > 1.D-10 )  then
                      If  ( str_to_upper(Ham_Imp_Kind) == "KONDO" )  then
                         I = -  Prefactor(nc) * nint(Imp_t(no,no1)/abs(Imp_t(no,no1)) ) 
                      else
                         I = - Prefactor(nc)
                      endif
                      nc1 = Invlist_f(1,no1)
                      if (  Prefactor(nc1) ==  0 )  then
                         Prefactor(nc1) = I
                      elseif ( Prefactor(nc1) /=  I  ) then
                         Particle_hole= .false.
                      endif
                   endif
                enddo
             endif
          enddo
          
          Do  nc = 1,  size(Prefactor,1)
             if  (Prefactor(nc)  == 0  )  then
                WRITE(error_unit,*) 'One  orbital  is not  linked to the  cluster.  Consider  adapting the  the parameter file', nc
                CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
             endif
          enddo
          
!!$          Do  i =  1, Latt_c%N
!!$             ix  = Latt_c%list(i,1) 
!!$             iy  = Latt_c%list(i,2)
!!$             nc  = invlist_c(i,1)
!!$             write(6,*) nc, ix,iy, Prefactor(nc) 
!!$          enddo
!!$          Do  no =  1, Latt_unit_f%Norb
!!$             nc  = invlist_f(1,no)
!!$             write(6,*) nc, no, Prefactor(nc) 
!!$          enddo

        end Subroutine Set_Prefactor


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
          Integer :: I,J
          complex (kind=kind(0.d0))  ::  ZZ
          Real    (kind=kind(0.d0))  :: X
          
          Do J = 1,Ndim
             Do I = 1,Ndim
                X =  real(Prefactor(I)*Prefactor(J), kind(0.d0))
                ZZ=cmplx(0.d0,0.d0,Kind(0.d0))
                if (I==J) ZZ=cmplx(1.d0,0.d0,Kind(0.d0))
                GR(I,J,nf_reconst) = ZZ - X*conjg(GR(J,I,nf_calc))
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
        Integer :: I,J
        real (kind=kind(0.d0)) :: X
        
        Do J = 1,NDIM
           Do I = 1,NDIM 
              X =  real(Prefactor(I)*Prefactor(J), kind(0.d0))
              G0T(I,J,nf_reconst) = -X*conjg(GT0(J,I,nf_calc))
              GT0(I,J,nf_reconst) = -X*conjg(G0T(J,I,nf_calc))
           enddo
         enddo
       end Subroutine GRT_reconstruction        
        
      end submodule ham_Kondo_impurities_smod
