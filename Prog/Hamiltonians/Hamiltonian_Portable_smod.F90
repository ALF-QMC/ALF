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

    submodule (Hamiltonian_main) ham_Portable_smod

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
      
      type, extends(ham_base) :: ham_Portable
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
        ! procedure, nopass :: ##PROCEDURE_NAME##  ! Some other procedure defined in ham_base
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Portable


      !#PARAMETERS START# VAR_Model_Generic
      !Integer              :: N_SUN        = 2        ! Number of colors
      !Integer              :: N_FL         = 1        ! Number of flavors
      real(Kind=Kind(0.d0)) :: Phi_X        = 0.d0     ! Twist along the L_1 direction, in units of the flux quanta
      real(Kind=Kind(0.d0)) :: Phi_Y        = 0.d0     ! Twist along the L_2 direction, in units of the flux quanta
      logical               :: Bulk         = .true.   ! Twist as a vector potential (.T.), or at the boundary (.F.)
      Integer               :: N_Phi        = 0        ! Total number of flux quanta traversing the lattice
      real(Kind=Kind(0.d0)) :: Dtau         = 0.1d0    ! Thereby Ltrot=Beta/dtau
      real(Kind=Kind(0.d0)) :: Beta         = 5.d0     ! Inverse temperature
      !#PARAMETERS END#

      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:)
      Integer, allocatable :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell
      logical  :: translation_symmetry

      integer :: n_spin

      ! variables for lattice
      integer  :: Norb
      Integer :: L1, L2
      Real (Kind=Kind(0.d0))  :: a1_p(3), a2_p(3), a3_p(3)
      real (kind=kind(0.d0)), allocatable :: Orb_a_p(:,:)

      !variables for hopping
      integer :: N_diag
      complex (kind=kind(0.d0)), allocatable :: Ham_t(:)
      integer, allocatable                   :: hop_list(:,:), hop_diag(:)


    contains
      
      module Subroutine Ham_Alloc_Portable
        allocate(ham_Portable::ham)
      end Subroutine Ham_Alloc_Portable

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Portable_read_write_parameters.F90"

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

           ! From dynamically generated file "Hamiltonian_Portable_read_write_parameters.F90"
           call read_parameters()
   
           !For now fix N_SUN and N_FL
           N_SUN  = 1
           N_FL   = 1
           N_spin = 2

           translation_symmetry = .true.
           Ltrot = nint(beta/dtau)
!           Thtrot = 0
!           if (Projector) Thtrot = nint(theta/dtau)
!           Ltrot = Ltrot+2*Thtrot
 
           ! read in information for lattice from file geometry.txt
           call read_latt

           ! read in information for hopping from file hoppings.txt
           call read_hop
           
           ! Setup the Bravais lattice
           call Ham_Latt
 
           ! Setup the hopping / single-particle part
           call Ham_Hop
! 
           ! Setup the interaction.
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
              Write(unit_info,*) 'Model is      :  Portable'
              Write(unit_info,*) '# unit cells  : ', Latt%N
              Write(unit_info,*) '# of orbitals : ', Latt_unit%Norb
!              Write(unit_info,*) 'Checkerboard  : ', Checkerboard
!              Write(unit_info,*) 'Symm. decomp  : ', Symm
!              if (Projector) then
!                 Write(unit_info,*) 'Projective version'
!                 Write(unit_info,*) 'Theta         : ', Theta
!                 Write(unit_info,*) 'Tau_max       : ', beta
!              else
                 Write(unit_info,*) 'Finite temperture version'
                 Write(unit_info,*) 'Beta          : ', Beta
!              endif
              Write(unit_info,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
              Write(unit_info,*) 'N_SUN         : ',   N_SUN
              Write(unit_info,*) 'N_FL          : ', N_FL
!              if (Projector) then
!                 Do nf = 1,N_FL
!                    Write(unit_info,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
!                    Write(unit_info,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
!                 enddo
!              endif
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

        Subroutine read_latt

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif

          implicit none
 
          integer                :: ierr, unit_latt, no, i
          Character (len=64)     :: file_latt



#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g


          If (Irank_g == 0) then
#endif
             File_latt = "geometry.txt"
#if defined(TEMPERING)
             write(File_latt,'(A,I0,A)') "Temp_",igroup,"/geometry.txt"
#endif
             Open(newunit=unit_latt, file=file_latt, status="old", action="read", iostat=ierr)
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'unable to open <geometry.txt>', ierr
                Call Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
             END IF

             read(unit_latt,*) L1, L2
             read(unit_latt,*) a1_p
             read(unit_latt,*) a2_p
             read(unit_latt,*) a3_p
             read(unit_latt,*) norb

             allocate(Orb_a_p(Norb,3))

             do no = 1, Norb
                read(unit_latt,*) i, Orb_a_p(no,:)
             enddo

             Close(unit_latt)
#ifdef MPI
          Endif

          CALL MPI_BCAST(L1       ,  1           ,MPI_INTEGER  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(L2       ,  1           ,MPI_INTEGER  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(a1_p     ,  size(a1_p,1),MPI_REAL8    ,0,Group_Comm,ierr)
          CALL MPI_BCAST(a2_p     ,  size(a2_p,1),MPI_REAL8    ,0,Group_Comm,ierr)
          CALL MPI_BCAST(a3_p     ,  size(a3_p,1),MPI_REAL8    ,0,Group_Comm,ierr)
          CALL MPI_BCAST(Norb     ,  1           ,MPI_INTEGER  ,0,Group_Comm,ierr)
      
          if (irank_g /= 0) allocate(Orb_a_p(Norb,3))
          CALL MPI_BCAST(Orb_a_p,  size(Orb_a_p) ,MPI_REAL8    ,0,Group_Comm,ierr)
#endif


        end subroutine read_latt

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the  Lattice
!--------------------------------------------------------------------
        Subroutine Ham_Latt


          Implicit none
          Real (Kind=Kind(0.d0))  :: d1_p(2), d2_p(2), L1_p(2), L2_p(2)
          integer :: no, ns, nc, i
          Real (Kind = Kind(0.d0) ) :: Zero = 1.0E-8


          latt_unit%Norb = N_spin*Norb
          allocate(latt_unit%Orb_pos_p(latt_unit%Norb,3))
          nc = 0
          ! numbering of orbitals
          ! latt_unit%Orb_pos_p(i       ,:) for i = 1, Norb : "spatial" orbitals for first  spin kind (up  )
          ! latt_unit%Orb_pos_p(Norb + i,:) for i = 1, Norb :                        second           (down)
          do ns = 1, N_spin
             do no = 1, Norb
                nc = nc + 1
                latt_unit%Orb_pos_p(nc,:) = Orb_a_p(no,1)*a1_p + Orb_a_p(no,2)*a2_p + (Orb_a_p(no,3)-dble(ns-1))*a3_p
             enddo
          enddo

          d1_p = a1_p(1:2)
          d2_p = a2_p(1:2)
          if (abs(a1_p(3)) > zero .or. abs(a2_p(3)) > zero ) then
             Write(error_unit,*) 'Unit cell vectors three dimensional, but ALF does not support three-dimensional lattices'
             CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
          endif

          L1_p    =  dble(L1)*d1_p
          L2_p    =  dble(L2)*d2_p
          Call Make_Lattice( L1_p, L2_p, d1_p,  d2_p, latt )

          Ndim = latt%N*latt_unit%Norb

          allocate ( list(ndim,2), Invlist(latt%N,latt_unit%Norb) )
          nc = 0
          do i = 1, latt%N
             do no = 1, latt_unit%Norb
                nc = nc + 1
                list(nc,1)    = I
                list(nc,2)    = no
                invlist(i,no) = nc
             enddo
         enddo


         Latt_unit%N_coord = 2

         do nc = 1, size(latt_unit%orb_pos_p,1)
         print *, "latt_unit%Orb_pos_p", latt_unit%Orb_pos_p(nc,:)
         enddo


        end Subroutine Ham_Latt


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the  Lattice
!--------------------------------------------------------------------

        Subroutine read_hop

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif

          implicit none
 
          integer                :: ierr, unit_hop, n_hop, nh, i
          integer                :: no1, no2, s1, s2, n1, n2
          Character (len=64)     :: file_hop
          real (kind=kind(0.d0)) :: x, y
          logical                :: diag


#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g


          If (Irank_g == 0) then
#endif
             File_hop = "hoppings.txt"
#if defined(TEMPERING)
             write(File_hop,'(A,I0,A)') "Temp_",igroup,"/hoppings.txt"
#endif
             Open(newunit=unit_hop, file=file_hop, status="old", action="read", iostat=ierr)
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'unable to open <hoppings.txt>', ierr
                Call Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
             END IF

             read(unit_hop,*) n_hop
             
             allocate( ham_t(n_hop), hop_list(n_hop,6) )
      
             do nh = 1, n_hop
                read(unit_hop,*) i, no1, s1, n1, n2, no2, s2, x, y
                hop_list(nh,1) = no1
                hop_list(nh,2) = s1
                hop_list(nh,3) = no2
                hop_list(nh,4) = s2
                hop_list(nh,5) = n1
                hop_list(nh,6) = n2
                ham_t(nh)      = cmplx( x, y, kind(0.d0))
             enddo

             Close(unit_hop)
#ifdef MPI
          Endif

          CALL MPI_BCAST(n_hop       ,  1             ,MPI_INTEGER  ,0,Group_Comm,ierr)

          if ( irank_g /= 0 ) allocate( ham_t(n_hop), hop_list(n_hop,6) )
          CALL MPI_BCAST(hop_list    ,  size(hop_list),MPI_INTEGER  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_t       ,  size(ham_t)   ,MPI_COMPLEX16,0,Group_Comm,ierr)
 
#endif

         ! if hop_diag(i) = 1 : hopping is     on-site (chemical potential)
         ! if hop_diag(i) = 0 : hopping is not on-site
         allocate( hop_diag(n_hop) )
         hop_diag = 0

         N_diag = 0
         do nh = 1, n_hop
            diag = hop_list(nh,1) == hop_list(nh,3) .and. hop_list(nh,2) == hop_list(nh,4) .and. &
                    &  hop_list(nh,5) == 0 .and. hop_list(nh,6) == 0
            if (diag) then
               hop_diag(nh) = 1
               N_diag = N_diag + 1
            endif
         enddo

         print *, irank_g, n_hop
         do nh = 1, n_hop
            print *, irank_g, hop_list(nh,:), hop_diag(nh)
            print *, irank_g, ham_t(nh)
         enddo

        end subroutine read_hop

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the Hopping
!--------------------------------------------------------------------
        Subroutine Ham_Hop

          Implicit none

          Integer :: nf , I, Ix, Iy, nh, nc, no
          real (kind=kind(0.d0)) :: ham_t_max, Zero = 1.0E-8

          allocate(Hopping_Matrix(N_Fl))

          ham_t_max = 0.d0
          do nh = 1, size(ham_t,1)
             if ( hop_diag(nh) == 0 .and. abs(dble (ham_t(nh))) > ham_t_max ) ham_t_max = abs(dble (ham_t(nh)))
             if ( hop_diag(nh) == 0 .and. abs(aimag(ham_t(nh))) > ham_t_max ) ham_t_max = abs(aimag(ham_t(nh)))
          enddo
          print *, "ham_t_max", ham_t_max

          do nf = 1, N_Fl
             Hopping_Matrix(N_Fl)%N_bonds = 0
             if ( abs(ham_t_max) > zero ) then
                hopping_matrix(nf)%N_bonds = size(hop_list,1) - N_diag
                allocate (hopping_matrix(nf)%List(hopping_matrix(nf)%N_bonds,4) )
                allocate (hopping_matrix(nf)%t(hopping_matrix(nf)%N_bonds) )

                allocate ( hopping_matrix(nf)%T_loc(Latt_unit%Norb) )
                hopping_matrix(nf)%T_loc = cmplx( 0.d0, 0.d0, kind(0.d0) )

                nc = 0
                do nh = 1, size(ham_t,1)
                   if (hop_diag(nh) == 0) then

                      nc = nc + 1
                      hopping_matrix(nf)%T(nc) = - ham_t(nh)
                      ! hopping_matrix(nf)%list(:,1) : orbital 1 = (orbital 1 of potential.txt + 1) + (spin of orbital 1)*(number of "spatial" orbitals)
                      ! hopping_matrix(nf)%list(:,2) : orbital 2 = (orbital 2 of potential.txt + 1) + (spin of orbital 2)*(number of "spatial" orbitals)
                      ! hopping_matrix(nf)%list(:,3) : shift of unit cell with vector a_1
                      ! hopping_matrix(nf)%list(:,4) : shift of unit cell with vector a_2
                      ! orbital I of potential.txt + 1 : orbitals in ALF start at 1 and in potential.txt at 0
                      hopping_matrix(nf)%list(nc,1) = (hop_list(nh,1)+1) + hop_list(nh,2)*Norb
                      hopping_matrix(nf)%list(nc,2) = (hop_list(nh,3)+1) + hop_list(nh,4)*Norb
                      hopping_matrix(nf)%list(nc,3) = hop_list(nh,5)
                      hopping_matrix(nf)%list(nc,4) = hop_list(nh,6)
               
                  else

                     no = (hop_list(nh,1)+1) + hop_list(nh,2)*Norb
                     hopping_matrix(nf)%t_loc(no) = ham_t(nh)

                  endif
                enddo

             endif

             hopping_matrix(nf)%N_Phi =  N_Phi
             hopping_matrix(nf)%Phi_X =  Phi_X
             hopping_matrix(nf)%Phi_Y =  Phi_Y
             hopping_matrix(nf)%Bulk  =  Bulk
             
             do nc = 1, size(hopping_matrix(nf)%T)
                print *, "list, T", nc, hopping_matrix(nf)%list(nc,:), hopping_matrix(nf)%T(nc)
             enddo

             do no = 1, size(hopping_matrix(nf)%t_loc)
                print *, "t_loc", no, hopping_matrix(nf)%t_loc(no)
             enddo

          enddo


          Call  Predefined_Hoppings_set_OPT(Hopping_Matrix,List,Invlist,Latt,  Latt_unit,  Dtau, .false.,  .false. , OP_T )




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

          Integer :: nf, I, nt
          Real (Kind=Kind(0.d0)) :: X

allocate(Op_V(1,N_FL))
do nf = 1,N_FL
   ! Fake hubbard interaction of weight 0.0 (last argument in the following call)
   Call Predefined_Int_U_SUN(  OP_V(1,nf), 1, N_SUN, DTAU, 0.0d0  )
enddo

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
               N = 1;   Filename = "Kin"
             case (2)
               N = 1;   Filename = "Pot"
             case (3)
               N = 1;   Filename = "Part"
             case (4)
               N = 1;   Filename = "Ener"
             case default
               Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
           enddo
 
!           ! Equal time correlators
!           Allocate ( Obs_eq(3) )
!           Do I = 1,Size(Obs_eq,1)
!             select case (I)
!             case (1)
!               Filename = "Green"
!             case (2)
!               Filename = "SpinZ"
!             case (3)
!               Filename = "Den"
!             case default
!               Write(6,*) ' Error in Alloc_obs '
!             end select
!             Nt = 1
!             Channel = '--'
!             Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
!           enddo
! 
!           If (Ltau == 1) then
!             ! Time-displaced correlators
!             Allocate ( Obs_tau(3) )
!             Do I = 1,Size(Obs_tau,1)
!               select case (I)
!               case (1)
!                 Channel = 'P' ; Filename = "Green"
!               case (2)
!                 Channel = 'PH'; Filename = "SpinZ"
!               case (3)
!                 Channel = 'PH'; Filename = "Den"
!               case default
!                 Write(6,*) ' Error in Alloc_obs '
!               end select
!               Nt = Ltrot+1-2*Thtrot
!               If(Projector) Channel = 'T0'
!               Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
!             enddo
!           endif

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
          Complex (Kind=Kind(0.d0)) :: ZP, ZS
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot
          Integer :: I, J, nf
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

          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo


          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Call Predefined_Hoppings_Compute_Kin(Hopping_Matrix,List,Invlist, Latt, Latt_unit, GRC, ZKin)
          Zkin = Zkin* dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf)
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS

          Obs_scal(4)%Obs_vec(1)  =    Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS


          ! Compute equal-time correlations

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

        end Subroutine OBSERT
        
    end submodule ham_Portable_smod
