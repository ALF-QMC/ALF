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
      Use Hamiltonian_Portable_input_mod

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
      !logical              :: Symm         = .true.   ! Whether symmetrization takes place
      !#PARAMETERS END#

      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:)
      Integer, allocatable :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell

      integer :: L1, L2

      !Variables for observables
      type (operator_matrix), allocatable :: obs_scal_ph(:), obs_corr_ph(:)
      type (obs_orbitals), allocatable :: orbitals_scal_ph(:), orbitals_corr_ph(:)

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

           Ltrot = nint(beta/dtau)
           Projector = .false.
!           Thtrot = 0
!           if (Projector) Thtrot = nint(theta/dtau)
!           Ltrot = Ltrot+2*Thtrot

           ! Setup the Bravais lattice
           call Ham_Latt

           ! Setup the hopping / single-particle part
           call Ham_Hop
 
           ! Setup the interaction.
           call Ham_V
 
           ! Setup the trival wave function, in case of a projector approach
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
              Write(unit_info,*) 'Symm. decomp  : ', Symm
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
        Subroutine Ham_Latt


          Implicit none
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          integer :: no, nc, i, norb
          Real (Kind=Kind(0.d0))  :: a_p(3,3)
          real (kind=kind(0.d0)), allocatable :: Orb_pos(:,:)

          ! read in information for lattice from file geometry.txt
          ! L1, L2, Norb, a_p, Orb_pos
          call read_latt(L1, L2, Norb, a_p, Orb_pos, Group_Comm)

          latt_unit%Norb = Norb
!          Latt_Unit%N_coord = ?
          a1_p(1) = a_p(1,1); a1_p(2) = a_p(1,2)
          a2_p(1) = a_p(2,1); a2_p(2) = a_p(2,2)
          allocate(latt_unit%Orb_pos_p(latt_unit%Norb,3))
          nc = 0
          ! numbering of orbitals:
          do no = 1, Norb
             nc = nc + 1
             latt_unit%Orb_pos_p(nc,:) = Orb_pos(no,1)*a_p(1,:) + Orb_pos(no,2)*a_p(2,:) + Orb_pos(no,3)*a_p(3,:)
          enddo
          L1_p    =  dble(L1) * a1_p
          L2_p    =  dble(L2) * a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, latt )

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

          Integer                :: nf , I, Ix, Iy, nh, nc, no, n_hop
          Integer                :: n, i1, no_i, nb, n_1, n_2, no_j, n1, n2
          integer                :: nj1, nj2, j, j1
          real (kind=kind(0.d0)) :: j_p(2), jp_p(2)
          real (kind=kind(0.d0)) :: ham_t_max, Zero = 1.0E-8
          integer                :: N_diag
          integer, allocatable   :: hop_diag(:)
          logical                :: diag
          Complex(Kind=Kind(0.d0))              :: Z
          type (operator_matrix), allocatable   :: hopping(:)

          ! read in information for hopping from file hoppings.txt
          allocate(hopping(N_Fl))
          call read_hop(hopping, Group_Comm)

          allocate(Hopping_Matrix(N_Fl))

          do nf = 1, N_Fl

             n_hop = size(hopping(nf)%g,1)

             ! if hop_diag(i) = 1 : hopping is     on-site (chemical potential)
             ! if hop_diag(i) = 0 : hopping is not on-site
             allocate( hop_diag(n_hop) )
             hop_diag = 0; N_diag = 0
             do nh = 1, n_hop
                diag = hopping(nf)%lattice(nh,1) == hopping(nf)%lattice(nh,4) .and. &
                     & hopping(nf)%lattice(nh,2) == 0 .and. hopping(nf)%lattice(nh,3) == 0
                if (diag) then
                   hop_diag(nh) = 1
                   N_diag = N_diag + 1
                endif
             enddo

             ham_t_max = 0.d0
             do nh = 1, n_hop
                if ( hop_diag(nh) == 0 .and. abs(hopping(nf)%g(nh)) > ham_t_max ) ham_t_max = abs(hopping(nf)%g(nh))
             enddo

             Hopping_Matrix(nf)%N_bonds = 0
             if ( abs(ham_t_max) > zero ) then
                hopping_matrix(nf)%N_bonds = n_hop - N_diag
                allocate (hopping_matrix(nf)%List(hopping_matrix(nf)%N_bonds,4) )
                allocate (hopping_matrix(nf)%t(hopping_matrix(nf)%N_bonds) )

                nc = 0
                do nh = 1, n_hop
                   if (hop_diag(nh) == 0) then

                      nc = nc + 1
                      hopping_matrix(nf)%T(nc) = - hopping(nf)%g(nh)
                      hopping_matrix(nf)%list(nc,1) = hopping(nf)%lattice(nh,1)+1
                      hopping_matrix(nf)%list(nc,2) = hopping(nf)%lattice(nh,4)+1
                      hopping_matrix(nf)%list(nc,3) = hopping(nf)%lattice(nh,2)
                      hopping_matrix(nf)%list(nc,4) = hopping(nf)%lattice(nh,3)
               
                   endif
                enddo
             endif

             allocate ( hopping_matrix(nf)%T_loc(Latt_unit%Norb) )
             hopping_matrix(nf)%T_loc = cmplx( 0.d0, 0.d0, kind(0.d0) )
             do nh = 1, n_hop
                if (hop_diag(nh) == 1) then
                   no = hopping(nf)%lattice(nh,1)+1
                   hopping_matrix(nf)%t_loc(no) = - hopping(nf)%g(nh)
                endif
             enddo

             hopping_matrix(nf)%N_Phi =  N_Phi
             hopping_matrix(nf)%Phi_X =  Phi_X
             hopping_matrix(nf)%Phi_Y =  Phi_Y
             hopping_matrix(nf)%Bulk  =  Bulk

             deallocate(hop_diag)

          enddo

!          Call  Predefined_Hoppings_set_OPT(Hopping_Matrix,List,Invlist,Latt,  Latt_unit,  Dtau, .false.,  symm , OP_T )

          select case (inquire_hop(Hopping_Matrix))
          case(0)  !  Zero
             allocate(Op_T(1,N_FL))
             do nf = 1,N_FL
                Call Op_make(Op_T(1,nf),1)
                Op_T(1,nf)%P(1)   = 1
                Op_T(1,nf)%O(1,1) = cmplx(0.d0,0.d0, kind(0.d0))
                Op_T(1,nf)%g      = 0.d0
                Op_T(1,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                Call Op_set(Op_T(1,nf))
             enddo
          case(1)  ! Diagonal
             allocate(Op_T(Ndim,N_FL))
             do nf = 1,N_FL
                do n = 1,ndim
                   Call Op_make(Op_T(n,nf),1)
                   Op_T(n,nf)%P(1)   = n
                   Op_T(n,nf)%O(1,1) =  Hopping_Matrix(nf)%T_Loc(list(n,2))
                   Op_T(n,nf)%g      = -Dtau
                   Op_T(n,nf)%alpha  =  cmplx(0.d0,0.d0, kind(0.D0))
                   Call Op_set(Op_T(n,nf))
                enddo
             enddo
          case default
             allocate(Op_T(1,N_FL))
             do nf = 1,N_FL
                !Write(6,*)
                Call Op_make(Op_T(1,nf),Ndim)   ! This is too restrictive for the Kondo type models. The hopping only occurs on one subsystem.
                N_Phi     = Hopping_Matrix(nf)%N_Phi
                Phi_X     = Hopping_Matrix(nf)%Phi_X
                Phi_Y     = Hopping_Matrix(nf)%Phi_Y
                Bulk      = Hopping_Matrix(nf)%Bulk
                !Write(6,*) N_Phi, Phi_X,Phi_Y, Bulk
                !Write(6,*) This(nf)%list
                DO I = 1, Latt%N
                   do Nb = 1, Hopping_Matrix(nf)%N_bonds
                      no_I = Hopping_Matrix(nf)%list(Nb,1)
                      no_J = Hopping_Matrix(nf)%list(Nb,2)
                      n_1  = Hopping_Matrix(nf)%list(Nb,3)
                      n_2  = Hopping_Matrix(nf)%list(Nb,4)
                      nj1  = latt%list(i,1) + n_1
                      nj2  = latt%list(i,2) + n_2
                      j_p  = dble(nj1)*latt%a1_p + dble(nj2)*latt%a2_p
                      N1 = 0; N2 = 0
                      Call npbc(jp_p, j_p, Latt%L1_p, Latt%L2_p,  N1, N2)
                      nj1 = nj1 - N1*L1
                      nj2 = nj2 - N2*L2
                      J    = Latt%invlist(nj1,nj2)
!                      print *, "1", n_1, latt%list(i,1), n1, nj1, j
!                      print *, "2", n_2, latt%list(i,2), n2, nj2, j
                      Z    = Generic_hopping(I,no_I, n_1, n_2, no_J, N_Phi, Phi_x,Phi_y, Bulk, Latt, Latt_Unit)
                      I1   = Invlist(I,no_I)
                      J1   = Invlist(J,no_J)
                      Op_T(1,nf)%O(I1,J1) = Hopping_Matrix(nf)%T(Nb)*Z
                      Op_T(1,nf)%O(J1,I1) = Conjg(Hopping_Matrix(nf)%T(Nb)*Z)
                   enddo
                   ! T(N_b=1..N_bonds)
                   ! List(N_b,1) = no_1
                   ! List(N_b,2) = no_2
                   ! List(N_b,3) = n_1
                   ! List(N_b,4) = n_2
                   ! H_[(i,no_1),(i + n_1 a_1 + n_2 a_2,no_2)] = T(N_b)
                   Do no_I = 1, Latt_Unit%Norb
                      I1   = Invlist(I,no_I)
                      Op_T(1,nf)%O(I1,I1) = Hopping_Matrix(nf)%T_Loc(no_I)
                   Enddo
                enddo
                Do I = 1,Ndim
                   Op_T(1,nf)%P(i) = i
                Enddo
                Op_T(1,nf)%g = -Dtau
                Op_T(1,nf)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                Call Op_set(Op_T(1,nf))
                !Do I = 1,Size(Op_T(1,nf)%E,1)
                !   Write(6,*) Op_T(1,nf)%E(I)
                !Enddo
             Enddo
          end select

        end Subroutine Ham_Hop

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> returns the number of interacting orbitals for every interaction term
!--------------------------------------------------------------------
  
        subroutine get_number_of_orbitals_per_interaction(this,orbitals)

           implicit none

           type (operator_matrix), intent(in) :: this(:,:)
           integer, intent(inout) :: orbitals(:,:)

           integer :: no, i, j1, j2, mk, n, no1, no2, nc, N_orbitals, x1, x2, nf
           integer, allocatable :: orbitals_tmp(:)

           i = 1

           do no = 1, size(this,1)
              do nf = 1, size(this,2)
                 mk = size(this(no,nf)%lattice,1 )
                 allocate( orbitals_tmp(2*mk))
                 orbitals_tmp = 0
                 nc   = 0
                 N_orbitals = 0

                 do n = 1, mk
                    nc  = nc + 1
                    j1  = latt%nnlist( i, this(no,nf)%lattice(n,1), this(no,nf)%lattice(n,2) )
                    no1 = (this(no,nf)%lattice(n,3)+1)
                    x1  = invlist(j1,no1)
                    if (.not. any(orbitals_tmp == x1)) N_orbitals = N_orbitals + 1
                    orbitals_tmp(nc) = x1

                    nc  = nc + 1
                    j2  = latt%nnlist( i, this(no,nf)%lattice(n,4), this(no,nf)%lattice(n,5) )
                    no2 = (this(no,nf)%lattice(n,6)+1)
                    x2  = invlist(j2,no2)
                    if (.not. any(orbitals_tmp == x2)) N_orbitals = N_orbitals + 1
                    orbitals_tmp(nc) = x2
                 enddo

                 orbitals(no,nf) = N_orbitals
                 deallocate ( orbitals_tmp )
              enddo
           enddo
      
        end subroutine get_number_of_orbitals_per_interaction

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> sets the matrix O and corresponding P for a given interaction term
!--------------------------------------------------------------------

        Subroutine set_op_v_matrix(i,op,this,orbital)

           implicit none

           integer, intent(in) :: i, orbital
           type (operator), intent(inout)      :: op
           type (operator_matrix), intent(in)  :: this

           integer, allocatable :: p(:)
           integer :: i1, no1, j1, x, i2, no2, j2, nc, mk
           integer :: x1, x2, n

           mk = size(this%lattice,1)
           allocate (p(orbital))
           p = 0
           nc = 0
           do n = 1, mk

              i1  = latt%nnlist( i, this%lattice(n,1), this%lattice(n,2) )
              no1 = (this%lattice(n,3)+1)
              j1  = invlist(i1,no1)
              x1  = -1
              do x = 1, nc
                 if (p(x) == j1) then
                    x1 = x
                    exit
                 endif
              enddo
              if (x1 == -1) then
                 nc = nc + 1
                 p(nc) = j1
                 x1 = nc
              endif

              i2  = latt%nnlist( i, this%lattice(n,4), this%lattice(n,5) )
              no2 = (this%lattice(n,6)+1)
              j2  = invlist(i2,no2)
              x2  = -1
              do x = 1, nc
                 if (p(x) == j2) then
                    x2 = x
                    exit
                 endif
              enddo
              if (x2 == -1) then
                 nc = nc + 1
                 p(nc) = j2
                 x2 = nc
              endif

              Op%O(x1,x2) =        this%g(n)
              Op%O(x2,x1) = CONJG( this%g(n) )

           enddo

           Op%P = p
           deallocate(p)

        end subroutine set_op_v_matrix

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

          Integer :: nf, I, nt, no, nc, i1, i2, no1, no2, n, j1, j2
          integer :: nc_o, x1, x2, x, N_ops
          type (operator_matrix),   allocatable   :: interaction(:,:)
          integer, allocatable :: orbitals(:,:)

          ! read in information for interaction from file potentials.txt
          call read_int(interaction, N_Fl, Group_Comm)
          N_ops = size(interaction,1)
          allocate (OP_V(N_ops*Latt%N,N_Fl), orbitals(N_ops,N_Fl))

          call get_number_of_orbitals_per_interaction(interaction,orbitals)

          nc = 0
          do i = 1, Latt%N
             do no = 1, n_ops
                nc = nc + 1
                do nf = 1, N_Fl
                   call Op_make(OP_V(nc,nf), orbitals(no,nf))
   
                   call set_op_v_matrix(i,op_v(nc,nf),interaction(no,nf),orbitals(no,nf))
                   Op_V(nc,nf)%g = sqrt(cmplx(-dtau*interaction(no,1)%u, 0.d0, kind(0.d0) ))
                   Op_V(nc,nf)%alpha = interaction(no,nf)%alpha
                   Op_V(nc,nf)%type  = 2
                   call Op_set( OP_V(nc,nf) )
            
                enddo
             enddo
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
          Integer    ::  i, N, Nt, N_obs_scal_ph, N_obs_corr_ph
          Character (len=64) ::  Filename
          Character (len=:), allocatable ::  Channel
          Character (len=64), allocatable ::  File(:)
          logical :: corr

           call read_obs_scal_ph(obs_scal_ph,file,Group_Comm)
           N_obs_scal_ph = size(obs_scal_ph,1)
           corr = .false.
           call set_obs_orbitals(obs_scal_ph,orbitals_scal_ph,corr)

           ! Scalar observables
           Allocate ( Obs_scal(N_obs_scal_ph+4) )
           Do I = 1, Size(Obs_scal,1)
             if (i <= N_obs_scal_ph   ) then
                N = 1;   Filename = file(i)
             elseif (i == N_obs_scal_ph+1 ) then
                N = 1;   Filename = "Kin"
             elseif (i == N_obs_scal_ph+2 ) then
                N = 1;   Filename = "Pot"
             elseif (i == N_obs_scal_ph+3 ) then
                N = 1;   Filename = "Part"
             elseif (i == N_obs_scal_ph+4 ) then
                N = 1;   Filename = "Ener"
             endif
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
           enddo

           deallocate(file)

           call read_obs_corr(obs_corr_ph,file,Group_Comm)
           N_obs_corr_ph = size(obs_corr_ph,1)
           corr = .true.
           call set_obs_orbitals(obs_corr_ph,orbitals_corr_ph,corr)
 
           ! Equal time correlators
           Allocate ( Obs_eq(1+N_obs_corr_ph) )
           Do I = 1,Size(Obs_eq,1)
             if (i <= N_obs_corr_ph) then
                Filename = file(i)
             elseif (i == N_obs_corr_ph+1) then
                Filename = "Green"
             endif
             Nt = 1
             Channel = '--'
             Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
           enddo
 
           If (Ltau == 1) then
             ! Time-displaced correlators
             Allocate ( Obs_tau(1+N_obs_corr_ph) )
             Do I = 1,Size(Obs_tau,1)
                if (i <= N_obs_corr_ph) then
                   Channel = 'PH'; Filename = file(i)
                elseif (i == N_obs_corr_ph+1) then
                   Channel = 'P' ; Filename = "Green"
                endif
               Nt = Ltrot+1-2*Thtrot
               If(Projector) Channel = 'T0'
               Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             enddo
           endif

           deallocate(file)

        End Subroutine Alloc_obs


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> creates list of interacting orbitals for each observable and lattice site
!--------------------------------------------------------------------

        Subroutine  set_obs_orbitals(obs,orbitals,corr)

           implicit none

           type (operator_matrix), intent(in)            :: obs(:)
           type (obs_orbitals), intent(out), allocatable :: orbitals(:)
           logical, intent(in) :: corr

           integer :: i1, no_i1, i2, no_i2, i, no, n, x1, x2, mk, m
           integer :: ns1, ns2, ns3, ns4, nf1, nf2, nf3, nf4

           allocate(orbitals(size(obs,1)))

           do no = 1, size(obs,1)

             mk = size(obs(no)%lattice,1)
             allocate(orbitals(no)%orb_list(latt%N,mk,2))

              do n = 1, mk
                 do i = 1, Latt%N

                    i1    = latt%nnlist(i,obs(no)%lattice(n,1),obs(no)%lattice(n,2))
                    no_i1 = obs(no)%lattice(n,3)+1
                    x1    = invlist(i1,no_i1)

                    i2    = latt%nnlist(i,obs(no)%lattice(n,4),obs(no)%lattice(n,5))
                    no_i2 = obs(no)%lattice(n,6)+1
                    x2    = invlist(i2,no_i2)

                    orbitals(no)%orb_list(i,n,1) = x1
                    orbitals(no)%orb_list(i,n,2) = x2

                 enddo
              enddo

              if (corr) then
                 allocate(orbitals(no)%nonzero_corr(mk,mk,2),orbitals(no)%nonzero_back(mk))
                 orbitals(no)%nonzero_corr = 0
                 orbitals(no)%nonzero_back = 0
                 do n = 1, mk
                    nf1 = obs_corr_ph(no)%flavor(n,1)
                    ns1 = obs_corr_ph(no)%color (n,1)
                    nf2 = obs_corr_ph(no)%flavor(n,2)
                    ns2 = obs_corr_ph(no)%color (n,2)

                    if (nf1 == nf2 .and. ns1 == ns2) orbitals(no)%nonzero_back(n) = 1
                    do m = 1, mk
                       nf3 = obs_corr_ph(no)%flavor(m,1)
                       ns3 = obs_corr_ph(no)%color (m,1)
                       nf4 = obs_corr_ph(no)%flavor(m,2)
                       ns4 = obs_corr_ph(no)%color (m,2)
                         
                       if ( (nf1 == nf2 .and. ns1 == ns2) .and. (nf3 == nf4 .and. ns3 == ns4) ) orbitals(no)%nonzero_corr(n,m,1) = 1
                       if ( (nf1 == nf4 .and. ns1 == ns4) .and. (nf2 == nf3 .and. ns2 == ns3) ) orbitals(no)%nonzero_corr(n,m,2) = 1
                    enddo
                 enddo
              endif

           enddo

        end subroutine set_obs_orbitals

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
          Complex (Kind=Kind(0.d0)) :: ZP, ZS, ZN, gn, gm
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Zlocal, ZZ, ZDen, Zeq
          Integer :: I, J, nf, i1, no_i, i2, no, n, no_1, no_2, j1, no_j, j2, imj, m
          integer :: x1, x2, y1, y2, no_i1, no_i2, no_j1, no_j2, nf1, nf2, ns1, ns2
          integer :: nf3, nf4, x3, x4, mk
          integer :: N_obs_scal_ph, N_obs_corr_ph
          ! Add local variables as needed

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

          ZS = ZS*Mc_step_weight
          ZN = cmplx(dble(N_SUN), 0.d0, kind(0.d0))
          
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

          N_obs_scal_ph = size(obs_scal_ph,1)
          N_obs_corr_ph   = size(obs_corr_ph,  1)

          do no = 1, N_obs_scal_ph
             Zlocal = cmplx(0.d0, 0.d0, kind(0.d0))
      
             mk = size(obs_scal_ph(no)%lattice,1)
             do n = 1, mk
                nf   = obs_scal_ph(no)%flavor(n,1)
                do i = 1, Latt%N
                   x1 = orbitals_scal_ph(no)%orb_list(i,n,1)
                   x2 = orbitals_scal_ph(no)%orb_list(i,n,2)

                   Zlocal = Zlocal +   obs_scal_ph(no)%g(n) *GRC(x1,x2,nf)
                enddo
             enddo
   
             Obs_scal(no)%Obs_vec(1)  =    Obs_scal(no)%Obs_vec(1) + Zlocal * dble(N_SUN) * ZP*ZS
          enddo

          !Additional observables like Kin, Pot, SpinZ, Green, ... for debugging
!          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
!          Call Predefined_Hoppings_Compute_Kin(Hopping_Matrix,List,Invlist, Latt, Latt_unit, GRC, ZKin)
!          Zkin = Zkin* dble(N_SUN)
!          Obs_scal(N_obs_scal_ph+1)%Obs_vec(1)  =    Obs_scal(N_obs_scal_ph+1)%Obs_vec(1) + Zkin *ZP* ZS

          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          Do I = 1, size(op_v,1)
             Zpot = Zpot + Predefined_Obs_V_Int(OP_V(i,:), GR, GRC, N_SUN )*(-Op_V(i,1)%g**2/dtau)
          Enddo
          Obs_scal(N_obs_scal_ph+2)%Obs_vec(1)  =  Obs_scal(N_obs_scal_ph+2)%Obs_vec(1) + Zpot * ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf)
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(N_obs_scal_ph+3)%Obs_vec(1)  =    Obs_scal(N_obs_scal_ph+3)%Obs_vec(1) + Zrho * ZP*ZS

          Obs_scal(N_obs_scal_ph+4)%Obs_vec(1)  =    Obs_scal(N_obs_scal_ph+4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS


          ! Compute equal-time correlations
          Do no = 1, N_obs_corr_ph
             Obs_eq(no)%N        = Obs_eq(no)%N + 1
             Obs_eq(no)%Ave_sign = Obs_eq(no)%Ave_sign + real(ZS,kind(0.d0))

             mk = size(obs_corr_ph(no)%lattice,1)
             do n = 1, mk
                nf1 = obs_corr_ph(no)%flavor(n,1)
                nf2 = obs_corr_ph(no)%flavor(n,2)
                gn = obs_corr_ph(no)%g(n)

                do I = 1, Latt%N
                   x1  = orbitals_corr_ph(no)%orb_list(i,n,1)
                   x2  = orbitals_corr_ph(no)%orb_list(i,n,2)

                   do m = 1, mk
                      nf3 = obs_corr_ph(no)%flavor(m,1)
                      nf4 = obs_corr_ph(no)%flavor(m,2)
                      gm = obs_corr_ph(no)%g(m)

                      do J = 1, Latt%N
                         imj = Latt%imj(I,J)
                         x3 = orbitals_corr_ph(no)%orb_list(j,m,1)
                         x4 = orbitals_corr_ph(no)%orb_list(j,m,2)

                         Zeq = cmplx(0.d0, 0.d0, kind(0.d0))
                         if (orbitals_corr_ph(no)%nonzero_corr(n,m,1) == 1) then
                            Zeq = Zeq + gn*gm * GRC(x1,x2,nf1)*GRC(x3,x4,nf3)
                         endif
                         if (orbitals_corr_ph(no)%nonzero_corr(n,m,2) == 1) then
                            Zeq = Zeq + gn*gm * GRC(x1,x4,nf1)*GR (x2,x3,nf2)
                         endif
                         Obs_eq(no)%Obs_latt(imj,1,1,1) = Obs_eq(no)%Obs_latt(imj,1,1,1) + Zeq*ZP*ZS
                      enddo
                   enddo
                   if (orbitals_corr_ph(no)%nonzero_back(n) == 1) then
                      Obs_eq(no)%Obs_latt0(1) = Obs_eq(no)%Obs_latt0(1) + gn * GRC(x1,x2,nf1) *ZP*ZS
                   endif
                enddo
             enddo
          enddo

          Call Predefined_Obs_eq_Green_measure  ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(N_obs_corr_ph+1) )

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
          Complex (Kind=Kind(0.d0)) :: ZP, ZS, Zden, ZZ, Ztau, ZN, dx1_x2, dx3_x4
          Complex (kind=Kind(0.d0)) :: gn, gm
          Integer :: N_obs_corr_ph, i1, i2, j1, j2, no_i, no_j, i, j, imj, no, ns
          Integer :: n, no_i1, x1, no_i2, x2, m, no_j1, y1, no_j2, y2, nf1, nf2
          Integer :: x3, x4, nf3, nf4, ns1, ns2, mk
          ! Add local variables as needed

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS * Mc_step_weight
          ZN = cmplx(dble(N_SUN), 0.d0, kind(0.d0))

          N_obs_corr_ph = size(obs_corr_ph,1)
   
          ! Compute observables
          Do no = 1, N_obs_corr_ph
             if (nt == 0) then
                Obs_tau(no)%N        = Obs_tau(no)%N + 1
                Obs_tau(no)%Ave_sign = Obs_tau(no)%Ave_sign + real(ZS,kind(0.d0))
             endif

             mk = size(obs_corr_ph(no)%lattice,1)
             do n = 1, mk
                nf1 = obs_corr_ph(no)%flavor(n,1)
                nf2 = obs_corr_ph(no)%flavor(n,2)
                gn  = obs_corr_ph(no)%g(n)
                do I = 1, Latt%N
                   x1 = orbitals_corr_ph(no)%orb_list(i,n,1)
                   x2 = orbitals_corr_ph(no)%orb_list(i,n,2)
                   if (x1 == x2) then
                      dx1_x2 = cmplx(1.d0, 0.d0, kind(0.d0))
                   else
                      dx1_x2 = cmplx(0.d0, 0.d0, kind(0.d0))
                   endif


                   do m = 1, mk
                      nf3 = obs_corr_ph(no)%flavor(m,1)
                      nf4 = obs_corr_ph(no)%flavor(m,2)
                      gm  = obs_corr_ph(no)%g(m)
                      do J = 1, Latt%N
                         imj = Latt%imj(I,J)
                         x3 = orbitals_corr_ph(no)%orb_list(j,m,1)
                         x4 = orbitals_corr_ph(no)%orb_list(j,m,2)
                         if (x3 == x4) then
                            dx3_x4 = cmplx(1.d0, 0.d0, kind(0.d0))
                         else
                            dx3_x4 = cmplx(0.d0, 0.d0, kind(0.d0))
                         endif

                         Ztau = cmplx(0.d0, 0.d0, kind(0.d0))
                         if (orbitals_corr_ph(no)%nonzero_corr(n,m,1) == 1) then
                            Ztau = Ztau + gn *gm * ((dx1_x2 - GTT(x2,x1,nf1))*(dx3_x4 - G00(x4,x3,nf3)))
                         endif
                         if (orbitals_corr_ph(no)%nonzero_corr(n,m,2) == 1) then
                            Ztau = Ztau - gn * gm * G0T(x4,x1,nf1)*GT0(x2,x3,nf2)
                         endif
                         Obs_tau(no)%Obs_latt(imj,nt+1,1,1) = Obs_tau(no)%Obs_latt(imj,nt+1,1,1) + Ztau*ZP*ZS

                      enddo
                   enddo
                   if (orbitals_corr_ph(no)%nonzero_back(n) == 1) then
                      Obs_tau(no)%Obs_latt0(1) = Obs_tau(no)%Obs_latt0(1) + gn*(dx1_x2 - GTT(x2,x1,nf1)) *ZP*ZS
                   endif
                enddo
             enddo
          enddo

          call Predefined_Obs_tau_Green_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, &
                  &   Obs_tau(N_obs_corr_ph+1) )

        end Subroutine OBSERT

        
    end submodule ham_Portable_smod
