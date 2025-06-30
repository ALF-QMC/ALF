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
      !#PARAMETERS END#

      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:)
      Integer, allocatable :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell

      integer :: n_spin, Norb

      !Variables for observables
      type (operator_matrix), allocatable :: obs_scal_ph(:), obs_eq_ph(:)

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

           Ltrot = nint(beta/dtau)
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
        Subroutine Ham_Latt


          Implicit none
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          integer :: no, ns, nc, i, L1, L2
          Real (Kind=Kind(0.d0))  :: a_p(3,3)
          real (kind=kind(0.d0)), allocatable :: Orb_pos(:,:)

          ! read in information for lattice from file geometry.txt
          ! L1, L2, Norb, a_p, Orb_pos
          call read_latt(L1, L2, Norb, a_p, Orb_pos, Group_Comm)

          latt_unit%Norb = N_spin*Norb
!          Latt_Unit%N_coord = ?
          a1_p(1) = a_p(1,1); a1_p(2) = a_p(1,2)
          a2_p(1) = a_p(2,1); a2_p(2) = a_p(2,2)
          allocate(latt_unit%Orb_pos_p(latt_unit%Norb,3))
          nc = 0
          ! numbering of orbitals:
          ! latt_unit%Orb_pos_p( (ns-1)*Norb + no       ,:) for no = 1, Norb : "spatial" orbitals for N_spin = ns
          do ns = 1, N_spin
             do no = 1, Norb
                nc = nc + 1
                latt_unit%Orb_pos_p(nc,:) = Orb_pos(no,1)*a_p(1,:) + Orb_pos(no,2)*a_p(2,:) + Orb_pos(no,3)*a_p(3,:)
             enddo
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
          real (kind=kind(0.d0)) :: ham_t_max, Zero = 1.0E-8
          integer                :: N_diag
          integer, allocatable   :: hop_diag(:)
          logical                :: diag
          type (operator_matrix)   :: hopping

          ! read in information for hopping from file hoppings.txt
          call read_hop(hopping, Group_Comm)

          n_hop = size(hopping%g,1)

          ! if hop_diag(i) = 1 : hopping is     on-site (chemical potential)
          ! if hop_diag(i) = 0 : hopping is not on-site
          allocate( hop_diag(n_hop) )
          hop_diag = 0; N_diag = 0
          do nh = 1, n_hop
             diag = hopping%list(nh,1) == hopping%list(nh,5) .and. &
                  & hopping%list(nh,2) == hopping%list(nh,6) .and. &
                  & hopping%list(nh,3) == 0 .and. hopping%list(nh,4) == 0
             if (diag) then
                hop_diag(nh) = 1
                N_diag = N_diag + 1
             endif
          enddo

          allocate(Hopping_Matrix(N_Fl))

          ham_t_max = 0.d0
          do nh = 1, n_hop
             if ( hop_diag(nh) == 0 .and. abs(hopping%g(nh)) > ham_t_max ) ham_t_max = abs(hopping%g(nh))
          enddo

          do nf = 1, N_Fl
             Hopping_Matrix(nf)%N_bonds = 0
             if ( abs(ham_t_max) > zero ) then
                hopping_matrix(nf)%N_bonds = n_hop - N_diag
                allocate (hopping_matrix(nf)%List(hopping_matrix(nf)%N_bonds,4) )
                allocate (hopping_matrix(nf)%t(hopping_matrix(nf)%N_bonds) )

                nc = 0
                do nh = 1, n_hop
                   if (hop_diag(nh) == 0) then

                      ! hopping_matrix(nf)%list(:,1) : orbital 1 = (orbital 1 of potentials.txt + 1) + (spin of orbital 1)*(number of "spatial" orbitals)
                      ! hopping_matrix(nf)%list(:,2) : orbital 2 = (orbital 2 of potentials.txt + 1) + (spin of orbital 2)*(number of "spatial" orbitals)
                      ! hopping_matrix(nf)%list(:,3) : shift of unit cell with vector a_1
                      ! hopping_matrix(nf)%list(:,4) : shift of unit cell with vector a_2
                      ! orbital I of potentials.txt + 1 : orbitals in ALF start at 1 and in potentials.txt at 0
                      nc = nc + 1
                      hopping_matrix(nf)%T(nc) = - hopping%g(nh)
                      hopping_matrix(nf)%list(nc,1) = (hopping%list(nh,1)+1) + hopping%list(nh,2)*Norb
                      hopping_matrix(nf)%list(nc,2) = (hopping%list(nh,5)+1) + hopping%list(nh,6)*Norb
                      hopping_matrix(nf)%list(nc,3) = hopping%list(nh,3)
                      hopping_matrix(nf)%list(nc,4) = hopping%list(nh,4)
               
                   endif
                enddo
             endif

             allocate ( hopping_matrix(nf)%T_loc(Latt_unit%Norb) )
             hopping_matrix(nf)%T_loc = cmplx( 0.d0, 0.d0, kind(0.d0) )
             do nh = 1, n_hop
                if (hop_diag(nh) == 1) then
                   no = (hopping%list(nh,1)+1) + hopping%list(nh,2)*Norb
                   hopping_matrix(nf)%t_loc(no) = - hopping%g(nh)
                endif
             enddo

             hopping_matrix(nf)%N_Phi =  N_Phi
             hopping_matrix(nf)%Phi_X =  Phi_X
             hopping_matrix(nf)%Phi_Y =  Phi_Y
             hopping_matrix(nf)%Bulk  =  Bulk

          enddo

          Call  Predefined_Hoppings_set_OPT(Hopping_Matrix,List,Invlist,Latt,  Latt_unit,  Dtau, .false.,  .false. , OP_T )

        end Subroutine Ham_Hop

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> returns the number of interacting orbitals for every interaction term
!--------------------------------------------------------------------
  
        subroutine get_number_of_orbitals_per_interaction(this)

           implicit none

           type (operator_matrix), allocatable, intent(inout) :: this(:)

           integer :: no, i, j1, j2, mk, n, no1, no2, nc, N_orbitals, x1, x2
           integer, allocatable :: orbitals_tmp(:)

           i = 1

           do no = 1, size(this,1)
              mk = size(this(no)%list,1 )
              allocate( orbitals_tmp(2*mk))
              orbitals_tmp = 0
              nc   = 0
              N_orbitals = 0

              do n = 1, mk
                 nc  = nc + 1
                 j1  = latt%nnlist( i, this(no)%list(n,1), this(no)%list(n,2) )
                 no1 = (this(no)%list(n,3)+1) + this(no)%list(n,4)*Norb
                 x1  = invlist(j1,no1)
                 if (.not. any(orbitals_tmp == x1)) N_orbitals = N_orbitals + 1
                 orbitals_tmp(nc) = x1

                 nc  = nc + 1
                 j2  = latt%nnlist( i, this(no)%list(n,5), this(no)%list(n,6) )
                 no2 = (this(no)%list(n,7)+1) + this(no)%list(n,8)*Norb
                 x2  = invlist(j2,no2)
                 if (.not. any(orbitals_tmp == x2)) N_orbitals = N_orbitals + 1
                 orbitals_tmp(nc) = x2
              enddo

              this(no)%N_orbitals = N_orbitals
              deallocate ( orbitals_tmp )
           enddo
      
        end subroutine get_number_of_orbitals_per_interaction

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> sets the matrix O and corresponding P for a given interaction term
!--------------------------------------------------------------------

        Subroutine set_op_v_matrix(i,op,this)

           implicit none

           integer, intent(in) :: i
           type (operator), intent(inout)      :: op
           type (operator_matrix), intent(in)  :: this

           integer, allocatable :: p(:)
           integer :: i1, no1, j1, x, i2, no2, j2, nc, mk
           integer :: x1, x2, n

           mk = size(this%list,1)
           allocate (p(this%N_orbitals))
           p = 0
           nc = 0
           do n = 1, mk

              i1  = latt%nnlist( i, this%list(n,1), this%list(n,2) )
              no1 = (this%list(n,3)+1) + this%list(n,4)*Norb
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

              i2  = latt%nnlist( i, this%list(n,5), this%list(n,6) )
              no2 = (this%list(n,7)+1) + this%list(n,8)*Norb
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
          type (operator_matrix),   allocatable   :: interaction(:)

          ! read in information for interaction from file potentials.txt
          call read_int(interaction, Group_Comm)
          call get_number_of_orbitals_per_interaction(interaction)

          N_ops = size(interaction,1)
          allocate (OP_V(N_ops*Latt%N,N_Fl))

          nc = 0
          do i = 1, Latt%N
             do no = 1, n_ops
                nc = nc + 1
                do nf = 1, N_Fl
                   call Op_make(OP_V(nc,nf), interaction(no)%N_orbitals)
   
                   call set_op_v_matrix(i,op_v(nc,nf),interaction(no))
                   Op_V(nc,nf)%g = sqrt(cmplx(-dtau*interaction(no)%u, 0.d0, kind(0.d0) ))
                   Op_V(nc,nf)%alpha = interaction(no)%alpha
                   Op_V(nc,nf)%type  = 2
                   call Op_set( OP_V(nc,nf) )
            
                enddo
             enddo
          enddo

        end Subroutine Ham_V

!--------------------------------------------------------------------

        Subroutine read_obs_eq

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif

          implicit none
 
          integer                :: ierr, unit_obs, mk, no, n, N_obs
          integer                :: i, no1, s1, i1, i2, j1, j2, no2, s2
          Character (len=64)     :: file_obs, file
          real (kind=kind(0.d0)) :: x, y


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
             File_obs = "obs_corr_ph.txt"
#if defined(TEMPERING)
             write(File_obs,'(A,I0,A)') "Temp_",igroup,"/obs_corr_ph.txt"
#endif
             Open(newunit=unit_obs, file=file_obs, status="old", action="read", iostat=ierr)
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'unable to open <obs_corr_ph.txt>', ierr
                Call Terminate_on_error(ERROR_FILE_NOT_FOUND,__FILE__,__LINE__)
             END IF
         
             !local correlation functions
             read(unit_obs,*) N_obs
             allocate(obs_eq_ph(N_obs))

             do no = 1, N_obs
                read(unit_obs,*) file, mk
                obs_eq_ph(no)%file = file
                allocate( obs_eq_ph(no)%list(mk,8), obs_eq_ph(no)%g(mk) )

                do n = 1, mk
                   read(unit_obs,*) i, i1, i2, no1, s1, j1, j2, no2, s2, x, y
                   obs_eq_ph(no)%list(n,1) = i1
                   obs_eq_ph(no)%list(n,2) = i2
                   obs_eq_ph(no)%list(n,3) = no1
                   obs_eq_ph(no)%list(n,4) = s1
                   obs_eq_ph(no)%list(n,5) = j1
                   obs_eq_ph(no)%list(n,6) = j2
                   obs_eq_ph(no)%list(n,7) = no2
                   obs_eq_ph(no)%list(n,8) = s2
                   obs_eq_ph(no)%g(n)      = cmplx( x, y, kind(0.d0) )
                enddo
             enddo

         

             Close(unit_obs)

#ifdef MPI
          Endif

          CALL MPI_BCAST(N_obs                     ,  1,MPI_INTEGER  ,0,Group_Comm,ierr)
          if (irank_g /= 0) allocate(obs_eq_ph(N_obs))
          do no = 1, N_obs
             CALL MPI_BCAST(obs_eq_ph(no)%file   , 64,MPI_CHARACTER,0,Group_Comm,ierr)
             if (irank_g == 0) mk = size(obs_eq_ph(no)%list,1)
             CALL MPI_BCAST(mk                   ,  1                       ,MPI_INTEGER  ,0,Group_Comm,ierr)
             if (irank_g /= 0) allocate( obs_eq_ph(no)%list(mk,8), obs_eq_ph(no)%g(mk) )
             CALL MPI_BCAST(obs_eq_ph(no)%list   ,  size(obs_eq_ph(no)%list),MPI_INTEGER  ,0,Group_Comm,ierr)
             CALL MPI_BCAST(obs_eq_ph(no)%g     ,  size(obs_eq_ph(no)%g)  ,MPI_COMPLEX16,0,Group_Comm,ierr)
          enddo


#endif

        end subroutine read_obs_eq

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
          Integer    ::  i, N, Nt, N_obs_scal_ph, N_obs_eq_ph
          Character (len=64) ::  Filename
          Character (len=:), allocatable ::  Channel

           call read_obs_scal_ph(obs_scal_ph,Group_Comm)
           N_obs_scal_ph = size(obs_scal_ph,1)

           ! Scalar observables
           Allocate ( Obs_scal(N_obs_scal_ph+4) )
           Do I = 1, Size(Obs_scal,1)
             if (i <= N_obs_scal_ph   ) then
                N = 1;   Filename = obs_scal_ph(i)%file
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

           call read_obs_eq
           N_obs_eq_ph = size(obs_eq_ph,1)
 
           ! Equal time correlators
           Allocate ( Obs_eq(3+N_obs_eq_ph) )
           Do I = 1,Size(Obs_eq,1)
             if (i <= N_obs_eq_ph) then
                Filename = obs_eq_ph(i)%File
             elseif (i == N_obs_eq_ph+1) then
               Filename = "Green"
             elseif (i == N_obs_eq_ph+2) then
               Filename = "SpinZ"
             elseif (i == N_obs_eq_ph+3) then
               Filename = "Den"
             endif
             Nt = 1
             Channel = '--'
             Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
           enddo
 
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
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Zlocal, ZZ, ZDen
          Integer :: I, J, nf, i1, no_i, i2, no, n, no_1, no_2, j1, no_j, j2, imj, m
          integer :: x1, x2, y1, y2, no_i1, no_i2, no_j1, no_j2, N_obs_scal_ph, N_obs_eq_ph
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

          N_obs_scal_ph = size(obs_scal_ph,1)
          N_obs_eq_ph   = size(obs_eq_ph,  1)

          do no = 1, N_obs_scal_ph
             Zlocal = cmplx(0.d0, 0.d0, kind(0.d0))
      
             do i = 1, Latt%N
                do n = 1, size(obs_scal_ph(no)%list,1)
                   no_1 = (obs_scal_ph(no)%list(n,1)+1) + obs_scal_ph(no)%list(n,2)*Norb
                   i1   = invlist(i,no_1)
               
                   j    = latt%nnlist(i,obs_scal_ph(no)%list(n,3),obs_scal_ph(no)%list(n,4))
                   no_2 = (obs_scal_ph(no)%list(n,5)+1) + obs_scal_ph(no)%list(n,6)*Norb
                   j1   = invlist(j,no_2)

                   Zlocal = Zlocal +       obs_scal_ph(no)%g(n) *GRC(i1,j1,1)
                   Zlocal = Zlocal + conjg(obs_scal_ph(no)%g(n))*GRC(j1,i1,1)
                enddo
             enddo
   
             Obs_scal(no)%Obs_vec(1)  =    Obs_scal(no)%Obs_vec(1) + Zlocal *ZP* ZS
          enddo


          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Call Predefined_Hoppings_Compute_Kin(Hopping_Matrix,List,Invlist, Latt, Latt_unit, GRC, ZKin)
          Zkin = Zkin* dble(N_SUN)
          Obs_scal(N_obs_scal_ph+1)%Obs_vec(1)  =    Obs_scal(N_obs_scal_ph+1)%Obs_vec(1) + Zkin *ZP* ZS

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
          Do no = 1, N_obs_eq_ph
             Obs_eq(no)%N        = Obs_eq(no)%N + 1
             Obs_eq(no)%Ave_sign = Obs_eq(no)%Ave_sign + real(ZS,kind(0.d0))


             do I = 1, Latt%N
                do n = 1, size(obs_eq_ph(no)%list,1)
                   i1    = latt%nnlist(i,obs_eq_ph(no)%list(n,1),obs_eq_ph(no)%list(n,2))
                   no_i1 = (obs_eq_ph(no)%list(n,3)+1) + obs_eq_ph(no)%list(n,4)*Norb
                   x1    = invlist(i1,no_i1)
                   i2    = latt%nnlist(i,obs_eq_ph(no)%list(n,5),obs_eq_ph(no)%list(n,6))
                   no_i2 = (obs_eq_ph(no)%list(n,7)+1) + obs_eq_ph(no)%list(n,8)*Norb
                   x2    = invlist(i2,no_i2)
                   do J = 1, Latt%N
                      imj = Latt%imj(I,J)
                      do m = 1, size(obs_eq_ph(no)%list,1)
                         j1    = latt%nnlist(j,obs_eq_ph(no)%list(m,1),obs_eq_ph(no)%list(m,2))
                         no_j1 = (obs_eq_ph(no)%list(m,3)+1) + obs_eq_ph(no)%list(m,4)*Norb
                         y1    = invlist(j1,no_j1)
                         j2    = latt%nnlist(j,obs_eq_ph(no)%list(m,5),obs_eq_ph(no)%list(m,6))
                         no_j2 = (obs_eq_ph(no)%list(m,7)+1) + obs_eq_ph(no)%list(m,8)*Norb
                         y2    = invlist(j2,no_j2)

                         Zlocal = cmplx(0.d0, 0.d0, kind(0.d0))
                         Zlocal = Zlocal +       obs_eq_ph(no)%g(n) *      obs_eq_ph(no)%g(m) * &
                                    &   (GRC(x1,x2,1)*GRC(y1,y2,1) + GRC(x1,y2,1)*GR(x2,y1,1))
                         Zlocal = Zlocal +       obs_eq_ph(no)%g(n) *conjg(obs_eq_ph(no)%g(m))* &
                                    &   (GRC(x1,x2,1)*GRC(y2,y1,1) + GRC(x1,y1,1)*GR(x2,y2,1))
                         Zlocal = Zlocal + conjg(obs_eq_ph(no)%g(n))*      obs_eq_ph(no)%g(m) * &
                                    &   (GRC(x2,x1,1)*GRC(y1,y2,1) + GRC(x2,y2,1)*GR(x1,y1,1))
                         Zlocal = Zlocal + conjg(obs_eq_ph(no)%g(n))*conjg(obs_eq_ph(no)%g(m))* &
                                    &   (GRC(x2,x1,1)*GRC(y2,y1,1) + GRC(x2,y1,1)*GR(x1,y2,1))
                         
                         Obs_eq(no)%Obs_latt(imj,1,1,1) = Obs_eq(no)%Obs_latt(imj,1,1,1) + Zlocal*ZP*ZS
                      enddo
                   enddo
                   Obs_eq(no)%Obs_latt0(1) = Obs_eq(no)%Obs_latt0(1) +  &
                           &     (obs_eq_ph(no)%g(n)*GRC(x1,x2,1) + conjg(obs_eq_ph(no)%g(n))*GRC(x2,x1,1))*ZP*ZS
                enddo
             enddo

          enddo



          Call Predefined_Obs_eq_Green_measure  ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(N_obs_eq_ph+1) )

          Obs_eq(N_obs_eq_ph+2)%N        = Obs_eq(N_obs_eq_ph+2)%N + 1
          Obs_eq(N_obs_eq_ph+2)%Ave_sign = Obs_eq(N_obs_eq_ph+2)%Ave_sign + real(ZS,kind(0.d0))
          Obs_eq(N_obs_eq_ph+3)%N        = Obs_eq(N_obs_eq_ph+3)%N + 1
          Obs_eq(N_obs_eq_ph+3)%Ave_sign = Obs_eq(N_obs_eq_ph+3)%Ave_sign + real(ZS,kind(0.d0))

          do I = 1, latt%N
             do no_i = 1, Norb
                i1 = invlist(i,no_i)
                i2 = invlist(i,no_i+Norb)
                do j = 1, Latt%N
                   imj  = latt%imj(I,J)
                   do no_j = 1, Norb
                      j1 = invlist(j,no_j)
                      j2 = invlist(j,no_j+Norb)
                      ZZ = (GRC(i1,i1,1)-GRC(i2,i2,1))*(GRC(j1,j1,1)-GRC(j2,j2,1)) + GRC(I1,j1,1)*GR(i1,j1,1) + GRC(I2,j2,1)*GR(i2,j2,1)
                      Obs_eq(N_obs_eq_ph+2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(N_obs_eq_ph+2)%Obs_Latt(imj,1,no_I,no_J) + ZZ*ZP*ZS
                      ZDen = (GRC(i1,i1,1)+GRC(I2,i2,1))*(GRC(j1,j1,1)+GRC(j2,j2,1)) + GRC(I1,j1,1)*GR(i1,j1,1) + GRC(I2,j2,1)*GR(i2,j2,1)
                      Obs_eq(N_obs_eq_ph+3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(N_obs_eq_ph+3)%Obs_Latt(imj,1,no_I,no_J) + ZDen*ZP*ZS
                   enddo
                enddo
                Obs_eq(N_obs_eq_ph+2)%Obs_Latt0(no_I) =  Obs_eq(N_obs_eq_ph+2)%Obs_Latt0(no_I) + (GRC(i1,i1,1)-GRC(i2,i2,1))*ZP*ZS
                Obs_eq(N_obs_eq_ph+3)%Obs_Latt0(no_I) =  Obs_eq(N_obs_eq_ph+3)%Obs_Latt0(no_I) + (GRC(i1,i1,1)+GRC(i2,i2,1))*ZP*ZS
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
