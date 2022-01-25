!  Copyright (C) 2021-2022 The ALF project
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
!> J. Schwab
!>
!> @brief
!> This module defines the  Hamiltonian and observables of several models
!> that correcond to Dirac systems undergoing a nematic quantum phase transition.
!>
!--------------------------------------------------------------------
   submodule (Hamiltonian_main) ham_Nematic_Dirac_mod
     
#ifdef MPI
      use mpi
#endif
      Use Operator_mod
      Use Observables
      Use Lattices_v3
      Use Random_Wrap

      Implicit none

      type, extends(ham_base) :: ham_Nematic_Dirac
      contains
         procedure, nopass :: Ham_Set
         procedure, nopass :: S0
         procedure, nopass :: Delta_S0_global
         procedure, nopass :: Alloc_obs
         procedure, nopass :: Obser
         procedure, nopass :: measure_hist
         procedure, nopass :: ObserT
         procedure, nopass :: Global_move
         procedure, nopass :: Hamiltonian_set_nsigma
#ifdef HDF5
         procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Nematic_Dirac
      
      integer, parameter :: dp=kind(0.d0)  ! double precision

      !#PARAMETERS START# VAR_Nematic_Dirac
      !Integer :: N_SUN = 2       ! SU(N) symmetry
      real(dp) :: dtau = 0.1d0    ! Imaginary time step size
      Character (len=64) :: Global_type = '' ! Type of global update. Possible values: 'Wolff', 'Geo', 'switch', 'flip'
      Integer  :: L1 = 4          ! Size of lattice in a1 direction
      Integer  :: L2 = 4          ! Size of lattice in a2 direction
      Integer  :: Model_vers = 1  ! Version of model. 1: C_2v model, 2: C_4v model
      real(dp) :: ham_t = 1.d0    ! Hopping amplitude of fermions
      real(dp) :: beta = 10.d0    ! Reciprocal temperature
      real(dp) :: Ham_h = 3.d0    ! Ising transverse field
      real(dp) :: Ham_J = 1.d0    ! Ferromagnetic Ising interaction
      real(dp) :: Ham_xi = 1.d0   ! Coupling strength Ising spins <-> fermions
      real(dp) :: Ham_xi2 = 0.d0  ! Static fermion hopping "distortion"
      real(dp) :: Ham_chem = 0.d0 ! Chemical potential
      real(dp) :: Global_J = 1.d0 ! J for proposing global updates
      real(dp) :: Global_h = 3.d0 ! h for proposing global updates
      real(dp) :: Phi_1 = 0.d0    ! Twisted boundary in a1 direction
      real(dp) :: Phi_2 = 0.d0    ! Twisted boundary in a2 direction
      Character (len=64) :: init_type = 'random' ! How to initialize Ising field. Possile values: 'random', 'up', 'down', 'updown'
      !#PARAMETERS END#

#ifdef MPI
      Integer        :: Isize, Irank, irank_g, isize_g, igroup
      Integer        :: STATUS(MPI_STATUS_SIZE)
#endif
      
     Integer                  :: N_coord, Norb
     Type (Lattice),   target :: Latt
     Type (Unit_cell), target :: Latt_unit
     Integer, allocatable     :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell


     real (Kind=Kind(0.d0)) :: Model_sign

     !>    Private variables for observing Z_x_ising
     Real (Kind=Kind(0.d0)) :: eq_x_ising, neq_x_ising
     Integer                         :: nBlub, nBlub2

     !>    Storage for the Ising action
     Real (Kind=Kind(0.d0)) :: DW_Ising_tau(-1:1), DW_Ising_Space(-1:1)
     Integer, allocatable   :: Ising_nnlist(:,:)

     !>    Variables for the Wolff cluster update
     Real (Kind=Kind(0.d0)) :: Wolff_addProb_space, Wolff_addProb_tau
     Integer                :: N_ising
     !>    Variables for the Geometric cluster update
     Real (Kind=Kind(0.d0)) :: Geo_AddProb_space, Geo_AddProb_tau
     Integer :: R_init(2)

     !>    Experimenting
     Integer :: n_global

   contains

   module subroutine Ham_Alloc_Nematic_Dirac()
     allocate(ham_Nematic_Dirac::ham)
   end subroutine Ham_Alloc_Nematic_Dirac

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Nematic_Dirac_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
      Subroutine Ham_Set()
         Implicit none
         Character (len=64) :: file_info
         integer            :: unit_info, ierr
         
#ifdef MPI
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
         call MPI_Comm_rank(Group_Comm, irank_g, ierr)
         call MPI_Comm_size(Group_Comm, isize_g, ierr)
         igroup           = irank/isize_g
#endif

           ! From dynamically generated file "Hamiltonian_Nematic_Dirac_read_write_parameters.F90"
           call read_parameters()

#ifdef MPI
           If (Irank_g == 0 ) then
#endif
#if defined(TEMPERING)
           write(file_info,'(A,I0,A)') "Temp_",igroup,"/info"
#else
           file_info = "info"
#endif
           OPEN(newunit=unit_info, file=file_info, status="unknown", position="append")
           Write(unit_info,*) '====================================='
           Write(unit_info,*) 'Model is            : ', 'Nematic_Dirac', Model_vers
           Write(unit_info,*) 'Global_type         : ', Global_type
           Write(unit_info,*) 'L1                  : ', L1
           Write(unit_info,*) 'L2                  : ', L2
           Write(unit_info,*) 'N_SUN               : ', N_SUN
           Write(unit_info,*) 'ham_t               : ', ham_t
           Write(unit_info,*) 'dtau                : ', dtau
           Write(unit_info,*) 'beta                : ', beta
           Write(unit_info,*) 'Ham_h               : ', Ham_h
           Write(unit_info,*) 'Ham_J               : ', Ham_J
           Write(unit_info,*) 'Ham_xi              : ', Ham_xi
           Write(unit_info,*) 'Ham_xi2             : ', Ham_xi2
           Write(unit_info,*) 'Ham_chem            : ', Ham_chem
           Write(unit_info,*) 'Global_J            : ', Global_J
           Write(unit_info,*) 'Global_h            : ', Global_h
           Write(unit_info,*) 'Phi_1               : ', Phi_1
           Write(unit_info,*) 'Phi_2               : ', Phi_2
           Write(unit_info,*) 'init_type           : ', init_type
           close(unit_info)
#ifdef MPI
          endif
#endif

          Call Ham_latt()

          N_FL = 1

          Call Ham_hop()
          Ltrot = nint(beta/dtau)
          Projector = .false.
          Thtrot = 0
          Symm = .false.

          Call Setup_Ising_action()
          call Ham_V()
      End Subroutine Ham_Set

!=============================================================================
        Subroutine Ham_Latt()
          Implicit none
          !Set the lattice

          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          Integer :: I, nc, no

          select case( Model_vers )
          case( 0, 1, 3 )
            Norb = 2
            N_coord   = 4
            a1_p(1) =  1.D0/sqrt(2.D0)  ; a1_p(2) =  1.D0/sqrt(2.D0)
            a2_p(1) =  1.D0/sqrt(2.D0)  ; a2_p(2) = -1.D0/sqrt(2.D0)
            L1_p    =  dble(L1)*a1_p
            L2_p    =  dble(L2)*a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
            If ( L1 == 1 .or. L2 == 1 ) then
              Write(6,*) ' One dimensional systems not implemented '
              Stop
            endif
          case( 2 )
            Norb = 2
            N_coord   = 4
            a1_p(1) = 1.D0  ; a1_p(2) = 0.D0
            a2_p(1) = 0.D0  ; a2_p(2) = 1.D0
            L1_p    =  dble(L1)*a1_p
            L2_p    =  dble(L2)*a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
            If ( L1 == 1 .or. L2 == 1 ) then
              Write(6,*) ' One dimensional systems not implemented '
              Stop
            endif
          case default
            Write(6,*) "Lattice not yet implemented!"
            Stop
          end select

          ! This is for the orbital structure.
          Ndim = Latt%N*Norb
          Allocate (List(Ndim,2), Invlist(Latt%N,Norb))
          nc = 0
          Do I = 1,Latt%N
             Do no = 1,Norb
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc
             Enddo
          Enddo

        end Subroutine Ham_Latt
        
        
        pure function twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, latt, L1, L2) result(f)
          Implicit none

          Integer      , intent(in) :: i, del_1, del_2, L1, L2
          real(dp)     , intent(in) :: phi_1, phi_2
          type(Lattice), intent(in) :: Latt
          complex(dp) :: f
          
          real(Kind=Kind(0.d0)), parameter :: pi = dacos(-1.d0)
          
          f = cmplx(1.d0,0.d0, kind(0.D0))
          if((latt%list(i,1)+del_1) > L1) then
            f = f * exp(cmplx(0.d0, 2.d0*pi*phi_1, kind(0.d0)))
          endif
          if((latt%list(i,1)+del_1) < 0) then
            f = f * exp(cmplx(0.d0, -2.d0*pi*phi_1, kind(0.d0)))
          endif
          if((latt%list(i,2)+del_2) > L2) then
            f = f * exp(cmplx(0.d0, 2.d0*pi*phi_2, kind(0.d0)))
          endif
          if((latt%list(i,2)+del_2) < 0) then
            f = f * exp(cmplx(0.d0, -2.d0*pi*phi_2, kind(0.d0)))
          endif
        end function twisted_boundary_f

!===================================================================================
        Subroutine Ham_hop()
          Implicit none

          Integer :: I, I1, J1, n, Ncheck, nc, nc1, del_1, del_2
          complex (Kind=Kind(0.d0)) :: t_temp

          Ncheck = 1
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
            Do nc = 1,Ncheck
              Call Op_make(Op_T(nc,n),Ndim)
              select case( Model_vers )
              case( 0 )
                DO I = 1, Latt%N
                  I1 = Invlist(I,1)
                  Do nc1 = 1,N_coord
                    select case (nc1)
                    case (1)
                      J1 = invlist(I,2)
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0) * (1 + ham_xi2)
                    case (2)
                      J1 = invlist(Latt%nnlist(I, 0,1),2)
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0) * (1 - ham_xi2)
                    case (3)
                      J1 = invlist(Latt%nnlist(I,-1,1),2)
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0) * (1 + ham_xi2)
                    case (4)
                      J1 = invlist(Latt%nnlist(I,-1,0),2)
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0) * (1 - ham_xi2)
                    case default
                      Write(6,*) ' Error in  Ham_Hop '
                    end select
                    Op_T(nc,n)%O(I1,J1) = t_temp
                    Op_T(nc,n)%O(J1,I1) = conjg(t_temp)
                  Enddo
                Enddo
              case( 1, 3 )
                DO I = 1, Latt%N
                  I1 = Invlist(I,1)
                  Do nc1 = 1,N_coord
                    select case (nc1)
                    case (1)
                      del_1 = 0
                      del_2 = 0
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0) * (1 + ham_xi2)
                    case (2)
                      del_1 = 0
                      del_2 = 1
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0) * (1 - ham_xi2)
                    case (3)
                      del_1 = -1
                      del_2 = 1
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0) * (1 - ham_xi2)
                    case (4)
                      del_1 = -1
                      del_2 = 0
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0) * (1 + ham_xi2)
                    case default
                      Write(6,*) ' Error in  Ham_Hop '
                    end select
                    t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                         & latt, L1, L2)
                    J1 = invlist(Latt%nnlist(I, del_1, del_2), 2)
                    Op_T(nc,n)%O(I1,J1) = t_temp
                    Op_T(nc,n)%O(J1,I1) = conjg(t_temp)
                  Enddo
                Enddo
              case( 2 )
                DO I = 1, Latt%N
                  I1 = Invlist(I,1)
                  Do nc1 = 1,N_coord
                    select case (nc1)
                    case (1)
                      del_1 = 1
                      del_2 = 0
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    case (2)
                      del_1 = -1
                      del_2 = 0
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    case (3)
                      del_1 = 0
                      del_2 = 1
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    case (4)
                      del_1 = 0
                      del_2 = -1
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    case default
                      Write(6,*) ' Error in  Ham_Hop '
                    end select
                    t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                         & latt, L1, L2)
                    J1 = invlist(Latt%nnlist(I, del_1, del_2), 2)
                    Op_T(nc,n)%O(I1,J1) = t_temp
                    Op_T(nc,n)%O(J1,I1) = conjg(t_temp)
                  Enddo
                Enddo
              case default
                Write(6,*) ' This hopping is not yet implemented '
                Stop
              end select

                Do I = 1,Ndim
                   Op_T(nc,n)%P(i) = i
                   Op_T(nc,n)%O(i,i) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                Enddo
                Op_T(nc,n)%g = -Dtau
                Op_T(nc,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                Call Op_set(Op_T(nc,n))
             enddo
          enddo
        end Subroutine Ham_hop

!===================================================================================
        Subroutine Ham_V()

          Implicit none

          Integer :: nf, I, nc1, del_1, del_2
          complex (Kind=Kind(0.d0)) :: t_temp

          select case( Model_vers )
          case( 0 )
            Allocate(Op_V(Latt%N,N_FL))
            do nf = 1,N_FL
              do I = 1, Latt%N
                call Op_make(Op_V(I,nf),5)
                Op_V(I,nf)%P(1) = Invlist(I,1)
                Do nc1 = 1,N_coord
                  select case (nc1)
                  case (1)
                    Op_V(I,nf)%P(nc1+1) = invlist(I,2)
                    t_temp = -cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                  case (2)
                    Op_V(I,nf)%P(nc1+1) = invlist(Latt%nnlist(I, 0,1),2)
                    t_temp =  cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                  case (3)
                    Op_V(I,nf)%P(nc1+1) = invlist(Latt%nnlist(I,-1,1),2)
                    t_temp = -cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                  case (4)
                    Op_V(I,nf)%P(nc1+1) = invlist(Latt%nnlist(I,-1,0),2)
                    t_temp =  cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                  case default
                    Write(6,*) ' Error in  Ham_V '
                  end select
                  Op_V(I,nf)%O(1   ,nc1+1) = t_temp
                  Op_V(I,nf)%O(nc1+1,1   ) = conjg(t_temp)
                Enddo
                Op_V(I,nf)%g      = cmplx(-dtau*ham_t*ham_xi,0.d0, kind(0.D0))
                Op_V(I,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                Op_V(I,nf)%type   = 1
                Call Op_set( Op_V(I,nf) )
              Enddo
            Enddo
          case( 1 )
            Allocate(Op_V(Latt%N,N_FL))
            do nf = 1,N_FL
              do I = 1, Latt%N
                call Op_make(Op_V(I,nf),5)
                Op_V(I,nf)%P(1) = Invlist(I,1)
                Do nc1 = 1,N_coord
                  select case (nc1)
                  case (1)
                    del_1 = 0
                    del_2 = 0
                    t_temp = -cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                  case (2)
                    del_1 = 0
                    del_2 = 1
                    t_temp =  cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                  case (3)
                    del_1 = -1
                    del_2 = 1
                    t_temp =  cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                  case (4)
                    del_1 = -1
                    del_2 = 0
                    t_temp = -cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                  case default
                    Write(6,*) ' Error in  Ham_V '
                  end select
                  t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                       & latt, L1, L2)
                  Op_V(I,nf)%P(nc1+1) = invlist(Latt%nnlist(I, del_1, del_2), 2)
                  Op_V(I,nf)%O(1   ,nc1+1) = t_temp
                  Op_V(I,nf)%O(nc1+1,1   ) = conjg(t_temp)
                Enddo
                Op_V(I,nf)%g      = cmplx(-dtau*ham_t*ham_xi,0.d0, kind(0.D0))
                Op_V(I,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                Op_V(I,nf)%type   = 1
                Call Op_set( Op_V(I,nf) )
              Enddo
            Enddo
          case( 3 )
            Allocate(Op_V(2*Latt%N,N_FL))
            do nf = 1,N_FL
              do I = 1, Latt%N
                call Op_make(Op_V(I         ,nf),3)
                call Op_make(Op_V(I+Latt%N,nf),3)
                Op_V(I       ,nf)%P(1) = Invlist(I,1)
                Op_V(I+Latt%N,nf)%P(1) = Invlist(I,1)
                Do nc1 = 1,4
                  select case (nc1)
                  case (1)
                    del_1 = 0
                    del_2 = 0
                    Op_V(I,nf)%P(2) = invlist(Latt%nnlist(I, del_1, del_2), 2)
                    t_temp = -cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                         & latt, L1, L2)
                    Op_V(I,nf)%O(1,2) = t_temp
                    Op_V(I,nf)%O(2,1) = conjg(t_temp)
                  case (2)
                    del_1 = 0
                    del_2 = 1
                    Op_V(I+Latt%N,nf)%P(2) = invlist(Latt%nnlist(I, del_1, del_2), 2)
                    t_temp =  cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                         & latt, L1, L2)
                    Op_V(I+Latt%N,nf)%O(1,2) = t_temp
                    Op_V(I+Latt%N,nf)%O(2,1) = conjg(t_temp)
                  case (3)
                    del_1 = -1
                    del_2 = 1
                    Op_V(I,nf)%P(3) = invlist(Latt%nnlist(I, del_1, del_2), 2)
                    t_temp =  cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                         & latt, L1, L2)
                    Op_V(I,nf)%O(1,3) = t_temp
                    Op_V(I,nf)%O(3,1) = conjg(t_temp)
                  case (4)
                    del_1 = -1
                    del_2 = 0
                    Op_V(I+Latt%N,nf)%P(3) = invlist(Latt%nnlist(I, del_1, del_2), 2)
                    t_temp = -cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    t_temp = t_temp * twisted_boundary_f(i, del_1, del_2, phi_1, phi_2, &
                                                         & latt, L1, L2)
                    Op_V(I+Latt%N,nf)%O(1,3) = t_temp
                    Op_V(I+Latt%N,nf)%O(3,1) = conjg(t_temp)
                  case default
                    Write(6,*) ' Error in  Ham_V '
                  end select
                Enddo
                Op_V(I,nf)%g      = cmplx(-dtau*ham_t*ham_xi,0.d0, kind(0.D0))
                Op_V(I,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                Op_V(I,nf)%type   = 1
                Call Op_set( Op_V(I,nf) )
                Op_V(I+Latt%N,nf)%g      = cmplx(-dtau*ham_t*ham_xi,0.d0, kind(0.D0))
                Op_V(I+Latt%N,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                Op_V(I+Latt%N,nf)%type   = 1
                Call Op_set( Op_V(I+Latt%N,nf) )
              Enddo
            Enddo
          case( 2 )
            Allocate(Op_V(Latt%N,N_FL))
            do nf = 1,N_FL
              do I = 1, Latt%N
                call Op_make(Op_V(I,nf),2)
                Op_V(I,nf)%P(1) = invlist(I,1)
                Op_V(I,nf)%P(2) = invlist(I,2)

                Op_V(I,nf)%O(1,2) = cmplx(0.d0, 1.d0, kind(0.D0))
                Op_V(I,nf)%O(2,1) = cmplx(0.d0,-1.d0, kind(0.D0))

                Op_V(I,nf)%g      = cmplx(-dtau*ham_t*ham_xi,0.d0, kind(0.D0))
                Op_V(I,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                Op_V(I,nf)%type   = 1
                Call Op_set( Op_V(I,nf) )
              Enddo
            Enddo
          case default
            Write(6,*) ' This interaction is not yet implemented '
            Stop
          end select
        end Subroutine Ham_V

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
          Real (Kind=Kind(0.d0)), Intent(In) :: Hs_new

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

!===================================================================================
        Subroutine  Hamiltonian_set_nsigma(Initial_field)

          ! The user can set the initial configuration

          Implicit none
         Real (Kind=Kind(0.d0)), allocatable, dimension(:,:), Intent(INOUT) :: Initial_field

          Integer :: I, nt, ierr

          !write(*,*) init_type
          Allocate( Initial_field(Size(OP_V,1), Ltrot) )

          select case (init_type)
          case ('random')
            deallocate( Initial_field )
          case ('up')
            Do nt = 1,Ltrot
              Do I = 1,Size(OP_V,1)
                Initial_field(I,nt)  = 1
              enddo
            enddo
          case ('down')
            Do nt = 1,Ltrot
              Do I = 1,Size(OP_V,1)
                Initial_field(I,nt)  = -1
              enddo
            enddo
          case ('updown')
            Do nt = 1,Ltrot
              Do I = 1,Size(OP_V,1)
                Initial_field(I,nt)  = 1
                If ( I > Latt%N ) Initial_field(I,nt)  = -1
              enddo
            enddo
#ifdef MPI
          case ('part_rand_up')
            If ( mod(Irank_g,2) == 0 ) then
              deallocate( Initial_field )
            else
              Initial_field(:,:)  = 1.d0
            endif
#endif
          case default
            Write(6,*) ' Error in  Hamiltonian_set_random_nsigma'
            STOP
          end select

        end Subroutine Hamiltonian_set_nsigma

!===================================================================================
        include "NematicDirac/Setup_Ising_action.f90"
!===================================================================================
        include "NematicDirac/Alloc_obs.f90"
!========================================================================
        include "NematicDirac/Global.f90"
!========================================================================
        include "NematicDirac/Obser.f90"
!=====================================================
        include "NematicDirac/ObserT.f90"
!==========================================================

   end submodule Ham_Nematic_Dirac_mod
