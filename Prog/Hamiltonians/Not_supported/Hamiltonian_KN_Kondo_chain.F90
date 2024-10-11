!  Copyright (C) 2016 - 2018 The ALF project
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
!> @param [public]  nsigma(:,:)
!> \verbatim Type(Fields)
!> Array containing all auxiliary fields. The first index runs through the operator sequence. The second
!> through the time slies.   \endverbatim
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

    Module Hamiltonian

      Use Operator_mod
      Use WaveFunction_mod
      Use Lattices_v3
      Use MyMats
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables
      Use Fields_mod
      Use LRC_Mod
      use iso_fortran_env, only: output_unit, error_unit


      Implicit none


      Type (Operator),     dimension(:,:), allocatable :: Op_V
      Type (Operator),     dimension(:,:), allocatable :: Op_T
      Type (WaveFunction), dimension(:),   allocatable :: WF_L
      Type (WaveFunction), dimension(:),   allocatable :: WF_R
      Type  (Fields)       :: nsigma
      Integer              :: Ndim
      Integer              :: N_FL
      Integer              :: N_SUN
      Integer              :: Ltrot
      Integer              :: Thtrot
      Logical              :: Projector
      Integer              :: Group_Comm
      Logical              :: Symm = .false.


      Type (Lattice),       private, target :: Latt
      Type (Lattice),       private :: Latt_unit
      Type (Unit_cell),     private, target :: Latt_unit_obs
      Integer,              private :: Norb, N_coord
      Integer,              private :: K,L
      Integer,              private :: Lx,Ly
      real (Kind=Kind(0.d0)),        private :: ham_T , ham_U, Ham_J, Ham_Jh
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta, Theta
      Integer, allocatable, private :: List(:,:), Invlist(:,:)          ! For orbital structure of Unit cell
      Character (len=64),   private :: Model, Lattice_Type
      Logical,              private :: K_space, Log_scale, Lin_scale
      real (Kind=Kind(0.d0)),        private :: W = 2.d0, Lambda
      real (Kind=Kind(0.d0)),        private, allocatable ::  eps(:), delta_eps(:), g(:)


!>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau



    contains

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

          integer            :: ierr, N_part, nf,i, L1, L2
          Character (len=64) :: file_info, file_para
          Real      (Kind=Kind(0.d0)) :: Delta
          logical            :: Checkerboard, Symm



          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model, N_SUN, Checkerboard, Symm


          NAMELIST /VAR_KN_Kondo/ L, K, N_SUN, ham_T,  Ham_U,  Ham_J, Ham_Jh,  Dtau, Beta, Theta, Projector, &
               &                  K_Space, Log_scale, Lin_scale, Lambda


#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
#endif

#ifdef MPI
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g

#endif
          ! Open files
#if defined(MPI)
          If (Irank_g == 0 ) then
#endif
             File_para = "parameters"
             File_info = "info"
#if defined(TEMPERING)
             write(File_para,'(A,I0,A)') "Temp_",igroup,"/parameters"
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif

             OPEN(UNIT=5,FILE=file_para,STATUS='old',ACTION='read',IOSTAT=ierr)
             OPEN(Unit = 50,file=file_info,status="unknown",position="append")
#ifdef MPI
          Endif
#endif

          ! Default values
          Projector = .false.
          K_Space   = .false.
          Log_scale = .false.
          Lambda    = 2.d0
          Theta     = 0.d0
          Thtrot    = 0
          N_FL      = 1
          L1        = 1
          L2        = 1
          N_SUN     = 2
          Checkerboard = .false.
          Symm         = .false.

#ifdef MPI
          If (Irank_g == 0 ) then
#endif
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'Ham_Set: unable to open <parameters>',ierr
                error stop 1
             END IF
             READ(5,NML=VAR_Lattice)
             READ(5,NML=VAR_KN_Kondo)
             Lx = L1
             Ly = L2

             Ltrot = nint(beta/dtau)
             if (Projector) Thtrot = nint(theta/dtau)
             Ltrot = Ltrot+2*Thtrot
             Write(50,*) '====================================='
             Write(50,*) 'Model is      :  Kondo Impurity  K Channels SU(N) '
             Write(50,*) 'Chain Length  : ', L
             Write(50,*) 'Channels  K   : ', K
             Write(50,*) 'Lx            : ', Lx
             Write(50,*) 'Ly            : ', Ly
             Write(50,*) 'SU(N)         : ', N_SUN
             Write(50,*) 'Flavors       : ', N_FL
             Write(50,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'Ham_U         : ', Ham_U
             Write(50,*) 'Ham_J  (Kondo): ', Ham_J
             Write(50,*) 'Ham_JH (Heis) : ', Ham_Jh
             If (Log_scale) then
                Write(50,*) 'Logscale       : '
                Write(50,*) 'Band Width     : ', W
                Write(50,*) 'Lambda         : ', Lambda
             Endif
             If (Lin_scale) then
                Write(50,*) 'Lin_scale      : '
                Write(50,*) 'Band Width     : ', W
                Write(50,*) 'Delta_E        : ', W/real(L,kind(0.d0))
             Endif
             if (Projector) then
                Write(50,*) 'Projective version'
                Write(50,*) 'Theta         : ', Theta
                Write(50,*) 'Tau_max       : ', beta
             else
                Write(50,*) 'Finite temperture version'
                Write(50,*) 'Beta          : ', Beta
             endif


#ifdef MPI
          Endif
          CALL MPI_BCAST(L           ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Lx          ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Ly          ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(K           ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(N_SUN       ,1  ,MPI_INTEGER,   0,Group_Comm,ierr)
          CALL MPI_BCAST(Ltrot       ,1  ,MPI_INTEGER,0,Group_Comm,ierr)
          CALL MPI_BCAST(Thtrot      ,1  ,MPI_INTEGER,0,Group_Comm,ierr)
          CALL MPI_BCAST(Projector   ,1  ,MPI_LOGICAL,0,Group_Comm,ierr)
          CALL MPI_BCAST(K_Space     ,1  ,MPI_LOGICAL,0,Group_Comm,ierr)
          CALL MPI_BCAST(Log_scale   ,1  ,MPI_LOGICAL,0,Group_Comm,ierr)
          CALL MPI_BCAST(Lin_scale   ,1  ,MPI_LOGICAL,0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_T       ,1  ,MPI_REAL8  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_J       ,1  ,MPI_REAL8  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_Jh      ,1  ,MPI_REAL8  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_U       ,1  ,MPI_REAL8  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(Dtau        ,1  ,MPI_REAL8  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(Lambda      ,1  ,MPI_REAL8  ,0,Group_Comm,ierr)
          CALL MPI_BCAST(Beta        ,1  ,MPI_REAL8  ,0,Group_Comm,ierr)
#endif



          Call Ham_Latt

          If ( Log_scale .or. Lin_scale ) Then
             !Here you can allocate the energies
             allocate ( eps(L), delta_eps(L), g(L) )
             Do I = 1,L
                g(i)  = real(L,kind(0.d0)) / ( W * 2.D0 )
              Enddo
             if (Log_scale) then
                Do I = 1,L/2
                   eps(I) =  - (lambda**( 1 -i) ) * W/2.d0
                Enddo
                Do I = L/2+1,L
                   eps(I) =   (lambda**(i-L)) * W/2.d0
                Enddo
             elseif (Lin_scale) then
                Delta = W/real(L,kind(0.d0))
                Do I = 1,L/2
                   eps(I) =  -W/2.d0   + real(i-1,kind(0.d0))*Delta
                Enddo
                Do I = L/2+1,L
                   eps(I) =   Delta*real( I - L/2,Kind(0.d0))
                Enddo
             else
                Write(error_unit,*) 'Ham_Set: No valid scale'
                error stop 1
             endif
             Do I = 1,L/2-1
                delta_eps(I) = abs(eps(I+1) - eps(I))
             Enddo
             delta_eps(L/2)   = abs(eps(L/2))
             delta_eps(L/2+1) = abs(eps(L/2))
             Do I = L/2+2,L
                delta_eps(I) = abs(eps(I-1) - eps(I))
             Enddo
#ifdef MPI
             If (Irank_g == 0 ) then
#endif
                Open (Unit=10,file="Discretizazion",status="unknown")
                Write(10,*) 'Disretization'
                Do i = 1,L
                   Write(10,*)  i, eps(i), delta_eps(i)
                enddo
                close(10)
#ifdef MPI
             endif
#endif
          Endif


          Call Ham_Hop


          if (Projector) then
             Call Ham_Trial

#ifdef MPI
             If (Irank_g == 0) then
#endif
                Do nf = 1,N_FL
                   Write(50,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                   Write(50,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
                enddo

#ifdef MPI
             Endif
#endif

          endif

#ifdef MPI
          If (Irank_g == 0 )  then
#endif
             close(50)
             Close(5)
#ifdef MPI
          endif
#endif

          call Ham_V


        end Subroutine Ham_Set
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the lattice
!--------------------------------------------------------------------
        Subroutine Ham_Latt

          Implicit none

          Integer ::  nc, I, n_L, n_K, no

          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2), x_p(2)



          ! Setup Latt_unit

          a1_p(1) =  1.0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0  ; a2_p(2) =  1.d0
          L2_p    =  real(L,Kind(0.d0))*a2_p
          L1_p    =  real(K,Kind(0.d0))*a1_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt_unit )
          Norb =  Latt_unit%N + 1

          ! Setup the lattice
          L1_p    =  Real(Lx,kind(0.d0))*a1_p
          L2_p    =  Real(Ly,Kind(0.d0))*a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
          If (  Latt%N == 1  ) then
             N_coord = 0 ! 0-dimensional
          elseif (  Latt%N   > 1 .and. Ly == 1  ) then
             N_coord = 1 ! 1-dimensional
          elseif (  Latt%N   > 1 .and. Ly > 1 .and. Lx > 1  ) then
             N_coord = 2 ! 2-dimensional
          else
             Write(error_unit,*) 'Ham_Latt: Error in Lx, Ly  (L1,L2) input'
             error stop 1
          endif
          !Write(6,*)  'N_coord is : ', N_coord
          Ndim = Latt%N*Norb

          Allocate (List(Ndim,2), Invlist(Latt%N,Norb))
          nc = 0
          Do I = 1,Latt%N
             ! Bath
             Do n_K = 1,K
                Do N_L = 1,L
                   x_p =  real(n_K,Kind(0.d0))*Latt_unit%a1_p + &
                        & real(n_L,Kind(0.d0))*Latt_unit%a2_p
                   no = Inv_R(x_p,Latt_unit)  ! Orbital
                   nc = nc + 1
                   List(nc,1)    = I ! Unit cell
                   List(nc,2)    = no
                   Invlist(I,no) = nc
                Enddo
             Enddo
             ! Impurity
             no = Latt_unit%N + 1
             nc = nc + 1
             List(nc,1)    = I ! Unit cell
             List(nc,2)    = no
             Invlist(I,no) = nc
          Enddo
          
          Latt_unit_obs%Norb = 1
          Latt_unit_obs%N_coord = 2
          allocate(Latt_unit_obs%Orb_pos_p(Latt_unit_obs%Norb,2))
          Latt_unit_obs%Orb_pos_p(1, :) = [0.d0, 0.d0]

        End Subroutine Ham_Latt

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hopping
!--------------------------------------------------------------------
        Subroutine Ham_Hop

          Implicit none

          Real (Kind=Kind(0.d0))  :: X_p(2), Delta_K
          Integer                 :: N_ch, I, Ip1, nc, i_0, no


          If  (Log_scale .or. Lin_scale) then
             allocate(Op_T(K*L*Latt%N,1))
             nc = 0
             Do i_0 = 1,Latt%N
                Do n_ch = 1,K
                   Do I = 1,L
                      nc = nc + 1
                      Call Op_make(Op_T(nc,1),1)
                      Op_T(nc,1)%O(1,1) = cmplx(eps(i)*delta_eps(i)*g(i), 0.d0, kind(0.D0))
                      X_p =  Latt_unit%a1_p* real(n_ch,kind(0.d0)) + Latt_unit%a2_p * real(i,kind(0.d0))
                      no  =  Inv_R(x_p,Latt_unit)
                      Op_T(nc,1)%P(1) =  Invlist(i_0,no)
                      !Write(6,*) i_0, no,  Op_T(nc,1)%P(1), Op_T(nc,1)%O(1,1)
                      Op_T(nc,1)%g      = -Dtau
                      Op_T(nc,1)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                      Call Op_set( Op_T(nc,1) )
                   Enddo
                   !Write(6,*)
                Enddo
             Enddo
          elseif (K_space) then
             Delta_K = 2.d0*acos(-1.d0)/real(L,Kind=Kind(0.d0))
             allocate(Op_T(K*L*Latt%N,1))
             nc = 0
             Do i_0 = 1,Latt%N
                Do n_ch = 1,K
                   Do I = 1,L
                      nc = nc + 1
                      Call Op_make(Op_T(nc,1),1)
                      Op_T(nc,1)%O(1,1) = cmplx( -2.d0*Ham_T*cos(Delta_K*real(I-1,Kind=Kind(0.d0))), &
                           &                     0.d0, kind(0.D0))
                      X_p =  Latt_unit%a1_p* real(n_ch,kind(0.d0)) + Latt_unit%a2_p * real(i,kind(0.d0))
                      no = Inv_R(x_p,Latt_unit)
                      Op_T(nc,1)%P(1) =  Invlist(i_0,no)
                      !Write(6,*) nc,  Op_T(nc,1)%P(1), Op_T(nc,1)%O(1,1)
                      Op_T(nc,1)%g      = -Dtau
                      Op_T(nc,1)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                      Call Op_set( Op_T(nc,1) )
                   Enddo
                   !Write(6,*)
                Enddo
             Enddo
          else
             allocate(Op_T(K*Latt%N,1))
             nc = 0
             do i_0 = 1,Latt%N
                Do n_ch = 1,K
                   nc = nc + 1
                   Call Op_make(Op_T(nc,1),L)
                   DO I = 1, L
                      Ip1 = I + 1
                      If (I == L ) Ip1 = 1
                      Op_T(n_ch,1)%O(I,Ip1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      Op_T(n_ch,1)%O(Ip1,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      X_p =  Latt_unit%a1_p* real(n_ch,kind(0.d0)) + Latt_unit%a2_p * real(I,kind(0.d0))
                      no  =  Inv_R(x_p,Latt_unit)
                      Op_T(n_ch,1)%P(I) =  Invlist(i_0,no)
                      !Write(6,*) i_0, n_ch, Op_T(n_ch,1)%P(I)
                   Enddo
                   !Write(6,*)
                   Op_T(n_ch,1)%g      = -Dtau
                   Op_T(n_ch,1)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                   Call Op_set( Op_T(n_ch,1) )
                Enddo
             Enddo
          endif
        End Subroutine Ham_Hop


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the interaction
!--------------------------------------------------------------------
        Subroutine Ham_V

          Implicit none

          Integer ::  n, nc, n_ch, I_c, I_f, I, i_0, i_1, I_f1
          Real (Kind=Kind(0.d0)) :: X, X_p(2)


          If (Log_scale .or. Lin_scale) then
             Allocate(Op_V((2*K + 1 + 2*N_coord)*Latt%N,1) )

             !  2*K*Latt%N        : for Kondo  (Kinetic and current terms) couplings with the  bath
             !  Latt%N            : for Hubbard on impurity sites
             !  2* N_coord*Latt%N : for Heisenberg coupling between the  impurity spins.

             !  You will still have to add the interaction between the
             !  the sites of the the Lattice Latt. There are Latt%N * N_coord*2  (if you take
             !  into account kinietic energy and current) of them.
             !  Note that for the SU(2) code, you could get rid of the current.

             nc = 0
             do i_0 = 1,Latt%N
                do n = 1,2*K
                   nc = nc + 1
                   Call Op_make(Op_V(nc,1),L+1)    ! Kondo
                enddo
                nc = nc + 1
                Call Op_make(Op_V (nc,1),1)        ! Hubbard
                do n = 1,2*N_coord
                   nc = nc + 1
                   Call Op_make(Op_V(nc,1),2  )    ! Kondo
                enddo

             enddo

             nc = 0
             do i_0 = 1, Latt%N
                I_f  = invlist(i_0, Latt_unit%N + 1)
                do n_ch = 1,K
                   ! Kinetic
                   nc = nc + 1
                   Op_V(nc,1)%P(L+1) = I_f
                   Do I = 1, L
                      X_p =  Latt_unit%a1_p* real(n_ch,kind(0.d0)) + Latt_unit%a2_p * real(i,kind(0.d0))
                      I_c  = Inv_R(x_p,Latt_unit)
                      Op_V(nc,1)%P(I ) = invlist(i_0,I_c)
                      !Write(6,*)  "Coupling: ", I_c, L
                      Op_V(nc,1)%O(I,L+1) = cmplx(g(i)*Delta_eps(i)  ,0.d0, kind(0.D0))
                      Op_V(nc,1)%O(L+1,I) = cmplx(g(i)*Delta_eps(i)  ,0.d0, kind(0.D0))
                   Enddo
                   Op_V(nc,1)%alpha  = cmplx(0.d0  ,0.d0, kind(0.d0) )
                   Op_V(nc,1)%g      = cmplx(sqrt(2.d0*Dtau*Ham_J/(4.d0*Real(L*N_SUN,Kind(0.d0)))),0.d0, kind(0.D0))
                   Op_V(nc,1)%type   = 2
                   Call Op_set( Op_V(nc,1) )
                   !Write(6,*) 'Log Kinetic non_zero:' , Op_V(nc,1)%N, Op_V(nc,1)%N_non_zero
                   ! Current
                   nc = nc + 1
                   Op_V(nc,1)%P(L+1) = I_f
                   Do I = 1, L
                      X_p =  Latt_unit%a1_p*real(n_ch,kind(0.d0)) +  Latt_unit%a2_p*real(I,kind(0.d0))
                      I_c  = Inv_R(x_p,Latt_unit)
                      !Write(6,*)  "Coupling: ", I_c, L
                      Op_V(nc,1)%P(I  ) = Invlist(i_0,I_c)
                      Op_V(nc,1)%O(I,L+1) = cmplx(0.d0  ,  g(i)*Delta_eps(i), kind(0.D0))
                      Op_V(nc,1)%O(L+1,I) = cmplx(0.d0  , -g(i)*Delta_eps(i), kind(0.D0))
                   Enddo
                   Op_V(nc,1)%alpha  = cmplx(0.d0  ,0.d0, kind(0.d0) )
                   Op_V(nc,1)%g      = cmplx(sqrt(2.d0*Dtau*Ham_J/(4.d0*Real(L*N_SUN,Kind(0.d0)))),0.d0, kind(0.D0))
                   Op_V(nc,1)%type   = 2
                   Call Op_set( Op_V(nc,1) )
                   !Write(6,*) 'Log Current non_zero:' , Op_V(nc,1)%N_non_zero
                Enddo
                ! Hubbard
                nc = nc + 1
                Op_V(nc,1)%P(1) =  I_f
                Op_V(nc,1)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                Op_V(nc,1)%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
                Op_V(nc,1)%g      = SQRT(CMPLX(-DTAU*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0)))
                Op_V(nc,1)%type   = 2
                Call Op_set( Op_V(nc,1) )

                ! Heisenberg
                do n = 1,N_coord
                   if (  n == 1 )  i_1 = latt%nnlist(i_0,1,0)
                   if (  n == 2 )  i_1 = latt%nnlist(i_0,0,1)
                   I_f1 = invlist(i_1, Latt_unit%N + 1)
                   !Kinetic
                   nc = nc + 1
                   Op_V(nc,1)%P(1) = I_f
                   Op_V(nc,1)%P(2) = I_f1
                   Op_V(nc,1)%O(1,2) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   Op_V(nc,1)%O(2,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   Op_V(nc,1)%alpha  = cmplx(0.d0  ,0.d0, kind(0.d0) )
                   Op_V(nc,1)%g      = sqrt(cmplx (Dtau*Ham_Jh/(4.d0*Real(N_SUN,Kind(0.d0))),0.d0, kind(0.D0) ) )
                   Op_V(nc,1)%type   = 2
                   Call Op_set( Op_V(nc,1) )
                   !Current
                   nc = nc + 1
                   Op_V(nc,1)%P(1) = I_f
                   Op_V(nc,1)%P(2) = I_f1
                   Op_V(nc,1)%O(1,2) = cmplx(0.d0 ,  1.d0, kind(0.D0))
                   Op_V(nc,1)%O(2,1) = cmplx(0.d0 , -1.d0, kind(0.D0))
                   Op_V(nc,1)%alpha  = cmplx(0.d0 ,  0.d0, kind(0.d0) )
                   Op_V(nc,1)%g      = sqrt(cmplx(Dtau*Ham_Jh/(4.d0*Real(N_SUN,Kind(0.d0))),0.d0, kind(0.D0) ) )
                   Op_V(nc,1)%type   = 2
                   Call Op_set( Op_V(nc,1) )
                enddo
             enddo
          elseif (K_space) then

             If (N_coord > 0 ) then
                write(error_unit,*) 'Ham_V: Heisenberg and K_space is not yet implemented'
                error stop 1
             endif
             Allocate(Op_V((2*K + 1 )*Latt%N,1))

             nc = 0
             do i_0 = 1,Latt%N
                do n = 1,2*K
                   nc = nc + 1
                   Call Op_make(Op_V(nc,1),L+1)    ! Kondo
                enddo
                nc = nc + 1
                Call Op_make(Op_V (nc,1),1)        ! Hubbard
             enddo


             nc = 0
             do i_0 = 1,Latt%N
                I_f  = invlist(i_0, Latt_unit%N + 1)
                do n_ch = 1,K
                   ! Kinetic
                   !Write(6,*) 'Kinetic non_zero:' , Op_V(nc,1)%N_non_zero
                   ! Current
                   nc = nc + 1
                   Op_V(nc,1)%P(L+1) = I_f
                   Do I = 1, L
                      X_p =  Latt_unit%a1_p*real(n_ch,kind(0.d0)) +  Latt_unit%a2_p*real(I,kind(0.d0))
                      I_c  = Inv_R(x_p,Latt_unit)
                      !Write(6,*)  "Coupling: ", I_c, L
                      Op_V(nc,1)%P(I  ) = Invlist(i_0,I_c)
                      Op_V(nc,1)%O(I,L+1) = cmplx(0.d0  ,  1.d0, kind(0.D0))
                      Op_V(nc,1)%O(L+1,I) = cmplx(0.d0  , -1.d0, kind(0.D0))
                   Enddo
                   Op_V(nc,1)%alpha  = cmplx(0.d0  ,0.d0, kind(0.d0) )
                   Op_V(nc,1)%g      = cmplx(sqrt(Dtau*Ham_J/(4.d0*Real(L*N_SUN,Kind(0.d0)))),0.d0, kind(0.D0))
                   Op_V(nc,1)%type   = 2
                   Call Op_set( Op_V(nc,1) )
                   !Write(6,*) 'Current non_zero:' , Op_V(nc,1)%N_non_zero
                Enddo
                nc = nc + 1
                Op_V(nc,1)%P(1) =  I_f
                Op_V(nc,1)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                Op_V(nc,1)%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
                Op_V(nc,1)%g      = SQRT(CMPLX(-DTAU*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0)))
                Op_V(nc,1)%type   = 2
                Call Op_set( Op_V(nc,1) )
             enddo
          else

             If (N_coord > 0 ) then
                write(error_unit,*) 'Ham_V: Heisenberg and default real space  is not yet implemented'
                error stop 1
             endif

             Allocate(Op_V((2*K + 1)*Latt%N,1))
             nc = 0
             do i_0 = 1,Latt%N
                do n = 1,2*K
                   nc = nc + 1
                   Call Op_make(Op_V(nc,1),2)    ! Kondo
                enddo
                nc = nc + 1
                Call Op_make(Op_V (nc,1),1)   ! Hubbard
             enddo
             nc = 0
             do i_0 = 1,Latt%N
                I_f  = invlist(i_0,Latt_unit%N + 1)
                do n_ch = 1,K
                   X_p =  Latt%a1_p*real(n_ch,kind(0.d0))
                   I_c  = Inv_R(x_p,Latt_unit)
                   !Write(6,*)  "Coupling: ", I_c
                   ! Kinetic
                   nc = nc + 1
                   Op_V(nc,1)%P(1) = invlist(i_0,I_c)
                   Op_V(nc,1)%P(2) = I_f
                   Op_V(nc,1)%O(1,2) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   Op_V(nc,1)%O(2,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   Op_V(nc,1)%alpha  = cmplx(0.d0  ,0.d0, kind(0.d0) )
                   Op_V(nc,1)%g      = cmplx(sqrt(Dtau*Ham_J/(4.d0*Real(N_SUN,Kind(0.d0)))),0.d0, kind(0.D0))
                   Op_V(nc,1)%type   = 2
                   Call Op_set( Op_V(nc,1) )
                   ! Current
                   nc = nc + 1
                   Op_V(nc,1)%P(1) = invlist(i_0,I_c)
                   Op_V(nc,1)%P(2) = I_f
                   Op_V(nc,1)%O(1,2) = cmplx(0.d0  ,  1.d0, kind(0.D0))
                   Op_V(nc,1)%O(2,1) = cmplx(0.d0  , -1.d0, kind(0.D0))
                   Op_V(nc,1)%alpha  = cmplx(0.d0  ,0.d0, kind(0.d0) )
                   Op_V(nc,1)%g      = cmplx(sqrt(Dtau*Ham_J/(4.d0*Real(N_SUN,Kind(0.d0)))),0.d0, kind(0.D0))
                   Op_V(nc,1)%type   = 2
                   Call Op_set( Op_V(nc,1) )
                Enddo
                nc = nc + 1
                Op_V(nc,1)%P(1) =  I_f
                Op_V(nc,1)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                Op_V(nc,1)%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
                Op_V(nc,1)%g      = SQRT(CMPLX(-DTAU*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0)))
                Op_V(nc,1)%type   = 2
                Call Op_set( Op_V(nc,1) )
             enddo
          Endif

        end Subroutine Ham_V

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the trial wave function
!--------------------------------------------------------------------

        Subroutine Ham_Trial

#if defined (MPI) || defined(TEMPERING)
          Use mpi
#endif
          Implicit none
          Integer  :: n, n_part

          Real (Kind=Kind(0.d0)), allocatable ::  Hop_trial(:,:),  U(:,:),  E(:)

          Integer                    :: I, I_f, I1_f,  nc, no, I1, I2, n_L, n_K, nf
          Real ( Kind = Kind(0.d0) ) :: Trial_T = 0.5,  x_p(2),  Trial_V = 0.5,  Trial_Jh = 0.1
          Complex (Kind=Kind(0.d0) ) :: Z_norm
#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup, Ierr
          Integer        :: STATUS(MPI_STATUS_SIZE)
#endif

#ifdef MPI
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g

#endif
          Allocate (Hop_trial(Ndim,Ndim), U(Ndim,Ndim), E(Ndim) )

          Hop_trial = 0.d0

          Do I = 1, Latt%N
             I_f = invlist(I,Norb)
             Do n_K = 1,K
                do n_L = 1,L
                   x_p =  real(n_K,Kind(0.d0))*Latt_unit%a1_p + &
                        & real(n_L,Kind(0.d0))*Latt_unit%a2_p
                   no = Inv_R(x_p,Latt_unit)
                   I1 = Invlist(I,   no)
                   Hop_trial(I1,I1)   =  delta_eps(n_L) * g(n_L) * eps(n_L)
                   Hop_trial(I_f,I1)  =  g(n_L) * delta_eps(n_L)*sqrt(2.d0/real(L,kind(0.d0)) )
                   Hop_trial(I1,I_f)  =  g(n_L) * delta_eps(n_L)*sqrt(2.d0/real(L,kind(0.d0)) )
                enddo
             enddo
          enddo

          Do I = 1, Latt%N
             I_f = invlist(I,Norb)
             Do n = 1,N_coord
                select case (n)
                case (1)
                   I1_f  = invlist( latt%nnlist(I,1,0), Norb)
                case (2)
                   I1_f  = invlist( latt%nnlist(I,0,1), Norb)
                case default
                   Write(error_unit,*) 'Error in Ham_Trial'
                   error stop 1
                end select
                Hop_trial(I_f , I1_f) = - Ham_Jh
                Hop_trial(I1_f, I_f ) = - Ham_Jh
             enddo
          enddo
          Call Diag(Hop_trial, U, E)

#if defined(MPI)
          If (Irank_g == 0 ) then
#endif
             Open (Unit=10,file="Discretizazion",status="unknown", position="Append")
             Write(10,*) 'K,L:', K,L
             Do nc = 1,Ndim
                Write(10,*) nc, E(nc)
             enddo
             Close(10)
#if defined(MPI)
          Endif
#endif

          N_part = Ndim/2  !  Half-filling.
          if (Abs(2 * N_part - Ndim) > 0 ) then
             Write(error_unit,*) 'Ham_Trial: The particle number is not integer. Ndim, Npart: ', Ndim, N_part
             error stop 1
          endif
          Allocate(WF_L(N_FL),WF_R(N_FL))
          do n=1,N_FL
             Call WF_alloc(WF_L(n),Ndim,N_part)
             Call WF_alloc(WF_R(n),Ndim,N_part)
          enddo

          Do nf = 1,N_FL
             do I2=1,N_part
                do I1=1,Ndim
                   WF_L(nf)%P(I1,I2)=U(I1,I2)
                   WF_R(nf)%P(I1,I2)=U(I1,I2)
                enddo
             enddo
             WF_L(nf)%Degen = E(N_part+1) - E(N_part)
             WF_R(nf)%Degen = E(N_part+1) - E(N_part)
          enddo

          Do nf = 1,N_FL
             Call WF_overlap(WF_L(nf), WF_R(nf), Z_norm)
             Write(6,*) " Z_norm ", Z_norm
          enddo

          Deallocate (Hop_trial, U, E )

        End Subroutine Ham_Trial
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


          ! Scalar observables
          Allocate ( Obs_scal(2) )
          Do I = 1,Size(Obs_scal,1)
             select case (I)
             case (1)
                N = 1;   Filename ="Part"
             case (2)
                N = 1;   Filename ="Double"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo


          ! Equal time correlators
          Allocate ( Obs_eq(2) )
          Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Filename ="Spin"
             case (2)
                Filename ="Den"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Nt = 1
             Channel = "--"
             Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_obs, Channel, dtau)
          enddo

          If (Ltau == 1) then
             ! Equal time correlators
             Allocate ( Obs_tau(2) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Channel = 'PH'; Filename ="Spin"
                case (2)
                   Channel = 'PH'; Filename ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = Ltrot+1-2*Thtrot
                If(Projector) Channel = 'T0'
                Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_obs, Channel, dtau)
             enddo
          endif
        End Subroutine Alloc_obs
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Global moves
!>
!> @details
!>  This routine generates a
!>  global update  and returns the propability T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)
!> @param [IN] nsigma_old,  Type(Fields)
!> \verbatim
!>  Old configuration. The new configuration is stored in nsigma.
!> \endverbatim
!> @param [OUT]  T0_Proposal_ratio Real
!> \verbatimam
!>  T0_Proposal_ratio  =  T0( sigma_new -> sigma_old ) /  T0( sigma_old -> sigma_new)
!> \endverbatim
!> @param [OUT]  Size_clust Real
!> \verbatim
!>  Size of cluster that will be flipped.
!> \endverbatim
!-------------------------------------------------------------------
        ! Functions for Global moves.  These move are not implemented in this example.
        Subroutine Global_move(T0_Proposal_ratio,nsigma_old,size_clust)

          Implicit none
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
          Type (Fields),  Intent(IN)  :: nsigma_old

          ! Local
          Integer :: N_op, N_tau, n1,n2, n

          T0_Proposal_ratio  = 1.d0
          size_clust         = 3.d0

          nsigma%f = nsigma_old%f
          nsigma%t = nsigma_old%t
          N_op  = size(nsigma_old%f,1)
          N_tau = size(nsigma_old%f,2)
          Do n  = 1, Nint(size_clust)
             n1 = nranf(N_op)
             n2 = nranf(N_tau)
             nsigma%f(n1,n2) = nsigma_old%flip(n1,n2)
          enddo

        End Subroutine Global_move
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
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none

          ! Arguments
          Type (Fields),  INTENT(IN) :: nsigma_old

          ! Local
          Integer :: I,n,n1,n2,n3,n4,nt,nt1, nc_F, nc_J, nc_h_p, nc_h_m


          Delta_S0_global = 1.d0


        end Function Delta_S0_global

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
        subroutine Obser(GR,Phase,Ntau)

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau

          !Local
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, ZZ, ZXY, Zdouble
          Integer :: I,J, imj, nf, dec, I1, J1, no_I, no_J,n, I_f,  i_0

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))


          Do nf = 1,N_FL
             Do I = 1,Ndim
                Do J = 1,Ndim
                   GRC(I, J, nf) = -GR(J, I, nf)
                Enddo
                GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
             Enddo
          Enddo
          ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >

          ! Compute scalar observables.
          Do I = 1,Size(Obs_scal,1)
             Obs_scal(I)%N         =  Obs_scal(I)%N + 1
             Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
          Enddo

          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf)
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zrho * ZP*ZS

          Zdouble = cmplx(0.d0,0.d0,Kind(0.d0))
          Z       = cmplx(dble(N_SUN), 0.d0, kind(0.D0))
          nf      = 1
          Do I_0     = 1,Latt%N
             I_f  =  invlist(I_0, Latt_unit%N + 1)
             Zdouble = Zdouble +  (Z* (  Grc(i_f,i_f,nf) - cmplx(0.5D0,0.d0,Kind(0.d0))) )**2 + &
                  &                Z * Grc(i_f,i_f,nf)*Gr(i_f,i_f,nf)
          Enddo
          Zdouble = Zdouble / Z
          Obs_scal(2)%Obs_vec(1)  =    Obs_scal(2)%Obs_vec(1) + Zdouble * ZP*ZS


          ! Compute spin-spin, Green, and den-den correlation functions  !  This is general N_SUN, and  N_FL = 1
          DO I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO
          Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
          Do I1 = 1,Latt%N
             I = invlist(I1, Latt_unit%N + 1)
             Do J1 = 1,Latt%N
                J    = invlist(J1, Latt_unit%N + 1)
                imj = Latt%imj(I1,J1)
                ! SpinT
                Obs_eq(1)%Obs_Latt(imj,1,1,1) =  Obs_eq(1)%Obs_Latt(imj,1,1,1) + &
                     &               Z * GRC(I,J,1) * GR(I,J,1) * ZP*ZS
                ! Den
                Obs_eq(2)%Obs_Latt(imj,1,1,1) =  Obs_eq(2)%Obs_Latt(imj,1,1,1)  +  &
                     &     (    GRC(I,I,1) * GRC(J,J,1) *Z     + &
                     &          GRC(I,J,1) * GR(I,J,1 )           &
                     &                                   ) * Z* ZP*ZS
             ENDDO
             Obs_eq(2)%Obs_Latt0(1) =  Obs_eq(2)%Obs_Latt0(1) +  Z * GRC(I,I,1) * ZP * ZS
          ENDDO


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
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE)
          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase

          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS
          Integer :: IMJ, I, J, I1, J1, no_I, no_J

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          If (NT == 0 ) then
             DO I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             ENDDO
          endif
          Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
          Do I1 = 1, Latt%N
             I    = Invlist(I1,Latt_unit%N +  1)
             no_I = 1
             Do J1 = 1,Latt%N
                J    =  Invlist(J1,Latt_unit%N +  1)
                no_J = 1
                imj = latt%imj(I1,J1)
                ! SpinZ
                Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                     &      - Z*G0T(J,I,1) * GT0(I,J,1) *ZP*ZS
                ! Den
                Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                     & + ( Z*Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I,I,1))*       &
                     &         (cmplx(1.d0,0.d0,kind(0.d0)) - G00(J,J,1))  -     &
                     &     Z * GT0(I,J,1)*G0T(J,I,1)                                ) * ZP * ZS
             Enddo
             Obs_tau(2)%Obs_Latt0(no_I) = Obs_tau(2)%Obs_Latt0(no_I) + &
                  &         Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I,I,1)) * ZP * ZS
          Enddo

        end Subroutine OBSERT
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Prints out the bins.  No need to change this routine.
!-------------------------------------------------------------------
        Subroutine  Pr_obs(LTAU)

          Implicit none

          Integer,  Intent(In) ::  Ltau

          !Local
          Integer :: I


          Do I = 1,Size(Obs_scal,1)
             Call  Print_bin_Vec(Obs_scal(I), Group_Comm)
          enddo
          Do I = 1,Size(Obs_eq,1)
             Call  Print_bin_Latt(Obs_eq(I), Group_Comm)
          enddo
          If (Ltau  == 1 ) then
             Do I = 1,Size(Obs_tau,1)
                Call  Print_bin_Latt(Obs_tau(I), Group_Comm)
             enddo
          endif

        end Subroutine Pr_obs

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Initializes observables to zero before each bins.  No need to change
!> this routine.
!-------------------------------------------------------------------
        Subroutine  Init_obs(Ltau)

          Implicit none
          Integer, Intent(In) :: Ltau

          ! Local
          Integer :: I

          Do I = 1,Size(Obs_scal,1)
             Call Obser_vec_Init(Obs_scal(I))
          Enddo

          Do I = 1,Size(Obs_eq,1)
             Call Obser_Latt_Init(Obs_eq(I))
          Enddo

          If (Ltau == 1) then
             Do I = 1,Size(Obs_tau,1)
                Call Obser_Latt_Init(Obs_tau(I))
             Enddo
          Endif

        end Subroutine Init_obs

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specify a global move on a given time slice tau.
!>
!> @details
!> @param[in] ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!> @param[out] T0_Proposal_ratio, Real
!> \verbatim
!>  T0_Proposal_ratio = T0( sigma_new -> sigma ) /  T0( sigma -> sigma_new)
!> \endverbatim
!> @param[out] S0_ratio, Real
!> \verbatim
!>  S0_ratio = e^( S_0(sigma_new) ) / e^( S_0(sigma) )
!> \endverbatim
!> @param[out] Flip_length  Integer
!> \verbatim
!>  Number of flips stored in the first  Flip_length entries of the array Flip_values.
!>  Has to be smaller than NDIM
!> \endverbatim
!> @param[out] Flip_list  Integer(Ndim)
!> \verbatim
!>  List of spins to be flipped: nsigma%f(Flip_list(1),ntau) ... nsigma%f(Flip_list(Flip_Length),ntau)
!>  Note that Ndim = size(Op_V,1)
!> \endverbatim
!> @param[out] Flip_value  Real(Ndim)
!> \verbatim
!>  Flip_value(:)= nsigma%flip(Flip_list(:),ntau)
!>  Note that Ndim = size(Op_V,1)
!> \endverbatim
!--------------------------------------------------------------------
        Subroutine Global_move_tau(T0_Proposal_ratio, S0_ratio, &
             &                     Flip_list, Flip_length,Flip_value,ntau)


          Implicit none
          Real (Kind = Kind(0.d0)),INTENT(OUT) :: T0_Proposal_ratio,  S0_ratio
          Integer                , INTENT(OUT) :: Flip_list(:)
          Real (Kind = Kind(0.d0)),INTENT(OUT) :: Flip_value(:)
          Integer, INTENT(OUT) :: Flip_length
          Integer, INTENT(IN)    :: ntau


          ! Local
          Integer :: n_op, n, ns
          Real (Kind=Kind(0.d0)) :: T0_proposal


        end Subroutine Global_move_tau

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

        Real (Kind=Kind(0.d0)), allocatable, dimension(:,:), Intent(OUT) :: Initial_field


      end Subroutine Hamiltonian_set_nsigma

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
      end Subroutine Overide_global_tau_sampling_parameters

      end Module Hamiltonian
