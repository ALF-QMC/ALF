      Module Hamiltonian

      Use Operator_mod
      Use WaveFunction_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables

      
      Implicit none


!>    Public variables. Have to be set by user 
      Type (Operator),     dimension(:,:), allocatable  :: Op_V
      Type (Operator),     dimension(:,:), allocatable  :: Op_T
      Type (WaveFunction), dimension(:),   allocatable  :: WF_L
      Type (WaveFunction), dimension(:),   allocatable  :: WF_R
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot, Thtrot ! Effective projection.  In units of Delta-tau.  
      Logical              :: Projector
!>    Defines MPI communicator 
      Integer              :: Group_Comm


!>    Privat variables 
      Type (Lattice),       private :: Latt 
      Integer,              private :: L1, L2
      real (Kind=Kind(0.d0)),        private :: ham_T , ham_U,  Ham_chem, Ham_h, Ham_J, Ham_xi, XB_X, Phi_X, Ham_tV
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta, Theta
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: One_dimensional, Checkerboard
      Integer,              private :: N_coord, Norb
      Integer, allocatable, private :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell
      



!>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau
      
!>    Storage for the Ising action
      Real (Kind=Kind(0.d0)),        private :: DW_Ising_tau(-1:1), DW_Ising_Space(-1:1)
      Integer, allocatable, private :: L_bond(:,:), L_bond_inv(:,:), Ising_nnlist(:,:)

      
    contains 


      Subroutine Ham_Set
#ifdef MPI
          Use mpi
#endif
          Implicit none

          integer :: ierr
          
          NAMELIST /VAR_Lattice/  L1, L2, Lattice_type, Model,  Checkerboard, N_SUN, Phi_X, XB_X

          NAMELIST /VAR_Hubbard/  ham_T, ham_chem, ham_U,  Dtau, Beta, Theta, Projector

          NAMELIST /VAR_Ising/    ham_T, ham_chem, ham_U, Ham_h, Ham_J, Ham_xi, Dtau, Beta, Theta, Projector

          NAMELIST /VAR_t_V/      ham_T, ham_chem, ham_tV, Dtau, Beta, Theta, Projector

#ifdef MPI
          Integer        :: Isize, Irank
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)

          If (Irank == 0 ) then
#endif
             N_SUN        = 1
             Checkerboard = .false.
             Phi_X        = 0.d0
             XB_X         = 1.d0
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(*,*) 'unable to open <parameters>',ierr
                STOP
             END IF
             READ(5,NML=VAR_lattice)
             CLOSE(5)
 
#ifdef MPI
          Endif
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(N_SUN       ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Phi_X       ,1  ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(XB_X        ,1  ,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Checkerboard,1  ,MPI_LOGICAL  , 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif
          Call Ham_latt

#ifdef MPI
          If (Irank == 0) then
#endif
             Open (Unit = 50,file="info",status="unknown",position="append")
             Write(50,*) '====================================='
             Write(50,*) 'Model is      : ', Model 
             Write(50,*) 'Lattice is    : ', Lattice_type
             Write(50,*) '# of orbitals : ', Ndim
             If (Lattice_type == "Square" .or. Lattice_type == "One_dimensional" ) then
                Write(50,*) 'X_boundary    : ', XB_X
                Write(50,*) 'Flux_X        : ', Phi_X
             Endif
             Write(50,*) 'Checkerboard  : ', Checkerboard
             Close(50)
#ifdef MPI
          Endif
#endif


          
#ifdef MPI
          If (Irank == 0 ) then
#endif
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             OPEN(Unit = 50,file="info",status="unknown",position="append")
#ifdef MPI
          Endif
#endif

          Projector = .false.
          Theta = 0.d0
          Thtrot = 0
          Select Case (Model)
          Case ("Hubbard_Mz")
             N_FL  = 2
             N_SUN = 1
#ifdef MPI
             If (Irank == 0 ) then
#endif
                READ(5,NML=VAR_Hubbard)
                Ltrot = nint(beta/dtau)
                if (Projector) Thtrot = nint(theta/dtau)
                Ltrot = Ltrot+2*Thtrot
                if (Projector) then
                  Write(50,*) 'Projective version'
                  Write(50,*) 'Theta         : ', Theta
                  Write(50,*) 'Tau_max       : ', beta
                else
                  Write(50,*) 'Finite temperture version'
                  Write(50,*) 'Beta          : ', Beta
                endif
                Write(50,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_U         : ', Ham_U
                Write(50,*) 'Ham_chem      : ', Ham_chem
#ifdef MPI
             Endif
#endif
#ifdef MPI
             CALL MPI_BCAST(Ltrot    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Thtrot   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Projector,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_U    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
          Case ("Hubbard_SU2")
             N_FL = 1
             N_SUN = 2
#ifdef MPI
             If (Irank == 0 ) then
#endif
                READ(5,NML=VAR_Hubbard)
                Ltrot = nint(beta/dtau)
                if (Projector) Thtrot = nint(theta/dtau)
                Ltrot = Ltrot+2*Thtrot
                if (Projector) then
                  Write(50,*) 'Projective version'
                  Write(50,*) 'Theta         : ', Theta
                  Write(50,*) 'Tau_max       : ', beta
                else
                  Write(50,*) 'Finite temperture version'
                  Write(50,*) 'Beta          : ', Beta
                endif
                Write(50,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_U         : ', Ham_U
                Write(50,*) 'Ham_chem      : ', Ham_chem
#ifdef MPI
             Endif
#endif
#ifdef MPI
             CALL MPI_BCAST(Ltrot    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Thtrot   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Projector,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_U    ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Dtau     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Beta     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
#endif
          Case ("Hubbard_SU2_Ising")
             N_FL = 1
             N_SUN = 2
#ifdef MPI
             If (Irank == 0 ) then
#endif
                READ(5,NML=VAR_Ising)
                Ltrot = nint(beta/dtau)
                if (Projector) Thtrot = nint(theta/dtau)
                Ltrot = Ltrot+2*Thtrot
                if (Projector) then
                  Write(50,*) 'Projective version'
                  Write(50,*) 'Theta         : ', Theta
                  Write(50,*) 'Tau_max       : ', beta
                else
                  Write(50,*) 'Finite temperture version'
                  Write(50,*) 'Beta          : ', Beta
                endif
                Write(50,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_U         : ', Ham_U
                Write(50,*) 'Ham_chem      : ', Ham_chem
                Write(50,*) 'Ham_xi        : ', Ham_xi
                Write(50,*) 'Ham_J         : ', Ham_J
                Write(50,*) 'Ham_h         : ', Ham_h
#ifdef MPI
             Endif
#endif
             If ( Lattice_type == "Honeycomb" .or. Lattice_type == "Pi_Flux" ) then 
                Write(6,*) "Hubbard_SU2_Ising is only implemented for a square lattice"
                Stop
             Endif
#ifdef MPI
             CALL MPI_BCAST(Ltrot    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Thtrot   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Projector,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_U    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Ham_xi   ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Ham_J    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Ham_h    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
             Call Setup_Ising_action
          Case ("t_V")
             N_FL  = 1
#ifdef MPI
             If (Irank == 0 ) then
#endif
                READ(5,NML=VAR_t_V )
                Ltrot = nint(beta/dtau)
                if (Projector) Thtrot = nint(theta/dtau)
                Ltrot = Ltrot+2*Thtrot
                if (Projector) then
                  Write(50,*) 'Projective version'
                  Write(50,*) 'Theta         : ', Theta
                  Write(50,*) 'Tau_max       : ', beta
                else
                  Write(50,*) 'Finite temperture version'
                  Write(50,*) 'Beta          : ', Beta
                endif
                Write(50,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
                Write(50,*) 'N_SUN         : ', N_SUN
                Write(50,*) 'N_FL          : ', N_FL
                Write(50,*) 't             : ', Ham_T
                Write(50,*) 'Ham_chem      : ', Ham_chem
                Write(50,*) 'Ham_tV         : ', Ham_tV
#ifdef MPI
             Endif
#endif
#ifdef MPI
             CALL MPI_BCAST(Ltrot    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Thtrot   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Projector,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(ham_tV   ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
          Case default 
             Write(6,*) "Model not yet implemented!"
             Stop
          end Select
          Call Ham_hop
          
          if (Projector) Call Ham_TrialWaveFunction


! #ifdef MPI
!           If (Irank == 0 )  then
! #endif
! #if defined(STAB1) 
!              Write(50,*) 'STAB1 is defined '
! #endif
! 
! #if defined(QRREF) 
!              Write(50,*) 'QRREF is defined '
! #endif
!              close(50)
!              Close(5)
! #ifdef MPI
!           endif
! #endif

          call Ham_V


        end Subroutine Ham_Set
!=============================================================================
        Subroutine Ham_Latt
          Implicit none
          !Set the lattice
          
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          Integer :: I, nc, no
          
          select case (Lattice_type)
          case("Square")
             If (L1==1 .or. L2==1 ) then
                Write(6,*) 'For one-dimensional lattices set : L2 = 1'
                stop
             endif
             Norb      = 1
             N_coord   = 2
             a1_p(1) =  1.0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.0  ; a2_p(2) =  1.d0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
          case("One_dimensional")
             If (L1 == 1 ) then 
                Write(6,*) ' For one dimensional systems set  L2 = 1 ' 
                Stop
             endif
             Norb      = 1
             N_coord   = 1
             N_coord   = 1
             a1_p(1) =  1.0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.0  ; a2_p(2) =  1.d0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
          case("Honeycomb")
             If (L1==1 .or. L2==1 ) then
                Write(6,*) 'For one-dimensional lattices set : L2 = 1'
                stop
             endif
             Norb    = 2
             N_coord = 3
             a1_p(1) =  1.D0   ; a1_p(2) =  0.d0
             a2_p(1) =  0.5D0  ; a2_p(2) =  sqrt(3.D0)/2.D0

             !del_p   =  (a2_p - 0.5*a1_p ) * 2.0/3.0
             L1_p    =  dble(L1) * a1_p
             L2_p    =  dble(L2) * a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
          case("Pi_Flux")
             If (L1==1 .or. L2==1 ) then
                Write(6,*) 'For one-dimensional lattices set : L2 = 1'
                stop
             endif
             Norb    = 2
             N_coord = 4
             a1_p(1) =  1.D0   ; a1_p(2) =   1.d0
             a2_p(1) =  1.D0   ; a2_p(2) =  -1.d0

             !del_p   =  (a2_p - 0.5*a1_p ) * 2.0/3.0
             L1_p    =  dble(L1) * (a1_p - a2_p)/2.d0
             L2_p    =  dble(L2) * (a1_p + a2_p)/2.d0
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
          case default 
             Write(6,*) "Lattice not yet implemented!"
             Stop
          end select
          !Call Print_latt(Latt)
          ! This is for the orbital structure.

          Ndim = Latt%N*Norb
          Allocate (List(Ndim,2), Invlist(Latt%N,Norb))
          nc = 0
          Do I = 1,Latt%N
             Do no = 1,Norb
                ! For the Honeycomb and pi-flux lattices no = 1,2 corresponds to the A,and B sublattice.
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc 
             Enddo
          Enddo

        end Subroutine Ham_Latt

!===================================================================================           
        Subroutine Ham_hop
          Implicit none

          !Setup the hopping
          !Per flavor, the  hopping is given by: 
          !  e^{-dtau H_t}  =    Prod_{n=1}^{Ncheck} e^{-dtau_n H_{n,t}}

          Integer :: I, I1, J1, I2, n, Ncheck,nc, nc1, no
          Complex (Kind=Kind(0.d0)) :: Z

          If ( .not. Checkerboard) then
             allocate(Op_T(1,N_FL))
             do n = 1,N_FL
                Call Op_make(Op_T(1,n),Ndim)
                nc = 1
                Select case (Lattice_type)
                Case ("Square")
                   Z  =  exp( cmplx(0.d0, 2.d0 * acos(-1.d0)*Phi_X/dble(L1), kind=kind(0.d0) ) )
                   DO I = 1, Latt%N
                      I1 = Latt%nnlist(I,1,0)
                      I2 = Latt%nnlist(I,0,1)
                      If ( Latt%list(I,1) == 0 ) then
                         Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*Z
                         Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*conjg(Z)
                      else
                         Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))*Z
                         Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))*conjg(Z)
                      endif
                      Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                   Enddo
                Case ("One_dimensional")
                   Z  =  exp( cmplx(0.d0, 2.d0 * acos(-1.d0)*Phi_X/dble(L1), kind=kind(0.d0) ) )
                   DO I = 1, Latt%N
                      I1 = Latt%nnlist(I,1,0)
                      If ( Latt%list(I,1) == 0 ) then
                         Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*Z
                         Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*conjg(Z)
                      else
                         Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))*Z
                         Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))*conjg(Z)
                      endif
                      Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                   Enddo
                Case ("Honeycomb")
                   DO I = 1, Latt%N
                      do no = 1,Norb
                         I1 = Invlist(I,no)
                         Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                      enddo
                      I1 = Invlist(I,1)
                      J1 = I1
                      Do nc1 = 1,N_coord
                         select case (nc1)
                         case (1)
                            J1 = invlist(I,2) 
                         case (2)
                            J1 = invlist(Latt%nnlist(I,1,-1),2) 
                         case (3)
                            J1 = invlist(Latt%nnlist(I,0,-1),2) 
                         case default
                            Write(6,*) ' Error in  Ham_Hop '  
                            Stop
                         end select
                         Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      Enddo
                   Enddo
                case("Pi_Flux")
                   DO I = 1, Latt%N
                      do no = 1,Norb
                         I1 = Invlist(I,no)
                         Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                      enddo
                      I1 = Invlist(I,1)
                      J1 = I1
                      Do nc1 = 1,N_coord
                         select case (nc1)
                         case (1)
                            J1 = invlist(I,2) 
                         case (2)
                            J1 = invlist(Latt%nnlist(I,0, 1),2) 
                         case (3)
                            J1 = invlist(Latt%nnlist(I,-1,1),2) 
                         case (4)
                            J1 = invlist(Latt%nnlist(I,-1,0),2) 
                         case default
                            Write(6,*) ' Error in  Ham_Hop '  
                            Stop
                         end select
                         if (nc1 == 1 ) then
                            Op_T(nc,n)%O(I1,J1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                            Op_T(nc,n)%O(J1,I1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                         Else
                            Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                            Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         endif
                      Enddo
                   Enddo
                case default 
                   Write(6,*) "Lattice not yet implemented!"
                   Stop
                end select
                Do I = 1,Ndim
                   Op_T(nc,n)%P(i) = i 
                Enddo
                if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                   Op_T(nc,n)%g = 0.d0
                else
                   Op_T(nc,n)%g = -Dtau
                endif
                Op_T(nc,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                Call Op_set(Op_T(nc,n))
                !Do I = 1,Size(Op_T(nc,n)%E,1)
                !   Write(6,*) Op_T(nc,n)%E(I)
                !Enddo
             Enddo
          else
             select case (Lattice_type)
             case ("Square")
                Allocate(Op_T(N_coord*Latt%N,N_FL))
                do n = 1,N_FL
                   do i  =  1, N_coord*Latt%N
                      call Op_make(Op_T(i,n),2)
                   enddo
                   nc = 0
                   do I = 1,Latt%N
                      I1 = I
                      do nc1 = 1,N_coord
                         nc = nc + 1
                         if (nc1 == 1 ) I2 = latt%nnlist(I,1,0) 
                         if (nc1 == 2 ) I2 = latt%nnlist(I,0,1)
                         Op_T(nc,n)%P(1) = I1
                         Op_T(nc,n)%P(2) = I2
                         Op_T(nc,n)%O(1,2) = cmplx(-Ham_T ,0.d0, kind(0.D0)) 
                         Op_T(nc,n)%O(2,1) = cmplx(-Ham_T ,0.d0, kind(0.D0))
                         Op_T(nc,n)%O(1,1) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0)) 
                         Op_T(nc,n)%O(2,2) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0))
                         if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                            Op_T(nc,n)%g = 0.d0
                         else
                            Op_T(nc,n)%g = -Dtau
                         endif
                         Op_T(nc,n)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                         Call Op_set( Op_T(nc,n) )
                      Enddo
                   Enddo
                Enddo
             case ("Pi_Flux")
                Allocate(Op_T(N_coord*Latt%N,N_FL))
                do n = 1,N_FL
                   !Write(6,*) 'N_coord, Latt%N ',  N_coord, Latt%N
                   do i  =  1, N_coord*Latt%N
                      call Op_make(Op_T(i,n),2)
                   enddo
                   nc = 0
                   do I = 1,Latt%N
                      I1 = Invlist(I,1)
                      do nc1 = 1,N_coord
                         nc = nc + 1
                         If (nc1 == 1 )  I2 = invlist(I,2)
                         If (nc1 == 2 )  I2 = invlist(Latt%nnlist(I,0, 1),2) 
                         If (nc1 == 3 )  I2 = invlist(Latt%nnlist(I,-1,1),2)
                         If (nc1 == 4 )  I2 = invlist(Latt%nnlist(I,-1,0),2) 
                         Op_T(nc,n)%P(1) = I1
                         Op_T(nc,n)%P(2) = I2
                         if (nc1 == 1 ) then
                            Op_T(nc,n)%O(1,2) = cmplx( Ham_T,    0.d0, kind(0.D0))
                            Op_T(nc,n)%O(2,1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                         Else
                            Op_T(nc,n)%O(1,2) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                            Op_T(nc,n)%O(2,1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         endif
                         Op_T(nc,n)%O(1,1) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0)) 
                         Op_T(nc,n)%O(2,2) = cmplx(-Ham_Chem/4.d0 ,0.d0, kind(0.D0))
                         if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
                            Op_T(nc,n)%g = 0.d0
                         else
                            Op_T(nc,n)%g = -Dtau
                         endif
                         Op_T(nc,n)%alpha  = cmplx( 0.d0, 0.d0, kind(0.D0) )
                         Call Op_set( Op_T(nc,n) )
                      Enddo
                   Enddo
                Enddo
             case default
                Write(6,*) "Checkeboard is not implemented for this lattice"
                Stop
             End select
          endif
          
        end Subroutine Ham_hop

!===================================================================================           
        Subroutine Ham_TrialWaveFunction
        
          Implicit none
          COMPLEX(Kind=Kind(0.d0)), allocatable :: H0(:,:), U(:,:)
          Real(Kind=Kind(0.d0)), allocatable :: En(:)
          
          COMPLEX(Kind=Kind(0.d0)) :: Z
          Integer :: n, n_part, i, I1, I2, J1, no, nc1
          
          N_part=Ndim/2
          
          Allocate(WF_L(N_FL),WF_R(N_FL))
          do n=1,N_FL
            Call WF_alloc(WF_L(n),Ndim,N_part)
            Call WF_alloc(WF_R(n),Ndim,N_part)
          enddo
          
          Allocate(H0(Ndim,Ndim),U(Ndim,Ndim),En(Ndim))
          do n = 1,N_FL
            H0=0.d0
            Select case (Lattice_type)
            Case ("Square")
                Z  =  exp( cmplx(0.d0, 2.d0 * acos(-1.d0)*Phi_X/dble(L1), kind=kind(0.d0) ) )
                DO I = 1, Latt%N
                  I1 = Latt%nnlist(I,1,0)
                  I2 = Latt%nnlist(I,0,1)
                  If ( Latt%list(I,1) == 0 ) then
                      H0(I,I1) = cmplx( Ham_T, 0.d0, kind(0.D0))*Z
                      H0(I1,I) = cmplx( Ham_T, 0.d0, kind(0.D0))*conjg(Z)
                  else
                      H0(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))*Z
                      H0(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))*conjg(Z)
                  endif
                  H0(I,I2) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                  H0(I2,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                  H0(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                Enddo
            Case ("One_dimensional")
                Z  =  exp( cmplx(0.d0, 2.d0 * acos(-1.d0)*Phi_X/dble(L1), kind=kind(0.d0) ) )
                DO I = 1, Latt%N
                  I1 = Latt%nnlist(I,1,0)
                  If ( Latt%list(I,1) == 0 ) then
                      H0(I,I1) = cmplx( Ham_T, 0.d0, kind(0.D0))*Z
                      H0(I1,I) = cmplx( Ham_T, 0.d0, kind(0.D0))*conjg(Z)
                  else
                      H0(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))*Z
                      H0(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))*conjg(Z)
                  endif
                  H0(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                Enddo
            Case ("Honeycomb")
                DO I = 1, Latt%N
                  do no = 1,Norb
                      I1 = Invlist(I,no)
                      H0(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                  enddo
                  I1 = Invlist(I,1)
                  J1 = I1
                  Do nc1 = 1,N_coord
                      select case (nc1)
                      case (1)
                        J1 = invlist(I,2) 
                      case (2)
                        J1 = invlist(Latt%nnlist(I,1,-1),2) 
                      case (3)
                        J1 = invlist(Latt%nnlist(I,0,-1),2) 
                      case default
                        Write(6,*) ' Error in  Ham_Hop '  
                        Stop
                      end select
                      If ( Latt%list(J1,1) == 0 .and. nc1==2 ) then
                        H0(I1,J1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                        H0(J1,I1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                      else
                        H0(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                        H0(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      endif
                  Enddo
                Enddo
!             case("Pi_Flux")
!                 DO I = 1, Latt%N
!                   do no = 1,Norb
!                       I1 = Invlist(I,no)
!                       H0(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
!                   enddo
!                   I1 = Invlist(I,1)
!                   J1 = I1
!                   Do nc1 = 1,N_coord
!                       select case (nc1)
!                       case (1)
!                         J1 = invlist(I,2) 
!                       case (2)
!                         J1 = invlist(Latt%nnlist(I,0, 1),2) 
!                       case (3)
!                         J1 = invlist(Latt%nnlist(I,-1,1),2) 
!                       case (4)
!                         J1 = invlist(Latt%nnlist(I,-1,0),2) 
!                       case default
!                         Write(6,*) ' Error in  Ham_Hop '  
!                         Stop
!                       end select
!                       if (nc1 == 1 ) then
!                         H0(I1,J1) = cmplx( Ham_T,    0.d0, kind(0.D0))
!                         H0(J1,I1) = cmplx( Ham_T,    0.d0, kind(0.D0))
!                       Else
!                         H0(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
!                         H0(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
!                       endif
!                   Enddo
!                 Enddo
            case default 
                Write(6,*) "Projector not yet implemented!"
                Stop
            end select
          
            Call Diag(H0,U,En)
            
            do I2=1,N_part
               do I1=1,Ndim
                  WF_L(n)%P(I1,I2)=U(I1,I2)
                  WF_R(n)%P(I1,I2)=U(I1,I2)
               enddo
            enddo
         enddo
          
        end Subroutine Ham_TrialWaveFunction

!===================================================================================           
        Subroutine Ham_V
          
          Implicit none 
          
          Integer :: nf, I, I1, I2,  nc, nc1,  J
          Real (Kind=Kind(0.d0)) :: X

          
          Select case (Model)
          Case ("Hubbard_SU2")  
             !Write(50,*) 'Model is ', Model
             Allocate(Op_V(Ndim,N_FL))
             do nf = 1,N_FL
                do i  = 1, Ndim
                   Call Op_make(Op_V(i,nf),1)
                enddo
             enddo
             Do nf = 1,N_FL
                nc = 0
                Do i = 1,Ndim
                   nc = nc + 1
                   Op_V(nc,nf)%P(1) = I
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                Enddo
             Enddo
          Case ("Hubbard_Mz") 
             Allocate(Op_V(Ndim,N_FL))
             do nf = 1,N_FL
                do i  = 1, Ndim
                   Call Op_make(Op_V(i,nf),1)
                enddo
             enddo
             Do nf = 1,N_FL
                nc = 0
                X = 1.d0
                if (nf == 2) X = -1.d0
                Do i = 1,Ndim
                   nc = nc + 1
                   Op_V(nc,nf)%P(1) = I
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
                   Op_V(nc,nf)%g      = X*SQRT(CMPLX(DTAU*ham_U/2.d0, 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                Enddo
             Enddo
          Case ("Hubbard_SU2_Ising") 
             Allocate(Op_V(3*Ndim,N_FL))
             do nf = 1,N_FL
                do i  =  1, N_coord*Ndim
                   call Op_make(Op_V(i,nf),2)
                enddo
                do i  = N_coord*Ndim +1 ,  N_coord*Ndim + Ndim ! For Hubbatd
                   Call Op_make(Op_V(i,nf),1)
                enddo
             enddo
             Do nc = 1,Ndim*N_coord   ! Runs over bonds.  Coordination number = 2.
                ! For the square lattice Ndim = Latt%N
                I1 = L_bond_inv(nc,1)
                I2 = I1
                if ( L_bond_inv(nc,2)  == 1 ) I2 = Latt%nnlist(I1,1,0)
                if ( L_bond_inv(nc,2)  == 2 ) I2 = Latt%nnlist(I1,0,1)
                Op_V(nc,1)%P(1) = I1
                Op_V(nc,1)%P(2) = I2
                Op_V(nc,1)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
                Op_V(nc,1)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
                Op_V(nc,1)%g = cmplx(-dtau*Ham_xi,0.D0,kind(0.D0))
                Op_V(nc,1)%alpha = cmplx(0d0,0.d0, kind(0.D0)) 
                Op_V(nc,1)%type =1
                Call Op_set( Op_V(nc,1) )
             Enddo
             
             Do i = 1,Ndim
                nc1 = N_coord*Ndim + i
                Op_V(nc1,1)%P(1)   = i
                Op_V(nc1,1)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                Op_V(nc1,1)%g      = sqrt(cmplx(-dtau*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0)))
                Op_V(nc1,1)%alpha  = cmplx(-0.5d0,0.d0, kind(0.d0))
                Op_V(nc1,1)%type   = 2
                Call Op_set( Op_V(nc1,1) )
             Enddo
          case ("t_V")
             Allocate(Op_V(N_coord*Latt%N,1))
             do i  =  1, N_coord*Latt%N
                call Op_make(Op_V(i,1),2)
             enddo
             select case (Lattice_type)
             case ("Square")
                !Write(6,*) "N_coord, Latt%N", N_coord, Latt%N, Dtau, ham_tV
                nc = 0
                do I = 1,Latt%N
                   I1 = I
                   do nc1 = 1,N_coord
                      nc = nc + 1
                      if (nc1 == 1 ) I2 = latt%nnlist(I,1,0) 
                      if (nc1 == 2 ) I2 = latt%nnlist(I,0,1)
                      Op_V(nc,1)%P(1) = I1
                      Op_V(nc,1)%P(2) = I2
                      Op_V(nc,1)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
                      Op_V(nc,1)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
                      Op_V(nc,1)%g      = SQRT(CMPLX(-DTAU*ham_tV, 0.D0, kind(0.D0))) 
                      Op_V(nc,1)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                      Op_V(nc,1)%type   = 2
                      Call Op_set( Op_V(nc,1) )
                   Enddo
                Enddo
             case ("Pi_Flux")
                nc = 0
                do I = 1,Latt%N
                   I1 = Invlist(I,1)
                   do nc1 = 1,N_coord
                      nc = nc + 1
                      If (nc1 == 1 )  I2 = invlist(I,2)
                      If (nc1 == 2 )  I2 = invlist(Latt%nnlist(I,0, 1),2) 
                      If (nc1 == 3 )  I2 = invlist(Latt%nnlist(I,-1,1),2)
                      If (nc1 == 4 )  I2 = invlist(Latt%nnlist(I,-1,0),2) 
                      Op_V(nc,1)%P(1) = I1
                      Op_V(nc,1)%P(2) = I2
                      Op_V(nc,1)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0)) 
                      Op_V(nc,1)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
                      Op_V(nc,1)%g     = SQRT(CMPLX(-DTAU*ham_tV, 0.D0, kind(0.D0))) 
                      Op_V(nc,1)%alpha = cmplx(0.d0, 0.d0, kind(0.D0))
                      Op_V(nc,1)%type  = 2
                      Call Op_set( Op_V(nc,1) )
                   Enddo
                Enddo
             case default
                Write(6,*) "Lattice for t-V  not implemented"
                Stop
             End select
          case default
             Write(6,*) "Model not yet implemented!"
             Stop
          end Select
          
        end Subroutine Ham_V

!===================================================================================           
        Real (Kind=Kind(0.d0)) function S0(n,nt)  
          Implicit none
          Integer, Intent(IN) :: n,nt
          Integer :: nt1,I
          S0 = 1.d0
          If ( Op_V(n,1)%type == 1 ) then 
             do i = 1,4
                S0 = S0*DW_Ising_space(nsigma(n,nt)*nsigma(Ising_nnlist(n,i),nt))
             enddo
             nt1 = nt +1 
             if (nt1 > Ltrot) nt1 = 1
             S0 = S0*DW_Ising_tau(nsigma(n,nt)*nsigma(n,nt1))
             nt1 = nt - 1 
             if (nt1 < 1  ) nt1 = Ltrot
             S0 = S0*DW_Ising_tau(nsigma(n,nt)*nsigma(n,nt1))
             If (S0 < 0.d0) Write(6,*) 'S0 : ', S0
          endif
          
        end function S0

!===================================================================================           
        Subroutine Setup_Ising_action
          
          ! This subroutine sets up lists and arrays so as to enable an 
          ! an efficient calculation of  S0(n,nt) 

          Integer :: nc, nth, n, n1, n2, n3, n4, I, I1, n_orientation
          Real (Kind=Kind(0.d0)) :: X_p(2)

          ! Setup list of bonds for the square lattice.
          Allocate (L_Bond(Latt%N,2),  L_bond_inv(Latt%N*N_coord,2) )
          
          nc = 0
          do nth = 1,2*N_coord  
             Do n1= 1, L1/2
                Do n2 = 1,L2
                   nc = nc + 1
                   n_orientation = 1
                   I1 = 1
                   If (nth == 1 ) then
                      X_p = dble(2*n1)*latt%a1_p + dble(n2)*latt%a2_p 
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 1
                   elseif (nth == 2) then
                      X_p = dble(2*n1)*latt%a1_p + dble(n2)*latt%a2_p  + latt%a1_p
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 1
                   elseif (nth == 3) then
                      X_p = dble(n2)*latt%a1_p + dble(2*n1)*latt%a2_p 
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 2
                   elseif (nth == 4) then
                      X_p = dble(n2)*latt%a1_p + dble(2*n1)*latt%a2_p  + latt%a2_p
                      I1 = Inv_R(X_p,Latt)
                      n_orientation = 2
                   endif
                   L_bond(I1,n_orientation) = nc
                   L_bond_inv(nc,1) = I1  
                   L_bond_inv(nc,2) = n_orientation 
                   ! The bond is given by  I1, I1 + a_(n_orientation).
                Enddo
             Enddo
          Enddo
          ! Setup the nearest neigbour lists for the Ising spins. 
          allocate(Ising_nnlist(2*Latt%N,4)) 
          do I  = 1,Latt%N
             n  = L_bond(I,1)
             n1 = L_bond(Latt%nnlist(I, 1, 0),2)
             n2 = L_bond(Latt%nnlist(I, 0, 0),2)
             n3 = L_bond(Latt%nnlist(I, 0,-1),2)
             n4 = L_bond(Latt%nnlist(I, 1,-1),2)
             Ising_nnlist(n,1) = n1
             Ising_nnlist(n,2) = n2
             Ising_nnlist(n,3) = n3
             Ising_nnlist(n,4) = n4
             n  = L_bond(I,2)
             n1 = L_bond(Latt%nnlist(I, 0, 1),1)
             n2 = L_bond(Latt%nnlist(I,-1, 1),1)
             n3 = L_bond(Latt%nnlist(I,-1, 0),1)
             n4 = L_bond(Latt%nnlist(I, 0, 0),1)
             Ising_nnlist(n,1) = n1
             Ising_nnlist(n,2) = n2
             Ising_nnlist(n,3) = n3
             Ising_nnlist(n,4) = n4
          enddo
          DW_Ising_tau  ( 1) = tanh(Dtau*Ham_h)
          DW_Ising_tau  (-1) = 1.D0/DW_Ising_tau(1)
          DW_Ising_Space( 1) = exp(-2.d0*Dtau*Ham_J) 
          DW_Ising_Space(-1) = exp( 2.d0*Dtau*Ham_J) 
          
        End Subroutine Setup_Ising_action
!===================================================================================           
        Subroutine  Alloc_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Ns,Nt,No
          Character (len=64) ::  Filename

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
                Ns = Latt%N;  No = Norb;  Filename ="Green"
             case (2)
                Ns = Latt%N;  No = Norb;  Filename ="SpinZ"
             case (3)
                Ns = Latt%N;  No = Norb;  Filename ="SpinXY"
             case (4)
                Ns = Latt%N;  No = Norb;  Filename ="Den"
             case (5)
                Ns = Latt%N;  No = Norb;  Filename ="SpinT"
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Nt = 1
             Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
          enddo

          If (Ltau == 1) then 
             ! Equal time correlators
             Allocate ( Obs_tau(4) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Ns = Latt%N; No = Norb;  Filename ="Green"
                case (2)
                   Ns = Latt%N; No = Norb;  Filename ="SpinZ"
                case (3)
                   Ns = Latt%N; No = Norb;  Filename ="SpinXY"
                case (4)
                   Ns = Latt%N; No = Norb;  Filename ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '  
                end select
                Nt = Ltrot+1-2*Thtrot
                Call Obser_Latt_make(Obs_tau(I),Ns,Nt,No,Filename)
             enddo
          endif
        end Subroutine Alloc_obs

        
!========================================================================
        ! Functions for Global moves.  These move are not implemented in this example.
        Subroutine Global_move(T0_Proposal_ratio,nsigma_old,size_clust)
          
          !>  The input is the field nsigma declared in this module. This routine generates a 
          !>  global update with  and returns the propability  
          !>  T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)  
          !>   
          Implicit none
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
          Integer, dimension(:,:),  allocatable, intent(in)  :: nsigma_old

          T0_Proposal_ratio  = 0.d0

        End Subroutine Global_move
!========================================================================
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none 
          
          !> Arguments
          Integer, dimension(:,:), allocatable, intent(IN) :: Nsigma_old

          Delta_S0_global = 0.d0

        end Function Delta_S0_global
!========================================================================
        Subroutine Obser(GR,Phase,Ntau)
          
          Implicit none
          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, ZZ, ZXY
          Integer :: I,J, imj, nf, dec, I1, J1, no_I, no_J,n
          
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
             

          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Do n  = 1,Size(Op_T,1)
             Do nf = 1,N_FL
                Do I = 1,Size(Op_T(n,nf)%O,1)
                   Do J = 1,Size(Op_T(n,nf)%O,2)
                      Zkin = Zkin +  Op_T(n,nf)%O(i, j)*Grc( Op_T(n,nf)%P(I), Op_T(n,nf)%P(J), nf )
                   ENddo
                Enddo
             Enddo
          Enddo
          Zkin = Zkin * dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          dec = 1
          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising" ) then
            dec = 1
          elseif (Model == "Hubbard_MZ") then
            dec = 2
          endif
          Do I = 1,Ndim
             ZPot = ZPot + Grc(i,i,1) * Grc(i,i, dec)
          Enddo
          Zpot = Zpot*ham_U
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


          
          ! Compute spin-spin, Green, and den-den correlation functions  !  This is general N_SUN, and  N_FL = 1
          DO I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising" .or. Model == "t_V"  ) then 
             Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                        &               Z * GRC(I1,J1,1) *  ZP*ZS 
                   ! SpinZ
                   Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
                        &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                   ! SpinXY
                   Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
                        &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                   ! SpinT
                   Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) + &
                        &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                   ! Den
                   Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J)  +  &
                        &     (    GRC(I1,I1,1) * GRC(J1,J1,1) *Z     + &
                        &          GRC(I1,J1,1) * GR(I1,J1,1 )           &
                        &                                   ) * Z* ZP*ZS
                ENDDO
                Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  Z * GRC(I1,I1,1) * ZP * ZS
             ENDDO
          elseif (Model == "Hubbard_Mz" ) Then
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                        &                          ( GRC(I1,J1,1) + GRC(I1,J1,2)) *  ZP*ZS 
                   ! SpinZ
                   ZZ =       GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                        &    (GRC(I1,I1,2) - GRC(I1,I1,1))*(GRC(J1,J1,2) - GRC(J1,J1,1))  
                   Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + ZZ* ZP*ZS
                   ! SpinXY
                   ! c^d_(i,u) c_(i,d) c^d_(j,d) c_(j,u)  +  c^d_(i,d) c_(i,u) c^d_(j,u) c_(j,d)
                   ZXY =  GRC(I1,J1,1) * GR(I1,J1,2) +  GRC(I1,J1,2) * GR(I1,J1,1) 
                   Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + ZXY* ZP*ZS
                   ! SpinT
                   Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                   !Den
                   Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) + &
                        & (   GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                        &   (GRC(I1,I1,2) + GRC(I1,I1,1))*(GRC(J1,J1,2) + GRC(J1,J1,1))     ) * ZP*ZS
                enddo
                Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  (GRC(I1,I1,1) + GRC(I1,I1,2)) * ZP*ZS
             enddo
          Endif
                

        end Subroutine Obser
!=====================================================
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
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising" .or. Model == "t_V"  ) then 
             Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & +  Z * GT0(I1,J1,1) * ZP* ZS
                   
                   ! SpinZ
                   Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                   
                   ! SpinXY
                   Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                   
                   ! Den
                   Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & + ( Z*Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1))*       &
                        &         (cmplx(1.d0,0.d0,kind(0.d0)) - G00(J1,J1,1))  -     &
                        &     Z * GT0(I1,J1,1)*G0T(J1,I1,1)                                ) * ZP * ZS
                Enddo
                Obs_tau(4)%Obs_Latt0(no_I) = Obs_tau(4)%Obs_Latt0(no_I) + &
                     &         Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1)) * ZP * ZS
             Enddo
          Elseif ( Model == "Hubbard_Mz"  ) then 
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   !Green
                   Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &   +   ( GT0(I1,J1,1) + GT0(I1,J1,2) ) * ZP* ZS

                   !SpinZ
                   Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                       & +  ( &
                       &    (GTT(I1,I1,1) -  GTT(I1,I1,2) ) * ( G00(J1,J1,1)  -  G00(J1,J1,2) )   &
                       &  - (G0T(J1,I1,1) * GT0(I1,J1,1)  +  G0T(J1,I1,2) * GT0(I1,J1,2) )    )*ZP*ZS

                   !SpinXY
                   Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &  - &
                        &   (G0T(J1,I1,1) * GT0(I1,J1,2)  +  G0T(J1,I1,2) * GT0(I1,J1,1))*ZP*ZS
                   !Den
                   Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & +  (                                        &  
                        &    (cmplx(2.D0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2) ) * &
                        &    (cmplx(2.D0,0.d0,kind(0.d0)) - G00(J1,J1,1) - G00(J1,J1,2) )   &
                        & -  ( G0T(J1,I1,1) * GT0(I1,J1,1) + G0T(J1,I1,2) * GT0(I1,J1,2) )  )*ZP*ZS     

                enddo
             
                Obs_tau(4)%Obs_Latt0(no_I) =  Obs_tau(4)%Obs_Latt0(no_I) + &
                     &       (cmplx(2.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2)) * ZP * ZS
             Enddo
          Endif
          
        end Subroutine OBSERT

!==========================================================        
        Subroutine  Pr_obs(LTAU)

          Implicit none

          Integer,  Intent(In) ::  Ltau
          
          !Local 
          Integer :: I


          Do I = 1,Size(Obs_scal,1)
             Call  Print_bin_Vec(Obs_scal(I),Group_Comm)
          enddo
          Do I = 1,Size(Obs_eq,1)
             Call  Print_bin_Latt(Obs_eq(I),Latt,dtau,Group_Comm)
          enddo
          If (Ltau  == 1 ) then
             Do I = 1,Size(Obs_tau,1)
                Call  Print_bin_Latt(Obs_tau(I),Latt,dtau,Group_Comm)
             enddo
          endif

        end Subroutine Pr_obs
!===================================================================================           
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
!---------------------------------------------------------------------
        Subroutine  Hamiltonian_set_random_nsigma
          
          ! The user can set the initial configuration
          
          Implicit none
          
          Integer :: I, nt
          
          Do nt = 1,Ltrot
             Do I = 1,Size(OP_V,1)
                nsigma(I,nt)  = 1
                if ( ranf_wrap()  > 0.5D0 ) nsigma(I,nt)  = -1
             enddo
          enddo
          
        end Subroutine Hamiltonian_set_random_nsigma

!---------------------------------------------------------------------
        Subroutine Global_move_tau(T0_Proposal_ratio, S0_ratio, &
             &                     Flip_list, Flip_length,Flip_value,ntau)

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
!> T0_proposal       = T0 ( sigma -> sigma_new )
!--------------------------------------------------------------------
          
          Implicit none 
          Real (Kind= kind(0.d0)), INTENT(INOUT) :: T0_Proposal_ratio,  S0_ratio
          Integer,    allocatable, INTENT(INOUT) :: Flip_list(:), Flip_value(:)
          Integer, INTENT(INOUT) :: Flip_length
          Integer, INTENT(IN)    :: ntau


          ! Local
          Integer :: n_op, n, ns
          Real (Kind=Kind(0.d0)) :: T0_proposal
          Flip_length = nranf(1)
          
          do n = 1,flip_length 
             n_op = nranf(size(OP_V,1))
             Flip_list(n)  = n_op
             ns            = nsigma(n_op,ntau)
             If ( OP_V(n_op,1)%type == 1 ) then 
                ns = nsigma(n_op,ntau)
                T0_Proposal       =  1.d0 - 1.d0/(1.d0+S0(n_op,ntau)) ! No move prob
                If ( T0_Proposal > Ranf_wrap() ) then
                   T0_Proposal_ratio =  1.d0 / S0(n_op,ntau)
                else
                   T0_Proposal_ratio = 0.d0
                endif
                S0_ratio          =  S0(n_op,ntau)
                Flip_value(n)     = - ns
             else
                Flip_value(n)     = NFLIPL(nsigma(n_op,ntau),nranf(3))
                T0_Proposal_ratio = 1.d0
                S0_ratio          = 1.d0
             endif
          Enddo


        end Subroutine Global_move_tau

      end Module Hamiltonian
