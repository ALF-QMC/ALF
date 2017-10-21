      Module Hamiltonian

      Use Operator_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables

      
      Implicit none


!>    Public variables. Have to be set by user 
      Type (Operator), dimension(:,:), allocatable  :: Op_V
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot
!>    Defines MPI communicator 
      Integer              :: Group_Comm

!>    Variables for updating scheme
!       Logical              :: Propose_S0, Global_moves
!       Integer              :: N_Global
      

!>    Private variables 
      Type (Lattice),       private :: Latt 
      Integer,              private :: L1, L2
      real (Kind=Kind(0.d0)),        private :: ham_T , ham_U,  Ham_chem, Ham_h, Ham_J, Ham_xi, Ham_xi2
      real (Kind=Kind(0.d0)),        private :: Dtau, Beta, Model_sign
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: One_dimensional
      Integer,              private :: N_coord, Norb, Model_vers
      Integer, allocatable, private :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell
      
!>    Private variables for observing Z_x_ising
      Real (Kind=Kind(0.d0)), private :: eq_x_ising, neq_x_ising
      Integer                         :: nBlub, nBlub2



!>    Private Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau
      
!>    Storage for the Ising action
      Real (Kind=Kind(0.d0)),        private :: DW_Ising_tau(-1:1), DW_Ising_Space(-1:1)
      Integer, allocatable, private :: L_bond(:,:), L_bond_inv(:,:), Ising_nnlist(:,:)

!>    Variables for the Wolff cluster update
      Real (Kind=Kind(0.d0)),        private :: addProb_space, addProb_tau
      Integer,                       private :: N_ising

!>    Variables for the Geometric cluster update
      Real (Kind=Kind(0.d0)),        private :: Geo_AddProb_space, Geo_AddProb_tau
!       Logical, allocatable, dimension(:,:), private :: Geo_cluster
      Integer, private :: R_init(2)
      
!>    Experimenting
      Integer, private :: i_sweep
      Complex (Kind=Kind(0.d0)) :: Z_m2, Z_m_0, Z_m2_0
      
    contains 


      Subroutine Ham_Set

#ifdef MPI
          Use mpi
#endif   
          Implicit none


          integer :: ierr
          Character (len=64) :: file1
          
          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model, Model_vers

          NAMELIST /VAR_Hubbard/  ham_T, ham_chem, ham_U, Dtau, Beta, N_SUN!, N_Global, Propose_S0

          NAMELIST /VAR_Ising/    Ham_h, Ham_J, Ham_xi, Ham_xi2

#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g
          !if ( irank_g == 0 )   write(6,*) "Mpi Test", igroup, isize_g
#endif
          

#ifdef MPI
          If (Irank == 0 ) then
#endif
             ! These are global variables
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(*,*) 'unable to open <parameters>',ierr
                STOP
             END IF
             Model_vers = 1
             READ(5,NML=VAR_lattice)
             CLOSE(5)
 
#ifdef MPI
          Endif
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Model_vers  ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
#endif
          if     ( Model_vers == 1 ) then
            Model_sign =  1.d0
          elseif ( Model_vers == 2 ) then
            Model_sign = -1.d0
          endif
          
          Call Ham_latt
          
          If ( Model == "NematicDirac") then
            N_FL = 1
            N_SUN = 2
          elseif ( Model == "NematicDirac3") then
            N_FL = 1
            N_SUN = 2
          elseif ( Model == "NematicDirac2") then
            N_FL = 1
            N_SUN = 2
          elseif ( Model == "yyhe") then
            N_FL = 2
            N_SUN = 1
          else
            Write(6,*) "Model not yet implemented!", Model
            Stop
          endif
          
!           N_Global= 0
!           Propose_S0 = .false.
#ifdef MPI
          If (Irank_g == 0 ) then
#endif
            File1 = "parameters"
#if defined(TEMPERING) 
            write(File1,'(A,I0,A)') "Temp_",igroup,"/parameters"
#endif
            ham_xi2 = 0
            OPEN(UNIT=5,FILE=File1,STATUS='old',ACTION='read',IOSTAT=ierr)
            READ(5,NML=VAR_Hubbard)
            If ( Model == "NematicDirac"  ) Read(5,NML=VAR_Ising)
            If ( Model == "NematicDirac2" ) Read(5,NML=VAR_Ising)
            If ( Model == "NematicDirac3" ) Read(5,NML=VAR_Ising)
            If ( Model == "yyhe"          ) Read(5,NML=VAR_Ising)
            CLOSE(5)
             
#ifdef MPI
          endif
          CALL MPI_BCAST(ham_T          ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_chem       ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_U          ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(Dtau           ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(Beta           ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(N_SUN          ,1,MPI_INTEGER,0,Group_Comm,ierr)
          CALL MPI_BCAST(Ham_xi   ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(Ham_J    ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(Ham_h    ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(Ham_xi2  ,1,MPI_REAL8,0,Group_Comm,ierr)
#endif

  
          Call Ham_hop
          Ltrot = nint(beta/dtau)

          If  ( Model == "NematicDirac"  ) Call Setup_Ising_action
          If  ( Model == "NematicDirac2" ) Call Setup_Ising_action
          If  ( Model == "NematicDirac3" ) Call Setup_Ising_action
          If  ( Model == "yyhe"          ) Call Setup_Ising_action
          
#if defined(TEMPERING)
           write(File1,'(A,I0,A)') "Temp_",igroup,"/info"
#else
           write(File1,'(A,I0)') "info"
#endif
          
#ifdef MPI
          If (Irank_g == 0) then
#endif
             Open (Unit = 50,file=File1,status="unknown",position="append")
             Write(50,*) '====================================='
             Write(50,*) 'Model is      : ', Model 
             Write(50,*) 'Version       : ', Model_vers 
             Write(50,*) 'Lattice is    : ', Lattice_type
             Write(50,*) '# of orbitals : ', Ndim
             Write(50,*) 'Beta          : ', Beta
             Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
             Write(50,*) 'N_SUN         : ', N_SUN
             Write(50,*) 'N_FL          : ', N_FL
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'Ham_U         : ', Ham_U
             Write(50,*) 'Ham_chem      : ', Ham_chem
             Write(50,*) 'Ham_xi        : ', Ham_xi
             Write(50,*) 'Ham_xi2       : ', Ham_xi2
             Write(50,*) 'Ham_J         : ', Ham_J
             Write(50,*) 'Ham_h         : ', Ham_h
#if defined(STAB1) 
             Write(50,*) 'STAB1 is defined '
#endif
#if defined(QRREF) 
             Write(50,*) 'QRREF is defined '
#endif
             close(50)
#ifdef MPI
          endif
#endif
          call Ham_V

        end Subroutine Ham_Set
!=============================================================================
        Subroutine Ham_Latt
          Implicit none
          !Set the lattice
          
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          Integer :: I, nc, no

          If ( Lattice_type =="BipartiteSquare" ) then
             Norb = 2
             N_coord   = 4
             One_dimensional = .false.
             a1_p(1) =  1.D0/sqrt(2.D0)  ; a1_p(2) =  1.D0/sqrt(2.D0)
             a2_p(1) =  1.D0/sqrt(2.D0)  ; a2_p(2) = -1.D0/sqrt(2.D0)
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             If ( L1 == 1 .or. L2 == 1 ) then 
                Write(6,*) ' One dimensional systems not implemented '
                Stop
             endif
          elseif ( (Lattice_type =="Square") .and. ( Model == "yyhe" ) ) then
             Norb = 2
             N_coord   = 4
             One_dimensional = .false.
             a1_p(1) = 1.D0  ; a1_p(2) = 0.D0
             a2_p(1) = 0.D0  ; a2_p(2) = 1.D0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             If ( L1 == 1 .or. L2 == 1 ) then 
                Write(6,*) ' One dimensional systems not implemented '
                Stop
             endif
          elseif ( (Lattice_type =="Square") .and. ( Model == "NematicDirac2" ) ) then
             Norb = 2
             N_coord   = 4
             One_dimensional = .false.
             a1_p(1) = 1.D0  ; a1_p(2) = 0.D0
             a2_p(1) = 0.D0  ; a2_p(2) = 1.D0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             If ( L1 == 1 .or. L2 == 1 ) then 
                Write(6,*) ' One dimensional systems not implemented '
                Stop
             endif
          else
             Write(6,*) "Lattice not yet implemented!"
             Stop
          endif
          
!           Write(6,*) "Lattice"
!           Do I =1,Latt%N
!             Write(6,*) Latt%list(i,1)*a1_p + Latt%list(i,2)*a2_p
!           Enddo
!           
!           Write(6,*) "BZ"
!           Do I =1,Latt%N
!             Write(6,*) Latt%listk(i,1)*Latt%b1_p + Latt%listk(i,2)*Latt%b2_p
!           Enddo


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

!===================================================================================           
        Subroutine Ham_hop
          Implicit none

          !Setup the hopping
          !Per flavor, the  hopping is given by: 
          !  e^{-dtau H_t}  =    Prod_{n=1}^{Ncheck} e^{-dtau_n H_{n,t}}

          Integer :: I, I1, J1, n, Ncheck, nc, nc1, no
          complex (Kind=Kind(0.d0)) :: t_temp

          Ncheck = 1
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
            Do nc = 1,Ncheck
              Call Op_make(Op_T(nc,n),Ndim)
              If ( Lattice_type =="BipartiteSquare" .and. (Model == "NematicDirac" .or. Model == "NematicDirac3") ) then
                DO I = 1, Latt%N
                  do no = 1,Norb
                    I1 = Invlist(I,no)
                    Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                  enddo
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
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0) * (1 + ham_xi2 * Model_sign)
                    case (4)
                      J1 = invlist(Latt%nnlist(I,-1,0),2) 
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0) * (1 - ham_xi2 * Model_sign)
                    case default
                      Write(6,*) ' Error in  Ham_Hop '  
                    end select
                    Op_T(nc,n)%O(I1,J1) = t_temp
                    Op_T(nc,n)%O(J1,I1) = conjg(t_temp)
                  Enddo
                Enddo
              ElseIf ( Lattice_type =="Square" .and. Model == "NematicDirac2" ) then
                DO I = 1, Latt%N
                  do no = 1,Norb
                    I1 = Invlist(I,no)
                    Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                  enddo
                  I1 = Invlist(I,1)
                  Do nc1 = 1,N_coord
                    select case (nc1)
                    case (1)
                      J1 = invlist(Latt%nnlist(I, 1, 0),2) 
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    case (2)
                      J1 = invlist(Latt%nnlist(I,-1, 0),2) 
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    case (3)
                      J1 = invlist(Latt%nnlist(I, 0, 1),2) 
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    case (4)
                      J1 = invlist(Latt%nnlist(I, 0,-1),2) 
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    case default
                      Write(6,*) ' Error in  Ham_Hop '  
                    end select
                    Op_T(nc,n)%O(I1,J1) = t_temp
                    Op_T(nc,n)%O(J1,I1) = conjg(t_temp)
                  Enddo
                Enddo
              elseif ( Lattice_type =="Square" .and. Model == "yyhe" ) then
                DO I = 1, Latt%N
                  do no = 1,Norb
                    I1 = Invlist(I,no)
                    Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                  enddo
                  I1 = Invlist(I,1)
                  Do nc1 = 1,N_coord
                    select case (nc1)
                    case (1)
                      J1 = invlist(I,2) 
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    case (2)
                      J1 = invlist(Latt%nnlist(I,-1, 0),2) 
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    case (3)
                      J1 = invlist(Latt%nnlist(I, 0,-1),2)
                      t_temp = -Ham_T * cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    case (4)
                      J1 = invlist(Latt%nnlist(I,-1,-1),2) 
                      t_temp = -Ham_T * cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    case default
                      Write(6,*) ' Error in  Ham_Hop '  
                    end select
                    If( n == 1 ) then
                      Op_T(nc,n)%O(I1,J1) = t_temp
                      Op_T(nc,n)%O(J1,I1) = conjg(t_temp)
                    else
                      Op_T(nc,n)%O(I1,J1) = conjg(t_temp)
                      Op_T(nc,n)%O(J1,I1) = t_temp
                    endif
                  Enddo
                Enddo
                
              else
                Write(6,*) ' This hopping is not yet implemented '
                Stop
              endif

                Do I = 1,Ndim
                   Op_T(nc,n)%P(i) = i 
                Enddo
                if ( abs(Ham_T) < 1.E-6 .and.  abs(Ham_chem) < 1.E-6 ) then 
                   Op_T(nc,n)%g = 0.d0
                else
                   Op_T(nc,n)%g = -Dtau
                endif
                Op_T(nc,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                Call Op_set(Op_T(nc,n)) 
             enddo
          enddo
        end Subroutine Ham_hop

!===================================================================================           
        Subroutine Ham_V
          
          Implicit none 
          
          Integer :: nf, I, nc1, I1
          complex (Kind=Kind(0.d0)) :: t_temp

          
          If  ( Lattice_type =="BipartiteSquare" .and. Model == "NematicDirac" ) then
            Allocate(Op_V(Latt%N,N_FL))
            do nf = 1,N_FL
              do I  =  1, Latt%N
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
                    t_temp = -cmplx(1, -1, kind(0.D0))/sqrt(2.D0) * Model_sign
                  case (4)
                    Op_V(I,nf)%P(nc1+1) = invlist(Latt%nnlist(I,-1,0),2) 
                    t_temp =  cmplx(1,  1, kind(0.D0))/sqrt(2.D0) * Model_sign
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
          elseIf  ( Lattice_type =="BipartiteSquare" .and. Model == "NematicDirac3" ) then
            Allocate(Op_V(2*Latt%N,N_FL))
            do nf = 1,N_FL
              do I  =  1, Latt%N
                call Op_make(Op_V(I       ,nf),3)
                call Op_make(Op_V(I+Latt%N,nf),3)
                Op_V(I       ,nf)%P(1) = Invlist(I,1)
                Op_V(I+Latt%N,nf)%P(1) = Invlist(I,1)
                Do nc1 = 1,4
                  select case (nc1)
                  case (1)
                    Op_V(I,nf)%P(2) = invlist(I,2) 
                    t_temp = -cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    Op_V(I,nf)%O(1,2) = t_temp
                    Op_V(I,nf)%O(2,1) = conjg(t_temp)
                  case (2)
                    Op_V(I+Latt%N,nf)%P(2) = invlist(Latt%nnlist(I, 0,1),2) 
                    t_temp =  cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
                    Op_V(I+Latt%N,nf)%O(1,2) = t_temp
                    Op_V(I+Latt%N,nf)%O(2,1) = conjg(t_temp)
                  case (3)
                    Op_V(I,nf)%P(3) = invlist(Latt%nnlist(I,-1,1),2) 
                    t_temp =  cmplx(1, -1, kind(0.D0))/sqrt(2.D0)
                    Op_V(I,nf)%O(1,3) = t_temp
                    Op_V(I,nf)%O(3,1) = conjg(t_temp)
                  case (4)
                    Op_V(I+Latt%N,nf)%P(3) = invlist(Latt%nnlist(I,-1,0),2) 
                    t_temp = -cmplx(1,  1, kind(0.D0))/sqrt(2.D0)
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
          elseIf  ( Lattice_type =="Square" .and. Model == "NematicDirac2" ) then
            Allocate(Op_V(Latt%N,N_FL))
            do nf = 1,N_FL
              do I  =  1, Latt%N
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
           elseIf  ( Lattice_type =="Square" .and. Model == "yyhe" ) then
            Allocate(Op_V(Latt%N,N_FL))
            do nf = 1,N_FL
              do I  =  1, Latt%N
                !no = 1
                  I1 = Invlist(I,1)
                  call Op_make(Op_V(I1,nf),4)
                  Op_V(I1,nf)%P(1) = invlist(I,1) 
                  Op_V(I1,nf)%P(2) = invlist(Latt%nnlist(I, 0,1),1) 
                  Op_V(I1,nf)%P(3) = invlist(I,2) 
                  Op_V(I1,nf)%P(4) = invlist(Latt%nnlist(I,-1,0),2) 
                  
                  Op_V(I1,nf)%O(1,2) = 1  
                  Op_V(I1,nf)%O(2,1) = 1  
                  Op_V(I1,nf)%O(3,4) = 1  
                  Op_V(I1,nf)%O(4,3) = 1  
                  
                  Op_V(I1,nf)%g      = cmplx(-dtau*ham_t*ham_xi,0.d0, kind(0.D0))
                  Op_V(I1,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                  Op_V(I1,nf)%type   = 1
                  Call Op_set( Op_V(I1,nf) )
                  
                !no = 2
                  I1 = Invlist(I,2)
                  call Op_make(Op_V(I1,nf),4)
                  Op_V(I1,nf)%P(1) = invlist(I,2) 
                  Op_V(I1,nf)%P(2) = invlist(Latt%nnlist(I, 0,1),2) 
                  Op_V(I1,nf)%P(3) = invlist(Latt%nnlist(I, 1,1),1) 
                  Op_V(I1,nf)%P(4) = invlist(Latt%nnlist(I, 0,1),1) 
                  
                  Op_V(I1,nf)%O(1,2) = -1  
                  Op_V(I1,nf)%O(2,1) = -1  
                  Op_V(I1,nf)%O(3,4) = -1  
                  Op_V(I1,nf)%O(4,3) = -1  
                  
                  Op_V(I1,nf)%g      = cmplx(-dtau*ham_t*ham_xi,0.d0, kind(0.D0))
                  Op_V(I1,nf)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
                  Op_V(I1,nf)%type   = 1
                  Call Op_set( Op_V(I1,nf) )
              Enddo
            Enddo
          else
             Write(6,*) ' This interaction is not yet implemented '
             Stop

          Endif
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
        include "NematicDirac/Setup_Ising_action.f"
!=================================================================================== 
        include "NematicDirac/Alloc_obs.f"
!========================================================================
        include "NematicDirac/Global.f"
!========================================================================
        include "NematicDirac/Obser.f"
!=====================================================
        include "NematicDirac/ObserT.f"
!==========================================================        
        Subroutine  Pr_obs(LTAU)

          Implicit none

          Integer,  Intent(In) ::  Ltau
          
          !Local 
          Integer :: I
          
!           Write(6,*) "Pr_obs"


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
          
!           Write(6,*) "Init_obs"
          nBlub = 0

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
          
          i_sweep = 1

        end Subroutine Init_obs

      end Module Hamiltonian
