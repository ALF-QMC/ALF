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
      use iso_fortran_env, only: output_unit, error_unit


      Implicit none

!>    Public variables. Have to be set by user
      Type (Operator), dimension(:,:), allocatable  :: Op_V
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Type (WaveFunction), dimension(:),   allocatable :: WF_L
      Type (WaveFunction), dimension(:),   allocatable :: WF_R
      Type (Fields)        :: nsigma
      Integer              :: Ndim
      Integer              :: N_FL
      Integer              :: N_SUN
      Integer              :: Ltrot
      Integer              :: Thtrot
      Logical              :: Projector
      Integer              :: Group_Comm
      Logical              :: Symm


!>    Privat variables
      Type (Lattice),        private :: Latt
      Integer,               private :: L1, L2
      real (Kind=Kind(0.d0)),private :: ham_T, Ham_chem, Ham_g, Ham_J,  Ham_K, Ham_h,  Ham_TZ2
      real (Kind=Kind(0.d0)),private :: Dtau, Beta
      Character (len=64),    private :: Model, Lattice_type
      Logical,               private :: One_dimensional
      Integer,               private :: N_coord, Norb
      Integer, allocatable,  private :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell



      !>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau

      !>    Storage for the Ising action
      Real (Kind=Kind(0.d0)), private :: DW_Ising_tau(-1:1), DW_Ising_Space(-1:1), DW_Ising_Flux(-1:1,-1:1)
      Real (Kind=Kind(0.d0)), private :: DW_Matter_tau(-1:1), DW_Ising_Matter(-1:1)
      Integer, allocatable  , private :: L_bond(:,:,:), L_bond_inv(:,:)


    contains


      Subroutine Ham_Set


#if defined (MPI) || defined(TEMPERING)
           use mpi
#endif
          Implicit none


          integer :: ierr
          Character (len=64) :: file1

          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model


          NAMELIST /VAR_Z2_Matter/ ham_T, Ham_chem, Ham_g, Ham_J,  Ham_K, Ham_h, &
               &                   Dtau, Beta, ham_TZ2


#ifdef MPI
          Integer        :: Isize, Irank, igroup, irank_g, isize_g
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
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'Ham_set: unable to open <parameters>',ierr
                error stop 1
             END IF
             READ(5,NML=VAR_lattice)
             CLOSE(5)
#ifdef MPI
          Endif
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif
          Call Ham_latt

          if ( Model == "Z2_Matter" ) then
             N_FL = 1
             N_SUN = 2
             If ( Lattice_type == "Honeycomb" ) then
                Write(error_unit,*) "Ham_set: Z2_Matter is only implemented for a square lattice"
                error stop 1
             Endif
          else
             Write(error_unit,*) "Ham_set: Model not yet implemented!"
             error stop 1
          endif


#if defined(TEMPERING)
          write(File1,'(A,I0,A)') "Temp_",igroup,"/parameters"
          OPEN(UNIT=5,File=file1,STATUS='old',ACTION='read',IOSTAT=ierr)
          ham_T = 0.d0; Ham_chem = 0.d0; Ham_g = 0.d0; Ham_J = 0.d0
          Ham_K = 0.d0; Ham_h = 0.d0
          READ(5,NML=VAR_Z2_Matter)
          CLOSE(5)
#else
#if defined(MPI)
          If (Irank == 0 ) then
#endif
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             ham_T = 0.d0; Ham_chem = 0.d0; Ham_g = 0.d0; Ham_J = 0.d0
             Ham_K = 0.d0; Ham_h = 0.d0
             READ(5,NML=VAR_Z2_Matter)
             CLOSE(5)
#ifdef MPI
          endif
          CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_TZ2  ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(ham_g    ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(Ham_J    ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(Ham_K    ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(Ham_h    ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,Group_Comm,ierr)
          CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,Group_Comm,ierr)
#endif
#endif

           Call Ham_hop
           Ltrot = nint(beta/dtau)

           If  ( Model == "Z2_Matter" )  Call Setup_Ising_action

#if defined(TEMPERING)
           write(File1,'(A,I0,A)') "Temp_",igroup,"/info"
#else
           File1 = "info"
#endif

#if defined(MPI) && !defined(TEMPERING)
           If (Irank == 0 ) then
#endif

              Open (Unit = 50,file=file1,status="unknown",position="append")
              Write(50,*) '====================================='
              Write(50,*) 'Model is      : ', Model
              Write(50,*) 'Lattice is    : ', Lattice_type
              Write(50,*) '# of orbitals : ', Ndim
              Write(50,*) 'Beta          : ', Beta
              Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
              Write(50,*) 'N_SUN         : ', N_SUN
              Write(50,*) 'N_FL          : ', N_FL
              Write(50,*) 't_fermion     : ', Ham_T
              Write(50,*) 't_Z2          : ', Ham_TZ2
              Write(50,*) 'J_Gauge_Z2    : ', Ham_J
              Write(50,*) 'h_Matter      : ', Ham_h
              Write(50,*) 'g_Z2          : ', Ham_g
              Write(50,*) 'Ham_chem      : ', Ham_chem
              Write(50,*) 'K_Gauge       : ', Ham_K
              close(50)
#if defined(MPI) && !defined(TEMPERING)
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

          If ( Lattice_type =="Square" ) then
             Norb = 1
             N_coord   = 2
             One_dimensional = .false.
             a1_p(1) =  1.0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.0  ; a2_p(2) =  1.d0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             If ( L1 == 1 .or. L2 == 1 ) then
                One_dimensional = .true.
                N_coord   = 1
                If (L1 == 1 ) then
                   Write(error_unit,*) 'Ham_Latt: For one dimensional systems set  L2 = 1 '
                   error stop 1
                endif
             endif
          else
             Write(error_unit,*) "Ham_Latt: Lattice not yet implemented!"
             error stop 1
          endif


          ! This is for the orbital structure.
          Ndim = Latt%N*Norb
          Allocate (List(Ndim,2), Invlist(Latt%N,Norb))
          nc = 0
          Do I = 1,Latt%N
             Do no = 1,Norb
                ! For the Honeycomb lattice no = 1,2 corresponds to the A,and B sublattices.
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

          Ncheck = 1
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
             Do nc = 1,Ncheck
                Call Op_make(Op_T(nc,n),Ndim)
                If ( Lattice_type =="Square" ) then
                   If (One_dimensional ) then
                      DO I = 1, Latt%N
                         I1 = Latt%nnlist(I, 1, 0)
                         Op_T(nc,n)%O(I,I1) = cmplx(0.d0, 0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I1,I) = cmplx(0.d0, 0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                      ENDDO
                   else
                      DO I = 1, Latt%N
                         I1 = Latt%nnlist(I,1,0)
                         I2 = Latt%nnlist(I,0,1)
                         Op_T(nc,n)%O(I,I1) = cmplx(0.d0,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I1,I) = cmplx(0.d0,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I,I2) = cmplx(0.d0,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I2,I) = cmplx(0.d0,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                      ENDDO
                   endif
                endif

                Do I = 1,Ndim
                   Op_T(nc,n)%P(i) = i
                Enddo
                if ( abs(Ham_chem) < 1.E-6 ) then
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

          Integer :: nf, I, I1, I2,  nc, nc1,  J, N_Bond_type
          Real (Kind=Kind(0.d0)) :: X

          If  ( Model == "Z2_Matter" ) then
             Allocate(Op_V(2*N_coord*Ndim +1, N_FL))
             do nf = 1,N_FL
                do i  =  1, 2*N_coord*Ndim
                   call Op_make(Op_V(i,nf),2)
                enddo
                call Op_make(Op_V(2*N_coord*Ndim + 1 ,nf),1 )
                ! This is for the reference spin on lattice site I = Latt%N
             enddo
             Do nc = 1,2*Ndim*N_coord   ! Runs over bonds.  Coordination number = 2.
                ! For the square lattice Ndim = Latt%N
                I1 = L_bond_inv(nc,1)
                if ( L_bond_inv(nc,2)  == 1 ) I2 = Latt%nnlist(I1,1,0)
                if ( L_bond_inv(nc,2)  == 2 ) I2 = Latt%nnlist(I1,0,1)
                Op_V(nc,1)%P(1) = I1
                Op_V(nc,1)%P(2) = I2
                Op_V(nc,1)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0))
                Op_V(nc,1)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
                If (L_bond_inv(nc,3) == 1 ) then
                   Op_V(nc,1)%g      = cmplx(dtau*ham_TZ2,0.D0,kind(0.D0))
                else
                   Op_V(nc,1)%g      = cmplx(dtau*Ham_t  ,0.D0,kind(0.D0))
                endif
                Op_V(nc,1)%alpha  = cmplx(0d0  ,0.d0, kind(0.D0))
                Op_V(nc,1)%type   = 1
                Call Op_set( Op_V(nc,1) )
             Enddo
             I1 = Latt%N
             nc = 2*Ndim*N_coord  +  1
             Op_V(nc,1)%P(1)   = I1
             Op_V(nc,1)%O(1,1) = cmplx(1.d0,0.d0, kind(0.D0))
             Op_V(nc,1)%g      = cmplx(0.d0,0.d0, kind(0.D0))
             Op_V(nc,1)%alpha  = cmplx(0.d0,0.d0, kind(0.D0))
             Op_V(nc,1)%type   = 1
             Call Op_set( Op_V(nc,1) )
          Else
             Write(error_unit,*) 'Ham_V: No model'
             error stop 1
          Endif
        end Subroutine Ham_V

!===================================================================================
        Real (Kind=Kind(0.d0)) function S0(n,nt,Hs_new)
          Implicit none
          Integer, Intent(IN) :: n,nt
          Real (Kind=Kind(0.d0)), Intent(In) :: Hs_new

          !Local
          Integer :: nt1,I, F1,F2,I1,I2,I3,  n_orientation, n_m

          !> Ratio for local spin-flip  of gauge field only.
          S0 = 1.d0
          If ( Op_V(n,1)%type == 1 ) then

             !L_bond_inv(nc,1) = I1
             !L_bond_inv(nc,2) = n_orientation
             !L_bond_inv(nc,3) = N_Bond_type

             If (L_bond_inv(n,3) == 1 ) then
                I              = L_bond_inv(n,1)
                n_orientation  = L_bond_inv(n,2)
                n_m            = L_bond(I,n_orientation,2)
                S0 = S0* DW_Ising_Matter(nsigma%i(n,nt)*nsigma%i(n_m,nt) )  ! Coupling to matter field.
                nt1 = nt +1
                if (nt1 > Ltrot) nt1 = 1
                S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt1))
                nt1 = nt - 1
                if (nt1 < 1  ) nt1 = Ltrot
                S0 = S0*DW_Ising_tau(nsigma%i(n,nt)*nsigma%i(n,nt1))
                ! Magnetic flux term
                I1 = L_bond_inv(n,1)
                if ( L_bond_inv(n,2) == 1 ) then
                   !     I2
                   !     I1 I3
                   I2 = Latt%nnlist(I1,0,1 )
                   I3 = Latt%nnlist(I1,1,0 )
                   F1 = nsigma%i(n,nt)*nsigma%i(L_bond(I1,2,1),nt)* nsigma%i(L_bond(I2,1,1),nt)*nsigma%i(L_bond(I3,2,1),nt)
                   !     I1
                   !     I2 I3
                   I2 = Latt%nnlist(I1,0,-1)
                   I3 = Latt%nnlist(I1,1,-1)
                   F2 = nsigma%i(n,nt)*nsigma%i(L_bond(I2,1,1),nt)* nsigma%i(L_bond(I2,2,1),nt)*nsigma%i(L_bond(I3,2,1),nt)
                else
                   !    I3
                   !    I2  I1
                   I2 = Latt%nnlist(I1,-1,0 )
                   I3 = Latt%nnlist(I1,-1,1 )
                   F1 = nsigma%i(n,nt)*nsigma%i(L_bond(I2,1,1),nt)* nsigma%i(L_bond(I2,2,1),nt)*nsigma%i(L_bond(I3,1,1),nt)
                   !    I2
                   !    I1  I3
                   I2 = Latt%nnlist(I1,0,1)
                   I3 = Latt%nnlist(I1,1,0)
                   F2 = nsigma%i(n,nt)*nsigma%i(L_bond(I1,1,1),nt)* nsigma%i(L_bond(I2,1,1),nt)*nsigma%i(L_bond(I3,2,1),nt)
                endif
                S0 = S0*DW_Ising_Flux(F1,F2)
             else
                Write(error_unit,*) 'Error in funtion S0'
                error stop 1
             endif

          endif

        end function S0

!===================================================================================
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
!> The move will be carried out with prbablity  T0 ( sigma -> sigma_new ).   If T0 ( sigma -> sigma_new ) > Ranf
!>  then T0_Proposal_ratio  will be initialized. Otherwise the latter quantity is set to zero.
!--------------------------------------------------------------------

          Implicit none
          Real (Kind= kind(0.d0)),INTENT(OUT) :: T0_Proposal_ratio, S0_ratio
          Integer                ,INTENT(OUT) :: Flip_list(:)
          Real (Kind= Kind(0.d0)),INTENT(out) :: Flip_value(:)
          Integer, INTENT(OUT) :: Flip_length
          Integer, INTENT(IN)    :: ntau


          !Local
          Integer                   ::  ns , nc, n_op, n_op1, ntau_p1, ntau_m1, I, n
          Integer, allocatable      ::  Isigma1(:),Isigma2(:),Isigma3(:)
          Real  (Kind = Kind(0.d0)) ::  S0_Matter, T0_Proposal

          ! Write(6,*) 'In GLob_move', m,direction,ntau, size(Flip_list,1), Size(Flip_value,1), Flip_list(1)
          ! Ising from n_op = 1,N_coord*Ndim
          ! Hubbard from n_op = N_coord*Ndim +1, Size(OP_V,1) = N_coord*Ndim +  Ndim
          ! Write(6,*) 'Global_move_tau ' , S0(Flip_list(1),ntau)

          Allocate (Isigma1(Latt%N), Isigma2(Latt%N), Isigma3(Latt%N) )

          I  =  nranf(Latt%N)
          Flip_length = 4
          S0_Matter = 1.d0
          do n = 1,4
             select case(n)
             case (1)
                n_op  = L_bond(I,1,2)
                n_op1 = L_bond(I,1,1)
                S0_Matter = S0_Matter*DW_Ising_Matter( nsigma%i(n_op,ntau) *  nsigma%i(n_op1,ntau) )
             case (2)
                n_op  = L_bond(I,2,2)
                n_op1 = L_bond(I,2,1)
                S0_Matter = S0_Matter*DW_Ising_Matter( nsigma%i(n_op,ntau) *  nsigma%i(n_op1,ntau) )
             case (3)
                n_op  = L_bond(latt%nnlist(I,-1,0),1,2)
                n_op1 = L_bond(latt%nnlist(I,-1,0),1,1)
                S0_Matter = S0_Matter*DW_Ising_Matter( nsigma%i(n_op,ntau) *  nsigma%i(n_op1,ntau) )
             case (4)
                n_op  = L_bond(latt%nnlist(I,0,-1),2,2)
                n_op1 = L_bond(latt%nnlist(I,0,-1),2,1)
                S0_Matter = S0_Matter*DW_Ising_Matter( nsigma%i(n_op,ntau) *  nsigma%i(n_op1,ntau) )
             case default
                Write(error_unit,*) 'Global_move_tau: Error'
                error stop 1
             end select
             Flip_list(n)  = n_op
             Flip_value(n) = nsigma%flip(n_op,ntau)
          enddo
          If ( I == Latt%N )   then
             Flip_length   = 5
             n             = 5
             n_op          = 2*N_coord*Ndim + 1
             Flip_list(n)  = n_op
             Flip_value(n) = nsigma%flip(n_op,ntau)
          endif

          ntau_p1 = ntau + 1
          if (ntau == Ltrot) ntau_p1 = 1
          ntau_m1 = ntau -1
          if (ntau == 1    ) ntau_m1 = Ltrot
          Call Hamiltonian_set_Z2_matter(Isigma1,ntau_p1)
          Call Hamiltonian_set_Z2_matter(Isigma2,ntau   )
          Call Hamiltonian_set_Z2_matter(Isigma3,ntau_m1)
          !  Check the dynamics and the ergodicity
          S0_Matter = S0_Matter*DW_Matter_tau ( Isigma1(I)*Isigma2(I) ) * DW_Matter_tau( Isigma2(I)*Isigma3(I) )
          T0_Proposal       =  1.d0 - 1.d0/(1.d0+S0_Matter)
          !  Move acceptance probability.
          If ( T0_Proposal > Ranf_wrap() )  then
             T0_Proposal_ratio =  1.d0 / S0_Matter
             !T0_Proposal       =  1.d0
             !T0_Proposal_ratio =  1.d0
          else
             T0_Proposal_ratio = 0.d0
          endif
          S0_ratio          =  S0_Matter

          Deallocate (Isigma1,Isigma2, Isigma3)
!!$          Flip_length    = 1
!!$          n_op = nranf(size(OP_V,1))
!!$          Flip_list(1)   = n_op
!!$          If ( OP_V(n_op,1)%type == 1 ) then
!!$             ns = nsigma(n_op,ntau)
!!$             T0_Proposal       =  1.d0 - 1.d0/(1.d0+S0(n_op,ntau)) ! No move prob
!!$             T0_Proposal_ratio =  1.d0 / S0(n_op,ntau)
!!$             S0_ratio          =  S0(n_op,ntau)
!!$             Flip_value(1)     = - ns
!!$          else
!!$             Flip_value(1)     = NFLIPL(nsigma(n_op,ntau),nranf(3))
!!$             T0_Proposal       = 1.d0
!!$             T0_Proposal_ratio = 1.d0
!!$             S0_ratio          = 1.d0
!!$          endif

        end Subroutine Global_move_tau
!===================================================================================
        Subroutine Global_move(T0_Proposal_ratio,nsigma_old,size_clust)
          !>  The input is the field nsigma declared in this module. This routine generates a
          !>  global update with  and returns the propability
          !>  T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)
          !>

          Implicit none
          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
          type (Fields),  Intent(IN)  :: nsigma_old
          !> nsigma_old contains a copy of nsigma upon entry


        End Subroutine Global_move
!===================================================================================
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none

          !> Arguments
          type (Fields),  Intent(IN)  :: nsigma_old
          !> Local
          Integer :: I,n,n1,n2,n3,n4,nt,nt1, nc_F, nc_J, nc_h_p, nc_h_m, n1_m, n4_m


          Delta_S0_global = 1.d0
          If ( Model == "Z2_Matter" ) then
             nc_F = 0
             nc_J = 0
             nc_h_p = 0
             nc_h_m = 0
             Do I = 1,Latt%N
                n1   = L_bond(I,1,1)
                n1_m = L_bond(I,1,2)
                n2   = L_bond(Latt%nnlist(I,1,0),2,1)
                n3   = L_bond(Latt%nnlist(I,0,1),1,1)
                n4   = L_bond(I,2,1)
                n4_m = L_bond(I,2,2)
                do nt = 1,Ltrot
                   nt1 = nt +1
                   if (nt == Ltrot) nt1 = 1
                   if (nsigma%i(n1,nt) == nsigma%i(n1,nt1) ) then
                      nc_h_p = nc_h_p + 1
                   else
                      nc_h_m = nc_h_m + 1
                   endif
                   if (nsigma_old%i(n1,nt) == nsigma_old%i(n1,nt1) ) then
                      nc_h_p = nc_h_p - 1
                   else
                      nc_h_m = nc_h_m - 1
                   endif

                   if (nsigma%i(n4,nt) == nsigma%i(n4,nt1) ) then
                      nc_h_p = nc_h_p + 1
                   else
                      nc_h_m = nc_h_m + 1
                   endif
                   if (nsigma_old%i(n4,nt) == nsigma_old%i(n4,nt1) ) then
                      nc_h_p = nc_h_p - 1
                   else
                      nc_h_m = nc_h_m - 1
                   endif

                   nc_F = nc_F + nsigma%i    (n1,nt)*nsigma%i    (n2,nt)*nsigma%i    (n3,nt)*nsigma%i    (n4,nt)  &
                        &      - nsigma_old%i(n1,nt)*nsigma_old%i(n2,nt)*nsigma_old%i(n3,nt)*nsigma_old%i(n4,nt)

                   nc_J = nc_J + nsigma%i(n1,nt)*nsigma%i(n1_m,nt) + &
                        &        nsigma%i(n4,nt)*nsigma%i(n4_m,nt) - &
                        &        nsigma_old%i(n1,nt)*nsigma_old%i(n1_m,nt) - &
                        &        nsigma_old%i(n4,nt)*nsigma_old%i(n4_m,nt)

                enddo
             enddo
             Delta_S0_global = ( sinh(Dtau*Ham_h)**nc_h_m ) * (cosh(Dtau*Ham_h)**nc_h_p) * &
                  &            exp( -Dtau*(Ham_K*real(nc_F,kind(0.d0)) + Ham_J*real(nc_J,kind(0.d0))))
          endif
        end Function Delta_S0_global
!===================================================================================
        Subroutine Setup_Ising_action

          ! This subroutine sets up lists and arrays so as to enable an
          ! an efficient calculation of  S0(n,nt)

          Integer :: nc, nth, n, n1, n2, n3, n4, I, I1, n_orientation, Ix, Iy, N_Bond_type
          Real (Kind=Kind(0.d0)) :: X_p(2)

          ! Setup list of bonds for the square lattice.
          Allocate ( L_Bond(Latt%N,2,2),  L_bond_inv(2*Latt%N*N_coord,3) )
          nc = 0
          Do N_Bond_type = 1,2  ! Two types of bonds for gauge Fields (N_Bond_Type = 1) and Matter fields (N_Bond_type = 2)
             DO I = 1,Latt%N
                Ix = Latt%list(I,1)
                Iy = Latt%list(I,2)
                if (mod(Ix + Iy,2) == 0 ) then
                   do n = 1,4
                      nc = nc + 1
                      select case (n)
                      case (1)
                         I1 = I                  ;  n_orientation  = 1
                      case (2)
                         I1 = I                  ;  n_orientation  = 2
                      case (3)
                         I1 = latt%nnlist(I,-1,0);  n_orientation  = 1
                      case (4)
                         I1 = latt%nnlist(I,0,-1);  n_orientation  = 2
                      case default
                         Write(6,*) ' Error in Setup_Ising_action '
                      end select
                      L_bond(I1,n_orientation,N_Bond_type) = nc
                      L_bond_inv(nc,1) = I1
                      L_bond_inv(nc,2) = n_orientation
                      L_bond_inv(nc,3) = N_Bond_type
                      ! The bond is given by  I1, I1 + a_(n_orientation).
                   enddo
                endif
             Enddo
          Enddo
          !Test
          !Do I = 1,Latt%N
          !   Write(6,*)
          !   Write(6,*) Latt%list(I,1), Latt%list(I,2), I
          !   Write(6,*) L_bond(I,1), L_bond(I,2), L_bond( latt%nnlist(I,-1,0),1), L_bond( latt%nnlist(I,0,-1),2)
          !Enddo
          DW_Ising_tau  ( 1) = tanh(Dtau*Ham_g)
          DW_Ising_tau  (-1) = 1.D0/DW_Ising_tau(1)
          DO n = -1,1,2
             do n1 = -1,1,2
                DW_Ising_Flux(n,n1) = exp( Dtau*Ham_K*(dble(n) +  dble(n1) ))/ exp(  -Dtau*Ham_K*(dble(n) +  dble(n1) ))
             enddo
          enddo
          DW_Ising_Matter( 1) = exp( 2.d0*Dtau*Ham_J)
          DW_Ising_Matter(-1) = exp(-2.d0*Dtau*Ham_J)
          DW_Matter_tau  ( 1) = tanh(Dtau*Ham_h)
          DW_Matter_tau  (-1) = 1.D0/DW_Matter_tau(1)

        End Subroutine Setup_Ising_action
!===================================================================================
        Subroutine  Alloc_obs(Ltau)

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Ns,Nt,No
          Character (len=64) ::  Filename

          ! Scalar observables
          Allocate ( Obs_scal(10) )
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
                N = 2;   Filename ="Flux"
             case (6)
                N = 2;   Filename ="X"
             case (7)
                N = L1/2
                If (L2 > L1 ) N = L2/2; Filename ="Wilson"
             case (8)
                N = L1/2 + L2/2;        Filename ="Vison"
             case (9)
                N = 1;        Filename ="KinZ2"
             case (10)
                N = 1;        Filename ="J"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo


          ! Equal time correlators
          Allocate ( Obs_eq(6) )
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
                Ns = Latt%N;  No = Norb;  Filename ="GreenZ2"
             case (6)
                Ns = Latt%N;  No = N_coord;  Filename ="Kin"
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
                Nt = Ltrot+1
                Call Obser_Latt_make(Obs_tau(I),Ns,Nt,No,Filename)
             enddo
          endif
        end Subroutine Alloc_obs

!========================================================================
        Subroutine Obser(GR,Phase,Ntau)


#if defined(Machine_Learning)
          use  mpi
#endif
          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau

          !Local
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, Z1, Z2, ZN
          Complex (Kind=Kind(0.d0)) :: ZQ, ZSTAR, ZQT, ZQTT
          Integer :: I,J, imj, nf, dec, I1, I2,I3,I4, J1, no_I, no_J,  iFlux_tot,  &
               &     no, no1, ntau1, L_Vison, L_Wilson, n, nx,ny
          Real (Kind=Kind(0.d0)) :: X_ave, X, XI1,XI2,XI3,XI4
          Integer,  allocatable  :: Isigma(:), Isigma1(:)
          Integer ::  IB_x, IB_y, Ix, Iy


#if defined(Machine_Learning)
          Character (len=64) :: File1
          Integer        :: Isize, Irank, IERR
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

          ZN =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))

          Allocate ( Isigma(Latt%N), Isigma1(Latt%N) )
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

          ntau1 = ntau + 1
          If (ntau == Ltrot)  ntau1 = 1
          Call Hamiltonian_set_Z2_matter(Isigma ,ntau )
          Call Hamiltonian_set_Z2_matter(Isigma1,ntau1)

          Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Latt%N
                I1 = latt%nnlist(I,1,0)
                I2 = latt%nnlist(I,0,1)
                Z1 = cmplx(real(Isigma(I)*Isigma(I1), kind(0.d0)), 0.d0,kind(0.d0))
                Z2 = cmplx(real(Isigma(I)*Isigma(I2), kind(0.d0)), 0.d0,kind(0.d0))
                Zkin = Zkin + Z1 * ( GRC(I,I1,1) + GRC(I1,I,1) ) + &
                     &        Z2 * ( GRC(I,I2,1) + GRC(I2,I,1) )
             Enddo
          Enddo
          Zkin = Zkin * dble(N_SUN)*Ham_t
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          X_ave = 0.d0
          Do I = 1,Latt%N
             X_ave = X_ave + DW_Ising_tau( Isigma(I)*Isigma1(I) )
          Enddo
          Zpot = cmplx(-X_ave*Ham_h,0.d0,kind(0.d0))
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


          iFlux_tot = 0
          Do I = 1, Ndim
             iFlux_tot = iFlux_tot + iFlux(I,Ntau,1)
          Enddo
          Obs_scal(5)%Obs_vec(1)  =   Obs_scal(5)%Obs_vec(1) + cmplx(dble(iFlux_tot),0.d0,kind(0.d0))*ZP*ZS
          iFlux_tot = 0
          Do I = 1, Ndim
             iFlux_tot = iFlux_tot + iFlux(I,Ntau,2)
          Enddo
          Obs_scal(5)%Obs_vec(2)  =   Obs_scal(5)%Obs_vec(2) + cmplx(dble(iFlux_tot),0.d0,kind(0.d0))*ZP*ZS


          ntau1 = ntau + 1
          If (ntau == Ltrot)  ntau1 = 1
          X_ave = 0.d0
          Do I = 1,Latt%N
             do no = 1,2
                X_ave = X_ave + DW_Ising_tau( nsigma%i(L_bond(I,no,1),ntau)*nsigma%i(L_bond(I,no,1),ntau1) )
             Enddo
          Enddo
          Obs_scal(6)%Obs_vec(1)  =  Obs_scal(6)%Obs_vec(1) + cmplx(X_ave,0.d0,kind(0.d0)) * ZP*ZS

          X_ave = 0.d0
          Do I = 1,Latt%N
             X_ave = X_ave + DW_Matter_tau( Isigma(I)*Isigma1(I) )
          Enddo
          Obs_scal(6)%Obs_vec(2)  =  Obs_scal(6)%Obs_vec(2) + cmplx(X_ave,0.d0,kind(0.d0)) * ZP*ZS


          Do I = 1,Latt%N
             do L_Wilson = 1,Size(Obs_scal(7)%Obs_vec,1)
                X  = 1.d0
                I1 = I
                I2 = I
                Do n = 1,L_Wilson
                   X = X *nsigma%i(L_bond(I1,1,1),ntau)
                   I1 = Latt%nnlist(I1,1,0)
                   X = X *nsigma%i(L_bond(I2,2,1),ntau)
                   I2 = Latt%nnlist(I2,0,1)
                enddo
                Do n = 1,L_Wilson
                   X = X *nsigma%i(L_bond(I1,2,1),ntau)
                   I1 = Latt%nnlist(I1,0,1)
                   X = X* nsigma%i(L_bond(I2,1,1),ntau)
                   I2 = Latt%nnlist(I2,1,0)
                enddo
                Obs_scal(7)%Obs_vec(L_Wilson) = Obs_scal(7)%Obs_vec(L_Wilson) + cmplx(X,0.d0,kind(0.d0)) *  ZP * ZS
             enddo
          Enddo

          Do I = 1, Latt%N
             I1 = I
             X  = 1.d0
             DO L_Vison = 1,L1/2
                X = X*DW_Ising_tau( nsigma%i(L_bond(I1 ,2,1),ntau)*nsigma%i(L_bond(I1 ,2,1),ntau1) )
                Obs_scal(8)%Obs_vec(L_Vison) =   Obs_scal(8)%Obs_vec(L_Vison) + cmplx(X,0.d0,kind(0.d0)) * ZP * ZS
                I1 = Latt%nnlist(I1,1,0)
             Enddo
             I1 = Latt%nnlist(I1,0,1)
             DO L_Vison = L1/2 + 1, L1/2 + L2/2
                X = X*DW_Ising_tau( nsigma%i(L_bond(I1,1,1),ntau)*nsigma%i(L_bond(I1 ,1,1),ntau1) )
                Obs_scal(8)%Obs_vec(L_Vison) =   Obs_scal(8)%Obs_vec(L_Vison) + cmplx(X,0.d0,kind(0.d0)) * ZP * ZS
                I1 = Latt%nnlist(I1,0,1)
             Enddo
          Enddo

          X  = 0.d0
          Z  = cmplx(0.d0,0.d0,kind(0.d0))
          Do I = 1, Latt%N
             I1 = latt%nnlist(I,1,0)
             I2 = latt%nnlist(I,0,1)
             X = X   +  Real( nsigma%i(L_bond(I,1,1),ntau)* Isigma(I)*Isigma(I1), kind=kind(0.d0) ) + &
                  &     Real( nsigma%i(L_bond(I,2,1),ntau)* Isigma(I)*Isigma(I2), kind=kind(0.d0) )
             Z = Z   + cmplx(Real( nsigma%i(L_bond(I,1,1),ntau),kind(0.d0)),0.d0,kind(0.d0)) * ( GRC(I,I1,1) + GRC(I1,I,1) )  &
                  &  + cmplx(Real( nsigma%i(L_bond(I,2,1),ntau),kind(0.d0)),0.d0,kind(0.d0)) * ( GRC(I,I2,1) + GRC(I2,I,1) )
          enddo
          Obs_scal (9)%Obs_vec(1) = Obs_scal(9 )%Obs_vec(1)  + Z*cmplx(Ham_TZ2*dble(N_SUN),0.d0,kind(0.d0))*ZP*ZS
          Obs_scal(10)%Obs_vec(1) = Obs_scal(10)%Obs_vec(1)  + cmplx(X*Ham_J  ,0.d0,kind(0.d0))*ZP*ZS


!!$          ZQ = cmplx(0.d0,0.d0,kind(0.d0))
!!$          DO I = 1,Latt%N
!!$             Z    =  ( cmplx(1.d0,0.d0,kind(0.d0)) - cmplx(2.d0,0.d0,kind(0.d0))*GRC(I,I,1) )**(N_SUN)
!!$             ZQ   = ZQ    + cmplx(DW_Ising_tau( Isigma(I)*Isigma1(I)),0.d0,kind(0.d0) ) * Z
!!$          Enddo
!!$          Obs_scal(6)%Obs_vec(1) = Obs_scal(6)%Obs_vec(1)  + ZQ * ZP*ZS


          ! Compute spin-spin, Green, and den-den correlation functions  !  This is general N_SUN, and  N_FL = 1
          DO I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO
          Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
          Do I1 = 1,Ndim
             I    = List(I1,1)
             no_I = List(I1,2)
             Do J1 = 1,Ndim
                J    = List(J1,1)
                no_J = List(J1,2)
                imj = latt%imj(I,J)
                ! Green_fermion
                Z1 = cmplx(real(Isigma(I1)*Isigma(J1), kind(0.d0)), 0.d0,kind(0.d0))

                Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                     &               Z * Z1*GRC(I1,J1,1) *  ZP*ZS

                ! Green_Z2
                Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) + &
                     &               Z * GRC(I1,J1,1) *  ZP*ZS

                ! SpinZ
                Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
                     &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                ! SpinXY
                Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
                     &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                ! Den
                Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J)  +  &
                     &     (    GRC(I1,I1,1) * GRC(J1,J1,1) *Z     + &
                     &          GRC(I1,J1,1) * GR(I1,J1,1 )          &
                     &                                     ) * Z* ZP*ZS




!!$             Z1 =   (cmplx(1.d0,0.d0,kind(0.d0)) - cmplx(2.d0,0.d0,kind(0.d0))*GRC(I1,I1,1)) *  &
!!$                  & (cmplx(1.d0,0.d0,kind(0.d0)) - cmplx<(2.d0,0.d0,kind(0.d0))*GRC(J1,J1,1)) +  &
!!$                  &  cmplx(4.d0,0.d0,kind(0.d0)) * GRC(I1,J1,1)*GR(I1,J1,1)
!!$             Z1 = Z1**(N_SUN)
!!$             ZQ = cmplx(DW_Ising_tau( Isigma(I1)*Isigma1(I1))*DW_Ising_tau( Isigma(J1)*Isigma1(J1)),0.d0,kind(0.d0) )*Z1
!!$             If (I1 == J1 )  ZQ = cmplx(1.d0,0.d0,kind(0.d0))
!!$             Obs_eq(6)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(6)%Obs_Latt(imj,1,no_I,no_J) + ZQ*ZP*ZS


             ENDDO
             Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  Z * GRC(I1,I1,1) * ZP * ZS
          ENDDO

          Do I = 1,Latt%N
             do no = 1,N_coord
                select  case(no)
                case(1)
                   I1 = Latt%nnlist(I,1,0)
                case(2)
                   I1 = Latt%nnlist(I,0,1)
                end select
                Do J = 1,Latt%N
                   imj = latt%imj(I,J)
                   do no1 = 1,N_coord
                      if ( no1 == 1) J1 = Latt%nnlist(J,1,0)
                      if ( no1 == 2) J1 = Latt%nnlist(J,0,1)

                      ! Kin-Kin
                      Z =  (  (GRC(I,I1,1) +  GRC(I1,I,1)) * (GRC(J,J1,1)  +  GRC(J1,J,1)) * ZN  + &
                           &   GRC(I ,J1,1)*GR(I1,J ,1) + &
                           &   GRC(I ,J ,1)*GR(I1,J1,1) + &
                           &   GRC(I1,J1,1)*GR(I ,J ,1) + &
                           &   GRC(I1,J ,1)*GR(I ,J1,1)   )* ZN * &
                           &   real(nsigma%i(L_bond(I ,no,1),ntau)*nsigma%i(L_bond(J ,no1,1),ntau), Kind(0.d0))

                      Obs_eq(6)%Obs_Latt(imj,1,no,no1)  = Obs_eq(6)%Obs_Latt(imj,1,no,no1)   +  Z * ZP*ZS


                   Enddo
                Enddo
                Z = (GRC(I,I1,1) +  GRC(I1,I,1))*ZN*Real(nsigma%i(L_bond(I ,no,1),ntau),kind(0.d0))
                Obs_eq(6)%Obs_Latt0(no) =  Obs_eq(6)%Obs_Latt0(no) +  Z * ZP*ZS
             enddo
          enddo

          Deallocate ( Isigma )

        end Subroutine Obser
!=====================================================
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE)
          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase

          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS
          Integer :: IMJ, I, J, I1, J1, no_I, no_J, NT1

          NT1 = NT
          If (NT == 0 ) NT1 = LTROT
          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          If (NT == 0 ) then
             DO I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             ENDDO
          endif
          If ( Model == "Z2_Matter" ) then
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

!===================================================================================

      Integer Function  iFlux(I,nt,nb_type)

        Implicit none

        Integer, INTENT(IN) :: I,nt, nb_type

        ! Local
        Integer :: n1,n2,n3,n4

        !   I3  I2
        !   I   I1
        n1  = L_bond(I,1,nb_type)
        n2  = L_bond(Latt%nnlist(I,1,0),2,nb_type)
        n3  = L_bond(Latt%nnlist(I,0,1),1,nb_type)
        n4  = L_bond(I,2,nb_type)
        iFlux =   nsigma%i(n1,nt)*nsigma%i(n2,nt)*nsigma%i(n3,nt)*nsigma%i(n4,nt)

      end Function iFlux

!===================================================================================
      Subroutine  Hamiltonian_set_nsigma(Initial_field)

        ! The user can set the initial configuration

        Implicit none

        Real (Kind=Kind(0.d0)), allocatable, dimension(:,:), Intent(OUT) :: Initial_field

        ! Local
        Integer :: I,nc, I1, nt, n_orientation
        Integer, allocatable::  Isigma(:), Isigma1(:)

        Allocate  (Initial_field(2*N_coord*Latt%N +1, Ltrot) )
        allocate  (Isigma(Latt%N), Isigma1(Latt%N) )
        Do nt = 1,Ltrot
           Do I = 1,Latt%N
              Isigma(I) = 1
              if ( ranf_wrap()  > 0.5D0 ) Isigma(I)  = -1
           enddo
           Initial_field(2*N_coord*Latt%N +1,nt) = real(Isigma(Latt%N), kind(0.d0))
           Do I = 1,Latt%N
              Do n_orientation = 1,2
                 nc = L_bond(I,n_orientation,2)
                 if (  n_orientation == 1 )  I1 = latt%nnlist(I,1,0)
                 if (  n_orientation == 2 )  I1 = latt%nnlist(I,0,1)
                 Initial_field(nc,nt) = real(Isigma(I)*Isigma(I1), kind(0.d0))
              enddo
           Enddo
           Initial_field(2*N_coord*Latt%N +1,nt)      = real(Isigma(Latt%N), kind(0.d0))
           Do I = 1, Latt%N
              if (mod( Latt%list(i,1) + latt%list(i,2), 2 ) == 0 ) then
                 Initial_field(L_bond(I,1,1),nt) =  1.d0
                 Initial_field(L_bond(I,2,1),nt) = -1.d0
              else
                 Initial_field(L_bond(I,1,1),nt) =  1.d0
                 Initial_field(L_bond(I,2,1),nt) =  1.d0
              endif
           Enddo


           do nc = 1,size(Initial_field,1)
              nsigma%f(nc,nt) = Initial_field(nc,nt)
           enddo
           Call Hamiltonian_set_Z2_matter(Isigma1,nt)
           Do nc = 1,Latt%N
              if ( Isigma(nc) .ne.  Isigma1(nc)  ) then
                 Write(error_unit,*) 'Error in Hamiltonian_set_Z2_matter'
                 error stop 1
              endif
           enddo

        enddo

        deallocate (Isigma, Isigma1)

        Call Print_fluxes

      end Subroutine Hamiltonian_set_nsigma
!===================================================================================
      Subroutine  Hamiltonian_set_Z2_matter(Isigma,nt)

        ! On input :  Link variables  nsigma(:,nt)
        ! On output:  The Z2_matter fields Isigma on the time slice.

        Implicit none

        Integer, Intent(IN)                  :: nt
        Integer, allocatable, INTENT(INOUT)  ::   Isigma(:)

        !Local
        Integer :: I, I1, nx, ny


        Isigma(Latt%N) = nsigma%i( 2*N_coord*Ndim + 1, nt )
        I = Latt%N
        do nx = 1,L1
           do ny = 1,L2
              I1 = latt%nnlist(I,0,1)
              Isigma(I1)  = Isigma(I)*nsigma%i(L_bond(I,2,2),nt)
              !Write(6,*) Latt%list(I,1), Latt%list(I,2), ' -> ', Latt%list(I1,1), Latt%list(I1,2)
              I = I1
           enddo
           I1          = latt%nnlist(I,1,0)
           Isigma(I1)  = Isigma(I)*nsigma%i(L_bond(I,1,2),nt)
           !Write(6,*) Latt%list(I,1), Latt%list(I,2), ' -> ', Latt%list(I1,1), Latt%list(I1,2)
           I = I1
        enddo

      end Subroutine Hamiltonian_set_Z2_matter

!===================================================================================
      Subroutine Hamiltonian_Print(Ntau)

        Integer, Intent(IN) :: Ntau

        Integer, allocatable :: Isigma(:)
        Integer :: I, Ix, Iy

        allocate (Isigma(Latt%N))

        Call Hamiltonian_set_Z2_matter(Isigma,ntau)

        Write(6,*)'-----'
        I = 1
        Do Iy = 1,L2
           Do Ix = 1,L1
              Write(6,"(I2,1x)", advance='no')  Isigma(I)
              I = Latt%nnlist(I,1,0)
           enddo
           Write(6,*)
           I = Latt%nnlist(I,0,1)
        enddo

        deallocate (Isigma)
      End Subroutine Hamiltonian_Print
!===================================================================================

      Subroutine Print_fluxes


#if defined (MPI) || defined(TEMPERING)
        use  mpi
#endif

        Implicit none


        ! Local
        Integer :: I,nt,ix, iy, n
        Character (len=64) :: File1


#ifdef MPI
        Integer        :: Isize, Irank, IERR, igroup, irank_g, isize_g
        Integer        :: STATUS(MPI_STATUS_SIZE)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        call MPI_Comm_rank(Group_Comm, irank_g, ierr)
        call MPI_Comm_size(Group_Comm, isize_g, ierr)
        igroup           = irank/isize_g
#endif

#if defined(TEMPERING)
        write(File1,'(A,I0,A)') "Temp_",igroup,"/Fluxes"
#else
        File1="Fluxes"
#endif

        Open (Unit=10,File=File1, status="unknown")
        Do nt = 1,Ltrot
           Do i  = 1,Ndim
              n = iFlux(I,nt,2)
              if (n == -1 ) then
                 ix = Latt%list(i,1)
                 iy = Latt%list(i,2)
                 Write(10,'(I4,2x,I4,2x,I4)')   IX, IY, NT
              endif
           Enddo
        enddo
        close(10)

      end Subroutine Print_fluxes
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

        Nt_sequential_start = 1
        Nt_sequential_end   = N_coord*Latt%N
        N_Global_tau        = Latt%N/4

      end Subroutine Overide_global_tau_sampling_parameters
!===================================================================================
!!$      Subroutine Test_Hamiltonian
!!$
!!$        Implicit none
!!$
!!$        Integer :: n,  nc, n_op, nt
!!$        Integer, allocatable :: nsigma_old(:,:)
!!$        Real (Kind=kind(0.d0)) :: X, X1, size_clust
!!$
!!$        n = size(Op_V,1)
!!$        allocate (nsigma_old(n,Ltrot))
!!$        do nc = 1,100
!!$           !nt  = nranf(Ltrot)
!!$           !n_op= nranf(n)
!!$           !if ( OP_V(n_op,1)%type == 1 ) then
!!$           !   X = S0(n_op,nt)
!!$           !   nsigma_old = nsigma
!!$           !   nsigma(n_op,nt) = -nsigma(n_op,nt)
!!$           !   X1 = Delta_S0_global(Nsigma_old)
!!$           !   Write(6,*) nc, X, X1
!!$           !endif
!!$           nsigma_old = nsigma
!!$           Call Global_move(X,nsigma_old,size_clust)
!!$           X1 = Delta_S0_global(Nsigma_old)
!!$           Write(6,*) nc, X, X1
!!$        enddo
!!$        deallocate (nsigma_old)
!!$
!!$        stop
!!$
!!$      end Subroutine Test_Hamiltonian

      end Module Hamiltonian
