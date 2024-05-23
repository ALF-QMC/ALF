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
      Type (WaveFunction), dimension(:),   allocatable  :: WF_L
      Type (WaveFunction), dimension(:),   allocatable  :: WF_R
      Type  (Fields)       :: nsigma
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot, Thtrot
      Logical              :: Projector
!>    Defines MPI communicator
      Integer              :: Group_Comm
      Logical              :: Symm =.false.


!>    Privat variables
      Type (Lattice),        private :: Latt
      Type (Unit_cell),      private :: Latt_unit
      Integer,               private :: L1, L2
      real (Kind=Kind(0.d0)),private :: ham_T , ham_U,  Ham_chem, Ham_h, Ham_J, Ham_xi, Ham_F
      real (Kind=Kind(0.d0)),private :: Dtau, Beta, Theta
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
      Integer, allocatable  , private :: L_bond(:,:), L_bond_inv(:,:), Ising_nnlist(:,:)


    contains


      Subroutine Ham_Set
#ifdef MPI
          Use mpi
#endif
          Implicit none


          integer :: ierr
          Character (len=64) :: file_info, file_para
          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model


          NAMELIST /VAR_Z2_Slave/  ham_T, ham_chem, ham_U, Dtau, Beta, &
               &                   Ham_h, Ham_J, Ham_xi, Ham_F


#ifdef MPI
          Integer        :: Isize, Irank, irank_g, isize_g, igroup
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
          call MPI_Comm_rank(Group_Comm, irank_g, ierr)
          call MPI_Comm_size(Group_Comm, isize_g, ierr)
          igroup           = irank/isize_g

#endif
          ! Global default values
          N_FL = 1; N_SUN = 2
          ham_T = 0.d0;  ham_U= 0.d0;   Ham_J= 0.d0;  Ham_F = 0.d0
#ifdef MPI
          If (Irank_g == 0 ) then
#endif
             File_para = "parameters"
             File_info = "info"
#ifdef TEMPERING
             write(File_para,'(A,I0,A)') "Temp_",igroup,"/parameters"
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif

             OPEN(UNIT=5,FILE=file_para,STATUS='old',ACTION='read',IOSTAT=ierr)
             IF (ierr /= 0) THEN
                WRITE(error_unit,*) 'Ham_set: unable to open <parameters>',ierr
                error stop 1
             END IF
             READ(5,NML=VAR_lattice)
             READ(5,NML=VAR_Z2_Slave)
             CLOSE(5)
#ifdef MPI
          Endif
          CALL MPI_BCAST(L1          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(L2          ,1  ,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(ham_T    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_chem ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_U    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Dtau     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Beta     ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Ham_xi   ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Ham_J    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Ham_h    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Ham_F    ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
          Call Ham_latt
          if ( Model .ne. "Z2_Slave" ) then
             Write(error_unit,*) "Ham_set: Model not yet implemented!"
             error stop 1
          endif
          If ( Lattice_type .ne. "Square" ) then
             Write(error_unit,*) "Ham_set: Z2_Slave is only implemented for a square lattice"
             error stop 1
          Endif

           Call Ham_hop
           Ltrot = nint(beta/dtau)
           Projector = .false.
           Theta = 0.d0
           Thtrot = 0

           Call Setup_Ising_action

#if defined(MPI)
           If (Irank_g == 0 ) then
#endif
              Open (Unit = 50,file=file_info,status="unknown",position="append")
              Write(50,*) '====================================='
              Write(50,*) 'Model is      : ', Model
              Write(50,*) 'Lattice is    : ', Lattice_type
              Write(50,*) '# of orbitals : ', Ndim
              Write(50,*) 'Beta          : ', Beta
              Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
              Write(50,*) 'N_SUN         : ', N_SUN
              Write(50,*) 'N_FL          : ', N_FL
              Write(50,*) 't             : ', Ham_T
              Write(50,*) 'Ham_U         : ', Ham_U
              Write(50,*) 'Ham_chem      : ', Ham_chem
              If ( Model == "Z2_Slave" ) then
                 Write(50,*) 'Ham_xi        : ', Ham_xi
                 Write(50,*) 'Ham_J         : ', Ham_J
                 Write(50,*) 'Ham_h         : ', Ham_h
                 Write(50,*) 'Ham_F         : ', Ham_F
              Endif
              close(50)
#if defined(MPI)
           endif
#endif
           call Ham_V

         end Subroutine Ham_Set

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the  Lattice
!--------------------------------------------------------------------
        Subroutine Ham_Latt

          Use Predefined_Lattices

          Implicit none
          ! Use predefined stuctures or set your own lattice.
          Call Predefined_Latt(Lattice_type, L1,L2,Ndim, List,Invlist,Latt,Latt_Unit)

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
                         Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                      ENDDO
                   else
                      DO I = 1, Latt%N
                         I1 = Latt%nnlist(I,1,0)
                         I2 = Latt%nnlist(I,0,1)
                         Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                      ENDDO
                   endif
                elseif ( Lattice_type=="Honeycomb" ) then
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
                         case (2)
                            J1 = invlist(Latt%nnlist(I,1,-1),2)
                         case (3)
                            J1 = invlist(Latt%nnlist(I,0,-1),2)
                         case default
                            Write(6,*) ' Error in  Ham_Hop '
                         end select
                         Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                         Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      Enddo
                   Enddo
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

          Integer :: nf, I, I1, I2,  nc, nc1,  J
          Real (Kind=Kind(0.d0)) :: X

          If  ( Model == "Z2_Slave" ) then
             Allocate(Op_V(N_coord*Ndim +1, N_FL))
             do nf = 1,N_FL
                do i  =  1, N_coord*Ndim
                   call Op_make(Op_V(i,nf),2)
                enddo
                call Op_make(Op_V(N_coord*Ndim + 1 ,nf),1 )
                ! This is for the reference spin on lattice site I = Latt%N
             enddo
             Do nc = 1,Ndim*N_coord   ! Runs over bonds.  Coordination number = 2.
                                      ! For the square lattice Ndim = Latt%N
                I1 = L_bond_inv(nc,1)
                if ( L_bond_inv(nc,2)  == 1 ) I2 = Latt%nnlist(I1,1,0)
                if ( L_bond_inv(nc,2)  == 2 ) I2 = Latt%nnlist(I1,0,1)
                Op_V(nc,1)%P(1) = I1
                Op_V(nc,1)%P(2) = I2
                Op_V(nc,1)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0))
                Op_V(nc,1)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
                Op_V(nc,1)%g      = cmplx(-dtau*Ham_xi,0.D0,kind(0.D0))
                Op_V(nc,1)%alpha  = cmplx(0d0  ,0.d0, kind(0.D0))
                Op_V(nc,1)%type   = 1
                Call Op_set( Op_V(nc,1) )
             Enddo
             I1 = Latt%N
             nc = Ndim*N_coord  +  1
             Op_V(nc,1)%P(1) = I1
             Op_V(nc,1)%O(1,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
             Op_V(nc,1)%g      = cmplx(0.d0,0.D0,kind(0.D0))
             Op_V(nc,1)%alpha  = cmplx(0d0,0.d0, kind(0.D0))
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
          !> Operator index
          Integer, Intent(IN) :: n
          !> Time slice
          Integer, Intent(IN) :: nt
          !> New local field on time slice nt and operator index n
          Real (Kind=Kind(0.d0)), Intent(In) :: Hs_new
          Integer :: nt1,I, F1,F2,I1,I2,I3

          !> Ratio for local spin-flip
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

             ! Magnetic flux term
             I1 = L_bond_inv(n,1)
             if ( L_bond_inv(n,2) == 1 ) then
                !     I2
                !     I1 I3
                I2 = Latt%nnlist(I1,0,1 )
                I3 = Latt%nnlist(I1,1,0 )
                F1 = nsigma%i(n,nt)*nsigma%i(L_bond(I1,2),nt)* nsigma%i(L_bond(I2,1),nt)*nsigma%i(L_bond(I3,2),nt)
                !     I1
                !     I2 I3
                I2 = Latt%nnlist(I1,0,-1)
                I3 = Latt%nnlist(I1,1,-1)
                F2 = nsigma%i(n,nt)*nsigma%i(L_bond(I2,1),nt)* nsigma%i(L_bond(I2,2),nt)*nsigma%i(L_bond(I3,2),nt)
             else
                !    I3
                !    I2  I1
                I2 = Latt%nnlist(I1,-1,0 )
                I3 = Latt%nnlist(I1,-1,1 )
                F1 = nsigma%i(n,nt)*nsigma%i(L_bond(I2,1),nt)* nsigma%i(L_bond(I2,2),nt)*nsigma%i(L_bond(I3,1),nt)
                !    I2
                !    I1  I3
                I2 = Latt%nnlist(I1,0,1)
                I3 = Latt%nnlist(I1,1,0)
                F2 = nsigma%i(n,nt)*nsigma%i(L_bond(I1,1),nt)* nsigma%i(L_bond(I2,1),nt)*nsigma%i(L_bond(I3,2),nt)
             endif
             S0 = S0*DW_Ising_Flux(F1,F2)
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
          Integer                   ::  ns , nc, n_op, ntau_p1, ntau_m1, I,n
          Integer, allocatable      ::  Isigma1(:),Isigma2(:),Isigma3(:)
          Real  (Kind = Kind(0.d0)) ::  S0_Matter, T0_Proposal

          ! Write(6,*) 'In GLob_move', m,direction,ntau, size(Flip_list,1), Size(Flip_value,1), Flip_list(1)
          ! Ising from n_op = 1,N_coord*Ndim
          ! Hubbard from n_op = N_coord*Ndim +1, Size(OP_V,1) = N_coord*Ndim +  Ndim
          ! Write(6,*) 'Global_move_tau ' , S0(Flip_list(1),ntau)

          Allocate (Isigma1(Latt%N), Isigma2(Latt%N), Isigma3(Latt%N) )

          I  =  nranf(Latt%N)
          Flip_length = 4
          do n = 1,4
             select case(n)
             case (1)
                n_op = L_bond(I,1)
             case (2)
                n_op = L_bond(I,2)
             case (3)
                n_op = L_bond(latt%nnlist(I,-1,0),1)
             case (4)
                n_op = L_bond(latt%nnlist(I,0,-1),2)
             case default
                Write(error_unit,*) 'Global_move_tau: Error'
                error stop 1
             end select
             Flip_list (n) = n_op
             Flip_value(n) = nsigma%flip(n_op,ntau)
          enddo
          If ( I == Latt%N )   then
             Flip_length   = 5
             n             = 5
             n_op          = N_coord*Ndim + 1
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
          S0_Matter = DW_Ising_tau ( Isigma1(I)*Isigma2(I) ) * DW_Ising_tau( Isigma2(I)*Isigma3(I) )
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

          !> Local
          Integer :: I,I1, nt,nt1, n, ns_old, n1,n2,n3,n4, dx, dy,dt,sx,sy,st
          Real (Kind=Kind(0.d0)) :: Weight, Ratio

          If (Model == "Z2_Slave" ) then

             I = nranf(Latt%N)
             n1  = L_bond(I,1)
             n2  = L_bond(I,2)
             n3  = L_bond(Latt%nnlist(I,-1,0),1)
             n4  = L_bond(Latt%nnlist(I,0,-1),2)
             do nt = 1,Ltrot
                nsigma%f(n1,nt) = nsigma%flip(n1,nt)
                nsigma%f(n2,nt) = nsigma%flip(n2,nt)
                nsigma%f(n3,nt) = nsigma%flip(n3,nt)
                nsigma%f(n4,nt) = nsigma%flip(n4,nt)
             enddo
          Endif

          Ratio = Delta_S0_global(Nsigma_old)
          Weight = 1.d0 - 1.d0/(1.d0+Ratio)
          If ( Weight < ranf_wrap() ) Then
             T0_Proposal_ratio = 0.d0
             nsigma%f  = nsigma_old%f
             nsigma%t  = nsigma_old%t
          else
             T0_Proposal_ratio = 1.d0/Ratio
          endif


        End Subroutine Global_move
!===================================================================================
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none

          !> Arguments
          type (Fields), Intent(IN)  :: nsigma_old
          !> Local

          Delta_S0_global = 1.d0

        end Function Delta_S0_global
!===================================================================================
        Subroutine Setup_Ising_action

          ! This subroutine sets up lists and arrays so as to enable an
          ! an efficient calculation of  S0(n,nt)

          Integer :: nc, nth, n, n1, n2, n3, n4, I, I1, n_orientation, Ix, Iy
          Real (Kind=Kind(0.d0)) :: X_p(2)

          ! Setup list of bonds for the square lattice.
          Allocate (L_Bond(Latt%N,2),  L_bond_inv(Latt%N*N_coord,2) )

          nc = 0
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
                   L_bond(I1,n_orientation) = nc
                   L_bond_inv(nc,1) = I1
                   L_bond_inv(nc,2) = n_orientation
                   ! The bond is given by  I1, I1 + a_(n_orientation).
                enddo
             endif
          Enddo
          !Test
          !Do I = 1,Latt%N
          !   Write(6,*)
          !   Write(6,*) Latt%list(I,1), Latt%list(I,2), I
          !   Write(6,*) L_bond(I,1), L_bond(I,2), L_bond( latt%nnlist(I,-1,0),1), L_bond( latt%nnlist(I,0,-1),2)
          !Enddo
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
          DO n = -1,1,2
             do n1 = -1,1,2
                DW_Ising_Flux(n,n1) = exp( Dtau*Ham_F*(dble(n) +  dble(n1) ))/ exp(  -Dtau*Ham_F*(dble(n) +  dble(n1) ))
             enddo
          enddo

        End Subroutine Setup_Ising_action
!===================================================================================
        Subroutine  Alloc_obs(Ltau)

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Nt
          Character (len=64) ::  Filename
          Character (len=:), allocatable ::  Channel

          ! Scalar observables
          Allocate ( Obs_scal(6) )
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
                N = 1;   Filename ="Flux"
             case (6)
                N = 1;   Filename ="Q"
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
                Filename = "Green"
             case (2)
                Filename = "SpinZ"
             case (3)
                Filename = "SpinXY"
             case (4)
                Filename = "Den"
             case (5)
                Filename = "Flux"
             case (6)
                Filename = "QQ"
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
             Nt = 1
             Channel = '--'
             Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
          enddo

          If (Ltau == 1) then
             ! Time displaced  correlators
             Allocate ( Obs_tau(6) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   Channel = 'P' ; Filename = "Green"
                case (2)
                   Channel = 'PH'; Filename = "SpinZ"
                case (3)
                   Channel = 'PH'; Filename = "SpinXY"
                case (4)
                   Channel = 'PH'; Filename = "Den"
                case (5)
                   Channel = 'PH'; Filename = "TauZ"
                case (6)
                   Channel = 'PH'
                   !Ns = Latt%N; No = 2   ;
                   Filename ="JJ" ! No = 2 for Jxx, Jyy, Jxy, Jyx
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = Ltrot+1-2*Thtrot
                If(Projector) Channel = 'T0'
                Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             enddo
          endif
        end Subroutine Alloc_obs

!========================================================================
        Subroutine Obser(GR,Phase,Ntau)
#if defined(Machine_Learning)
          Use mpi
#endif
          Implicit none


          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          !Local
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS, Z1, Z2
          Complex (Kind=Kind(0.d0)) :: ZQ, ZSTAR, ZQT, ZQTT
          Integer :: I,J, imj, nf, dec, I1, I2,I3,I4, J1, no_I, no_J,  iFlux_tot,  &
               &     no, ntau1, L_Vison, L_Wilson, n, nx,ny
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
          Zkin = Zkin * dble(N_SUN)* dble(Ham_Xi)
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
             iFlux_tot = iFlux_tot + iFlux(I,Ntau)
          Enddo
          Obs_scal(5)%Obs_vec(1)  =   Obs_scal(5)%Obs_vec(1) + cmplx(dble(iFlux_tot),0.d0,kind(0.d0))*ZP*ZS

          ZQ = cmplx(0.d0,0.d0,kind(0.d0))
          DO I = 1,Latt%N
             Z    =  ( cmplx(1.d0,0.d0,kind(0.d0)) - cmplx(2.d0,0.d0,kind(0.d0))*GRC(I,I,1) )**(N_SUN)
             ZQ   = ZQ    + cmplx(DW_Ising_tau( Isigma(I)*Isigma1(I)),0.d0,kind(0.d0) ) * Z
          Enddo
          Obs_scal(6)%Obs_vec(1) = Obs_scal(6)%Obs_vec(1)  + ZQ * ZP*ZS


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
                ! Green
                Z1 = cmplx(real(Isigma(I1)*Isigma(J1), kind(0.d0)), 0.d0,kind(0.d0))

                Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                     &               Z * Z1*GRC(I1,J1,1) *  ZP*ZS
                ! SpinZ
                Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
                     &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                ! SpinXY
                Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
                     &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
                ! Den
                Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J)  +  &
                     &     (    GRC(I1,I1,1) * GRC(J1,J1,1) *Z     + &
                     &          GRC(I1,J1,1) * GR(I1,J1,1 )           &
                     &                                   ) * Z* ZP*ZS
                !Flux
                Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) + &
                     &         cmplx(real(iFlux(I1,Ntau)*iFlux(J1,Ntau),kind(0.d0)),0.d0,kind(0.d0))*ZP*ZS

                Z1 =   (cmplx(1.d0,0.d0,kind(0.d0)) - cmplx(2.d0,0.d0,kind(0.d0))*GRC(I1,I1,1)) *  &
                     & (cmplx(1.d0,0.d0,kind(0.d0)) - cmplx(2.d0,0.d0,kind(0.d0))*GRC(J1,J1,1)) +  &
                     &  cmplx(4.d0,0.d0,kind(0.d0)) * GRC(I1,J1,1)*GR(I1,J1,1)
                Z1 = Z1**(N_SUN)
                ZQ = cmplx(DW_Ising_tau( Isigma(I1)*Isigma1(I1))*DW_Ising_tau( Isigma(J1)*Isigma1(J1)),0.d0,kind(0.d0) )*Z1
                If (I1 == J1 )  ZQ = cmplx(1.d0,0.d0,kind(0.d0))
                Obs_eq(6)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(6)%Obs_Latt(imj,1,no_I,no_J) + ZQ*ZP*ZS

             ENDDO
             Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  Z * GRC(I1,I1,1) * ZP * ZS
             Obs_eq(5)%Obs_Latt0(no_I) =  Obs_eq(5)%Obs_Latt0(no_I) +  &
                  &   cmplx(real(iFlux(I1,Ntau),kind(0.d0)),0.d0,kind(0.d0)) * ZP * ZS
          ENDDO

          Deallocate ( Isigma, Isigma1 )

        end Subroutine Obser
!=====================================================
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE)
          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase

          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS, Z1
          Integer :: IMJ, I, J, I1, J1, no_I, no_J, NT1, no
          Integer,  allocatable  :: Isigma(:), IsigmaT(:)
          Complex (Kind=Kind(0.d0)), allocatable ::  J_tmp0(:,:), J_tmpT(:,:)

          Allocate ( Isigma(Latt%N), IsigmaT(Latt%N), J_tmp0(Latt%N,2), J_tmpT(Latt%N,2) )
          NT1 = NT
          If (NT == 0 ) NT1 = LTROT
          Call Hamiltonian_set_Z2_matter(Isigma ,Ltrot )
          Call Hamiltonian_set_Z2_matter(IsigmaT,NT1)

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          If (NT == 0 ) then
             DO I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             ENDDO
          endif
          If ( Model == "Z2_Slave"  ) then
             Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   !  GT0 =  < T c_I1(nt) c^dag_J1(0) >
                   Z1 = cmplx(real( IsigmaT(I1)*Isigma(J1), kind(0.d0) ), 0.d0,kind(0.d0))
                   Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & +  Z * GT0(I1,J1,1) * Z1* ZP* ZS

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
                   ! Iz Iz correlations
                   Obs_tau(5)%Obs_Latt(imj,nt+1,no_I,no_J) = Obs_tau(5)%Obs_Latt(imj,nt+1,no_I,no_J) + &
                        &  cmplx( real( IsigmaT(I1)*Isigma(J1), kind(0.d0) ) , 0.d0,kind(0.d0)) * ZP * ZS
                Enddo
                Obs_tau(4)%Obs_Latt0(no_I) = Obs_tau(4)%Obs_Latt0(no_I) + &
                     &         Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1)) * ZP * ZS
             Enddo
             ! Current-Current correlations
             Do I = 1,Latt%N
                do no = 1,2
                   if (no == 1)  I1 = Latt%nnlist(I,1,0)
                   if (no == 2)  I1 = Latt%nnlist(I,0,1)
                   J_tmp0(I,no)  = cmplx(0.d0,-1.d0,Kind(0.d0)) * Z * ( G00(I1,I,1) - G00(I,I1,1) )
                   J_tmpT(I,no)  = cmplx(0.d0,-1.d0,Kind(0.d0)) * Z * ( GTT(I1,I,1) - GTT(I,I1,1) )
                enddo
             Enddo
             do no_I = 1,2
                Do no_J = 1,2
                   DO I = 1,Latt%N
                      if (no_I == 1 ) I1 = Latt%nnlist(I,1,0)
                      if (no_I == 2 ) I1 = Latt%nnlist(I,0,1)
                      Do J = 1,Latt%N
                         if (no_J == 1 ) J1 = Latt%nnlist(J,1,0)
                         if (no_J == 2 ) J1 = Latt%nnlist(J,0,1)
                         imj = latt%imj(I,J)
                         Z1 = cmplx(real( IsigmaT(I)*IsigmaT(I1) * Isigma(J)*Isigma(J1), kind(0.d0) ), 0.d0,kind(0.d0))
                         Obs_tau(6)%Obs_Latt(imj,nt+1,no_I,no_J) = Obs_tau(6)%Obs_Latt(imj,nt+1,no_I,no_J)  +                     &
                              &                              (                                                                    &
                              &                                    J_tmpT(I,no_I)* J_tmp0(J,no_J)           +                     &
                              &                                 Z*( + G0T(J1,I,1) * GT0(I1,J ,1)  + G0T(J ,I1,1) * GT0(I,J1,1)    &
                              &                                     - G0T(J ,I,1) * GT0(I1,J1,1)  - G0T(J1,I1,1) * GT0(I,J, 1) )  &
                              &                              ) * Z1 * ZP * ZS
                      Enddo
                   Enddo
                Enddo
             Enddo
          Endif

          Deallocate ( Isigma, IsigmaT, J_tmp0, J_tmpT )

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
             Call  Print_bin_Latt(Obs_eq(I),Group_Comm)
          enddo
          If (Ltau  == 1 ) then
             Do I = 1,Size(Obs_tau,1)
                Call  Print_bin_Latt(Obs_tau(I),Group_Comm)
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

      Integer Function  iFlux(I,nt)

        Implicit none

        Integer, INTENT(IN) :: I,nt

        ! Local
        Integer :: n1,n2,n3,n4

        !   I3  I2
        !   I   I1
        n1  = L_bond(I,1)
        n2  = L_bond(Latt%nnlist(I,1,0),2)
        n3  = L_bond(Latt%nnlist(I,0,1),1)
        n4  = L_bond(I,2)
        iFlux =   nsigma%i(n1,nt)*nsigma%i(n2,nt)*nsigma%i(n3,nt)*nsigma%i(n4,nt)

      end Function iFlux

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

        ! The user can set the initial configuration

        Implicit none

        Real (Kind=Kind(0.d0)), allocatable, dimension(:,:), Intent(OUT) :: Initial_field


        Integer :: I,nc, I1, nt
        Integer, allocatable::  Isigma(:), Isigma1(:)
        Allocate  (Initial_field(N_coord*Latt%N +1, Ltrot) )

        allocate (Isigma(Latt%N), Isigma1(Latt%N) )
        Do nt = 1,Ltrot
           Do I = 1,Latt%N
              Isigma(I) = 1
              if ( ranf_wrap()  > 0.5D0 ) Isigma(I)  = -1
           enddo
           Initial_field(N_coord*Latt%N +1,nt) = real(Isigma(Latt%N), kind(0.d0))
           do nc = 1,N_coord*Latt%N
              I = L_bond_inv(nc,1)
              if (  L_bond_inv(nc,2) == 1 )  I1 = latt%nnlist(I,1,0)
              if (  L_bond_inv(nc,2) == 2 )  I1 = latt%nnlist(I,0,1)
              Initial_field(nc,nt) = real(Isigma(I)*Isigma(I1), kind(0.d0))
           enddo
           nsigma%f = Initial_field
           Call Hamiltonian_set_Z2_matter(Isigma1,nt)
           Do nc = 1,Latt%N
              if ( Isigma(nc) .ne.  Isigma1(nc)  ) then
                 Write(error_unit,*) 'Error in Hamiltonian_set_Z2_matter.Hamiltonian_set_nsigma'
                 error stop 1
              endif
           enddo
           !Write(6,*)
           !Do nc = 1,Latt%N
           !   Write(6,*) Isigma(nc), Isigma1(nc)
           !enddo
        enddo
        deallocate (Isigma, Isigma1)

        !Call Print_fluxes

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

        Isigma(Latt%N) = nsigma%i( N_coord*Ndim + 1, nt )
        I = Latt%N
        do nx = 1,L1
           do ny = 1,L2
              I1 = latt%nnlist(I,0,1)
              Isigma(I1)  = Isigma(I)*nsigma%i(L_bond(I,2),nt)
              !Write(6,*) Latt%list(I,1), Latt%list(I,2), ' -> ', Latt%list(I1,1), Latt%list(I1,2)
              I = I1
           enddo
           I1          = latt%nnlist(I,1,0)
           Isigma(I1)  = Isigma(I)*nsigma%i(L_bond(I,1),nt)
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
        Use mpi
#endif
        Implicit none

        ! Local
        Integer :: I,nt,ix, iy, n
        Character (len=64) :: File1

#ifdef MPI
        Integer        :: Isize, Irank, IERR
        Integer        :: STATUS(MPI_STATUS_SIZE)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif

#if defined(TEMPERING)
        write(File1,'(A,I0,A)') "Temp_",Irank,"/Fluxes"
#else
        File1="Fluxes"
#endif
        Open (Unit=10,File=File1, status="unknown")
        Do nt = 1,Ltrot
           Do i  = 1,Ndim
              n = iFlux(I,nt)
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
      end Subroutine Overide_global_tau_sampling_parameters




!--------------------------------------------------------------------
!       Subroutine Test_Hamiltonian
!
!         Implicit none
!
!         Integer :: n,  nc, n_op, nt
!         Integer, allocatable :: nsigma_old(:,:)
!         Real (Kind=kind(0.d0)) :: X, X1, size_clust
!
!         n = size(Op_V,1)
!         allocate (nsigma_old(n,Ltrot))
!         do nc = 1,100
!            !nt  = nranf(Ltrot)
!            !n_op= nranf(n)
!            !if ( OP_V(n_op,1)%type == 1 ) then
!            !   X = S0(n_op,nt)
!            !   nsigma_old = nsigma
!            !   nsigma(n_op,nt) = -nsigma(n_op,nt)
!            !   X1 = Delta_S0_global(Nsigma_old)
!            !   Write(6,*) nc, X, X1
!            !endif
!            nsigma_old = nsigma
!            Call Global_move(X,nsigma_old,size_clust)
!            X1 = Delta_S0_global(Nsigma_old)
!            Write(6,*) nc, X, X1
!         enddo
!         deallocate (nsigma_old)
!
!
!       end Subroutine Test_Hamiltonian

      end Module Hamiltonian
