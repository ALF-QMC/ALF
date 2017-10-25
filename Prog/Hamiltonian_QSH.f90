    Module Hamiltonian

      Use Operator_mod
      Use Lattices_v3 
      Use MyMats 
      Use Random_Wrap
      Use Files_mod
      Use Matrix
      Use Observables      

      Type (Operator), dimension(:,:), allocatable  :: Op_V
      Type (Operator), dimension(:,:), allocatable  :: Op_T
      Integer, allocatable :: nsigma(:,:)
      Integer              :: Ndim,  N_FL,  N_SUN,  Ltrot
!>    Defines MPI communicator 
      Integer              :: Group_Comm


      ! What is below is  private 
      Type (Lattice),       private :: Latt
      Integer,              private :: L1, L2, n_Lambda
      real (Kind=8),        private :: ham_T , ham_U,  Ham_chem, Lambda, Ham_lamb
      real (Kind=8),        private :: Dtau, Beta, Phi_x
      Character (len=64),   private :: Model, Lattice_type
      Logical,              private :: One_dimensional
      Integer,              private :: N_coord, Norb 
      Integer, allocatable, private :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell


!>    Privat Observables
      Type (Obser_Vec ),  private, dimension(:), allocatable ::   Obs_scal
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_eq
      Type (Obser_Latt),  private, dimension(:), allocatable ::   Obs_tau
      
    contains 


      Subroutine Ham_Set

          Implicit none

#ifdef MPI
          include 'mpif.h'
#endif   

          integer :: ierr

          
          NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model, Phi_x

          NAMELIST /VAR_Hubbard/  ham_T, ham_chem, ham_U, Dtau, Beta, Lambda, n_Lambda

          NAMELIST /VAR_QSH/ ham_T, ham_chem, Ham_lamb, Dtau, Beta

#ifdef MPI
          Integer        :: Isize, Irank
          Integer        :: STATUS(MPI_STATUS_SIZE)
          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
          CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
          
          


#ifdef MPI
          If (Irank == 0 ) then
#endif
             Phi_x = 0.d0
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
          CALL MPI_BCAST(Phi_x       ,1  ,MPI_REAL8  ,   0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Model       ,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
          CALL MPI_BCAST(Lattice_type,64 ,MPI_CHARACTER, 0,MPI_COMM_WORLD,IERR)
#endif
          Call Ham_latt

          If ( Model == "Hubbard_Mz") then
             N_FL = 2
             N_SUN = 1
          elseif  ( Model == "Hubbard_SU2" ) then
             N_FL = 1
             N_SUN = 2
          elseif  ( (Model == "QSH_SU2") .or. (Model == "QSH_Z2") &
          &   .or. (Model == "Hubbard_SSC") .or. (Model == "QSH_U1") ) then
             N_FL = 1
             N_SUN = 1

          else
             Write(6,*) "Model not yet implemented!"
             Stop
          endif
#ifdef MPI
          If (Irank == 0 ) then
#endif
             OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
             If (( Model == "QSH_SU2").or. (Model == "QSH_Z2") .or. (Model == "QSH_U1")) then
             READ(5,NML=VAR_QSH)
             else
             READ(5,NML=VAR_Hubbard)
             endif
             CLOSE(5)
#ifdef MPI
          endif

          If ((Model == "QSH_SU2").or.(Model == "QSH_Z2").or.(Model == "QSH_U1")) then
          CALL MPI_BCAST(ham_T    ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_chem ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_lamb ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Dtau     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Beta     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          else
          CALL MPI_BCAST(ham_T    ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_chem ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ham_U    ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Dtau     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Beta     ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(Lambda   ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(n_Lambda ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          endif

#endif
          Call Ham_hop

          Ltrot = nint(beta/dtau)
#ifdef MPI
          If (Irank == 0) then
#endif
             Open (Unit = 50,file="info",status="unknown",position="append")

             if ((Model == "QSH_SU2").or.(Model == "QSH_Z2").or.(Model == "QSH_U1")) then
             Write(50,*) '=============QSH========='
             Write(50,*) 'Model is      : ', Model
             Write(50,*) 'L1,L2         : ', L1,L2
             Write(50,*) 'Beta          : ', Beta
             Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
             Write(50,*) 'N_SUN         : ', N_SUN
             Write(50,*) 'N_FL          : ', N_FL
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'Ham_chem      : ', Ham_chem
             Write(50,*) 'Ham_lamb      : ', Ham_lamb
             else
             Write(50,*) '=============Canonical========='
             Write(50,*) 'Model is      : ', Model 
             Write(50,*) 'L1,L2         : ', L1,L2
             Write(50,*) 'Phi_x         : ', Phi_x
             Write(50,*) 'Beta          : ', Beta
             Write(50,*) 'dtau,Ltrot    : ', dtau,Ltrot
             Write(50,*) 'N_SUN         : ', N_SUN
             Write(50,*) 'N_FL          : ', N_FL
             Write(50,*) 't             : ', Ham_T
             Write(50,*) 'Ham_U         : ', Ham_U
             Write(50,*) 'Ham_chem      : ', Ham_chem
             Write(50,*) 'Lambda        : ', Lambda
             Write(50,*) 'n_Lambda      : ', n_Lambda
             endif
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
          Real (Kind=8)  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
          Integer :: I, nc, no


          If ( Lattice_type =="Square" ) then
             a1_p(1) =  1.0  ; a1_p(2) =  0.d0
             a2_p(1) =  0.0  ; a2_p(2) =  1.d0
             L1_p    =  dble(L1)*a1_p
             L2_p    =  dble(L2)*a2_p
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
             !Write(6,*)  'Lattice: ', Ndim
             One_dimensional = .false.
             N_coord   = 2
             Norb = 1
             If ( L1 == 1 .or. L2 == 1 ) then 
                One_dimensional = .true.
                N_coord   = 1
                If (L1 == 1 ) then 
                   Write(6,*) ' For one dimensional systems set  L2 = 1 ' 
                   Stop
                endif
             endif

             elseif ( Lattice_type=="Honeycomb" ) then
             if ((Model.eq."QSH_SU2").or.(Model.eq."QSH_Z2") .or.  &
             & ( Model.eq."QSH_U1" ) .or.(Model.eq."Hubbard_SSC")) then
              Norb    = 4
             else 
              Norb    = 2
             endif
              N_coord = 3
              a1_p(1) =  1.D0   ; a1_p(2) =  0.d0
              a2_p(1) =  0.5D0  ; a2_p(2) =  sqrt(3.D0)/2.D0
 
              !del_p   =  (a2_p - 0.5*a1_p ) * 2.0/3.0
              L1_p    =  dble(L1) * a1_p
              L2_p    =  dble(L2) * a2_p
              Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )

          else
             Write(6,*) "Lattice not yet implemented!"
             Stop
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

          Integer :: I, I1, I2, n, Ncheck, nc, no, nc1, J1
          Complex  (Kind=8) :: Z
          Real     (Kind=8) :: Pi

          Ncheck = 1
          Pi = acos(-1.d0)
          Z = exp( cmplx(0.d0,Phi_x*2.d0*pi/real(L1,kind(0.d0)),kind(0.d0)))
          allocate(Op_T(Ncheck,N_FL))
          do n = 1,N_FL
             Do nc = 1,Ncheck
                Call Op_make(Op_T(nc,n),Ndim)
                If (One_dimensional ) then 
                   DO I = 1, Latt%N
                      I1 = Latt%nnlist(I, 1, 0)
                      Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))*Z
                      Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))*Conjg(Z)
                      Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                   ENDDO
                elseif ( Lattice_type=="Squre" ) then
                   DO I = 1, Latt%N
                      I1 = Latt%nnlist(I,1,0)
                      I2 = Latt%nnlist(I,0,1)
                      Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))*Z
                      Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))*Conjg(Z)
                      Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                   ENDDO

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

                        if ((Model.eq."QSH_SU2").or.(Model.eq."QSH_Z2").or. &
                      & (Model.eq."QSH_U1") .or. (Model.eq."Hubbard_SSC")) then
                         I1 = Invlist(I,3)
                         Do nc1 = 1,N_coord
                            select case (nc1)
                            case (1)
                               J1 = invlist(I,4)
                            case (2)
                               J1 = invlist(Latt%nnlist(I,1,-1),4)
                            case (3)
                               J1 = invlist(Latt%nnlist(I,0,-1),4)
                            case default
                               Write(6,*) ' Error in  Ham_Hop '
                            end select
                            Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,  0.d0,  kind(0.D0))
                            Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,  0.d0,  kind(0.D0))
                          Enddo
                         endif
                     Enddo
                endif

                
                Do I = 1, Ndim
                   Op_T(nc,n)%P(i) = i 
                Enddo
                if ( abs(Ham_T) < 1.E-6 .and.  abs(Ham_chem) < 1.E-6 ) then 
                   Op_T(nc,n)%g = 0.d0
                else
                   Op_T(nc,n)%g = -Dtau
                endif
                Op_T(nc,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
                !Write(6,*) 'In Ham_hop', Ham_T
                Call Op_set(Op_T(nc,n)) 
                !Write(6,*) 'In Ham_hop 1'
                !Do I = 1,Latt%N
                !   Write(6,*) Op_T(n)%E(i)
                !enddo
                !Call Op_exp( cmplx(-Dtau,0.d0), Op_T(n), Exp_T   (:,:,n) )
                !Call Op_exp( cmplx( Dtau,0.d0), Op_T(n), Exp_T_M1(:,:,n) )
             enddo
          enddo
        end Subroutine Ham_hop

!===================================================================================           
        
        Subroutine Ham_V
          
          Implicit none 
          
          Integer :: nf, I, nc, nth, nz, nb, J1, J2, I1, I2, I3, I4, no
          Integer :: I1u,I1d,I2u,I2d,I3u,I3d,J1u,J1d,J2u,J2d,J3u,J3d,Bond(6,4)
          Integer :: I_up, I_do, J_up, J_do
          Real (Kind=8) :: X,Xp(6)

          
          If (Model == "Hubbard_SU2")  then
             !Write(50,*) 'Model is ', Model
!            Allocate(Op_V(Latt%N+n_Lambda,N_FL))
             Allocate(Op_V(Ndim+n_Lambda,N_FL))
             do nf = 1,N_FL
!                do i  = 1, Latt%N
                do i  = 1, Ndim
                   Call Op_make(Op_V(i,nf),1)
                enddo
                do nth = 1,n_lambda
                   Call Op_make(Op_V(Latt%N+nth,nf),Latt%N)
                enddo
             enddo
             Do nf = 1,N_FL
                nc = 0
!                Do i = 1,Latt%N
                Do i = 1, Ndim
                   nc = nc + 1
                   Op_V(nc,nf)%P(1) = I
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*ham_U/(DBLE(N_SUN)), 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(-0.5d0,0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                   ! The operator reads:  
                   !  g*s*( c^{dagger} O c  + alpha ))
                   ! with s the HS field.
                Enddo
                Do nth = 1,n_Lambda
                   nc = nc + 1
                   Do  I = 1,latt%N
                      Op_V(nc,nf)%P(I) = I
                      Op_V(nc,nf)%O(I,I) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   enddo
                   Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*Lambda/dble(n_lambda), 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(-DBLE(Latt%N)*0.5d0,0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                Enddo
             Enddo
          Elseif (Model == "QSH0_SU2") then
          print*, 'N_FL=', N_FL
             Allocate(Op_V(Latt%N*18,N_FL))
             do nf = 1,N_FL
               do i = 1, Latt%N*18    
                  Call Op_make(Op_V(i,nf),4)  
               enddo
             enddo
             nc = 0
             nf = 1
             Do I  = 1,Latt%N
              do nb = 1, 3    
              do nz = 1, 3 
              do no = 1, 2
              
               nc = nc + 1
               I1 = invlist(I,1+no/2)
               Op_V(nc,nf)%P(1) = I1
               I2 = invlist(I,3+no/2) 
               Op_V(nc,nf)%P(3) = I2               

               select case (nb)
               case (1)
               J1 = invlist(Latt%nnlist(I,1,0), 1+no/2)
               Op_V(nc,nf)%P(2) = J1              
               J2 = invlist(Latt%nnlist(I,1,0), 3+no/2)
               Op_V(nc,nf)%P(4) = J2
               X = 1.d0 * (-1) ** no
               case (2)
               J1 = invlist(Latt%nnlist(I,0,1), 1+no/2)
               Op_V(nc,nf)%P(2) = J1              
               J2 = invlist(Latt%nnlist(I,0,1), 3+no/2)
               Op_V(nc,nf)%P(4) = J2         
               X = -1.d0 * (-1) ** no
               case (3)
               J1 = invlist(Latt%nnlist(I,-1,1),1+no/2)
               Op_V(nc,nf)%P(2) = J1    
               J2 = invlist(Latt%nnlist(I,-1,1),3+no/2)
               Op_V(nc,nf)%P(4) = J2
               X = 1.d0 * (-1) ** no
               end select
                
               select case (nz)
               case (1)
               Op_V(nc,nf)%O(1,2) =  cmplx( 0, X, kind(0.D0))
               Op_V(nc,nf)%O(2,1) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(3,4) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(4,3) =  cmplx( 0, X, kind(0.D0))

               case (2)
               Op_V(nc,nf)%O(1,4) =  cmplx( 0, X, kind(0.D0))
               Op_V(nc,nf)%O(4,1) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(2,3) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(3,2) =  cmplx( 0, X, kind(0.D0))

               case (3)
               Op_V(nc,nf)%O(1,4) =  cmplx(-X, 0, kind(0.D0))
               Op_V(nc,nf)%O(4,1) =  cmplx(-X, 0, kind(0.D0))
               Op_V(nc,nf)%O(2,3) =  cmplx( X, 0, kind(0.D0))
               Op_V(nc,nf)%O(3,2) =  cmplx( X, 0, kind(0.D0))
               end select
 

               Op_V(nc,nf)%g = SQRT(CMPLX(DTAU*ham_lamb/dble(N_SUN), 0.D0, kind(0.D0)))              
               Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
               Op_V(nc,nf)%type   = 2
               Call Op_set( Op_V(nc,nf) )

              enddo
              enddo
              enddo
             Enddo

             Elseif (Model == "QSH0_U1") then
             Allocate(Op_V(Latt%N*12,N_FL))
             do nf = 1,N_FL
               do i = 1, Latt%N*12    
                  Call Op_make(Op_V(i,nf),4)  
               enddo
             enddo
             nc = 0
             nf = 1
             Do I  = 1,Latt%N
              do nb = 1, 3    
              do nz = 1, 2 
              do no = 1, 2
              
               nc = nc + 1
               I1 = invlist(I,1+no/2)
               Op_V(nc,nf)%P(1) = I1
               I2 = invlist(I,3+no/2) 
               Op_V(nc,nf)%P(3) = I2               

               select case (nb)
               case (1)
               J1 = invlist(Latt%nnlist(I,1,0), 1+no/2)
               Op_V(nc,nf)%P(2) = J1              
               J2 = invlist(Latt%nnlist(I,1,0), 3+no/2)
               Op_V(nc,nf)%P(4) = J2
               X = 1.d0 * (-1) ** no
               case (2)
               J1 = invlist(Latt%nnlist(I,0,1), 1+no/2)
               Op_V(nc,nf)%P(2) = J1              
               J2 = invlist(Latt%nnlist(I,0,1), 3+no/2)
               Op_V(nc,nf)%P(4) = J2         
               X = -1.d0 * (-1) ** no
               case (3)
               J1 = invlist(Latt%nnlist(I,-1,1),1+no/2)
               Op_V(nc,nf)%P(2) = J1    
               J2 = invlist(Latt%nnlist(I,-1,1),3+no/2)
               Op_V(nc,nf)%P(4) = J2
               X = 1.d0 * (-1) ** no
               end select
                
               select case (nz)
               case (1)
               Op_V(nc,nf)%O(1,2) =  cmplx( 0, X, kind(0.D0))
               Op_V(nc,nf)%O(2,1) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(3,4) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(4,3) =  cmplx( 0, X, kind(0.D0))

               case (2)
               Op_V(nc,nf)%O(1,4) =  cmplx( 0, X, kind(0.D0))
               Op_V(nc,nf)%O(4,1) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(2,3) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(3,2) =  cmplx( 0, X, kind(0.D0))

               end select
 

               Op_V(nc,nf)%g = SQRT(CMPLX(DTAU*ham_lamb/dble(N_SUN), 0.D0, kind(0.D0)))              
               Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
               Op_V(nc,nf)%type   = 2
               Call Op_set( Op_V(nc,nf) )

              enddo
              enddo
              enddo
             Enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
          Elseif (Model == "QSH_Z2") then

             Allocate(Op_V(Latt%N,N_FL))
             do nf = 1,N_FL
               do i = 1, Latt%N
                  Call Op_make(Op_V(i,nf),12)
               enddo
             enddo
             nc = 0
             nf = 1


             Do I  = 1,Latt%N

              nc = nc + 1
              do nb = 1, 3
              do no = 1, 2

               select case (nb)
               case (1)
               J1 = invlist(Latt%nnlist(I,1,0), 1+no/2)
               Op_V(nc,nf)%P(1+no/2) = J1
               J2 = invlist(Latt%nnlist(I,1,0), 3+no/2)
               Op_V(nc,nf)%P(3+no/2) = J2
               case (2)
               select case (no)
               case (1)
               J1 = invlist(Latt%nnlist(I,0,1), 1+no/2)
               Op_V(nc,nf)%P(5+no/2) = J1
               J2 = invlist(Latt%nnlist(I,0,1), 3+no/2)
               Op_V(nc,nf)%P(7+no/2) = J2
               case (2)
               J1 = invlist(Latt%nnlist(I,1,-1), 1+no/2)
               Op_V(nc,nf)%P(5+no/2) = J1
               J2 = invlist(Latt%nnlist(I,1,-1), 3+no/2)
               Op_V(nc,nf)%P(7+no/2) = J2
               end select
               case (3)
               J1 = invlist(I,1+no/2)
               Op_V(nc,nf)%P(9+no/2) = J1
               J2 = invlist(I,3+no/2)
               Op_V(nc,nf)%P(11+no/2)= J2
               end select

               enddo
               enddo
               do nb = 1, 3
               do no = 1, 2
               X = (-1.d0) ** no

               I_up = 1 + (nb-1) * 4 + no/2
               I_do = 3 + (nb-1) * 4 + no/2
               if (nb.lt.3) then
               J_up = 1 + nb * 4 + no/2
               J_do = 3 + nb * 4 + no/2
               else
               J_up = 1 + no/2
               J_do = 3 + no/2
               endif

               Op_V(nc,nf)%O(I_up, J_up) =  cmplx( 0, X, kind(0.D0))
               Op_V(nc,nf)%O(J_up, I_up) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(I_do, J_do) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(J_do, I_do) =  cmplx( 0, X, kind(0.D0))
               enddo
               enddo

               Op_V(nc,nf)%g = SQRT(CMPLX( DTAU*ham_lamb, 0.D0, kind(0.D0)))
               Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
               Op_V(nc,nf)%type   = 2
               Call Op_set( Op_V(nc,nf) )

             Enddo    
!------------------------------------------------------------------------------------1
          Elseif (Model == "QSH_U1") then
          print*, 'N_FL=', N_FL
             Allocate(Op_V(Latt%N*2,N_FL))
             do nf = 1,N_FL
               do i = 1, Latt%N*2
                  Call Op_make(Op_V(i,nf), 12)
               enddo
             enddo
             print*, 'finish'
             nc = 0
             nf = 1
             Do I  = 1,Latt%N

              do nz = 1, 2
              nc = nc + 1
              do nb = 1, 3
              do no = 1, 2

               print*, nz, nb, no
               select case (nb)
               case (1)
               J1 = invlist(Latt%nnlist(I,1,0), 1+no/2)
               Op_V(nc,nf)%P(1+no/2) = J1
               J2 = invlist(Latt%nnlist(I,1,0), 3+no/2)
               Op_V(nc,nf)%P(3+no/2) = J2
               case (2)
               select case (no)
               case (1)
               J1 = invlist(Latt%nnlist(I,0,1), 1+no/2)
               Op_V(nc,nf)%P(5+no/2) = J1
               J2 = invlist(Latt%nnlist(I,0,1), 3+no/2)
               Op_V(nc,nf)%P(7+no/2) = J2
               case (2)
               J1 = invlist(Latt%nnlist(I,1,-1), 1+no/2)
               Op_V(nc,nf)%P(5+no/2) = J1
               J2 = invlist(Latt%nnlist(I,1,-1), 3+no/2)
               Op_V(nc,nf)%P(7+no/2) = J2
               end select
               case (3)
               J1 = invlist(I,1+no/2)
               Op_V(nc,nf)%P(9+no/2) = J1
               J2 = invlist(I,3+no/2)
               Op_V(nc,nf)%P(11+no/2)= J2
               end select
               enddo
               enddo

               do nb = 1, 3
               do no = 1, 2
               X = (-1.d0) ** no

               I_up = 1 + (nb-1) * 4 + no/2
               I_do = 3 + (nb-1) * 4 + no/2
               if (nb.lt.3) then
               J_up = 1 + nb * 4 + no/2
               J_do = 3 + nb * 4 + no/2
               else
               J_up = 1 + no/2
               J_do = 3 + no/2
               endif

               select case (nz)
               case (1)
               Op_V(nc,nf)%O(I_up, J_up) =  cmplx( 0, X, kind(0.D0))
               Op_V(nc,nf)%O(J_up, I_up) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(I_do, J_do) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(J_do, I_do) =  cmplx( 0, X, kind(0.D0))
               case (2)
               Op_V(nc,nf)%O(I_up, J_do) =  cmplx( 0, X, kind(0.D0))
               Op_V(nc,nf)%O(J_do, I_up) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(J_up, I_do) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(I_do, J_up) =  cmplx( 0, X, kind(0.D0))
               end select
               enddo
               enddo

               Op_V(nc,nf)%g = SQRT(CMPLX(DTAU*ham_lamb/dble(N_SUN), 0.D0, kind(0.D0)))
               Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
               Op_V(nc,nf)%type   = 2
               Call Op_set( Op_V(nc,nf) )

              enddo
             Enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
          Elseif (Model == "QSH_SU2") then
          print*, 'N_FL=', N_FL
             Allocate(Op_V(Latt%N*3,N_FL))
             print*, 'allocate'
             do nf = 1,N_FL
               do i = 1, Latt%N*3
                  Call Op_make(Op_V(i,nf), 12)
               enddo
             enddo
             print*, 'finish'
             nc = 0
             nf = 1
             Do I  = 1,Latt%N

              do nz = 1, 3
              nc = nc + 1
              do nb = 1, 3
              do no = 1, 2

               print*, nz, nb, no
               select case (nb)
               case (1)
               J1 = invlist(Latt%nnlist(I,1,0), 1+no/2)
               Op_V(nc,nf)%P(1+no/2) = J1
               J2 = invlist(Latt%nnlist(I,1,0), 3+no/2)
               Op_V(nc,nf)%P(3+no/2) = J2
               case (2)
               select case (no)
               case (1)
               J1 = invlist(Latt%nnlist(I,0,1), 1+no/2)
               Op_V(nc,nf)%P(5+no/2) = J1
               J2 = invlist(Latt%nnlist(I,0,1), 3+no/2)
               Op_V(nc,nf)%P(7+no/2) = J2
               case (2)
               J1 = invlist(Latt%nnlist(I,1,-1), 1+no/2)
               Op_V(nc,nf)%P(5+no/2) = J1
               J2 = invlist(Latt%nnlist(I,1,-1), 3+no/2)
               Op_V(nc,nf)%P(7+no/2) = J2
               end select     
               case (3)
               J1 = invlist(I,1+no/2)
               Op_V(nc,nf)%P(9+no/2) = J1
               J2 = invlist(I,3+no/2)
               Op_V(nc,nf)%P(11+no/2)= J2
               end select

               enddo
               enddo

               do nb = 1, 3
               do no = 1, 2
               X = (-1.d0) ** no

               I_up = 1 + (nb-1) * 4 + no/2
               I_do = 3 + (nb-1) * 4 + no/2
               if (nb.lt.3) then
               J_up = 1 + nb * 4 + no/2
               J_do = 3 + nb * 4 + no/2
               else
               J_up = 1 + no/2
               J_do = 3 + no/2
               endif
               select case (nz)
               case (1)
               Op_V(nc,nf)%O(I_up, J_up) =  cmplx( 0, X, kind(0.D0))
               Op_V(nc,nf)%O(J_up, I_up) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(I_do, J_do) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(J_do, I_do) =  cmplx( 0, X, kind(0.D0))
               case (2)
               Op_V(nc,nf)%O(I_up, J_do) =  cmplx( 0, X, kind(0.D0))
               Op_V(nc,nf)%O(J_do, I_up) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(J_up, I_do) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(I_do, J_up) =  cmplx( 0, X, kind(0.D0))
               case (3)
               Op_V(nc,nf)%O(I_up, J_do) =  cmplx(-X, 0, kind(0.D0))
               Op_V(nc,nf)%O(J_do, I_up) =  cmplx(-X, 0, kind(0.D0))   
               Op_V(nc,nf)%O(J_up, I_do) =  cmplx( X, 0, kind(0.D0))
               Op_V(nc,nf)%O(I_do, J_up) =  cmplx( X, 0, kind(0.D0))
               end select
               enddo
               enddo
               Op_V(nc,nf)%g = SQRT(CMPLX(DTAU*ham_lamb/dble(N_SUN), 0.D0, kind(0.D0)))
               Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
               Op_V(nc,nf)%type   = 2
               Call Op_set( Op_V(nc,nf) )
              enddo
             Enddo   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
          Elseif (Model == "QSHlyh_Z2") then
!-----------------------------------------------------------------------------------------------------!
             Allocate(Op_V(Latt%N,N_FL))
             do nf = 1,N_FL
               do i = 1, Latt%N
                  Call Op_make(Op_V(i,nf),12)
               enddo
             enddo
               Bond(1,1)=1
               Bond(1,2)=2
               Bond(1,3)=3
               Bond(1,4)=4

               Bond(2,1)=1
               Bond(2,2)=2
               Bond(2,3)=5
               Bond(2,4)=6

               Bond(3,1)=3
               Bond(3,2)=4
               Bond(3,3)=5
               Bond(3,4)=6

               Bond(4,1)=7
               Bond(4,2)=8
               Bond(4,3)=9
               Bond(4,4)=10

               Bond(5,1)=7
               Bond(5,2)=8
               Bond(5,3)=11
               Bond(5,4)=12

               Bond(6,1)=9
               Bond(6,2)=10
               Bond(6,3)=11
               Bond(6,4)=12

               Xp(1)=-1.d0
               Xp(2)=1.d0
               Xp(3)=-1.d0
               Xp(4)=1.d0
               Xp(5)=1.d0
               Xp(6)=-1.d0
!--------------------------------------------------!
             nc=0
             nf=1
             do i = 1,Latt%N
                nc=nc+1
                print*,nc
               I1u = invlist(I,1)
               Op_V(nc,nf)%P(1) = I1u
               I1d = invlist(I,3)
               Op_V(nc,nf)%P(2) = I1d
!--------------------------------------------------!
               I2u = invlist(Latt%nnlist(I,1,0), 1)
               Op_V(nc,nf)%P(3) = I2u
               I2d = invlist(Latt%nnlist(I,1,0), 3)
               Op_V(nc,nf)%P(4) = I2d
!--------------------------------------------------!
               I3u = invlist(Latt%nnlist(I,0,1), 1)
               Op_V(nc,nf)%P(5) = I3u
               I3d = invlist(Latt%nnlist(I,0,1), 3)
               Op_V(nc,nf)%P(6) = I3d
!==================================================!
               J1u = invlist(I,2)
               Op_V(nc,nf)%P(7) = J1u
               J1d = invlist(I,4)
               Op_V(nc,nf)%P(8) = J1d
!--------------------------------------------------!
               J2u = invlist(Latt%nnlist(I,1,0), 2)
               Op_V(nc,nf)%P(9) = J2u
               J2d = invlist(Latt%nnlist(I,1,0), 4)
               Op_V(nc,nf)%P(10) = J2d
!--------------------------------------------------!
               J3u = invlist(Latt%nnlist(I,1,-1), 2)
               Op_V(nc,nf)%P(11) = J3u
               J3d = invlist(Latt%nnlist(I,1,-1), 4)
               Op_V(nc,nf)%P(12) = J3d
!--------------------------------------------------!
              do nb=1,6
                 Op_V(nc,nf)%O(Bond(nb,1),Bond(nb,3)) =  cmplx( 0, Xp(nb), kind(0.D0))
                 Op_V(nc,nf)%O(Bond(nb,3),Bond(nb,1)) =  cmplx( 0,-Xp(nb), kind(0.D0))
                 Op_V(nc,nf)%O(Bond(nb,2),Bond(nb,4)) =  cmplx( 0,-Xp(nb), kind(0.D0))
                 Op_V(nc,nf)%O(Bond(nb,4),Bond(nb,2)) =  cmplx( 0, Xp(nb), kind(0.D0))

!                 Op_V(nc,nf)%O(Bond(nb,1),Bond(nb,4)) =  cmplx( 0, Xp(nb), kind(0.D0))
!                 Op_V(nc,nf)%O(Bond(nb,2),Bond(nb,3)) =  cmplx( 0, Xp(nb), kind(0.D0))
!                 Op_V(nc,nf)%O(Bond(nb,3),Bond(nb,2)) =  cmplx( 0,-Xp(nb), kind(0.D0))
!                 Op_V(nc,nf)%O(Bond(nb,4),Bond(nb,1)) =  cmplx( 0,-Xp(nb), kind(0.D0))

!                 Op_V(nc,nf)%O(Bond(nb,1),Bond(nb,4)) =  cmplx( -Xp(nb),0, kind(0.D0))
!                 Op_V(nc,nf)%O(Bond(nb,2),Bond(nb,3)) =  cmplx(  Xp(nb),0, kind(0.D0))
!                 Op_V(nc,nf)%O(Bond(nb,3),Bond(nb,2)) =  cmplx(  Xp(nb),0, kind(0.D0))
!                 Op_V(nc,nf)%O(Bond(nb,4),Bond(nb,1)) =  cmplx( -Xp(nb),0, kind(0.D0))

              end do
               Op_V(nc,nf)%g = SQRT(CMPLX( DTAU*ham_lamb, 0.D0, kind(0.D0)))
               Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
               Op_V(nc,nf)%type   = 2
               Call Op_set( Op_V(nc,nf) )
               end do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
          Elseif (Model == "QSH_SU2lyh") then
!-----------------------------------------------------------------------------------------------------!
             Allocate(Op_V(Latt%N,N_FL*3))
             do nf = 1,N_FL
               do i = 1, Latt%N*3
                  Call Op_make(Op_V(i,nf),12)
               enddo
             enddo
               Bond(1,1)=1
               Bond(1,2)=2
               Bond(1,3)=3
               Bond(1,4)=4

               Bond(2,1)=1
               Bond(2,2)=2
               Bond(2,3)=5
               Bond(2,4)=6

               Bond(3,1)=3
               Bond(3,2)=4
               Bond(3,3)=5
               Bond(3,4)=6

               Bond(4,1)=7
               Bond(4,2)=8
               Bond(4,3)=9
               Bond(4,4)=10

               Bond(5,1)=7
               Bond(5,2)=8
               Bond(5,3)=11
               Bond(5,4)=12

               Bond(6,1)=9
               Bond(6,2)=10
               Bond(6,3)=11
               Bond(6,4)=12

               Xp(1)=-1.d0
               Xp(2)=1.d0
               Xp(3)=-1.d0
               Xp(4)=1.d0
               Xp(5)=1.d0
               Xp(6)=-1.d0
!-------------------------------------------------!
             nc=0
             nf=1
             do i = 1,Latt%N
                do nz=1,3
                nc=nc+1
               I1u = invlist(I,1)
               Op_V(nc,nf)%P(1) = I1u
               I1d = invlist(I,3)
               Op_V(nc,nf)%P(2) = I1d
!--------------------------------------------------!
               I2u = invlist(Latt%nnlist(I,1,0), 1)
               Op_V(nc,nf)%P(3) = I2u
               I2d = invlist(Latt%nnlist(I,1,0), 3)
               Op_V(nc,nf)%P(4) = I2d
!--------------------------------------------------!
               I3u = invlist(Latt%nnlist(I,0,1), 1)
               Op_V(nc,nf)%P(5) = I3u
               I3d = invlist(Latt%nnlist(I,0,1), 3)
               Op_V(nc,nf)%P(6) = I3d
!==================================================!
               J1u = invlist(I,2)
               Op_V(nc,nf)%P(7) = J1u
               J1d = invlist(I,4)
               Op_V(nc,nf)%P(8) = J1d
!--------------------------------------------------!
               J2u = invlist(Latt%nnlist(I,1,0), 2)
               Op_V(nc,nf)%P(9) = J2u
               J2d = invlist(Latt%nnlist(I,1,0), 4)
               Op_V(nc,nf)%P(10) = J2d
!--------------------------------------------------!
               J3u = invlist(Latt%nnlist(I,1,-1), 2)
               Op_V(nc,nf)%P(11) = J3u
               J3d = invlist(Latt%nnlist(I,1,-1), 4)
               Op_V(nc,nf)%P(12) = J3d
!--------------------------------------------------!
               do nb=1,6
                 select case (nz)
                 case (1)
                 Op_V(nc,nf)%O(Bond(nb,1),Bond(nb,3)) =  cmplx( 0, Xp(nb), kind(0.D0))
                 Op_V(nc,nf)%O(Bond(nb,3),Bond(nb,1)) =  cmplx( 0,-Xp(nb), kind(0.D0))
                 Op_V(nc,nf)%O(Bond(nb,2),Bond(nb,4)) =  cmplx( 0,-Xp(nb), kind(0.D0))
                 Op_V(nc,nf)%O(Bond(nb,4),Bond(nb,2)) =  cmplx( 0, Xp(nb), kind(0.D0))
                 case (2)
                 Op_V(nc,nf)%O(Bond(nb,1),Bond(nb,4)) =  cmplx( 0, Xp(nb), kind(0.D0))
                 Op_V(nc,nf)%O(Bond(nb,2),Bond(nb,3)) =  cmplx( 0, Xp(nb), kind(0.D0))
                 Op_V(nc,nf)%O(Bond(nb,3),Bond(nb,2)) =  cmplx( 0,-Xp(nb), kind(0.D0))
                 Op_V(nc,nf)%O(Bond(nb,4),Bond(nb,1)) =  cmplx( 0,-Xp(nb), kind(0.D0))
                 case (3)
                 Op_V(nc,nf)%O(Bond(nb,1),Bond(nb,4)) =  cmplx( -Xp(nb),0, kind(0.D0))
                 Op_V(nc,nf)%O(Bond(nb,2),Bond(nb,3)) =  cmplx(  Xp(nb),0, kind(0.D0))
                 Op_V(nc,nf)%O(Bond(nb,3),Bond(nb,2)) =  cmplx(  Xp(nb),0, kind(0.D0))
                 Op_V(nc,nf)%O(Bond(nb,4),Bond(nb,1)) =  cmplx( -Xp(nb),0, kind(0.D0))
                 end select
               end do

               Op_V(nc,nf)%g = SQRT(CMPLX( DTAU*ham_lamb, 0.D0, kind(0.D0)))
               Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
               Op_V(nc,nf)%type   = 2
               Call Op_set( Op_V(nc,nf) )
               end do
               end do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
          Elseif (Model == "QSH0_Z2") then

             Allocate(Op_V(Latt%N*6,N_FL))
             do nf = 1,N_FL
               do i = 1, Latt%N*6  
                  Call Op_make(Op_V(i,nf),4)  
               enddo
             enddo
             nc = 0
             nf = 1
             Do I  = 1,Latt%N
              do nb = 1, 3    
              do no = 1, 2
              
               nc = nc + 1
               I1 = invlist(I,1+no/2)
               Op_V(nc,nf)%P(1) = I1
               I2 = invlist(I,3+no/2) 
               Op_V(nc,nf)%P(3) = I2               

               select case (nb)
               case (1)
               J1 = invlist(Latt%nnlist(I,1,0), 1+no/2)
               Op_V(nc,nf)%P(2) = J1              
               J2 = invlist(Latt%nnlist(I,1,0), 3+no/2)
               Op_V(nc,nf)%P(4) = J2
               X = 1.d0 * (-1) ** no

               case (2)
               J1 = invlist(Latt%nnlist(I,0,1), 1+no/2)
               Op_V(nc,nf)%P(2) = J1              
               J2 = invlist(Latt%nnlist(I,0,1), 3+no/2)
               Op_V(nc,nf)%P(4) = J2         
               X = -1.d0 * (-1) ** no

               case (3)
               J1 = invlist(Latt%nnlist(I,-1,1),1+no/2)
               Op_V(nc,nf)%P(2) = J1    
               J2 = invlist(Latt%nnlist(I,-1,1),3+no/2)
               Op_V(nc,nf)%P(4) = J2
               X = 1.d0 * (-1) ** no
               end select
                
               Op_V(nc,nf)%O(1,2) =  cmplx( 0, X, kind(0.D0))
               Op_V(nc,nf)%O(2,1) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(3,4) =  cmplx( 0,-X, kind(0.D0))
               Op_V(nc,nf)%O(4,3) =  cmplx( 0, X, kind(0.D0))

!               Op_V(nc,nf)%O(1,4) =  cmplx( 0, X, kind(0.D0))
!               Op_V(nc,nf)%O(4,1) =  cmplx( 0,-X, kind(0.D0))
!               Op_V(nc,nf)%O(2,3) =  cmplx( 0,-X, kind(0.D0))
!               Op_V(nc,nf)%O(3,2) =  cmplx( 0, X, kind(0.D0))

!               Op_V(nc,nf)%O(1,4) =  cmplx(-X, 0, kind(0.D0))
!               Op_V(nc,nf)%O(4,1) =  cmplx(-X, 0, kind(0.D0))
!               Op_V(nc,nf)%O(2,3) =  cmplx( X, 0, kind(0.D0))
!               Op_V(nc,nf)%O(3,2) =  cmplx( X, 0, kind(0.D0))


               Op_V(nc,nf)%g = SQRT(CMPLX( DTAU*ham_lamb, 0.D0, kind(0.D0)))
               Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
               Op_V(nc,nf)%type   = 2
               Call Op_set( Op_V(nc,nf) )

              enddo
              enddo
             Enddo

          Elseif (Model == "Hubbard_SSC") then
             Allocate(Op_V(Latt%N*2,N_FL))
             do nf = 1,N_FL
                do i  = 1, Latt%N*2
                   Call Op_make(Op_V(i,nf),2)
                enddo
             enddo
             Do nf = 1, N_FL
                nc = 0
                Do i = 1, Latt%N
                Do no = 1, 2
                   nc = nc +1
                   Op_V(nc,nf)%P(1) = invlist(I, 1 + no/2)
                   Op_V(nc,nf)%P(2) = invlist(I, 3 + no/2)
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0)) 
                   Op_V(nc,nf)%O(2,2) = cmplx(1.d0, 0.d0, kind(0.D0))
                   Op_V(nc,nf)%g    = SQRT(CMPLX(-DTAU*ham_U/2.d0, 0.D0, kind(0.D0)))
                   Op_V(nc,nf)%alpha  = cmplx(-1.0d0, 0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                Enddo
                Enddo
              Enddo

          Elseif (Model == "Hubbard_Mz")  then
             !Write(50,*) 'Model is ', Model
             Allocate(Op_V(Latt%N + n_lambda ,N_FL))
             do nf = 1,N_FL
                do i  = 1, Latt%N
                   Call Op_make(Op_V(i,nf),1)
                enddo
                do nth = 1,n_lambda
                   Call Op_make(Op_V(Latt%N+nth,nf),Latt%N)
                enddo
             enddo
             Do nf = 1,N_FL
                nc = 0
                X = 1.d0
                if (nf == 2) X = -1.d0
                Do i = 1,Latt%N
                   nc = nc + 1
                   Op_V(nc,nf)%P(1) = I
                   Op_V(nc,nf)%O(1,1) = cmplx(1.d0, 0.d0, kind(0.D0))
                   Op_V(nc,nf)%g      = X*SQRT(CMPLX(DTAU*ham_U/2.d0, 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                   ! The operator reads:  
                   ! g*s*( c^{dagger} O c   - alpha ))
                   ! with s the HS field.
                   ! Write(6,*) nc,nf, Op_V(nc,nf)%g 
                Enddo
                Do nth = 1,n_Lambda
                   nc = nc + 1
                   Do  I = 1,latt%N
                      Op_V(nc,nf)%P(I) = I
                      Op_V(nc,nf)%O(I,I) = cmplx(1.d0  ,0.d0, kind(0.D0))
                   enddo
                   Op_V(nc,nf)%g      = SQRT(CMPLX(-DTAU*Lambda/dble(n_lambda), 0.D0, kind(0.D0))) 
                   Op_V(nc,nf)%alpha  = cmplx(-DBLE(Latt%N)*0.5d0,0.d0, kind(0.D0))
                   Op_V(nc,nf)%type   = 2
                   Call Op_set( Op_V(nc,nf) )
                enddo
             Enddo
          Endif
        end Subroutine Ham_V

!===================================================================================           
        Real (Kind=8) function S0(n,nt)  
          Implicit none
          Integer, Intent(IN) :: n,nt 
          S0 = 1.d0
          
        end function S0

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

          If ((Model == "QSH_SU2").or.(Model == "QSH_Z2").or.(Model == "Hubbard_SSC").or.(Model == "QSH_U1")) then
          Allocate ( Obs_eq(6) )
          Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Ns = Latt%N;  No = Norb;  Filename ="Green"
             case (2)
                Ns = Latt%N;  No = 6;     Filename ="Bond_QSH"
             case (3)
                Ns = Latt%N;  No = 2;     Filename ="SSC"
             case (4)
                Ns = Latt%N;  No = 1;     Filename ="Bond_tr"       
             case (5)
                Ns = Latt%N;  No = 2;     Filename ="SpinZ"     
             case (6)
                Ns = Latt%N;  No = 2;     Filename ="Den"
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Nt = 1
             Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
          enddo

          Else
          Allocate ( Obs_eq(4) )
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
             case default
                Write(6,*) ' Error in Alloc_obs '  
             end select
             Nt = 1
             Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
          enddo
          Endif


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
          
          Implicit none
          
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer, INTENT(IN)          :: Ntau
          
          !Local 
          Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z, ZP,ZS,Z1,Z2,Z3,Z4
          Integer :: I,J, imj, nf, dec, I1, J1, no_I, no_J, Sz
          Integer :: norb_I, norb_J, I1_u, I1_d, I2_u, I2_d      
          Integer :: J1_u, J1_d, J2_u, J2_d, nb_I, nb_J, X_I, X_J


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
          Do nf = 1,N_FL
             Do J = 1,Ndim
                Zkin = Zkin + sum(Op_T(1,nf)%O(:, j)*Grc(:, j, nf))
             ENddo
          Enddo
          Zkin = Zkin * dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS

          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising" ) then
            dec = 1
          else
            dec = 2
          endif
          Do I = 1,Ndim
!             ZPot = ZPot + Grc(i,i,1) * Grc(i,i, dec)
          Enddo
!          Zpot = Zpot*ham_U
!          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS

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
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising"  ) then 
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
                   ! Den
                   Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J)  +  &
                        &     (    GRC(I1,I1,1) * GRC(J1,J1,1) *Z     + &
                        &          GRC(I1,J1,1) * GR(I1,J1,1 )           &
                        &                                   ) * Z* ZP*ZS
                ENDDO
                Obs_eq(4)%Obs_Latt0(no_I) =  Obs_eq(4)%Obs_Latt0(no_I) +  Z * GRC(I1,I1,1) * ZP * ZS
             ENDDO

          elseif ((Model == "QSH_SU2").or.(Model == "QSH_Z2").or.(Model == "Hubbard_SSC").or.(Model == "QSH_U1")) Then
            ! Green
            Do I = 1, Latt%N
               Do norb_I = 1, Norb
                  I1 = invlist(I, norb_I)             
                  Do J = 1, Latt%N
                  Do norb_J = 1, Norb
                  J1 = invlist(J, norb_J)    
                  imj = latt%imj(I,J)               
                  Obs_eq(1)%Obs_Latt(imj,1,norb_I,norb_J) = Obs_eq(1)%Obs_Latt(imj,1,norb_I,norb_J) + &
                      &       GRC(I1,J1,1) *  ZP*ZS
                  Enddo
                  Enddo
               Enddo           
            Enddo

            ! SSC 
            Do I = 1, Latt%N
               Do norb_I = 1, 2
                  I1_u = invlist(I, norb_I)
                  I1_d = invlist(I, norb_I + 2)
                  Do J = 1, Latt%N
                  Do norb_J = 1, 2
                  J1_u = invlist(J, norb_J)
                  J1_d = invlist(J, norb_J + 2)
                  imj = latt%imj(I,J)

                  Z = GRC(J1_u,I1_u,1)*GRC(J1_d,I1_d,1) + GRC(J1_u,I1_d,1)*GR(I1_u,J1_d,1) &
                  & + GRC(I1_d,J1_d,1)*GRC(I1_u,J1_u,1) + GRC(I1_d,J1_u,1)*GR(J1_d,I1_u,1)
                  Obs_eq(3)%Obs_Latt(imj,1,norb_I,norb_J) = Obs_eq(3)%Obs_Latt(imj,1,norb_I,norb_J) + &
                      &  Z * ZP * ZS

                  Z = GRC(I1_u,J1_u,1) * GR(I1_u,J1_u,1) + GRC(I1_d,J1_d,1) * GR(I1_d,J1_d,1) &
                  & - GRC(I1_u,J1_d,1) * GR(I1_u,J1_d,1) - GRC(I1_d,J1_u,1) * GR(I1_d,J1_u,1) &
                  & +(GRC(I1_u,I1_u,1) -GRC(I1_d,I1_d,1))*(GRC(J1_u,J1_u,1) -GRC(J1_d,J1_d,1))
                  Obs_eq(5)%Obs_Latt(imj,1,norb_I,norb_J) = Obs_eq(5)%Obs_Latt(imj,1,norb_I,norb_J) + &
                      &  Z * ZP * ZS

                  Z = GRC(I1_u,J1_u,1) * GR(I1_u,J1_u,1) + GRC(I1_d,J1_d,1) * GR(I1_d,J1_d,1) &
                  & + GRC(I1_u,J1_d,1) * GR(I1_u,J1_d,1) + GRC(I1_d,J1_u,1) * GR(I1_d,J1_u,1) &
                  & +(GRC(I1_u,I1_u,1) +GRC(I1_d,I1_d,1))*(GRC(J1_u,J1_u,1) +GRC(J1_d,J1_d,1))
                  Obs_eq(6)%Obs_Latt(imj,1,norb_I,norb_J) = Obs_eq(6)%Obs_Latt(imj,1,norb_I,norb_J) + &
                      &  Z * ZP * ZS
                 
                  Enddo
                  Enddo
                  Z = GRC(I1_u,I1_u,1) + GRC(I1_d,I1_d,1)
                  Obs_eq(6)%Obs_Latt0(norb_I) = Obs_eq(6)%Obs_Latt0(norb_I) + Z * ZP * ZS
               Enddo
            Enddo

            ! QSH correlation
            Do  I = 1, Latt%N
              norb_I  = 0
              do nb_I = 1,3
              do no_I = 1,2
                norb_I = norb_I + 1
                I1_u = invlist(I,1+no_I/2)
                I1_d = invlist(I,3+no_I/2)                                  
                select case (nb_I)
                case (1)
                I2_u = invlist(Latt%nnlist(I,1,0), 1+no_I/2)                            
                I2_d = invlist(Latt%nnlist(I,1,0), 3+no_I/2)               
                X_I = 1.d0 * (-1) ** no_I
                case (2)
                I2_u = invlist(Latt%nnlist(I,0,1), 1+no_I/2)
                I2_d = invlist(Latt%nnlist(I,0,1), 3+no_I/2)
                X_I =-1.d0 * (-1) ** no_I
                case (3)
                I2_u = invlist(Latt%nnlist(I,-1,1),1+no_I/2)
                I2_d = invlist(Latt%nnlist(I,-1,1),3+no_I/2)
                X_I = 1.d0 * (-1) ** no_I
                end select

                Do J= 1, Latt%N
                   norb_J  = 0 
                   do nb_J = 1,3
                   do no_J = 1,2
                    norb_J = norb_J + 1
                    J1_u = invlist(J,1+no_J/2)
                    J1_d = invlist(J,3+no_J/2)                                  
                    select case (nb_J)
                    case (1)
                    J2_u = invlist(Latt%nnlist(J,1,0), 1+no_J/2)                            
                    J2_d = invlist(Latt%nnlist(J,1,0), 3+no_J/2)               
                    X_J = 1.d0 * (-1) ** no_J
                    case (2)
                    J2_u = invlist(Latt%nnlist(J,0,1), 1+no_J/2)
                    J2_d = invlist(Latt%nnlist(J,0,1), 3+no_J/2)
                    X_J =-1.d0 * (-1) ** no_J
                    case (3)
                    J2_u = invlist(Latt%nnlist(J,-1,1),1+no_J/2)
                    J2_d = invlist(Latt%nnlist(J,-1,1),3+no_J/2)
                    X_J = 1.d0 * (-1) ** no_J
                    end select

                    imj = latt%imj(I,J) 
                    Sz = 1
                    Z1 =  GRC(I1_u,I2_u,1) * GRC(J1_u,J2_u,1) + GRC(I1_d,I2_d,1) * GRC(J1_d,J2_d,1) &
                    &   + GRC(I1_u,J2_u,1) * GR( I2_u,J1_u,1) + GRC(I1_d,J2_d,1) * GR( I2_d,J1_d,1) &
                  &+Sz*(- GRC(I1_u,I2_u,1) * GRC(J1_d,J2_d,1) - GRC(I1_d,I2_d,1) * GRC(J1_u,J2_u,1) &
                    &   - GRC(I1_u,J2_d,1) * GR( I2_u,J1_d,1) - GRC(I1_d,J2_u,1) * GR( I2_d,J1_u,1))

                    Z2 =  GRC(I2_u,I1_u,1) * GRC(J1_u,J2_u,1) + GRC(I2_d,I1_d,1) * GRC(J1_d,J2_d,1) &
                    &   + GRC(I2_u,J2_u,1) * GR( I1_u,J1_u,1) + GRC(I2_d,J2_d,1) * GR( I1_d,J1_d,1) &
                  &+Sz*(- GRC(I2_u,I1_u,1) * GRC(J1_d,J2_d,1) - GRC(I2_d,I1_d,1) * GRC(J1_u,J2_u,1) &
                    &   - GRC(I2_u,J2_d,1) * GR( I1_u,J1_d,1) - GRC(I2_d,J2_u,1) * GR( I1_d,J1_u,1))

                    Z3 =  GRC(I1_u,I2_u,1) * GRC(J2_u,J1_u,1) + GRC(I1_d,I2_d,1) * GRC(J2_d,J1_d,1) &
                    &   + GRC(I1_u,J1_u,1) * GR( I2_u,J2_u,1) + GRC(I1_d,J1_d,1) * GR( I2_d,J2_d,1) &
                  &+Sz*(- GRC(I1_u,I2_u,1) * GRC(J2_d,J1_d,1) - GRC(I1_d,I2_d,1) * GRC(J2_u,J1_u,1) &
                    &   - GRC(I1_u,J1_d,1) * GR( I2_u,J2_d,1) - GRC(I1_d,J1_u,1) * GR( I2_d,J2_u,1))

                    Z4 =  GRC(I2_u,I1_u,1) * GRC(J2_u,J1_u,1) + GRC(I2_d,I1_d,1) * GRC(J2_d,J1_d,1) &
                    &   + GRC(I2_u,J1_u,1) * GR( I1_u,J2_u,1) + GRC(I2_d,J1_d,1) * GR( I1_d,J2_d,1) &
                  &+Sz*(- GRC(I2_u,I1_u,1) * GRC(J2_d,J1_d,1) - GRC(I2_d,I1_d,1) * GRC(J2_u,J1_u,1) &
                    &   - GRC(I2_u,J1_d,1) * GR( I1_u,J2_d,1) - GRC(I2_d,J1_u,1) * GR( I1_d,J2_u,1))

                    Z = - Z1 + Z2 + Z3 - Z4
                    Obs_eq(2)%Obs_Latt(imj,1,norb_I,norb_J) = Obs_eq(2)%Obs_Latt(imj,1,norb_I,norb_J) &
                    & + X_I * X_J * Z * ZP * ZS

                    if (norb_I.eq.norb_J) then
                    Obs_eq(4)%Obs_Latt(imj,1,1,1) = Obs_eq(4)%Obs_Latt(imj,1,1,1) &
                    & + X_I * X_J * Z * ZP * ZS
                    endif

                   enddo
                   enddo
                Enddo
                Obs_eq(2)%Obs_Latt0(norb_I) = Obs_eq(2)%Obs_Latt0(norb_I) &
                & + GRC(I1_u,I2_u,1) - GRC(I1_d,I2_d,1)
              enddo
              enddo
            Enddo



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
                   Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + &
                        & (   GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                        &    (GRC(I1,I1,2) - GRC(I1,I1,1))*(GRC(J1,J1,2) - GRC(J1,J1,1))    ) * ZP*ZS
                   ! SpinXY
                   ! c^d_(i,u) c_(i,d) c^d_(j,d) c_(j,u)  +  c^d_(i,d) c_(i,u) c^d_(j,u) c_(j,d)
                   Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + &
                     & (   GRC(I1,J1,1) * GR(I1,J1,2) +  GRC(I1,J1,2) * GR(I1,J1,1)    ) * ZP*ZS
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
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising"  ) then 
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

!========================================================================

!=====================================================
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
!====================================================================================


!---------------------------------------------------------------------
        Real (Kind=kind(0.d0)) Function Delta_S0_global(Nsigma_old)

          !>  This function computes the ratio:
          !e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
          Implicit none 
          
          !> Arguments
          Integer, dimension(:,:), allocatable, intent(IN) :: Nsigma_old
        end Function Delta_S0_global
!---------------------------------------------------------------------





!---------------------------------------------------------------------
        Subroutine  Hamiltonian_set_random_nsigma
          
          ! The user can set the initial configuration
          
          Implicit none
          
          Integer :: I, nt
          
          print*, Size(OP_V,1)
          Do nt = 1,Ltrot
             Do I = 1,Size(OP_V,1)
                nsigma(I,nt)  = 1
                if ( ranf_wrap()  > 0.5D0 ) nsigma(I,nt)  =-1
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

        end Subroutine Global_move_tau



      end Module Hamiltonian
