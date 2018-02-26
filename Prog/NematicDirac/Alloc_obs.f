        Subroutine  Alloc_obs(Ltau) 

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N, Ns,Nt,No
          Character (len=64) ::  Filename
          
          nBlub2 = 100
          
          IF ( (Lattice_type =="BipartiteSquare" .and. Model == "NematicDirac")  .or. &
             & (Lattice_type =="Square"          .and. Model == "NematicDirac2") .or. &
             & (Lattice_type =="BipartiteSquare" .and. Model == "NematicDirac3")  ) then

          ! Scalar observables
          Allocate ( Obs_scal(6) )
          Do I = 1,Size(Obs_scal,1)
            select case (I)
            case (1)
              N = 1;   Filename ="Kin"
            case (2)
              N = 1;   Filename ="ising_z"
            case (3)
              N = 1;   Filename ="F_by_xi"
            case (4)
              N = 1;   Filename ="ising_x"
              eq_x_ising  = tanh(Dtau*Ham_h)
              neq_x_ising = 1/tanh(Dtau*Ham_h)
            case (5)
              N = 3;   Filename ="m"
            case (6)
              N = 2; Filename ="chi2"
!             case (7)
!               N = 1; Filename ="m_tau"
!             case (8)
!               N = 1; Filename ="m2_tau"
!             case (9)
!               N = 500; Filename ="m_auto"
!             case (10)
!               N = 500; Filename ="m2_auto"
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
              Ns = Latt%N;  No = 1;        Filename ="IsingX"
              Nt = 1
            case (2)
              Ns = Latt%N;  No = 1;        Filename ="IsingZ"
              Nt = 1
            case (3)
              Ns = Latt%N;  No = Norb;     Filename ="Green"
              Nt = 1
            case (4)
              Ns = Latt%N;  No = 1;     Filename ="IsingZT"
              Nt = Ltrot+1
            case (5)
              Ns = Latt%N;  No = 1;     Filename ="IsingXT"
              Nt = Ltrot+1
            case (6)
              Ns = Latt%N;  No = Norb;     Filename ="SpinZ"
              Nt = 1
            case default
              Write(6,*) ' Error in Alloc_obs '  
            end select
!             Nt = 1
            Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
          enddo

          If (Ltau == 1) then 
             ! Equal time correlators
             Allocate ( Obs_tau(1) )
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
          
          elseIF ( Lattice_type =="Square" .and. Model == "yyhe" ) then
          
          ! Scalar observables
          Allocate ( Obs_scal(5) )
          Do I = 1,Size(Obs_scal,1)
            select case (I)
            case (1)
              N = 1;   Filename ="Kin"
            case (2)
              N = 1;   Filename ="ising_z"
            case (3)
              N = 1;   Filename ="F_by_xi"
            case (4)
              N = 1;   Filename ="ising_x"
              eq_x_ising  = tanh(Dtau*Ham_h)
              neq_x_ising = 1/tanh(Dtau*Ham_h)
            case (5)
              N = 3;   Filename ="m"
            case default
              Write(6,*) ' Error in Alloc_obs '  
            end select
            Call Obser_Vec_make(Obs_scal(I),N,Filename)
          enddo


          ! Equal time correlators
          Allocate ( Obs_eq(3) )
          Do I = 1,Size(Obs_eq,1)
            select case (I)
            case (1)
              Ns = Latt%N;  No = 2;        Filename ="IsingX"
              Nt = 1
            case (2)
              Ns = Latt%N;  No = 2;        Filename ="IsingZ"
              Nt = 1
            case (3)
              Ns = Latt%N;  No = Norb;     Filename ="Green"
              Nt = 1
            case default
              Write(6,*) ' Error in Alloc_obs '  
            end select
!             Nt = 1
            Call Obser_Latt_make(Obs_eq(I),Ns,Nt,No,Filename)
          enddo
          
          endif

        end Subroutine Alloc_obs
