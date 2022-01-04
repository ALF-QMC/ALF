        Subroutine  Alloc_obs(Ltau)

          Implicit none
          Integer, Intent(In) :: Ltau
          Integer    ::  i, N,Nt, No
          Character (len=64) ::  name
          Character (len=2)  ::  Channel

          nBlub2 = 100

          ! Scalar observables
          If ( Model_vers == 3 ) then
            Allocate ( Obs_scal(11) )
          else
            Allocate ( Obs_scal(9) )
          endif
          Do I = 1,Size(Obs_scal,1)
            select case (I)
            case (1)
              N = 1;   name ="Part"
            case (2)
              N = 1;   name ="ising_z"
            case (3)
              N = 3;   name ="Kin_Pot_E"
            case (4)
              N = 1;   name ="ising_x"
              eq_x_ising  = tanh(Dtau*Ham_h)
              neq_x_ising = 1/tanh(Dtau*Ham_h)
            case (5)
              N = 3;   name ="m"
            case (6)
              N = 2;   name ="chi2"

            case (7)
              N = 1;   name ="ising_z_alt"
            case (8)
              N = 1;   name ="ising_x_alt"
            case (9)
              N = 3;   name ="m_alt"

            case (10)
              N = 1; name ="d_Z"
            case (11)
              N = 1; name ="d_X"
            case default
              Write(6,*) ' Error in Alloc_obs '
            end select
            Call Obs_scal(I)%make(N,name)
          enddo

          ! Histogram observables
          If ( Model_vers == 3 ) then
            Allocate ( Obs_hist(5) )
          else
            Allocate ( Obs_hist(3) )
          endif
          Do I = 1,Size(Obs_hist,1)
            select case (I)
            case (1)
              N = 400; name = "X"
            case (2)
              N = 400; name = "m"
            case (3)
              N = 400; name = "B"
            case (4)
              N = 400; name = "d_Z"
            case (5)
              N = 400; name = "d_X"
            ! case (5)
            !   N = 1; name ="m12"
            ! case (6)
            !   N = 1; name ="X12"
            case default
              Write(6,*) ' Error in Alloc_obs '
            end select
            Call Obs_hist(I)%make(N, 1.5d0, -0.5d0, name)
          enddo


          ! Equal time correlators
          Allocate ( Obs_eq(5) )
          Do I = 1,Size(Obs_eq,1)
            select case (I)
            case (1)
              if ( Model_vers == 3 ) then
                No = 2
              else
                No = 1
              endif
              name ="IsingX"
              Nt = 1
            case (2)
              if ( Model_vers == 3 ) then
                No = 2
              else
                No = 1
              endif
              name ="IsingZ"
              Nt = 1
            case (3)
              No = Norb; name ="Green"
              Nt = 1
            case (4)
              No = 1;     name ="IsingZT"
              Nt = Ltrot+1
            case (5)
              No = 1; name ="IsingXT"
              Nt = Ltrot+1
            case default
              Write(6,*) ' Error in Alloc_obs '
            end select
            !Nt = 1
            Channel = '--'
            Call Obs_eq(I)%make_norb(Nt, No, name, Latt, Channel, dtau)
          enddo

          If (Ltau == 1) then
             ! Equal time correlators
             Allocate ( Obs_tau(1) )
             Do I = 1,Size(Obs_tau,1)
                select case (I)
                case (1)
                   No = Norb;  name ="Green"
                case (2)
                   No = Norb;  name ="GreenUp"
                case (3)
                   No = Norb;  name ="GreenDow"
                case (4)
                   No = Norb;  name ="Den"
                case default
                   Write(6,*) ' Error in Alloc_obs '
                end select
                Nt = Ltrot+1-2*Thtrot
                Channel = '--'
                Call Obs_tau(I)%make_norb(Nt, No, name, Latt, Channel, dtau)
             enddo
          endif

        end Subroutine Alloc_obs
