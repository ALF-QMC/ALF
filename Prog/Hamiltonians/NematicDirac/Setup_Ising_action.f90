        Subroutine Setup_Ising_action()

          ! This subroutine sets up lists and arrays so as to enable an
          ! an efficient calculation of  S0(n,nt)

          Implicit none
          Integer :: I
          n_global = 0

          select case( Model_vers )
          case( 0,1 )
            allocate(Ising_nnlist(Latt%N,4))
            N_ising = Latt%N
            do I = 1, Latt%N
               Ising_nnlist(I,1) = Latt%nnlist(I, 1, 0)
               Ising_nnlist(I,2) = Latt%nnlist(I, 0, 1)
               Ising_nnlist(I,3) = Latt%nnlist(I,-1, 0)
               Ising_nnlist(I,4) = Latt%nnlist(I, 0,-1)
            enddo
          case( 3 )
            allocate(Ising_nnlist(2*Latt%N,4))
            N_ising = Latt%N
            do I = 1, Latt%N
               Ising_nnlist(I,1) = Latt%nnlist(I, 1, 0)
               Ising_nnlist(I,2) = Latt%nnlist(I, 0, 1)
               Ising_nnlist(I,3) = Latt%nnlist(I,-1, 0)
               Ising_nnlist(I,4) = Latt%nnlist(I, 0,-1)
               Ising_nnlist(I+Latt%N,1) = Latt%nnlist(I, 1, 0) + Latt%N
               Ising_nnlist(I+Latt%N,2) = Latt%nnlist(I, 0, 1) + Latt%N
               Ising_nnlist(I+Latt%N,3) = Latt%nnlist(I,-1, 0) + Latt%N
               Ising_nnlist(I+Latt%N,4) = Latt%nnlist(I, 0,-1) + Latt%N
            enddo
          case( 2 )
            allocate(Ising_nnlist(Latt%N,4))
            N_ising = Latt%N
            do I = 1, Latt%N
               Ising_nnlist(I,1) = Latt%nnlist(I, 1, 0)
               Ising_nnlist(I,2) = Latt%nnlist(I, 0, 1)
               Ising_nnlist(I,3) = Latt%nnlist(I,-1, 0)
               Ising_nnlist(I,4) = Latt%nnlist(I, 0,-1)
            enddo
          case default
            Write(6,*) ' Error in Setup_Ising_action '
            Stop
          end select

          !print*, Global_J, Global_h

          DW_Ising_tau  ( 1) = tanh(Dtau*Ham_h)
          DW_Ising_tau  (-1) = 1.D0/DW_Ising_tau(1)
          DW_Ising_Space( 1) = exp(-2.d0*Dtau*Ham_J)
          DW_Ising_Space(-1) = exp( 2.d0*Dtau*Ham_J)

          Wolff_addProb_space = 1 - exp(-2.d0*Dtau*Global_J)
          Wolff_addProb_tau   = 1 - tanh(Dtau*Global_h)

          Geo_addProb_space = 1 - exp(-4.d0*Dtau*Ham_J)
          Geo_addProb_tau   = 1 - tanh(Dtau*Ham_h)**2
        End Subroutine Setup_Ising_action
