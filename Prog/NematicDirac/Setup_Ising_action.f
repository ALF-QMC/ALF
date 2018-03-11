        Subroutine Setup_Ising_action
          
          ! This subroutine sets up lists and arrays so as to enable an 
          ! an efficient calculation of  S0(n,nt) 

          Integer :: nc, nth, n, n1, n2, n3, n4, I, I1, n_orientation
          Real (Kind=Kind(0.d0)) :: X_p(2)
        
          If  ( Model == "Hubbard_SU2_Ising" ) then
            ! Setup list of bonds for the square lattice.
            Allocate (L_Bond(Latt%N,2),  L_bond_inv(Latt%N*N_coord,2) )
          
            nc = 0
            do nth = 1,2*N_coord  
               Do n1= 1, L1/2
                  Do n2 = 1,L2
                     nc = nc + 1
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
          elseIF ( Lattice_type =="BipartiteSquare" .and. Model == "NematicDirac" ) then
            allocate(Ising_nnlist(Latt%N,4))
            N_ising = Latt%N
            do I = 1, Latt%N
               Ising_nnlist(I,1) = Latt%nnlist(I, 1, 0)
               Ising_nnlist(I,2) = Latt%nnlist(I, 0, 1)
               Ising_nnlist(I,3) = Latt%nnlist(I,-1, 0)
               Ising_nnlist(I,4) = Latt%nnlist(I, 0,-1)
            enddo
          elseIF ( Lattice_type =="BipartiteSquare" .and. Model == "NematicDirac3" ) then
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
          elseIF ( Lattice_type =="Square" .and. Model == "NematicDirac2" ) then
            allocate(Ising_nnlist(Latt%N,4))
            N_ising = Latt%N
            do I = 1, Latt%N
               Ising_nnlist(I,1) = Latt%nnlist(I, 1, 0)
               Ising_nnlist(I,2) = Latt%nnlist(I, 0, 1)
               Ising_nnlist(I,3) = Latt%nnlist(I,-1, 0)
               Ising_nnlist(I,4) = Latt%nnlist(I, 0,-1)
            enddo
          elseIF ( Lattice_type =="Square" .and. Model == "yyhe" ) then
            allocate(Ising_nnlist(Ndim,4))
            N_ising = Ndim
            do I = 1, Latt%N
              !no=1
              I1 = Invlist(I,1)
              Ising_nnlist(I1,1) = Invlist(I,2)
              Ising_nnlist(I1,2) = Invlist(Latt%nnlist(I, 0,-1),2) 
              Ising_nnlist(I1,3) = Invlist(Latt%nnlist(I,-1,-1),2)
              Ising_nnlist(I1,4) = Invlist(Latt%nnlist(I,-1, 0),2)
              
              !no=2
              I1 = Invlist(I,2)
              Ising_nnlist(I1,1) = Invlist(Latt%nnlist(I, 1, 1),1)
              Ising_nnlist(I1,2) = Invlist(Latt%nnlist(I, 1, 0),1) 
              Ising_nnlist(I1,3) = Invlist(I,1)
              Ising_nnlist(I1,4) = Invlist(Latt%nnlist(I, 0, 1),1)
            enddo
          else
            Write(6,*) ' Error in Setup_Ising_action '
            Stop
          endif
          
          DW_Ising_tau  ( 1) = tanh(Dtau*Ham_h)
          DW_Ising_tau  (-1) = 1.D0/DW_Ising_tau(1)
          DW_Ising_Space( 1) = exp(-2.d0*Dtau*Ham_J) 
          DW_Ising_Space(-1) = exp( 2.d0*Dtau*Ham_J) 
          
          addProb_space = 1 - exp(-2.d0*Dtau*Ham_J)
          addProb_tau   = 1 - tanh(Dtau*Ham_h)
          Geo_addProb_space = 1 - exp(-4.d0*Dtau*Ham_J)
          Geo_addProb_tau   = 1 - tanh(Dtau*Ham_h)**2
        End Subroutine Setup_Ising_action
