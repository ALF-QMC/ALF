!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module provides a set of predefined trial wave functions.
!>
!
!--------------------------------------------------------------------

    Module Predefined_Trial

      Use runtime_error_mod
      Use Lattices_v3
      Use Operator_mod
      Use WaveFunction_mod
      Use MyMats
      Use Predefined_Hoppings
      Use Random_wrap
      use iso_fortran_env, only: output_unit, error_unit

      Implicit none

    contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!>    Sets the trial wave function corresponding to the solution of the non-interacting
!>    tight binding Hamiltonian on the given lattice. Twisted boundary conditions (Phi_X=0.01)
!>    are implemented so as to generate a non-degenerate trial wave functions.
!> @param [in]  Lattice_type
!>    Character(64)
!> \verbatim
!>    Square,  Honeycomb, Pi_Flux
!> \endverbatim
!> @param [in]  Latt_unit
!>    Type(Unit_cell)
!> \verbatim
!>     Contains number of orbitals per unit cell and positions, as well as coordination number
!> \endverbatim
!> @param [in]  Ndim
!>    Integer
!> \verbatim
!>     Number of orbitals
!> \endverbatim
!> @param [in]  List, Invlist
!>    Integer(:,:)
!> \verbatim
!>    List(I=1.. Ndim,1)    =   Unit cell of site I
!>    List(I=1.. Ndim,2)    =   Orbital index  of site I
!>    Invlist(Unit_cell,Orbital) = site I
!> \endverbatim
!> @param [in]    Latt
!>    Type(Lattice)
!> \verbatim
!>    The Lattice
!> \endverbatim
!> @param [in]  N_part
!>    Integer
!> \verbatim
!>    Particle number for each flavor
!> \endverbatim
!> @param [in]  N_FL
!>    Integer
!> \verbatim
!>    Flavor
!> \endverbatim
!> @param [out]  WF_L, WF_R
!>    Type(Wavefunction)(N_FL)
!> \verbatim
!>    Wavefunction
!>    Also sets the degeneracy:  E(N_part + 1) - E(N_part). Energy eigenvalues are ordered in ascending order.
!> \endverbatim
!>
!------------------------------------------------------------------
      subroutine Predefined_TrialWaveFunction(Lattice_type, Ndim, List, Invlist, Latt, Latt_unit, &
           &                                  N_part, alpha, N_FL, N_slat, WF_L, WF_R)

        Implicit none
        Character (len=64), Intent(IN)                :: Lattice_type
        Integer, Intent(IN)                           :: Ndim, N_FL, N_slat, N_part
        Integer, Intent(IN), Dimension(:,:)           :: List, Invlist
        real(Kind=kind(0.d0)), intent(in)             :: alpha
        Type(Lattice),   Intent(in)                   :: Latt
        Type(Unit_cell), Intent(in)                   :: Latt_Unit
        Type(WaveFunction), Intent(out), Dimension(:,:), allocatable :: WF_L, WF_R

        Type(Operator),  dimension(:,:), allocatable  :: OP_tmp
        Type (Hopping_Matrix_type), allocatable       :: Hopping_Matrix_tmp(:)
        Real (Kind=Kind(0.d0))                        :: Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, Dimer, mass
        Real (Kind=Kind(0.d0))                        :: ham_U, sgn_i, sgn_updn, rmu
        Logical                                       :: Checkerboard, Symm, Kekule_Trial, hatree_fock, l_cmplx_trial

        Type (Lattice)                                :: Latt_Kekule
        Real (Kind=Kind(0.d0))  :: A1_p(2), A2_p(2), L1_p(2), L2_p(2), x_p(2),x1_p(2), hop(3), del_p(2)
        Real (Kind=Kind(0.d0))  :: delta = 0.01, Ham_T1, Ham_T2, Ham_Tperp

        Integer :: N, nf, I, I1, I2, nc, nc1, IK_u, I_u, J1, lp, J, N_Phi, ns, no, k1
        Logical :: Test=.false. ,  Bulk =.true.
        Complex (Kind=Kind(0.d0)) :: Z_norm

        Real (Kind=Kind(0.d0) ), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
               & ham_tx_vec(:), ham_ty_vec(:), Ham_T2_vec(:), Ham_lambda_vec(:) 
        Integer, allocatable ::   N_Phi_vec(:)

        Allocate(WF_L(N_FL,N_slat),WF_R(N_FL,N_slat))
        do ns=1,N_slat
        do n=1,N_FL
           Call WF_alloc(WF_L(n,ns),Ndim,N_part)
           Call WF_alloc(WF_R(n,ns),Ndim,N_part)
        enddo
        enddo

        Allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
             &    ham_tx_vec(N_FL), ham_ty_vec(N_FL), Ham_lambda_vec(N_FL), N_Phi_vec(N_FL) )

        Checkerboard  = .false.
        Kekule_Trial  = .false.
        hatree_fock   = .false.
        l_cmplx_trial = .false.
        Symm          = .false.

        N_Phi       = 0
        Phi_X       = 0.d0
        Phi_Y       = 0.d0
        Bulk        = .false.
        Ham_T       = 1.d0
        Ham_T2      = 0.d0
        Ham_Tperp   = 0.d0
        Ham_Chem    = 0.d0
        Dtau        = 1.d0

        N_Phi_vec      = N_Phi
        Phi_X_vec      = Phi_X
        Phi_Y_vec      = Phi_Y
        Ham_T_vec      = Ham_T
        Ham_Tperp_vec  = Ham_Tperp
        Ham_T2_vec     = Ham_T2
        Ham_Chem_vec   = Ham_Chem
        Ham_lambda_vec = 0.d0

        Select case (Lattice_type)

        case ("Honeycomb")
           if ( l_cmplx_trial ) then 

                Ham_T_vec      = 1.d0
                Ham_lambda_vec = 0.d0
                Ham_Chem_vec   = 0.d0
                Phi_X_vec      = 0.01
                Phi_Y_vec      = 0.02
                !Phi_X_vec(1)   = 0.01
                !Phi_X_vec(2)   =-0.01
                !Phi_Y_vec(1)   = 0.02
                !Phi_Y_vec(2)   =-0.02
                call  Set_Default_hopping_parameters_honeycomb(Hopping_Matrix_tmp, Ham_T_vec, Ham_Lambda_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                     &                                         Bulk,  N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )

           else

                Ham_T = 1.d0
                allocate(Op_Tmp(1,N_FL))
                do n = 1,N_FL
                   call Op_make(Op_Tmp(1,n),Ndim)
                enddo
                do I = 1,Latt%N
                   I1 = Invlist(I,1)
                   do nc1 = 1,Latt_unit%N_coord
                      select case (nc1)
                      case (1)
                         J1 = invlist(I,2)
                      case (2)
                         J1 = invlist(Latt%nnlist(I,1,-1),2)
                      case (3)
                         J1 = invlist(Latt%nnlist(I,0,-1),2)
                      case default
                         Write(error_unit,*) 'Error in  Predefined_TrialWaveFunction'
                         error stop 1
                      end select
                      delta = 0.005d0*ranf_wrap()
                      delta = 0.d0
                      do n = 1,N_FL
                          Op_Tmp(1,n)%O(I1,J1) = cmplx(-Ham_T-delta,0.d0, kind(0.D0))
                          Op_Tmp(1,n)%O(J1,I1) = cmplx(-Ham_T-delta,0.d0, kind(0.D0))
                      enddo
                   enddo
                   do n = 1,N_FL
                      rmu = 0.d0
                      rmu = 0.002d0*(ranf_wrap()-0.5d0)
                      Op_Tmp(1,n)%O(I1,I1) = cmplx(rmu,0.d0, kind(0.D0))
                   enddo
                enddo
                do n = 1,N_FL
                   do I = 1,Ndim
                      op_tmp(1,n)%P(i) = i
                   enddo
                   op_Tmp(1,n)%g    = cmplx(1.d0,0.d0,kind(0.d0))
                   op_Tmp(1,n)%alpha= cmplx(0.d0,0.d0,kind(0.D0))
                   call op_set(Op_Tmp(1,n))
                enddo
        
           endif

        Case ("Square")
           if ( hatree_fock ) then
              Ham_T       = 1.d0
              Ham_T2      =-0.3d0
              Ham_U       = 8.d0
              Allocate(Op_Tmp(1,N_FL))
              do n = 1,N_FL
                 Call Op_make(Op_Tmp(1,n),Ndim)
                 Do I = 1,Latt%N
                    I1 = Invlist(I,1)
                    Do nc1 = 1,4
                       select case (nc1)
                       case (1)
                          J1 = invlist(Latt%nnlist(I, 1, 0),1)
                          Op_Tmp(1,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                          Op_Tmp(1,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                       case (2)
                          J1 = invlist(Latt%nnlist(I, 0, 1),1)
                          Op_Tmp(1,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                          Op_Tmp(1,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                       case (3)
                          J1 = invlist(Latt%nnlist(I, 1, 1),1)
                          Op_Tmp(1,n)%O(I1,J1) = cmplx(-Ham_T2,    0.d0, kind(0.D0))
                          Op_Tmp(1,n)%O(J1,I1) = cmplx(-Ham_T2,    0.d0, kind(0.D0))
                       case (4)
                          J1 = invlist(Latt%nnlist(I,-1, 1),1)
                          Op_Tmp(1,n)%O(I1,J1) = cmplx(-Ham_T2,    0.d0, kind(0.D0))
                          Op_Tmp(1,n)%O(J1,I1) = cmplx(-Ham_T2,    0.d0, kind(0.D0))
                       case default
                          Write(error_unit,*) 'Error in  Predefined_TrialWaveFunction'
                          CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
                       end select
                    Enddo
                    sgn_i = 1.d0
                    if (  mod(Latt%List(I,1) + Latt%List(I,2),2) .eq. 0  ) sgn_i = -1.d0
                    sgn_updn = 1.d0
                    if ( n .eq. 2 ) sgn_updn = -1.d0
                    mass = -2.d0*ham_U*0.97*sgn_i*sgn_updn
                    Op_Tmp(1,n)%O(I1,I1) = cmplx(mass,    0.d0, kind(0.D0))
                 enddo
                 do I = 1,Ndim
                    Op_Tmp(1,n)%P(i) = i
                 Enddo
                 Op_Tmp(1,n)%g    = cmplx(1.d0, 0.d0,kind(0.d0))
                 Op_Tmp(1,n)%alpha= cmplx(0.d0,0.d0, kind(0.D0))
                 Call Op_set(Op_Tmp(1,n))
              Enddo
           else
               Ham_T_vec    = 1.d0
               !Ham_T2_vec   = -0.3d0
               Phi_X_vec    = 0.01
               Call  Set_Default_hopping_parameters_square(Hopping_Matrix_tmp,Ham_T_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
                    &                                      Bulk,  N_Phi_vec, N_FL, &
                    &                                      List, Invlist, Latt, Latt_unit )
           endif
        !Case ("N_leg_ladder")
        !   Ham_T_vec     = 1.d0
        !   Ham_Tperp_vec = 1.d0
        !   Phi_X_vec     = 0.01
        !   Call  Set_Default_hopping_parameters_n_leg_ladder(Hopping_Matrix_tmp, Ham_T_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, &
        !        &                                            Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, &
        !        &                                            List, Invlist, Latt, Latt_unit )
        Case ("Bilayer_square")
           Ham_T_vec     = 1.d0
           Ham_T2_vec    = 0.d0
           Ham_Tperp_vec = 1.d0
           Phi_X_vec     = 0.00
           Call  Set_Default_hopping_parameters_Bilayer_square(Hopping_Matrix_tmp,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
                &                                              Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
                &                                              List, Invlist, Latt, Latt_unit )
        Case ("Bilayer_honeycomb")
           Ham_T_vec     = 1.d0
           Ham_T2_vec    = 0.d0
           Ham_Tperp_vec = 1.d0
           Phi_X_vec     = 0.00
           Call  Set_Default_hopping_parameters_Bilayer_honeycomb(Hopping_Matrix_tmp,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, &
                &                                                 Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
                &                                                 List, Invlist, Latt, Latt_unit )

        case ('N_leg_ladder')
             Allocate(op_tmp(1,n_fl))
             do n = 1,n_fl
                call op_make(op_tmp(1,n),ndim)
                do i = 1,ndim
                   op_tmp(1,n)%P(i) = i
                enddo
                op_tmp(1,n)%g    = cmplx(1.d0,0.d0,kind(0.d0))
                op_tmp(1,n)%alpha= cmplx(0.d0,0.d0,kind(0.D0))
             enddo

             do i = 1,latt%n
                do no = 1,Latt_unit%Norb -1
                   delta = 0.005d0*ranf_wrap()
                   ham_t = 1.d0+delta

                   i1 = invlist(i,no)
                   j1 = invlist(latt%nnlist(i,1,0),no)
                   k1 = invlist(i,no+1)
                   
                   ham_tx_vec(1) = ham_t;
                   ham_tx_vec(2) = ham_t*alpha;
                   ham_ty_vec(1) = ham_t*alpha;
                   ham_ty_vec(2) = ham_t;
                   
                   do n = 1, n_fl
                      op_tmp(1,n)%o(i1,j1) = cmplx(-ham_tx_vec(n),0.d0, kind(0.D0))
                      op_tmp(1,n)%o(j1,i1) = cmplx(-ham_tx_vec(n),0.d0, kind(0.D0))
                      op_tmp(1,n)%o(i1,k1) = cmplx(-ham_ty_vec(n),0.d0, kind(0.D0))
                      op_tmp(1,n)%o(k1,i1) = cmplx(-ham_ty_vec(n),0.d0, kind(0.D0))
                   enddo
                enddo
                delta = 0.005d0*ranf_wrap()
                ham_t = 1.d0+delta
                no = Latt_unit%Norb 
                i1 = invlist(i,no)
                j1 = invlist(latt%nnlist(i,1,0),no)
                ham_tx_vec(1) = ham_t;
                ham_tx_vec(2) = ham_t*alpha;
                do n = 1, n_fl
                   op_tmp(1,n)%o(i1,j1) = cmplx(-ham_tx_vec(n),0.d0, kind(0.D0))
                   op_tmp(1,n)%o(j1,i1) = cmplx(-ham_tx_vec(n),0.d0, kind(0.D0))
                enddo
             enddo
             
             do n = 1, n_fl
                Call op_set(op_tmp(1,n))
             enddo

        case default
           Write(error_unit,*) 'No predefined trial wave function for this lattice.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        end Select


        If ( l_cmplx_trial .and. ( .not. hatree_fock ) .and. ( lattice_type .ne. "N_leg_ladder") )   &
             &     call  Predefined_Hoppings_set_OPT(Hopping_Matrix_tmp,List,Invlist,Latt,  Latt_unit,  Dtau, Checkerboard, Symm, OP_tmp )

        Do nf = 1,N_FL
           Call Diag(Op_tmp(1,nf)%O,Op_tmp(1,nf)%U,Op_tmp(1,nf)%E)
           do I2=1,N_part
              do I1=1,Ndim
                 WF_L(nf,1)%P(I1,I2)=Op_tmp(1,nf)%U(I1,I2)
                 WF_R(nf,1)%P(I1,I2)=Op_tmp(1,nf)%U(I1,I2)
              enddo
           enddo
           WF_L(nf,1)%Degen = Op_tmp(1,nf)%E(N_part+1) - Op_tmp(1,nf)%E(N_part)
           WF_R(nf,1)%Degen = Op_tmp(1,nf)%E(N_part+1) - Op_tmp(1,nf)%E(N_part)
        enddo

        Do nf = 1,N_FL
           Call WF_overlap(WF_L(nf,1), WF_R(nf,1), Z_norm)
           WF_L(nf,1)%P(:,:) = WF_L(nf,1)%P(:,:)/sqrt(dble(N_slat)) 
           WF_R(nf,1)%P(:,:) = WF_R(nf,1)%P(:,:)/sqrt(dble(N_slat))  
        enddo


        Do nf = 1, N_FL
        do ns = 2, n_slat
           WF_L(nf,ns)%Degen  = WF_L(nf,1)%Degen
           WF_R(nf,ns)%Degen  = WF_R(nf,1)%Degen
           WF_L(nf,ns)%P(:,:) = WF_L(nf,1)%P(:,:)
           WF_R(nf,ns)%P(:,:) = WF_R(nf,1)%P(:,:)
        enddo
        enddo

        Do nf = 1,N_FL
           Call Op_clear(OP_tmp(1,nf),Ndim)
        enddo
        Deallocate (OP_tmp)
        Call Predefined_hoppings_clear(Hopping_Matrix_tmp)

        Deallocate (Ham_T_vec, Ham_Tperp_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec,  N_Phi_vec )
        deallocate (ham_tx_vec, ham_ty_vec)

      end Subroutine Predefined_TrialWaveFunction

     end Module Predefined_Trial
