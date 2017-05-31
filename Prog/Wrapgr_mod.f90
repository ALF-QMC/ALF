!  Copyright (C) 2016, 2017 The ALF project
! 
!  This file is part of the ALF project.
! 
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
!     
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!     
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

Module Wrapgr_mod

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Handles global updates on a single time slice
!
!--------------------------------------------------------------------

  Use Hamiltonian
  Use MyMats 
  Use Operator_mod
  Use Control
  Use Random_Wrap

  Implicit none

  
  !> Privat 
  Complex (Kind=Kind(0.d0)),  private, allocatable ::  GR_ST(:,:,:)
  
Contains
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Allocate, Deallocate space
!--------------------------------------------------------------------
  Subroutine Wrapgr_alloc
    Implicit none
    Allocate (GR_ST(Ndim,Ndim,N_FL) )
  end Subroutine Wrapgr_alloc
  
  Subroutine  Wrapgr_dealloc
    Implicit none
    deallocate ( GR_ST )
  end Subroutine Wrapgr_dealloc

!--------------------------------------------------------------------
  SUBROUTINE WRAPGRUP(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Given the green function matrix GR at time NTAU  the routine   propagates 
!> it to time  NTAU + 1 and carries  out an update of the fields at time NTAU+1
!> NTAU: [0:LTROT-1]
!
!--------------------------------------------------------------------
 
    Use Hamiltonian
    Use Hop_mod
    Implicit none
    
    Interface
       Subroutine Upgrade2(GR,N_op,NT,PHASE,Op_dim,ns_new, Prev_Ratiotot, S0_ratio, T0_proposal_ratio, Acc,mode) 
         Use Hamiltonian
         Use Random_wrap
         Use Control
         Implicit none  
         Complex (Kind=Kind(0.d0)) :: GR(Ndim,Ndim, N_FL)
         Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot
         Integer, INTENT(IN)       :: N_op, Nt, Op_dim
         Complex (Kind=Kind(0.d0)) :: Phase
         Integer                   :: ns_new
         Real    (Kind=Kind(0.d0)) :: S0_ratio, T0_Proposal_ratio
         Character (Len=64)        :: Mode
         Logical                   :: Acc
       end Subroutine Upgrade2
    end Interface
    
    ! Arguments
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), allocatable ::  GR(:,:,:)
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) ::  PHASE
    INTEGER, INTENT(IN) :: NTAU
    LOGICAL, INTENT(IN) :: Propose_S0
    INTEGER, INTENT(IN) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau

    !Local 
    Integer :: nf, N_Type, NTAU1,n, m, Flip_value(1)
    Complex (Kind=Kind(0.d0)) :: Mat_TMP(Ndim,Ndim),  Prev_Ratiotot
    Real    (Kind=Kind(0.d0)) :: X, T0_proposal,  T0_Proposal_ratio,  S0_ratio
    Character (Len=64)        :: Mode
    Logical                   :: Acc
    
    ! Wrap up, upgrade ntau1.  with B^{1}(tau1) 
    NTAU1 = NTAU + 1
    Do nf = 1,N_FL
       CALL HOP_MOD_mmthr( GR(:,:,nf),  MAT_TMP,nf)
       CALL HOP_MOD_mmthl_m1(MAT_TMP,GR(:,:,nf), nf )
       !CALL MMULT ( MAT_TMP,    Exp_T(:,:,nf), GR(:,:,nf)        )
       !CALL MMULT ( GR(:,:,nf), MAT_TMP      , Exp_T_M1(:,:,nf)  )
    Enddo
    Do n = Nt_sequential_start,Nt_sequential_end
       ! Write(6,*) 'Hi'
       Do nf = 1, N_FL
          X = Phi(nsigma(n,ntau1),Op_V(n,nf)%type)
          N_type = 1
          Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),X,Ndim,N_Type)
       enddo
       nf = 1
       T0_proposal       = 1.5D0
       T0_Proposal_ratio = 1.D0
       S0_ratio          = S0(n,ntau1)
       if ( Propose_S0 ) then
          If ( Op_V(n,nf)%type == 1)  then
             T0_proposal       = 1.d0 - 1.d0/(1.d0+S0_ratio)
             T0_Proposal_ratio = 1.d0/S0_ratio
          endif
       Endif
       If ( T0_proposal > ranf_wrap() ) Then
          if  (Op_V(n,nf)%type == 1 ) Flip_value(1)  = -nsigma(n,ntau1)
          if  (Op_V(n,nf)%type == 2 ) Flip_value(1)  =  NFLIPL(nsigma(n,ntau1),nranf(3))
          mode = "Final"
          Prev_Ratiotot = cmplx(1.d0,0.d0,kind(0.d0))
          Call Upgrade2(GR,n,ntau1,PHASE,Op_V(n,nf)%N_non_Zero,Flip_value(1), Prev_Ratiotot, S0_ratio,T0_Proposal_ratio, Acc, mode ) 
       Endif
       do nf = 1,N_FL
          N_type =  2
          Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),X,Ndim,N_Type)
       enddo
    Enddo

    If ( N_Global_tau > 0 ) then 
       m         = Nt_sequential_end
       !if ( Nt_sequential_start >  Nt_sequential_end ) m = Nt_sequential_start 
       Call Wrapgr_Random_update(GR,m,ntau1, PHASE, N_Global_tau )
       Call Wrapgr_PlaceGR(GR,m, Size(OP_V,1), ntau1)
    Endif

  END SUBROUTINE WRAPGRUP


!--------------------------------------------------------------------    
  SUBROUTINE WRAPGRDO(GR,NTAU,PHASE,Propose_S0,Nt_sequential_start, Nt_sequential_end, N_Global_tau)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Given the green function matrix GR at time NTAU  the routine   carries out an 
!> update of the fields at time NTAU and  propagates  it to time NTAU-1
!> NTAU: [LTROT:1]
!
!--------------------------------------------------------------------    
    Use Hamiltonian
    Use MyMats
    Use Hop_mod
    Implicit None
    
    Interface
       Subroutine Upgrade2(GR,N_op,NT,PHASE,Op_dim,ns_new, Prev_Ratiotot, S0_ratio, T0_proposal_ratio, Acc, mode) 
         Use Hamiltonian
         Use Random_wrap
         Use Control
         Implicit none  
         Complex (Kind=Kind(0.d0)) :: GR(Ndim,Ndim, N_FL)
         Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot
         Integer, INTENT(IN)       :: N_op, Nt, Op_dim
         Complex (Kind=Kind(0.d0)) :: Phase
         Integer                   :: ns_new
         Real    (Kind=Kind(0.d0)) :: S0_ratio, T0_Proposal_ratio
         Character (Len=64)        :: Mode
         Logical                   :: Acc
       end Subroutine Upgrade2
    End Interface
    
    ! Given GREEN at time NTAU => GREEN at time NTAU - 1,
    ! Upgrade NTAU  [LTROT:1]
    
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), allocatable :: GR(:,:,:)
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: PHASE
    Integer, INTENT(IN) :: NTAU
    LOGICAL, INTENT(IN) :: Propose_S0
    INTEGER, INTENT(IN) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau
    
    ! Local
    Integer :: nf, N_Type, n, m, Flip_value(1)
    Complex (Kind=Kind(0.d0)) :: Mat_TMP(Ndim,Ndim),  Prev_Ratiotot
    Real    (Kind=Kind(0.d0)) :: spin, T0_proposal,  T0_Proposal_ratio,  S0_ratio
    Character (Len=64)        :: Mode
    Logical                   :: Acc


    If ( N_Global_tau > 0 ) then 
       m         = Size(OP_V,1)
       !Write(6,*) 'Call Ran_up ', m,ntau
       Call Wrapgr_Random_update(GR,m,ntau, PHASE, N_Global_tau )
       Call Wrapgr_PlaceGR(GR,m, Nt_sequential_end, ntau)
    Endif

    
    Do n =  Nt_sequential_end, Nt_sequential_start, -1
       N_type = 2
       nf = 1
       spin = Phi(nsigma(n,ntau),Op_V(n,nf)%type)
       do nf = 1,N_FL
          Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), spin, Ndim, N_Type)
       enddo
       !Write(6,*) 'Upgrade : ', ntau,n 
       nf = 1
       T0_proposal       = 1.5D0
       T0_Proposal_ratio = 1.D0
       S0_ratio          = S0(n,ntau)
       if ( Propose_S0 ) then
          If ( Op_V(n,nf)%type == 1)  then
             T0_proposal       = 1.d0 - 1.d0/(1.d0+S0_ratio)
             T0_Proposal_ratio = 1.d0/S0_ratio
          endif
       Endif
       If ( T0_proposal > ranf_wrap() ) Then
          if  (Op_V(n,nf)%type == 1 ) Flip_value(1)  = -nsigma(n,ntau)
          if  (Op_V(n,nf)%type == 2 ) Flip_value(1)  =  NFLIPL(nsigma(n,ntau),nranf(3))
          mode = "Final"
          Prev_Ratiotot = cmplx(1.d0,0.d0,kind(0.d0))
          Call Upgrade2(GR,n,ntau,PHASE,Op_V(n,nf)%N_non_Zero,Flip_value(1), Prev_Ratiotot, S0_ratio,T0_Proposal_ratio, Acc, mode ) 
       Endif
     
       !Call Upgrade(GR,n,ntau,PHASE,Op_V(n,1)%N_non_zero) 
       ! The spin has changed after the upgrade!
       nf = 1
       spin = Phi(nsigma(n,ntau),Op_V(n,nf)%type)
       N_type = 1
       do nf = 1,N_FL
          Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), spin, Ndim, N_Type )
       enddo
    enddo
    DO nf = 1,N_FL
       Call Hop_mod_mmthl   (GR(:,:,nf), MAT_TMP, nf)
       Call Hop_mod_mmthr_m1(MAT_TMP, GR(:,:,nf), nf)
       !CALL MMULT(MAT_TMP   , GR(:,:,nf)      , Exp_T(:,:,nf) )
       !CALL MMULT(GR(:,:,nf), Exp_T_M1(:,:,nf), MAT_TMP       )
    enddo
    
  end SUBROUTINE WRAPGRDO
  

!--------------------------------------------------------------------
  Subroutine  Wrapgr_PlaceGR(GR,m,m1,ntau)
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> The Green function on a given time slice reads
!> G(tau) = (  1 + B(tau) B(tau-1,0)  B(beta,tau))^(-1) with B(tau) =   U_n e^(d_n) U_n^(dag) .... U_1 e^(V_1) U_1^(dag)  e^(-dtau H_t)
!> On input you have 
!> G(tau,m)  = [ 1 + U_m e^(d_m) U_m^(dag) U_m^(dag) ... U_1 e^(V_1) U_1^(dag) e^(-dtau H_t) B(tau-1,0) 
!>                   B(Beta,tau)  U_n e^(d_n) U_n^(dag) ...U_(m+1) e^(d_(m+1)) U_(m+1)^(dag) U_(m+1) ] 
!> On output you have 
!> G(tau,  m1 )  
!> 
!> Note that m,m1 in [0,n]
!--------------------------------------------------------------------
    
    Implicit none

    !Arguments
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: GR
    Integer, INTENT(IN) :: m,m1, ntau

    !Local 
    Integer :: n, nf, N_Type
    Real (Kind=Kind(0.d0)) :: X

    If (m == m1)  then 
       return
    elseif  ( m1 > m  ) then
       !Write(6,*) "Wrapup from ",  m + 1, "to",  m1, " on tau=",  ntau
       Do n = m+1,m1
          Do nf = 1, N_FL
             X = Phi(nsigma(n,ntau),Op_V(n,nf)%type)
             N_type = 1
             Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),X,Ndim,N_Type)
          enddo
          do nf = 1,N_FL
             N_type =  2
             Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),X,Ndim,N_Type)
          enddo
       Enddo
    elseif  (m1 < m ) then
       !Write(6,*) "Wrapdo from ",  m, "to",  m1 + 1 
       Do n =  m, m1+1 ,-1 
          N_type = 2
          nf = 1
          X = Phi(nsigma(n,ntau),Op_V(n,nf)%type)
          do nf = 1,N_FL
             Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), X, Ndim, N_Type)
          enddo
          !Write(6,*) 'Upgrade : ', ntau,n 
          nf = 1
          X = Phi(nsigma(n,ntau),Op_V(n,nf)%type)
          N_type = 1
          do nf = 1,N_FL
             Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), X, Ndim, N_Type )
          enddo
       enddo
    endif
    
  end Subroutine Wrapgr_PlaceGR



!--------------------------------------------------------------------
  Subroutine  Wrapgr_Random_update( GR,m,ntau, PHASE, N_Global_tau )
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> On input: 
!> GR(tau,m) as defined in  Global_tau_mod_PlaceGR and the direction of updating scheme
!> direction=u --> You are visiting the time slices from tau = 1  to tau =Ltrot
!> direction=d --> You are visiting the time slices from tau = Ltrot   to tau = 1
!> 
!> The routine calls  
!> Global_move_tau(T0_Proposal_ratio, Flip_list, Flip_length,ntau,m,direction)
!> in the Hamiltonian module and then carries out the update  
!> 
!> On output
!> 
!> Flip_length==1  
!>        Green function is on  GR(tau,Flip_list(1) +1 )  if direction = u 
!>        Green function is on  GR(tau,Flip_list(1) -1 )  if direction = d
!>        This is valid if the move has or has not been accepted. 
!>
!> Flip_length > 1 
!>        Let m_min = min(Flip_list), m_max = max(Flip_list)
!>        direction = u -->  On output Green on m_max is accepted. Green is on m_min if not accepted. 
!>        direction = d -->  On output Green on m_min if accepted. Green is on m_max if not accepted.
!--------------------------------------------------------------------
        
    Implicit none
    
    Interface 
       Subroutine Upgrade2(GR,N_op,NT,PHASE,Op_dim,ns_new, Prev_Ratiotot,S0_ratio, T0_Proposal_ratio, Acc, Mode)
         Use Hamiltonian
         Use Random_wrap
         Use Control
         Implicit none  
         Complex (Kind=Kind(0.d0)) :: GR(Ndim,Ndim, N_FL)
         Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot
         Integer, INTENT(IN)       :: N_op, Nt, Op_dim
         Complex (Kind=Kind(0.d0)) :: Phase
         Integer                   :: ns_new
         Real    (Kind=Kind(0.d0)) :: S0_ratio, T0_Proposal_ratio
         Logical, INTENT(INOUT)    :: Acc
         Character (Len=64)        :: Mode
       end Subroutine Upgrade2
    end Interface

    ! Arguments 
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: GR
    Integer,           INTENT(INOUT) :: m
    Integer,           INTENT(IN)    :: ntau, N_Global_tau 
    Complex  (Kind=Kind(0.d0)), INTENT(INOUT) :: PHASE
    


    ! Space for local varaibles
    Integer                   :: n, Flip_length, nf, N_Type, ng_c, Flip_count
    Real    (Kind=Kind(0.d0)) :: T0_Proposal_ratio, T0_proposal, X,S0_ratio
    COMPLEX (Kind=Kind(0.d0)) :: Prev_Ratiotot 
    Logical                   :: Acc
    Character (Len=64)        :: Mode
    Integer,      allocatable :: Flip_list(:), Flip_value(:), Flip_value_st(:)
    Real    (Kind=Kind(0.d0)) :: Zero = 10D-8

    Allocate ( Flip_list(Size(Op_V,1)), Flip_value(Size(Op_V,1)), Flip_value_st(Size(Op_V,1)) )

    Do ng_c = 1,N_Global_tau
       ! New configuration
       Call Global_move_tau(T0_Proposal_ratio, S0_ratio,  Flip_list, Flip_length,Flip_value,ntau )
       !Write(6,*)  "Calling global move",  m, Flip_list(1), nsigma(Flip_list(1),ntau),Flip_value(1)
       If ( T0_Proposal_ratio  >  Zero )  Then
          ! Order the list
          Call wrapgr_sort(Flip_length,Flip_list,Flip_value)
          If ( Flip_length > 1 ) then
             Do Flip_count = 1, Flip_length-1 
                Flip_value_st(Flip_count)  = nsigma( Flip_list(Flip_count), ntau  )
             Enddo
          endif
          Prev_Ratiotot = cmplx(1.d0,0.d0,kind(0.d0))
          !Write(6,*) "-----", Flip_length
          do Flip_count = 1,Flip_length
             n = Flip_list(Flip_count)
             !Write(6,*)  "PlaceGR",  m, n-1,ntau
             Call Wrapgr_PlaceGR(GR,m, n-1, ntau)
             !Write(6,*)  "Back from PlaceGR",  m, n-1,ntau
             If ( Flip_count == 1 .and. Flip_length > 1 ) GR_st = Gr
             Do nf = 1, N_FL
                X = Phi(nsigma(n,ntau),Op_V(n,nf)%type)
                N_type = 1
                Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),X,Ndim,N_Type)
             enddo
             nf = 1
             If (Flip_count <  Flip_length)  then 
                mode = "Intermediate"
                Call Upgrade2(GR,n,ntau,PHASE,Op_V(n,nf)%N_non_Zero,Flip_value(Flip_count), &
                     &        Prev_Ratiotot, S0_ratio, T0_Proposal_ratio, Acc, mode ) 
                do nf = 1,N_FL
                   N_type =  2
                   Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),X,Ndim,N_Type)
                enddo
             else
                !Write(6,*)  "Call Up mode final", n,ntau
                mode = "Final"
                Call Upgrade2(GR,n,ntau,PHASE,Op_V(n,nf)%N_non_Zero,Flip_value(Flip_count), &
                     &        Prev_Ratiotot, S0_ratio, T0_Proposal_ratio, Acc, mode ) 
                !Write(6,*)  "Back from up mode final", n,ntau
                !Write(6,*)  "Acceptance", Acc
                do nf = 1,N_FL
                   N_type =  2
                   Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),X,Ndim,N_Type)
                enddo
             endif
             m = n
          enddo
          If ( .not. Acc .and. Flip_length > 1 ) then
             Gr = Gr_st
             m = Flip_list(1) - 1
             Do Flip_count = 1, Flip_length-1 
                nsigma( Flip_list(Flip_count), ntau  ) = Flip_value_st(Flip_count)  
             Enddo
          Endif
          !If (Acc) Call Hamiltonian_Print(Ntau)
       endif
    Enddo
    
    Deallocate ( Flip_list, Flip_value, Flip_value_st )
    

  end Subroutine Wrapgr_Random_update

!----------------------------------------------------------------------------
  subroutine Wrapgr_sort(Flip_length,Flip_list,Flip_value)

    ! Arguments
    Integer, INTENT(IN) :: Flip_length
    Integer, INTENT(INOUT), allocatable :: Flip_list(:), Flip_value(:)
    
    ! Local
    integer :: swaps            ! number of swaps made in one pass
    integer :: nc               ! loop variable
    integer :: temp             ! temporary holder for making swap
    if ( Flip_length == 1 ) return 
    
    !Write(6,*) 'Before sort'
    !DO nc = 1,Flip_length
    !   Write(6,*) Flip_list(nc),  Flip_value(nc)
    !Enddo

    do 
       swaps      = 0           ! Initially, we've made no swaps
       do nc = 1, (Flip_length - 1)
          if ( Flip_list(nc)  >  Flip_list(nc+1) ) then
             temp              = Flip_list(nc  ) 
             Flip_list(nc)     = Flip_list(nc+1) 
             Flip_list(nc+1)   = temp
             temp              = Flip_value(nc  ) 
             Flip_value(nc)    = Flip_value(nc+1) 
             Flip_value(nc+1)  = temp
             swaps             = swaps + 1
          end if
       end do
       if ( swaps == 0 ) exit ! do count swaps
    end do

    !Write(6,*) 'After sort'
    !DO nc = 1,Flip_length
    !   Write(6,*) Flip_list(nc),  Flip_value(nc)
    !Enddo

  end subroutine Wrapgr_sort
  
!----------------------------------------------------------------------------
  
  Subroutine Wrapgr_Test(Gr,ntau)
    
    ! Arguments 
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: GR
    Integer :: ntau
    
    ! Local 
    Integer :: m, m1, N_size, nf, nth
    Real (Kind=kind(0.d0)) :: Xmean, Xmax

    !Input is the Green function on time slice tau 
    N_size = size(OP_V,1)
    GR_ST = GR
    m = N_Size
    m1 = 0
    DO nth = 1,10
       call Wrapgr_PlaceGR(GR,m,m1,ntau)        
       Write(6,*) m, m1
       m = m1
       m1 =  nranf(N_Size)-1
    enddo
    call Wrapgr_PlaceGR(GR,m,N_Size,ntau)        
    Write(6,*) m, N_size
    
    Xmax = 0.d0
    Do nf = 1,N_FL
       Call COMPARE(GR_st(:,:,nf),GR(:,:,nf),XMAX,XMEAN)
    Enddo
    Write(6,*)  'Compare Global_tau_mod_Test ', Xmax
    Deallocate ( GR_ST )
    
    
  End Subroutine Wrapgr_Test

end Module Wrapgr_mod
