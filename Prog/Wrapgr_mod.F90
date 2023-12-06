!  Copyright (C) 2016 - 2022 The ALF project
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
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
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

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Handles global updates on a single time slice
!
!--------------------------------------------------------------------

Module Wrapgr_mod


  Use Hamiltonian_main
  Use MyMats 
  Use Operator_mod
  Use Control
  Use Random_Wrap
  Use Fields_mod
  Use Hamiltonian_main
  Use Hop_mod
  use upgrade_mod

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
    Implicit none
    
    ! Arguments
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), allocatable ::  GR(:,:,:)
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) ::  PHASE
    INTEGER, INTENT(IN) :: NTAU
    LOGICAL, INTENT(IN) :: Propose_S0
    INTEGER, INTENT(IN) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau

    !Local 
    Integer :: nf, nf_eff, N_Type, NTAU1,n, m
    Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot, HS_Field, HS_New
    Real    (Kind=Kind(0.d0)) :: T0_proposal,  T0_Proposal_ratio,  S0_ratio
    Character (Len=64)        :: Mode
    Logical                   :: Acc, toggle1
    
    ! Wrap up, upgrade ntau1.  with B^{1}(tau1) 
    NTAU1 = NTAU + 1
    Do nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       CALL HOP_MOD_mmthr   (GR(:,:,nf), nf, ntau1 )
       CALL HOP_MOD_mmthl_m1(GR(:,:,nf), nf, ntau1 )
    Enddo
    Do n = Nt_sequential_start,Nt_sequential_end
       Do nf_eff = 1, N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          HS_Field =  nsigma%f(n,ntau1) 
          N_type = 1
          Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau1)
       enddo
       nf = 1
       T0_proposal       = 1.5D0
       T0_Proposal_ratio = 1.D0
       Hs_New =   nsigma%flip(n,ntau1) 
       S0_ratio          = ham%S0(n,ntau1, Hs_New )
       if ( Propose_S0 ) then
          If ( Op_V(n,nf)%type == 1)  then
             T0_proposal       = 1.d0 - 1.d0/(1.d0+S0_ratio)
             T0_Proposal_ratio = 1.d0/S0_ratio
          endif
       Endif
       If ( T0_proposal > ranf_wrap() ) Then
          !Write(6,*) 'Hi', n, Op_V(n,nf)%type, T0_Proposal_ratio, S0_ratio  
          mode = "Final"
          Prev_Ratiotot = cmplx(1.d0,0.d0,kind(0.d0))
          Call Upgrade2(GR,n,ntau1,PHASE,HS_new, Prev_Ratiotot, S0_ratio,T0_Proposal_ratio, Acc, mode ) 
       else
          toggle1 = .false.
          Call Control_upgrade_eff(toggle1)
       Endif
       do nf_eff = 1,N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          N_type =  2
          Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau1)
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
    Implicit None
    
    ! Given GREEN at time NTAU => GREEN at time NTAU - 1,
    ! Upgrade NTAU  [LTROT:1]
    
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), allocatable :: GR(:,:,:)
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT) :: PHASE
    Integer, INTENT(IN) :: NTAU
    LOGICAL, INTENT(IN) :: Propose_S0
    INTEGER, INTENT(IN) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau
    
    ! Local
    Integer :: nf, nf_eff, N_Type, n, m
    Complex (Kind=Kind(0.d0)) :: Prev_Ratiotot, HS_Field, HS_New
    Real    (Kind=Kind(0.d0)) :: T0_proposal,  T0_Proposal_ratio,  S0_ratio
    Character (Len=64)        :: Mode
    Logical                   :: Acc, toggle1


    If ( N_Global_tau > 0 ) then 
       m         = Size(OP_V,1)
       !Write(6,*) 'Call Ran_up ', m,ntau
       Call Wrapgr_Random_update(GR,m,ntau, PHASE, N_Global_tau )
       Call Wrapgr_PlaceGR(GR,m, Nt_sequential_end, ntau)
    Endif

    
    Do n =  Nt_sequential_end, Nt_sequential_start, -1
       N_type = 2
       nf = 1
       HS_Field = nsigma%f(n,ntau) 
       do nf_eff = 1,N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), HS_Field, Ndim, N_Type,ntau)
       enddo
       !Write(6,*) 'Upgrade : ', ntau,n 
       nf = 1
       T0_proposal       = 1.5D0
       T0_Proposal_ratio = 1.D0
       HS_new            = nsigma%flip(n,ntau) 
       S0_ratio          = ham%S0(n,ntau,HS_new)
       if ( Propose_S0 ) then
          If ( Op_V(n,nf)%type == 1)  then
             T0_proposal       = 1.d0 - 1.d0/(1.d0+S0_ratio)
             T0_Proposal_ratio = 1.d0/S0_ratio
          endif
       Endif
       If ( T0_proposal > ranf_wrap() ) Then
          mode = "Final"
          Prev_Ratiotot = cmplx(1.d0,0.d0,kind(0.d0))
          Call Upgrade2(GR,n,ntau,PHASE,HS_new, Prev_Ratiotot, S0_ratio,T0_Proposal_ratio, Acc, mode ) 
       else
          toggle1 = .false.
          Call Control_upgrade_eff(toggle1)
       Endif
     
       !Call Upgrade(GR,n,ntau,PHASE,Op_V(n,1)%N_non_zero) 
       ! The spin has changed after the upgrade!
       nf = 1
       HS_Field = nsigma%f(n,ntau)  
       N_type = 1
       do nf_eff = 1,N_FL_eff
          nf=Calc_Fl_map(nf_eff)
          Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), HS_Field, Ndim, N_Type, ntau )
       enddo
    enddo
    DO nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       Call Hop_mod_mmthl   (GR(:,:,nf), nf,ntau)
       Call Hop_mod_mmthr_m1(GR(:,:,nf), nf,ntau)
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
    Integer :: n, nf, nf_eff, N_Type 
    Complex (Kind=Kind(0.d0)) :: HS_Field

    If (m == m1)  then 
       return
    elseif  ( m1 > m  ) then
       !Write(6,*) "Wrapup from ",  m + 1, "to",  m1, " on tau=",  ntau
       Do n = m+1,m1
          Do nf_eff = 1, N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             HS_Field = nsigma%f(n,ntau)
             N_type = 1
             Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
          enddo
          do nf_eff = 1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             N_type =  2
             Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
          enddo
       Enddo
    elseif  (m1 < m ) then
       !Write(6,*) "Wrapdo from ",  m, "to",  m1 + 1 
       Do n =  m, m1+1 ,-1 
          N_type = 2
          nf = 1
          HS_Field = nsigma%f(n,ntau)
          do nf_eff = 1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), HS_Field, Ndim, N_Type,ntau)
          enddo
          !Write(6,*) 'Upgrade : ', ntau,n 
          nf = 1
          HS_Field= nsigma%f(n,ntau)
          N_type = 1
          do nf_eff = 1,N_FL_eff
             nf=Calc_Fl_map(nf_eff)
             Call Op_Wrapdo( Gr(:,:,nf), Op_V(n,nf), HS_Field, Ndim, N_Type, ntau )
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

    ! Arguments 
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: GR
    Integer,           INTENT(INOUT) :: m
    Integer,           INTENT(IN)    :: ntau, N_Global_tau 
    Complex  (Kind=Kind(0.d0)), INTENT(INOUT) :: PHASE
    


    ! Space for local variables
    Integer                   :: n, Flip_length, nf, nf_eff, N_Type, ng_c, Flip_count
    Real    (Kind=Kind(0.d0)) :: T0_Proposal_ratio, T0_proposal,S0_ratio
    COMPLEX (Kind=Kind(0.d0)) :: Prev_Ratiotot, HS_Field, HS_New 
    Logical                   :: Acc
    Character (Len=64)        :: Mode
    Integer,      allocatable :: Flip_list(:)
    Complex (Kind=Kind(0.d0)), allocatable :: Flip_value(:), Flip_value_st(:)
    Real    (Kind=Kind(0.d0)) :: Zero = 10D-8

    Allocate ( Flip_list(Size(Op_V,1)), Flip_value(Size(Op_V,1)), Flip_value_st(Size(Op_V,1)) )

    Do ng_c = 1,N_Global_tau
       ! New configuration
       Call ham%Global_move_tau(T0_Proposal_ratio, S0_ratio,  Flip_list, Flip_length,Flip_value,ntau )
       !Write(6,*)  "Calling global move",  m, Flip_list(1), nsigma(Flip_list(1),ntau),Flip_value(1)
       If ( T0_Proposal_ratio  >  Zero )  Then
          ! Order the list
          Call wrapgr_sort(Flip_length,Flip_list,Flip_value)
          If ( Flip_length > 1 ) then
             Do Flip_count = 1, Flip_length-1 
                Flip_value_st(Flip_count)  = nsigma%f( Flip_list(Flip_count), ntau  )
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
             Do nf_eff = 1, N_FL_eff
                nf=Calc_Fl_map(nf_eff)
                HS_Field = nsigma%f(n,ntau)
                N_type = 1
                Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
             enddo
             nf = 1
             If (Flip_count <  Flip_length)  then 
                mode = "Intermediate"
                HS_new = Flip_value(Flip_count)
                Call Upgrade2(GR,n,ntau,PHASE, HS_new , &
                     &        Prev_Ratiotot, S0_ratio, T0_Proposal_ratio, Acc, mode ) 
                do nf_eff = 1,N_FL_eff
                   nf=Calc_Fl_map(nf_eff)
                   N_type =  2
                   Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
                enddo
             else
                !Write(6,*)  "Call Up mode final", n,ntau
                mode = "Final"
                HS_new = Flip_value(Flip_count)
                Call Upgrade2(GR,n,ntau,PHASE,HS_new, &
                     &        Prev_Ratiotot, S0_ratio, T0_Proposal_ratio, Acc, mode ) 
                !Write(6,*)  "Back from up mode final", n,ntau
                !Write(6,*)  "Acceptance", Acc
                do nf_eff = 1,N_FL_eff
                   nf=Calc_Fl_map(nf_eff)
                   N_type =  2
                   Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),HS_Field,Ndim,N_Type,ntau)
                enddo
             endif
             m = n
          enddo
          If ( .not. Acc .and. Flip_length > 1 ) then
             Gr = Gr_st
             m = Flip_list(1) - 1
             Do Flip_count = 1, Flip_length-1 
                nsigma%f( Flip_list(Flip_count), ntau  ) = Flip_value_st(Flip_count)
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
    Integer, INTENT(INOUT), allocatable :: Flip_list(:)
    Complex   (Kind=Kind(0.d0)), INTENT(INOUT), allocatable :: Flip_value(:)
    
    ! Local
    integer :: swaps            ! number of swaps made in one pass
    integer :: nc               ! loop variable
    integer :: temp, n          ! temporary holder for making swap
    Complex (Kind=Kind(0.d0))      :: X

    
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
             X                   = Flip_value(nc   ) 
             Flip_value(nc)    = Flip_value(nc+1 ) 
             Flip_value(nc+1)  = X
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
    Integer :: m, m1, N_size, nf, nf_eff, nth
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
    Do nf_eff = 1,N_FL_eff
       nf=Calc_Fl_map(nf_eff)
       Call COMPARE(GR_st(:,:,nf),GR(:,:,nf),XMAX,XMEAN)
    Enddo
    Write(6,*)  'Compare Global_tau_mod_Test ', Xmax
    Deallocate ( GR_ST )
    
    
  End Subroutine Wrapgr_Test

end Module Wrapgr_mod
