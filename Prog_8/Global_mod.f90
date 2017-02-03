!  Copyright (C) 2016 The ALF project
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

Module Global_mod

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Handles global updates.
!> 
!
!--------------------------------------------------------------------

  Use Hamiltonian
  Use MyMats 
  Use Operator_mod
  Use Control

  Implicit none

  
Contains

  Subroutine Global_Updates(Phase,GR,UR,DR,VR, UL,DL,VL,Stab_nt, UST, VST, DST)

    Implicit none

    Interface
       SUBROUTINE WRAPUL(NTAU1, NTAU, UL ,DL, VL)
         Use Hamiltonian
         Implicit none
         COMPLEX (Kind=Kind(0.d0)) :: UL(Ndim,Ndim,N_FL), VL(Ndim,Ndim,N_FL)
         COMPLEX (Kind=Kind(0.d0)) :: DL(Ndim,N_FL)
         Integer :: NTAU1, NTAU
       END SUBROUTINE WRAPUL
       SUBROUTINE CGR(PHASE,NVAR, GRUP, URUP,DRUP,VRUP, ULUP,DLUP,VLUP)
         Use UDV_Wrap_mod
         Implicit None
         COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(In)    :: URUP, VRUP, ULUP, VLUP
         COMPLEX(Kind=Kind(0.d0)), Dimension(:)  , Intent(In)    :: DLUP, DRUP
         COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(Inout) :: GRUP
         COMPLEX(Kind=Kind(0.d0)) :: PHASE
         INTEGER         :: NVAR
       END SUBROUTINE CGR
    end Interface
    
!>  Arguments
    COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT)                   :: Phase
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:)  , INTENT(INOUT), allocatable :: DL, DR
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: UL, VL, UR, VR
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(INOUT), allocatable :: GR
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:)  , INTENT(INOUT), allocatable :: DST
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:,:), INTENT(INOUT), allocatable :: UST,  VST
    INTEGER, dimension(:),     INTENT   (IN), allocatable      :: Stab_nt
    
!>  Local variables.
    Integer :: NST, NSTM, NF, NT, NT1, NVAR,N, N1,N2, I, NC
    Integer, Dimension(:,:),  allocatable :: nsigma_old
    Real    (Kind=Kind(0.d0)) :: T0_Proposal_ratio, Weight
    Complex (Kind=Kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Z, Ratiotot, Phase_old, Phase_new
    Complex (Kind=Kind(0.d0)), allocatable :: Det_vec_old(:,:), Det_vec_new(:,:), Phase_Det_new(:), Phase_Det_old(:)
    Logical :: TOGGLE, L_Test

 

!>  On entry and on exit the left storage is full, and the Green function is on time slice 0 and the phase is set.

    n1 = size(nsigma,1)
    n2 = size(nsigma,2)
    NSTM = Size(UST,3)
    Allocate ( nsigma_old(n1,n2) )
    Allocate ( Det_vec_old(NDIM,N_FL), Det_vec_new(NDIM,N_FL) ) 
    Allocate ( Phase_Det_new(N_FL), Phase_Det_old(N_FL) )

    
    L_test = .false.
    ! Set old weight. 
    Phase_old =cmplx(1.d0,0.d0,kind(0.d0))
    do nf = 1,N_Fl
       Call Compute_Fermion_Det(Z,Det_Vec_old(:,nf),UL(:,:,nf),DL(:,nf),VL(:,:,nf))
       Phase_det_old(nf) = Z
       Phase_old = Phase_old*Z
    Enddo
    call Op_phase(Phase_old,OP_V,Nsigma,N_SUN) 
    If (L_test) then
       Write(6,*) 'Testing global : ',  Phase_old, Phase
    Endif


    If (L_test) then 
       ! Testing    
       Do nf = 1,N_FL
          CALL INITD(UR(:,:,nf),Z_ONE)
          CALL INITD(VR(:,:,nf),Z_ONE)
          DR(:,nf) = Z_ONE
       Enddo
       NVAR = 1
       Phase = cmplx(1.d0,0.d0,kind(0.d0))
       do nf = 1,N_Fl
          CALL CGR(Z, NVAR, GR(:,:,nf), UR(:,:,nf),DR(:,nf),VR(:,:,nf),  UL(:,:,nf),DL(:,nf),VL(:,:,nf)  )
          Phase = Phase*Z
       Enddo
       call Op_phase(Phase,OP_V,Nsigma,N_SUN) 
       Do Nf = 1,N_FL
          Call DET_C_LU(GR(:,:,nf),Det_vec_new(:,nf),Ndim)
          Z = Phase_det_old(nf)
          DO I = 1,Ndim
             Z = Z*Det_vec_new(I,nf)*Det_vec_old(I,nf)
          Enddo
          Write(6,*) 'Testing weight: ', Z
       Enddo
    Endif

    !> Store old configuration
    nsigma_old = nsigma 
    !> Phase_old, Phase_det_old and Det_vec_old  are all set. 
    NC = 0
    Do n = 1,N_Global
       !> Draw a new spin configuration. This is provides by the user in the Hamiltonian module
       !> Note that nsigma is a variable in the module Hamiltonian
       Call Global_move(T0_Proposal_ratio,nsigma_old)
       If (T0_Proposal_ratio > 1.D-24) then
          NC = NC + 1
          !> Compute the new Green function
          DO nf = 1,N_FL
             CALL INITD(UL(:,:,Nf),Z_ONE)
             DL(:,nf) = Z_ONE
             CALL INITD(VL(:,:,nf),Z_ONE)
          ENDDO
          DO NST = NSTM-1,1,-1
             NT1 = Stab_nt(NST+1)
             NT  = Stab_nt(NST  )
             !Write(6,*) NT1,NT, NST
             CALL WRAPUL(NT1,NT,UL,DL, VL)
             Do nf = 1,N_FL
                UST(:,:,NST,nf) = UL(:,:,nf)
                VST(:,:,NST,nf) = VL(:,:,nf)
                DST(:  ,NST,nf) = DL(:  ,nf)
             ENDDO
          ENDDO
          NT1 = stab_nt(1)
          CALL WRAPUL(NT1,0, UL ,DL, VL)
          !You could now compute the det directly here.
          Phase_new = cmplx(1.d0,0.d0,kind(0.d0))
          do nf = 1,N_Fl
             Call Compute_Fermion_Det(Z,Det_Vec_new(:,nf),UL(:,:,nf),DL(:,nf),VL(:,:,nf))
             Phase_det_new(nf) = Z
             Phase_new = Phase_new*Z
          Enddo
          call Op_phase(Phase_new,OP_V,Nsigma,N_SUN) 
          
          Ratiotot = Compute_Ratio_Global(Phase_Det_old, Phase_Det_new, &
               &                               Det_vec_old, Det_vec_new, nsigma_old, T0_Proposal_ratio) 
          
          !Write(6,*) 'Ratio_global: ', Ratiotot
          
          Weight = abs(  real(Phase_old * Ratiotot, kind=Kind(0.d0))/real(Phase_old,kind=Kind(0.d0)) )
          
          Z = Phase_old * Ratiotot/ABS(Ratiotot)
          Call Control_PrecisionP_Glob(Z,Phase_new)
          !Write(6,*) Z, Phase_new
          
          
          TOGGLE = .false. 
          if ( Weight > ranf_wrap() )  Then
             TOGGLE = .true.
             Phase_old     = Phase_new
             Phase_det_old = Phase_det_new
             nsigma_old    = nsigma
             Det_vec_old   = Det_vec_new
          else
             nsigma = nsigma_old
          endif
          Call Control_upgrade_Glob(TOGGLE)
       endif
    Enddo

    If (NC > 0 ) then
       If (.not.TOGGLE) then
          DO nf = 1,N_FL
             CALL INITD(UL(:,:,Nf),Z_ONE)
             DL(:,nf) = Z_ONE
             CALL INITD(VL(:,:,nf),Z_ONE)
          ENDDO
          DO NST = NSTM-1,1,-1
             NT1 = Stab_nt(NST+1)
             NT  = Stab_nt(NST  )
             !Write(6,*) NT1,NT, NST
             CALL WRAPUL(NT1,NT,UL,DL, VL)
             Do nf = 1,N_FL
                UST(:,:,NST,nf) = UL(:,:,nf)
                VST(:,:,NST,nf) = VL(:,:,nf)
                DST(:  ,NST,nf) = DL(:  ,nf)
             ENDDO
          ENDDO
          NT1 = stab_nt(1)
          CALL WRAPUL(NT1,0, UL ,DL, VL)
       Endif
       !Compute the Green functions so as to provide correct starting point for the sequential updates.
       NVAR = 1
       Phase =cmplx(1.d0,0.d0,kind(0.d0))
       do nf = 1,N_Fl
          CALL CGR(Z, NVAR, GR(:,:,nf), UR(:,:,nf),DR(:,nf),VR(:,:,nf),  UL(:,:,nf),DL(:,nf),VL(:,:,nf)  )
          Phase = Phase*Z
       Enddo
       call Op_phase(Phase,OP_V,Nsigma,N_SUN) 
    endif
    
    Deallocate ( nsigma_old)
    Deallocate ( Det_vec_old  , Det_vec_new  ) 
    Deallocate ( Phase_Det_new, Phase_Det_old )


  End Subroutine Global_Updates



!--------------------------------------------------------------------
  Complex (Kind=Kind(0.d0)) Function Compute_Ratio_Global(Phase_Det_old, Phase_Det_new, &
       &                    Det_vec_old, Det_vec_new, nsigma_old, T0_Proposal_ratio)
!--------------------------------------------------------------------
!> @author
!> Fakher Assaad 
!>
!> @brief 
!> This fucntion computes ratio of weights  T0(nsigma--> nsigma_old) W(nsigma)/ 
!>                                          T0(nsigma_old--> nsigma) W(nsigma_old) =
!>                                          T0_Proposal_ratio   W(nsigma)/  W(nsigma_old)
!> 
!> Note that the new configuration, nsigma, is contained in the Hamiltonian moddule
!> The fermionic determinant stems from the routine Compute_Fermion_Det. 
!> 
!--------------------------------------------------------------------

    
    Implicit none

    !> Arguments
    Complex (Kind=Kind(0.d0)), allocatable, INTENT(IN) :: Phase_Det_old(:), Phase_Det_new(:), &
         &                                                Det_vec_old(:,:), Det_vec_new(:,:)
    Real    (Kind=Kind(0.d0)) :: T0_proposal_ratio 
    Integer,    allocatable :: nsigma_old(:,:)

    !> Local 
    Integer                                :: Nf, i, nt
    Complex (Kind=Kind(0.d0)) :: Z, Z1
    Real    (Kind=Kind(0.d0)) :: X


    X = 1.d0
    Do nf = 1,N_Fl
       DO I = 1,Ndim
          X= X*real(Det_vec_new(I,nf),kind(0.d0)) / Real(Det_vec_old(I,nf),kind(0.d0) )
       enddo
    enddo
    Z = cmplx(X,0.d0,kind(0.d0))
    Do nf = 1,N_FL
       Z = Z*Phase_Det_new(nf)/Phase_Det_old(nf)
    enddo
    Z = Z**N_SUN 

    Do I = 1,Size(Op_V,1)
       If (Op_V(i,1)%type == 2) then 
          X = 0.d0
          Do nt = 1,Ltrot
             if ( nsigma(i,nt) /= nsigma_old(i,nt) )  then 
                Z = Z * cmplx( Gaml(nsigma(i,nt),2)/Gaml(nsigma_old(i,nt),2),0.d0,kind(0.d0) ) 
                X = X + Phi(nsigma(i,nt),2) - Phi(nsigma_old(i,nt),2)
             endif
          Enddo
          Do nf = 1,N_FL
             Z = Z * exp(cmplx( X*Real(N_SUN,Kind(0.d0)), 0.d0,kind(0.d0)) * Op_V(i,nf)%g * Op_V(i,nf)%alpha )
          Enddo
       endif   
    Enddo
    Z =  Z * cmplx( Delta_S0_global(Nsigma_old),0.d0,kind(0.d0) )
    Z =  Z * cmplx( T0_Proposal_ratio, 0.d0,kind(0.d0))

    Compute_Ratio_Global = Z


  end Function Compute_Ratio_Global



!--------------------------------------------------------------------
  subroutine Compute_Fermion_Det(Phase,Det_Vec,UL,DL,VL)
!--------------------------------------------------------------------
!> @author 
!> Fakher Assaad 
!>
!> @brief 
!> Computes det( 1 +  VL*DL*UL)   = Phase * Det_vec(1)*..*Det_vec(Ndim)
!> Note that Phase is a unit complex number and the Det_vec contains only 
!> positive, implying real, numbers.
!>  
!--------------------------------------------------------------------

    Use  UDV_Wrap_mod

    Implicit none
    
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:), Intent(IN)   ::  UL, VL
    COMPLEX (Kind=Kind(0.d0)), Dimension(:),   Intent(In)   ::  DL
    Complex (Kind=Kind(0.d0)), Dimension(:),   Intent(OUT)  ::  Det_Vec
    Complex (Kind=Kind(0.d0)) :: Phase

    !> Local variables
    Integer ::  N_size, NCON, J
    COMPLEX (Kind=Kind(0.d0)) :: alpha,beta, Z, Z1
    COMPLEX (Kind=Kind(0.d0)), Dimension(:,:), Allocatable ::  TP, U, V
    COMPLEX (Kind=Kind(0.d0)), Dimension(:), Allocatable :: D

    N_size = SIZE(DL,1)
    NCON  = 0
    alpha = cmplx(1.d0,0.d0,kind(0.d0))
    beta  = cmplx(0.d0,0.d0,kind(0.d0))
    Allocate (TP(N_Size,N_Size), U(N_size,N_Size), D(N_size), V(N_size,N_size) )
    TP = CT(UL)
    DO J = 1,N_size
       TP(:,J) = TP(:,J) +  VL(:,J)*DL(J)
    ENDDO
    Call  UDV_WRAP_Pivot(TP,U,D,V,NCON,N_size,N_Size)
    Z  = DET_C(V, N_size) ! Det destroys its argument
    Call MMULT(TP,UL,U)
    Z1 = Det_C(TP,N_size) 
    Phase   = Z*Z1/ABS(Z*Z1)
    Det_vec = D
    Det_vec(1) = Det_vec(1)*ABS(Z*Z1)

    Deallocate ( TP,U,D,V )

  end subroutine Compute_Fermion_Det

end Module Global_mod
