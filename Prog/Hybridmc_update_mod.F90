!  Copyright (C) 2016 - 2019 The ALF project
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


      Module Hybrid_mod
        
        Use Hamiltonian
        Use UDV_State_mod
        Use Control
        Use Hop_mod
        Use Global_mod

        Implicit none
        
        real(Kind=Kind(0.d0))    , allocatable, private ::  pfield(:,:)         !momentum field
        real(Kind=Kind(0.d0))    , allocatable, private ::  Forces_fer(:,:)     !Fermion MD force
        real(Kind=Kind(0.d0))    , allocatable, private ::  Forces_bos(:,:)     !Boson   MD force
        Complex (Kind=Kind(0.d0)), allocatable, private ::  xfield_it_tmp(:), xfield_iw_tmp(:)

      Contains
!--------------------------------------------------------------------
!> @author 
!> Zihong Liu
!
!> @brief 
!> Fermion force calculation
!> At start of MD move
!>   On input GR is on the first time slice and  the storage is full with
!> left propagations.   Udvr  and Udvl are on time slice 1.
!>   On output. The  field configuration is  updated.  GR, Udvr,  Udvl and Udvst are as on input but with the
!> updated configuration.  
!> Equal time measurments as well as time displaced ones  is projector is true are also carried out. 
!> 
!--------------------------------------------------------------------


      subroutine Hybrid_MD_update(Phase, GR, GR_Tilde, Test, udvr, udvl, Stab_nt, udvst, LOBS_ST, LOBS_EN)

        Use UDV_State_mod
        Implicit none

        Interface
           SUBROUTINE WRAPUR(NTAU, NTAU1, UDVR)
             Use Hamiltonian
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVR
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUR
           SUBROUTINE WRAPUL(NTAU1, NTAU, udvl)
             Use Hamiltonian
             Use UDV_State_mod
             Implicit none
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: udvl
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUL
           SUBROUTINE CGR(PHASE,NVAR, GRUP, udvr, udvl)
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), INTENT(IN) :: udvl, udvr
             COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(Inout) :: GRUP
             COMPLEX(Kind=Kind(0.d0)) :: PHASE
             INTEGER         :: NVAR
           END SUBROUTINE CGR
        end Interface

        !  Arguments
        COMPLEX (Kind=Kind(0.d0)),                                INTENT(INOUT) :: Phase
        CLASS   (UDV_State), DIMENSION(:), ALLOCATABLE,           INTENT(INOUT) :: udvl, udvr
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:),  allocatable,INTENT(INOUT) :: GR, GR_Tilde
        CLASS(UDV_State), Dimension(:,:), ALLOCATABLE,            INTENT(INOUT) :: udvst
        INTEGER, dimension(:),   allocatable,                     INTENT(IN)    :: Stab_nt
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:),    allocatable,INTENT(INOUT) :: TEST
        Integer, intent(in) :: LOBS_ST, LOBS_EN


        !  Local variables.
        Integer :: NST, NSTM, NF, NT, NTAU, NTAU1, NT1, NVAR,N, N1,N2, I, NC, N_part,j, N_Type
        Real    (Kind=Kind(0.d0)) :: T0_Proposal_ratio, Weight
        Complex (Kind=Kind(0.d0)) :: Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0)), Z, Z1, Ratiotot, Phase_old, Phase_new
        Complex (Kind=Kind(0.d0)), allocatable :: Det_vec_test(:,:), Phase_Det_new(:), Phase_Det_old(:)
        Real    (Kind=Kind(0.d0)), allocatable :: Det_vec_old(:,:), Det_vec_new(:,:)
        Real    (Kind=Kind(0.d0)), allocatable :: pfield_old(:,:)
        Type   (Fields)   :: nsigma_old

        Complex (Kind=Kind(0.d0)) :: Ratio(2)
        Logical :: TOGGLE, L_Test
        Real    (Kind=Kind(0.d0)) :: size_clust
        Real    (Kind=Kind(0.d0)) :: ratio_2_test
        Character (Len=64)  :: storage
        Real    (Kind=Kind(0.d0)) :: spin

        n1 = size(nsigma%f,1)
        n2 = size(nsigma%f,2)
        NSTM = Size(udvst, 1)
        call nsigma_old%make(n1, n2)
        Allocate ( pfield_old(n1, n2) )

        Allocate ( Det_vec_old(NDIM,N_FL), Det_vec_new(NDIM,N_FL), Det_vec_test(NDIM,N_FL) )
        Allocate ( Phase_Det_new(N_FL), Phase_Det_old(N_FL) )

        ! initial momentum field and R field drawn from Gaussian distribution
        call Hybrid_pfield_initial

        ! Set up the old fermion determinant
        storage = "Empty"
        Call Compute_Fermion_Det(Phase_det_old,Det_Vec_old, udvl, udvst, Stab_nt, storage)
        Phase_old = cmplx(1.d0,0.d0,kind=kind(0.d0))
        Do nf = 1,N_Fl
           Phase_old = Phase_old*Phase_det_old(nf)
        Enddo
        Call Op_phase(Phase_old,OP_V,Nsigma,N_SUN)

        ! Store old configuration
        nsigma_old%f = nsigma%f
        nsigma_old%t = nsigma%t
        ! Store initial pfield
        pfield_old = pfield

        ! First calculate the fermion force and boson force
        call Hybrid_refresh_Greensfunction(Phase, GR, udvr, udvl, Stab_nt, udvst)
        call Hybrid_cal_force_fer(Phase, GR, Test, udvr, udvl, Stab_nt, udvst)
        call Hybrid_cal_force_bos
        ! Start of Molecular Dynamics
        do i = 1, mdstep
            call Hybrid_md_splitting(Phase, GR, Test, udvr, udvl, Stab_nt, udvst)
        enddo

        ! calculate the ratio
        storage = "Empty"
        Call Compute_Fermion_Det(Phase_det_new,Det_Vec_new, udvl, udvst, Stab_nt, storage)

        Phase_new = cmplx(1.d0,0.d0,kind=kind(0.d0))
        Do nf = 1,N_Fl
           Phase_new = Phase_new*Phase_det_new(nf)
        Enddo
        Call Op_phase(Phase_new,OP_V,Nsigma,N_SUN)
        Ratiotot = Hybrid_Compute_Ratio(Phase_Det_old, Phase_Det_new, &
                                Det_vec_old, Det_vec_new, nsigma_old, pfield_old, Ratio)
        
        Weight = abs(  real( Phase_old * Ratiotot, kind=Kind(0.d0))/real(Phase_old,kind=Kind(0.d0)) )

        Z = Phase_old * Ratiotot/ABS(Ratiotot)
        Call Control_PrecisionP_Glob(Z,Phase_new)

        ! Decide whether the MD updates are accepted
        TOGGLE = .false.
        if ( Weight > ranf_wrap() )  Then
           TOGGLE = .true.
           Phase_old     = Phase_new
           Phase_det_old = Phase_det_new
           nsigma_old%t  = nsigma%t
           nsigma_old%f  = nsigma%f
           Det_vec_old   = Det_vec_new
        else
           nsigma%t = nsigma_old%t
           nsigma%f = nsigma_old%f
        endif
        Call Control_upgrade_Hybrid(TOGGLE)

        !If (.not. TOGGLE) then
        !   call Hybrid_refresh_Greensfunction(Phase, GR, udvr, udvl, Stab_nt, udvst)
        !Endif
        call Hybrid_refresh_Greensfunction(Phase, GR, udvr, udvl, Stab_nt, udvst)

        ! Measure
        NST = 1
        DO NTAU = 0, LTROT-1
           NTAU1 = NTAU + 1
           
           Do nf = 1,N_FL
              CALL HOP_MOD_mmthr   (GR(:,:,nf), nf )
              CALL HOP_MOD_mmthl_m1(GR(:,:,nf), nf )
           Enddo
           Do n = 1, size(OP_V,1) 
              Do nf = 1, N_FL
                 spin = nsigma%phi(n,ntau1) ! Phi(nsigma(n,ntau1),Op_V(n,nf)%type)
                 N_type = 1
                 Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),spin,Ndim,N_Type)
                 N_type = 2
                 Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),spin,Ndim,N_Type)
              enddo
           enddo

           If (NTAU1 == Stab_nt(NST) ) then 
              NT1 = Stab_nt(NST-1)
              CALL WRAPUR(NT1, NTAU1, udvr)
              Z = cmplx(1.d0, 0.d0, kind(0.D0))
              Do nf = 1, N_FL
                 ! Read from storage left propagation from LTROT to  NTAU1
                 udvl(nf) = udvst(NST, nf)
                 ! Write in storage right prop from 1 to NTAU1
                 udvst(NST, nf) = udvr(nf)
                 NVAR = 1
                 IF (NTAU1 .GT. LTROT/2) NVAR = 2
                 TEST(:,:) = GR(:,:,nf)
                 CALL CGR(Z1, NVAR, GR(:,:,nf), UDVR(nf), UDVL(nf))
                 Z = Z*Z1
                 Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
              ENDDO
              call Op_phase(Z,OP_V,Nsigma,N_SUN) 
              Call Control_PrecisionP(Z,Phase)
              Phase = Z
              NST = NST + 1
           ENDIF
       
           IF (NTAU1.GE. LOBS_ST .AND. NTAU1.LE. LOBS_EN ) THEN
              If (Symm) then
                 Call Hop_mod_Symm(GR_Tilde,GR)
                 CALL Obser( GR_Tilde, PHASE, Ntau1 )
              else
                 CALL Obser( GR, PHASE, Ntau1 )
              endif
           endif

        enddo

        ! Deallocate the tmp array
        call nsigma_old%clear
        Deallocate ( Det_vec_old  , Det_vec_new, Det_vec_test  )
        Deallocate ( Phase_Det_new, Phase_Det_old )
        Deallocate ( pfield_old )

      end subroutine Hybrid_MD_update

      subroutine Hybrid_cal_force_bos
        Implicit none
        
        ! External function
        complex(Kind=Kind(0.d0)), external :: zdotu, zdotc

        ! local
        Integer :: NSTM, n, nf, NST, NTAU, nt, nt1, Ntau1, NVAR, N_Type, I, J

        Forces_bos = 0.d0
        do i = 1, size(OP_V,1) 
            xfield_it_tmp = dcmplx(nsigma%f(i,:))
            call onedimension_fft( ltrot, xfield_it_tmp )
            call zscal(ltrot, dcmplx(1.d0/dble(ltrot),0.d0), xfield_it_tmp, 1 )
            call vzmul(ltrot, xfield_it_tmp, w_coeffi, xfield_iw_tmp)
            call onedimension_invfft( ltrot, xfield_iw_tmp )
            do j = 1, ltrot
               Forces_bos(i, j) = -2.d0*dble(xfield_iw_tmp(j)) + nsigma%phi(i,j)
            enddo
        enddo

      endsubroutine Hybrid_cal_force_bos

      SUBROUTINE Hybrid_cal_force_fer(Phase, GR, Test, udvr, udvl, Stab_nt, udvst)
        Implicit none
        
        Interface
           SUBROUTINE WRAPUR(NTAU, NTAU1, UDVR)
             Use Hamiltonian
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVR
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUR
           SUBROUTINE WRAPUL(NTAU1, NTAU, UDVL)
             Use Hamiltonian
             Use UDV_State_mod
             Implicit none
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVL
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUL
           SUBROUTINE CGR(PHASE,NVAR, GRUP, udvr, udvl)
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), INTENT(IN) :: UDVL, UDVR
             COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(Inout) :: GRUP
             COMPLEX(Kind=Kind(0.d0)) :: PHASE
             INTEGER         :: NVAR
           END SUBROUTINE CGR
        end Interface
        
        CLASS(UDV_State), intent(inout), allocatable, dimension(:  ) :: udvl, udvr
        CLASS(UDV_State), intent(inout), allocatable, dimension(:,:) :: udvst
        Complex (Kind=Kind(0.d0)), intent(inout) :: Phase
        Complex (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:)   :: Test
        COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:,:) :: GR
        Integer, intent(in),  dimension(:), allocatable :: Stab_nt


        !Local
        Integer :: NSTM, n, nf, NST, NTAU, nt, nt1, Ntau1, NVAR, N_Type, I, J
        Complex (Kind=Kind(0.d0)) :: Z, Z1
        Real    (Kind=Kind(0.d0)) :: spin
        
        NSTM = Size(udvst, 1)
        
        Forces_fer = cmplx(0.d0,0.d0,Kind(0.d0))
        do nf = 1,N_FL
           if (Projector) then
              CALL udvr(nf)%reset('r',WF_R(nf)%P)
           else
              CALL udvr(nf)%reset('r')
           endif
        Enddo
        NST = 1
        DO NTAU = 0, LTROT-1
           NTAU1 = NTAU + 1
           
           Do nf = 1,N_FL
              CALL HOP_MOD_mmthr   (GR(:,:,nf), nf )
              CALL HOP_MOD_mmthl_m1(GR(:,:,nf), nf )
           Enddo
           Do n = 1, size(OP_V,1) 
              Do nf = 1, N_FL
                 spin = nsigma%phi(n,ntau1) ! Phi(nsigma(n,ntau1),Op_V(n,nf)%type)
                 N_type = 1
                 Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),spin,Ndim,N_Type)
                 N_type = 2
                 Call Op_Wrapup(Gr(:,:,nf),Op_V(n,nf),spin,Ndim,N_Type)
              enddo
              ! Compute forces here.
              if (OP_V(n,1)%type == 3 ) then
                 Do nf = 1, N_Fl
                    Z = cmplx(0.d0,0.d0,Kind(0.d0))
                    do I = 1,size(OP_V(n,nf)%P,1)
                       do J = 1,size(OP_V(n,nf)%P,1)
                          Z1 =  cmplx(0.d0,0.d0,Kind(0.d0))
                          if ( I == J ) Z1 = cmplx(1.d0,0.d0,Kind(0.d0))
                          Z  = Z + Op_V(n,nf)%O(I,J) * ( Z1 - Gr(Op_V(n,nf)%P(J),Op_V(n,nf)%P(I), nf) )
                       Enddo
                    Enddo
                    Forces_fer(n,ntau1) = Forces_fer(n,ntau1) - dble(Op_V(n,nf)%g*Z*cmplx(real(N_SUN,Kind(0.d0)),0.d0,Kind(0.d0)) )
                 Enddo
              endif
           enddo

           If (NTAU1 == Stab_nt(NST) ) then 
              NT1 = Stab_nt(NST-1)
              CALL WRAPUR(NT1, NTAU1, udvr)
              Z = cmplx(1.d0, 0.d0, kind(0.D0))
              Do nf = 1, N_FL
                 ! Read from storage left propagation from LTROT to  NTAU1
                 udvl(nf) = udvst(NST, nf)
                 ! Write in storage right prop from 1 to NTAU1
                 udvst(NST, nf) = udvr(nf)
                 NVAR = 1
                 IF (NTAU1 .GT. LTROT/2) NVAR = 2
                 TEST(:,:) = GR(:,:,nf)
                 CALL CGR(Z1, NVAR, GR(:,:,nf), UDVR(nf), UDVL(nf))
                 Z = Z*Z1
                 Call Control_PrecisionG(GR(:,:,nf),Test,Ndim)
              ENDDO
              call Op_phase(Z,OP_V,Nsigma,N_SUN) 
              Call Control_PrecisionP(Z,Phase)
              Phase = Z
              NST = NST + 1
           ENDIF
           
        enddo

        Call Control_Force_fer(Forces_fer,Group_Comm)

      END SUBROUTINE Hybrid_cal_force_fer

      subroutine Hybrid_md_splitting(Phase, GR, Test, udvr, udvl, Stab_nt, udvst)
        Implicit none

        !  Arguments
        COMPLEX (Kind=Kind(0.d0)),                                INTENT(INOUT) :: Phase
        CLASS   (UDV_State), DIMENSION(:), ALLOCATABLE,           INTENT(INOUT) :: udvl, udvr
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:,:),  allocatable,INTENT(INOUT) :: GR
        CLASS(UDV_State), Dimension(:,:), ALLOCATABLE,            INTENT(INOUT) :: udvst
        INTEGER, dimension(:),   allocatable,                     INTENT(IN)    :: Stab_nt
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:),    allocatable,INTENT(INOUT) :: TEST

        ! local variables
        integer :: nf, i, splitting_M

        !!! leapfrog with splitting Hamiltonian
        splitting_M = 10

        ! firstly update the momentum using force fermion
        pfield = pfield + dt/2.d0*Forces_fer

        ! M inner leapfrog dynamics with effective stepsize \epsilon /M
        do i =1, splitting_M
            pfield = pfield + (dt/dble(splitting_M))/2.d0*Forces_bos

            nsigma%f = nsigma%f + (dt/dble(splitting_M))*pfield

            ! calculate the force boson
            call Hybrid_cal_force_bos
            pfield = pfield + (dt/dble(splitting_M))/2.d0*Forces_bos
        enddo

        ! lastly update the momentum using force fermion (prepare the force fermion before update)
        call Hybrid_refresh_Greensfunction(Phase, GR, udvr, udvl, Stab_nt, udvst)
        call Hybrid_cal_force_fer(Phase, GR, Test, udvr, udvl, Stab_nt, udvst)
        pfield = pfield + dt/2.d0*Forces_fer

      endsubroutine Hybrid_md_splitting

      
      SUBROUTINE Hybrid_pfield_initial
        Implicit none

        integer :: i, nt, Nr, Ntau

        Nr = size(nsigma%f,1)
        Nt = size(nsigma%f,2)

        Do i = 1,Nr
           Do ntau = 1,Nt
               pfield(i, ntau) = rang_wrap()
           Enddo
        Enddo
        
      end SUBROUTINE Hybrid_pfield_initial

      Complex (Kind=Kind(0.d0)) Function  Hybrid_Compute_Ratio(Phase_Det_old, Phase_Det_new, &
                                Det_vec_old, Det_vec_new, nsigma_old, pfield_old, Ratio)

        Implicit none

        ! Arguments
        Complex (Kind=Kind(0.d0)), allocatable, INTENT(IN) :: Phase_Det_old(:), Phase_Det_new(:)
        REAL    (Kind=Kind(0.d0)), allocatable, INTENT(IN) :: Det_vec_old(:,:), Det_vec_new(:,:)
        REAL    (Kind=Kind(0.d0)), allocatable, INTENT(IN) :: pfield_old(:,:)
        Type    (Fields),             INTENT(IN)  :: nsigma_old
        Complex (Kind=Kind(0.d0)),    INTENT(out) :: Ratio(2)

        ! Local
        Integer  :: Nf, i, nt
        Complex (Kind=Kind(0.d0)) :: Z, Z1, z_ene_tmp, z_ene_tmp2
        Real    (Kind=Kind(0.d0)) :: X, Ratio_2, ener_p_old, ener_p_new

        ! External function
        complex(Kind=Kind(0.d0)), external :: zdotu, zdotc

        Ratio = cmplx(0.d0,0.d0,kind(0.d0))
        Ratio_2 = 0.d0
        Do nf = 1,N_Fl
           DO I = 1,Ndim
              Ratio_2 = Ratio_2 +  Det_vec_new(I,nf) - Det_vec_old(I,nf)
           enddo
        enddo
        Ratio(1) = cmplx(1.d0,0.d0,kind(0.d0))
        Do nf = 1,N_FL
           Ratio(1) = Ratio(1) *  Phase_Det_new(nf)/Phase_Det_old(nf)
        enddo
        Ratio(1) = Ratio(1)**N_SUN
        Ratio_2 = real(N_SUN,kind(0.d0))*Ratio_2

        Do I = 1,Size(Op_V,1)
           X = 0.d0
           Do nt = 1,Ltrot
              Ratio(1) = Ratio(1) * cmplx( nsigma%Gama(i,nt)/nsigma_old%Gama(i,nt),0.d0,kind(0.d0) )
              Ratio(1) = Ratio(1) * exp(cmplx( 0.5d0*(pfield_old(i,nt)**2-pfield(i,nt)**2),0.d0,kind(0.d0) ))
              X = X + nsigma%Phi(i,nt) - nsigma_old%Phi(i,nt)
           Enddo
           Do nf = 1,N_FL
              Ratio(1) = Ratio(1) * exp(cmplx( X*Real(N_SUN,Kind(0.d0)), 0.d0,kind(0.d0)) * Op_V(i,nf)%g * Op_V(i,nf)%alpha )
           Enddo
           xfield_it_tmp = dcmplx(nsigma%f(i,:))
           call onedimension_fft( ltrot, xfield_it_tmp )
           call zscal(ltrot, dcmplx(sqrt(1.d0/dble(ltrot)),0.d0), xfield_it_tmp, 1 )
           call vzmul(ltrot, xfield_it_tmp, w_coeffi, xfield_iw_tmp)
           z_ene_tmp = zdotc(ltrot, xfield_it_tmp, 1, xfield_iw_tmp, 1)
           xfield_it_tmp = dcmplx(nsigma_old%f(i,:))
           call onedimension_fft( ltrot, xfield_it_tmp )
           call zscal(ltrot, dcmplx(sqrt(1.d0/dble(ltrot)),0.d0), xfield_it_tmp, 1 )
           call vzmul(ltrot, xfield_it_tmp, w_coeffi, xfield_iw_tmp)
           z_ene_tmp2 = zdotc(ltrot, xfield_it_tmp, 1, xfield_iw_tmp, 1)
           Ratio(1) = Ratio(1) * exp(-(z_ene_tmp-z_ene_tmp2))
        Enddo

        Ratio(2) = Ratio_2
        Hybrid_Compute_Ratio = Ratio(1)*exp(Ratio(2))

      end function Hybrid_Compute_Ratio

      SUBROUTINE Hybrid_setup
        Implicit none

        Integer :: Nr,Nt
        
        Nr = size(nsigma%f,1)
        Nt = size(nsigma%f,2)
        Allocate ( pfield    (Nr,Nt) )
        Allocate ( Forces_fer(Nr,Nt) )
        Allocate ( Forces_bos(Nr,Nt) )
        Allocate ( xfield_it_tmp(Nt) )
        Allocate ( xfield_iw_tmp(Nt) )
      end SUBROUTINE Hybrid_setup


      SUBROUTINE Hybrid_clear
        Implicit none
        Deallocate ( pfield     )
        Deallocate ( Forces_fer )
        Deallocate ( Forces_bos )
        Deallocate ( xfield_it_tmp )
        Deallocate ( xfield_iw_tmp )
      end SUBROUTINE Hybrid_clear

      SUBROUTINE Hybrid_refresh_Greensfunction(Phase, GR, udvr, udvl, Stab_nt, udvst)
        Implicit none
        
        Interface
           SUBROUTINE WRAPUR(NTAU, NTAU1, UDVR)
             Use Hamiltonian
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVR
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUR
           SUBROUTINE WRAPUL(NTAU1, NTAU, UDVL)
             Use Hamiltonian
             Use UDV_State_mod
             Implicit none
             CLASS(UDV_State), intent(inout), allocatable, dimension(:) :: UDVL
             Integer :: NTAU1, NTAU
           END SUBROUTINE WRAPUL
           SUBROUTINE CGR(PHASE,NVAR, GRUP, udvr, udvl)
             Use UDV_Wrap_mod
             Use UDV_State_mod
             Implicit None
             CLASS(UDV_State), INTENT(IN) :: UDVL, UDVR
             COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(Inout) :: GRUP
             COMPLEX(Kind=Kind(0.d0)) :: PHASE
             INTEGER         :: NVAR
           END SUBROUTINE CGR
        end Interface
        
        CLASS(UDV_State), intent(inout), allocatable, dimension(:  ) :: udvl, udvr
        CLASS(UDV_State), intent(inout), allocatable, dimension(:,:) :: udvst
        Complex (Kind=Kind(0.d0)), intent(inout) :: Phase
        COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:,:,:) :: GR
        Integer, intent(in),  dimension(:), allocatable :: Stab_nt


        !Local
        Integer :: NSTM, n, nf, NST, NTAU, nt, nt1, Ntau1, NVAR, N_Type, I, J
        Complex (Kind=Kind(0.d0)) :: Z, Z1
        Real    (Kind=Kind(0.d0)) :: spin
        
        NSTM = Size(Stab_nt,1) - 1 ! is the same with NSTM = Size(udvst, 1)

        Do nf = 1,N_FL
           if (Projector) then
              CALL udvl(nf)%reset('l',WF_L(nf)%P)
              CALL udvst(NSTM, nf)%reset('l',WF_L(nf)%P)
           else
              CALL udvl(nf)%reset('l')
              CALL udvst(NSTM, nf)%reset('l')
           endif
        ENDDO
        
        DO NST = NSTM-1,1,-1
           NT1 = Stab_nt(NST+1)
           NT  = Stab_nt(NST  )
           !Write(6,*)'Hi', NT1,NT, NST
           CALL WRAPUL(NT1, NT, UDVL)
           Do nf = 1,N_FL
              UDVST(NST, nf) = UDVL(nf)
           ENDDO
        ENDDO
        NT1 = stab_nt(1)
        CALL WRAPUL(NT1, 0, UDVL)
        
        do nf = 1,N_FL
           if (Projector) then
              CALL udvr(nf)%reset('r',WF_R(nf)%P)
           else
              CALL udvr(nf)%reset('r')
           endif
        ENDDO
        
        NVAR = 1
        Phase = cmplx(1.d0, 0.d0, kind(0.D0))
        do nf = 1,N_Fl
           CALL CGR(Z, NVAR, GR(:,:,nf), UDVR(nf), UDVL(nf))
           Phase = Phase*Z
        Enddo
        call Op_phase(Phase,OP_V,Nsigma,N_SUN)
        
      END SUBROUTINE Hybrid_refresh_Greensfunction

    end Module Hybrid_mod
