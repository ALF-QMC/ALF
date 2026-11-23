!  Copyright (C) 2026 The ALF project
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

Module qmc_walker_mod

  !--------------------------------------------------------------------
  !> @author
  !> ALF-project
  !
  !> @brief
  !> State of one Monte Carlo walker: Green function, phase, the current
  !> left/right UDV propagations and the storage of UDV states at the
  !> stabilisation time slices.
  !>
  !> All large per-walker arrays live in this one type, so that their
  !> allocation (and, in a GPU build, their mapping to the device) happens
  !> in one place.
  !>
  !> @details
  !> - GR(Ndim,Ndim,N_FL): equal-time Green function; GR_Tilde and Test are
  !>   work arrays of the same shape (symmetrised G, precision check).
  !> - Phase, Phase_array(N_FL): total phase and phase per flavour.
  !> - udvl(N_FL_eff), udvr(N_FL_eff): current left and right propagations.
  !> - Storage: udvst(NSTM, N_FL_eff) holds the UDV states at the time slices
  !>   Stab_nt(0:NSTM), with Stab_nt(0) = 0 and Stab_nt(NSTM) = Ltrot.
  !>
  !> The UDV components are ALLOCATABLE, as the routines that receive them
  !> (WRAPUL, WRAPUR, Global_Updates, Langevin_HMC, Tau_m, Tau_p, ...)
  !> declare their dummy arguments ALLOCATABLE.
  !--------------------------------------------------------------------

  Use UDV_State_mod, only: UDV_State
  Use Hamiltonian_main, only: Ndim, N_FL, N_FL_eff, Calc_Fl_map, Projector, WF_L, WF_R

  Implicit none
  private

  Type, public :: qmc_walker_t
     ! Green function and work arrays of the same shape
     Complex (Kind=Kind(0.d0)), Dimension(:,:,:), Allocatable :: GR, GR_Tilde
     Complex (Kind=Kind(0.d0)), Dimension(:,:),   Allocatable :: Test
     ! Phase: total and per flavour
     Complex (Kind=Kind(0.d0))                                :: Phase
     Complex (Kind=Kind(0.d0)), Dimension(:),     Allocatable :: Phase_array
     ! Current left and right propagations
     Class (UDV_State),         Dimension(:),     Allocatable :: udvl, udvr
     ! Storage of UDV states at the stabilisation time slices
     Integer                                                  :: NSTM
     Integer,                   Dimension(:),     Allocatable :: Stab_nt
     Class (UDV_State),         Dimension(:,:),   Allocatable :: udvst
  end Type qmc_walker_t

  public :: Walker_layout, Walker_alloc, Walker_dealloc

  contains

    !--------------------------------------------------------------------
    !> @brief
    !> Sets the storage layout: NSTM stabilisation intervals of Nwrap time
    !> slices (the last one may be shorter), Stab_nt(0:NSTM) with
    !> Stab_nt(0) = 0 and Stab_nt(NSTM) = Ltrot.
    !--------------------------------------------------------------------
    Subroutine Walker_layout(walker, Ltrot, Nwrap)
      Type (qmc_walker_t), intent(inout) :: walker
      Integer,             intent(in)    :: Ltrot, Nwrap
      Integer :: n

      If ( mod(Ltrot,Nwrap) == 0  ) then
         walker%NSTM = Ltrot/Nwrap
      else
         walker%NSTM = Ltrot/Nwrap + 1
      endif
      allocate ( walker%Stab_nt(0:walker%NSTM) )
      walker%Stab_nt(0) = 0
      do n = 1,walker%NSTM -1
         walker%Stab_nt(n) = Nwrap*n
      enddo

      walker%Stab_nt(walker%NSTM) = Ltrot
    end Subroutine Walker_layout

    !--------------------------------------------------------------------
    !> @brief
    !> Allocates the walker arrays for the model in Hamiltonian_main and
    !> sets the UDV states to the identity (or to the trial wave functions
    !> in the projector code). Walker_layout must have been called.
    !--------------------------------------------------------------------
    Subroutine Walker_alloc(walker)
      Type (qmc_walker_t), intent(inout) :: walker
      Integer :: nf, nf_eff, n

      allocate(walker%Phase_array(N_FL))
      Allocate ( walker%Test(Ndim,Ndim), walker%GR(NDIM,NDIM,N_FL), walker%GR_Tilde(NDIM,NDIM,N_FL)  )
      ALLOCATE(walker%udvl(N_FL_eff), walker%udvr(N_FL_eff), walker%udvst(walker%NSTM, N_FL_eff))
      do nf_eff = 1, N_FL_eff
         nf=Calc_Fl_map(nf_eff)
         do n = 1, walker%NSTM
            CALL walker%udvst(n,nf_eff)%alloc(ndim)
         ENDDO
         if (Projector) then
            CALL walker%udvl(nf_eff)%init(ndim,'l',WF_L(nf)%P)
            CALL walker%udvr(nf_eff)%init(ndim,'r',WF_R(nf)%P)
            CALL walker%udvst(walker%NSTM, nf_eff)%reset('l',WF_L(nf)%P)
         else
            CALL walker%udvl(nf_eff)%init(ndim,'l')
            CALL walker%udvr(nf_eff)%init(ndim,'r')
            CALL walker%udvst(walker%NSTM, nf_eff)%reset('l')
         endif
      enddo
    end Subroutine Walker_alloc

    !--------------------------------------------------------------------
    !> @brief
    !> Releases everything Walker_layout and Walker_alloc allocated.
    !--------------------------------------------------------------------
    Subroutine Walker_dealloc(walker)
      Type (qmc_walker_t), intent(inout) :: walker
      Integer :: nf_eff, n

      DO nf_eff = 1, N_FL_eff
         CALL walker%udvl(nf_eff)%dealloc
         CALL walker%udvr(nf_eff)%dealloc
         do n = 1, walker%NSTM
            CALL walker%udvst(n,nf_eff)%dealloc
         ENDDO
      ENDDO
      DEALLOCATE(walker%udvl, walker%udvr, walker%udvst)
      DEALLOCATE(walker%GR, walker%TEST, walker%Stab_nt,walker%GR_Tilde)
      deallocate(walker%Phase_array)
    end Subroutine Walker_dealloc

end Module qmc_walker_mod
