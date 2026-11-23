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

Module qmc_driver_mod

  !--------------------------------------------------------------------
  !> @author
  !> ALF-project
  !
  !> @brief
  !> Steps of the Monte Carlo run that main calls in order: the global
  !> continuous-field updates (Langevin, HMC, MALA) of one walker.
  !>
  !> @details
  !> These updates are not ported to the GPU at first; a GPU build can guard
  !> them here in one place.
  !--------------------------------------------------------------------

  Use qmc_walker_mod, only: qmc_walker_t
  Use Hamiltonian_main, only: Projector
  Use Langevin_HMC_mod, only: Langevin_HMC, Langevin_HMC_Reset_storage, Metropolis_Langevin
  Use Tau_m_mod, only: Tau_m
  Use Tau_p_mod, only: Tau_p
  Use QMC_runtime_var, only: get_LOBS_ST, get_LOBS_EN, get_Ltau, get_Nwrap, &
       &                     get_N_HMC_sweeps, get_sequential, get_Global_MALA_moves

  Implicit none
  private

  public :: Langevin_update, HMC_update

  contains

    !--------------------------------------------------------------------
    !> @brief
    !> One Langevin update of all fields, with the equal-time measurements
    !> done inside the update and, for Ltau = 1, the time-displaced ones.
    !--------------------------------------------------------------------
    Subroutine Langevin_update(walker)
      Type (qmc_walker_t), intent(inout) :: walker
      Integer :: NST

      !  Carry out a Langevin update and calculate equal time observables.
      Call Langevin_HMC%update(walker%Phase, walker%GR, walker%GR_Tilde, walker%Test, walker%udvr, walker%udvl, walker%Stab_nt, walker%udvst, &
           &                   get_LOBS_ST(), get_LOBS_EN(), get_Ltau())

      IF ( get_Ltau() == 1 ) then
         If (Projector) then 
            NST = 0 
            Call Tau_p ( walker%udvl, walker%udvr, walker%udvst, walker%GR, walker%PHASE, walker%NSTM, walker%STAB_NT, NST, get_LOBS_ST(), get_LOBS_EN())
            call Langevin_HMC%set_L_Forces(.true.)
         else
            Call Tau_m( walker%udvst, walker%GR, walker%PHASE, walker%NSTM, get_Nwrap(), walker%STAB_NT, get_LOBS_ST(), get_LOBS_EN() )
            call Langevin_HMC%set_L_Forces(.true.)
         endif
      endif
    end Subroutine Langevin_update

    !--------------------------------------------------------------------
    !> @brief
    !> get_N_HMC_sweeps() HMC trajectories, then the measurements (time-displaced
    !> for Ltau = 1) unless the sequential sweep follows.
    !--------------------------------------------------------------------
    Subroutine HMC_update(walker)
      Type (qmc_walker_t), intent(inout) :: walker
      Integer :: NST, n

      ! A preceding MALA sweep may have changed both the forces and determinants.
      if ( get_sequential() .or. get_Global_MALA_moves() ) call Langevin_HMC%set_L_Forces(.False.)
      Do n=1, get_N_HMC_sweeps()
          !  Carry out a Langevin update and calculate equal time observables.
          Call Langevin_HMC%update(walker%Phase, walker%GR, walker%GR_Tilde, walker%Test, walker%udvr, walker%udvl, walker%Stab_nt, walker%udvst, &
               &                   get_LOBS_ST(), get_LOBS_EN(), get_Ltau())
          if (n /= get_N_HMC_sweeps()) then
             Call Langevin_HMC%calc_Forces(walker%Phase, walker%GR, walker%GR_Tilde, walker%Test, walker%udvr, walker%udvl, walker%Stab_nt, walker%udvst,&
                  &  get_LOBS_ST(), get_LOBS_EN(), .True. )
             Call Langevin_HMC_Reset_storage(walker%Phase, walker%GR, walker%udvr, walker%udvl, walker%Stab_nt, walker%udvst)
             call Langevin_HMC%set_L_Forces(.true.)
          endif
      enddo

      !Do time-displaced measurements if needed, else set Calc_Obser_eq=.True. for the very first leapfrog ONLY
      If ( .not. get_sequential() ) then
         IF ( get_Ltau() == 1 ) then
            If (Projector) then 
               NST = 0 
               Call Tau_p ( walker%udvl, walker%udvr, walker%udvst, walker%GR, walker%PHASE, walker%NSTM, walker%STAB_NT, NST, get_LOBS_ST(), get_LOBS_EN())
            else
               Call Tau_m( walker%udvst, walker%GR, walker%PHASE, walker%NSTM, get_Nwrap(), walker%STAB_NT, get_LOBS_ST(), get_LOBS_EN() )
            endif
         else
            Call Langevin_HMC%calc_Forces(walker%Phase, walker%GR, walker%GR_Tilde, walker%Test, walker%udvr, walker%udvl, walker%Stab_nt, walker%udvst,&
            &  get_LOBS_ST(), get_LOBS_EN(), .True. )
            Call Langevin_HMC_Reset_storage(walker%Phase, walker%GR, walker%udvr, walker%udvl, walker%Stab_nt, walker%udvst)
         endif
         call Langevin_HMC%set_L_Forces(.true.)
      endif
    end Subroutine HMC_update

end Module qmc_driver_mod
