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

! git.h defines GIT, GIT_COMMIT_HASH and GIT_BRANCH for the info output
#include "git.h"

Module qmc_setup_mod

  !--------------------------------------------------------------------
  !> @author
  !> ALF-project
  !
  !> @brief
  !> Setup steps of a run that main calls once before the bin loop (model,
  !> fields, update schemes, HDF5 file) and the output to the info file.
  !--------------------------------------------------------------------

  Use Hamiltonian_main, only: ham, Alloc_Ham, Validate_Ham_Variables, Setup_flavor_map, &
       &                      N_FL, N_FL_eff, Ltrot, Thtrot, Op_T, Op_V, nsigma, WF_L, WF_R, &
       &                      Calc_Fl, Calc_Fl_map, Projector, Group_Comm, leap_frog_bulk
  Use QMC_runtime_var, only: apply_CPU_MAX_policy, apply_Langevin_overrides, &
       &                     check_compatibility_reconstruct_greens_function, &
       &                     check_langevin_schemes_and_variables, check_MALA_variables, &
       &                     check_update_schemes_compatibility, get_Amplitude, get_CPU_MAX, &
       &                     get_Delta_t_Langevin_HMC, get_Delta_t_MALA_global, &
       &                     get_Delta_t_MALA_global_tau, get_Delta_t_MALA_sequential, &
       &                     get_Global_MALA_moves, get_Global_moves, get_Global_tau_MALA_moves, &
       &                     get_Global_tau_moves, get_ham_name, get_HMC, get_Langevin, &
       &                     get_Leapfrog_Steps, get_LOBS_EN, get_LOBS_ST, get_Ltau, get_Max_Force, &
       &                     get_MAX_Force_MALA_global, get_Max_Force_MALA_global_tau, &
       &                     get_Max_Force_MALA_sequential, get_N_Global, get_N_Global_MALA_sweeps, &
       &                     get_N_Global_tau, get_N_Global_tau_MALA, get_N_HMC_sweeps, get_NBin, &
       &                     get_NSweep, get_Nt_sequential_end, get_Nt_sequential_start, get_Nwrap, &
       &                     get_Propose_S0, get_sequential, get_Sequential_MALA, &
       &                     get_Tempering_calc_det, set_default_global_tau_sampling, &
       &                     set_default_values_measuring_interval, set_N_Global_tau, &
       &                     set_N_Global_tau_MALA, set_Nt_sequential_end, set_Nt_sequential_start
#if defined(TEMPERING)
  Use QMC_runtime_var, only: get_N_exchange_steps, get_N_Tempering_frequency
#endif
  Use runtime_error_mod, only: Terminate_on_error, ERROR_HAMILTONIAN
  Use Fields_mod, only: Fields_init
  Use Hop_mod, only: Hop_mod_init
  Use Wrapgr_mod, only: Wrapgr_alloc
  Use set_random, only: Set_Random_number_Generator
  Use Langevin_HMC_mod, only: Langevin_HMC, Metropolis_Langevin
  Use Natural_Constants, only: Eps_small
  Use mpi_shared_memory, only: use_mpi_shm
  Use comm_mod, only: is_root, comm_size
  Use iso8601_datetime_mod, only: iso8601_datetime
  Use iso_fortran_env, only: output_unit, error_unit
#ifdef HDF5
  Use hdf5
  Use h5lt
#endif

  Implicit none
  private

  public :: Setup_model

  contains


    !--------------------------------------------------------------------
    !> @brief
    !> Allocates and sets the Hamiltonian, checks the time-dependent couplings
    !> and the trial wave functions, sets the flavour map and the default
    !> measuring interval. Returns N_op, the number of interaction operators.
    !--------------------------------------------------------------------
    Subroutine Setup_model(N_op)
      Integer, intent(out) :: N_op
      Integer :: n, i, nf, nf_eff

      Call Fields_init(get_Amplitude())
      Call Alloc_Ham(get_ham_name())
      leap_frog_bulk = .false.
      Call ham%Ham_set()
      Call Validate_Ham_Variables()
      If ( allocated(Calc_Fl) ) then
         if ( any(.not. Calc_Fl) ) call check_compatibility_reconstruct_greens_function()
      endif
      ! Test  if  user  has  specified  correct  array  size  for time dependent Hamiltonians
      N_op = Size(OP_V,1)
      do n = 1, N_op
          do i = 1, N_Fl
              if (Op_V(n,i)%get_g_t_alloc()) then
                  if (size(Op_V(n,i)%g_t,1) /= Ltrot) then
                      write(error_unit,*) "Array size of time-dependent coupling Op_V%g_t has to be Ltrot!"
                      CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
                  endif
              endif
          enddo
      enddo
      do n = 1, size(Op_T,1)
          do i = 1, N_Fl
              if (Op_T(n,i)%get_g_t_alloc()) then
                  if (size(Op_T(n,i)%g_t,1) /= Ltrot) then
                      write(error_unit,*) "Array size of time-dependent coupling Op_T%g_t has to be Ltrot!"
                      CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
                  endif
              endif
          enddo
      enddo


      ! Flavours to compute: Calc_Fl (default all), N_FL_eff, Calc_Fl_map, reconstruction_needed
      Call Setup_flavor_map()

      if(Projector) then
         do nf_eff=1,N_fl_eff
            nf=Calc_Fl_map(nf_eff)
            if (.not. allocated(WF_R(nf)%P) .or. .not. allocated(WF_L(nf)%P)) then
               write(error_unit,*) "Projector is selected but there are no trial wave functions!"
               CALL Terminate_on_error(ERROR_HAMILTONIAN,__FILE__,__LINE__)
            endif
         enddo
      endif

      Call set_default_values_measuring_interval(Thtrot, Ltrot, Projector)
    end Subroutine Setup_model

end Module qmc_setup_mod
