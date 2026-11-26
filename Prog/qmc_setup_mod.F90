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
  Use QMC_runtime_var, only: apply_Langevin_overrides, &
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

#if defined(HDF5)
  public :: Init_hdf5_file
#endif
  public :: Write_info_end, Write_info_start, Write_info_parameters, Setup_update_schemes, Setup_fields, Setup_model

  contains


    !--------------------------------------------------------------------
    !> @brief
    !> Allocates and sets the Hamiltonian, checks the time-dependent couplings
    !> and the trial wave functions, sets the flavour map, the default
    !> measuring interval and the kinetic propagators (Hop_mod_init).
    !> Returns N_op, the number of interaction operators.
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

      ! Kinetic propagators exp(-Dtau T)
      Call Hop_mod_init
    end Subroutine Setup_model


    !--------------------------------------------------------------------
    !> @brief
    !> Auxiliary fields: types and flip protocols, random number generator
    !> and initial configuration.
    !--------------------------------------------------------------------
    Subroutine Setup_fields(N_op)
      Integer, intent(in) :: N_op
      Integer :: n, Seed_in
      Character (len=64) :: file_seeds
      Complex (Kind=Kind(0.d0)), allocatable, dimension(:,:) :: Initial_field

      call nsigma%make(N_op, Ltrot)
      Do n = 1,N_op
         nsigma%t(n)              = OP_V(n,1)%type
         nsigma%flip_protocol(n)  = OP_V(n,1)%flip_protocol
      Enddo

      File_seeds="seeds"
      Call Set_Random_number_Generator(File_seeds,Seed_in)
      !Write(6,*) Seed_in

      Call ham%Hamiltonian_set_nsigma(Initial_field)
      if (allocated(Initial_field)) then
         Call nsigma%in(Group_Comm,Initial_field)
         deallocate(Initial_field)
      else
         Call nsigma%in(Group_Comm)
      endif
    end Subroutine Setup_fields


    !--------------------------------------------------------------------
    !> @brief
    !> Sets up the update schemes: global-tau sampling parameters (and the
    !> Wrapgr work arrays), Langevin/HMC and global MALA; checks their
    !> compatibility.
    !--------------------------------------------------------------------
    Subroutine Setup_update_schemes(N_op)
      Integer, intent(in) :: N_op
      Integer :: n
      Integer :: tmp_Nt_sequential_start, tmp_Nt_sequential_end, tmp_N_Global_tau, tmp_N_Global_tau_MALA

      ! Global-tau sampling parameters. This comes before the Langevin
      ! overrides below, which switch Global_tau_moves off.
      If ( .not. ( get_Global_tau_moves() .or. get_Global_tau_MALA_moves() ) )  then
         ! This  corresponds to the default updating scheme
         call set_default_global_tau_sampling(Size(OP_V,1))
      else
         !  Gives the possibility to set parameters in the Hamiltonian file
         tmp_Nt_sequential_start = get_Nt_sequential_start()
         tmp_Nt_sequential_end   = get_Nt_sequential_end()
         tmp_N_Global_tau        = get_N_Global_tau()
         tmp_N_Global_tau_MALA   = get_N_Global_tau_MALA()
         Call ham%Overide_global_tau_sampling_parameters(tmp_Nt_sequential_start,tmp_Nt_sequential_end,tmp_N_Global_tau,tmp_N_Global_tau_MALA)
         call set_Nt_sequential_start(tmp_Nt_sequential_start)
         call set_Nt_sequential_end(tmp_Nt_sequential_end)
         call set_N_Global_tau(tmp_N_Global_tau)
         call set_N_Global_tau_MALA(tmp_N_Global_tau_MALA)
      endif
      If (get_N_Global_tau() > 0 .or. get_N_Global_tau_MALA() > 0 .or. get_Sequential_MALA() ) then
         Call Wrapgr_alloc
      endif

    !   Sequential = .true.
      !TODO: check if sequential is done if some fields are discrete (Warning or error termination?)
      if ( get_Langevin() .or.  get_HMC() ) then
         if ( get_Langevin() ) then
              if (is_root(Group_Comm)) then
                  Call check_langevin_schemes_and_variables()
#if defined(TEMPERING)
                  if ( get_N_exchange_steps() > 0 ) then
                      write(output_unit,*) "Langevin mode does not allow tempering updates."
                      write(output_unit,*) "Overwriting N_exchange_steps to 0."
                  endif
#endif
              endif
             call apply_Langevin_overrides()
         endif
         Call Langevin_HMC%make(get_Langevin(), get_HMC() , get_Global_MALA_moves(), get_Delta_t_Langevin_HMC(), get_Max_Force(), get_Leapfrog_Steps())
      else
         Call Langevin_HMC%set_Update_scheme(get_Langevin(), get_HMC(), .false. )
      endif
      if (get_Global_MALA_moves()) then
         Call Metropolis_Langevin%make(.False., .False., get_Global_MALA_moves(), get_Delta_t_MALA_global(), get_MAX_Force_MALA_global(), get_Leapfrog_Steps())
      else
         Call Metropolis_Langevin%set_Update_scheme(.False., .False., get_Global_MALA_moves() )
      endif 
      Call check_MALA_variables()
      Call check_update_schemes_compatibility()

      if ( get_Sequential_MALA() .or. get_Global_tau_MALA_moves() ) then
       Do n = 1,N_op
        if ( nsigma%t(n) /= 3 ) then
           write(output_unit,*)
           WRITE(output_unit,*) 'Warning:    Not all fields are of type 3.'
           WRITE(output_unit,*) 'Fields that are not of type 3 will not be updated with MALA updates.'
           write(output_unit,*)
         exit
        endif
       enddo
      endif
    end Subroutine Setup_update_schemes


    !--------------------------------------------------------------------
    !> @brief
    !> Writes the run parameters to the info file (group root only).
    !> chunk_size_gb is read only in MPI builds.
    !--------------------------------------------------------------------
    Subroutine Write_info_parameters(file_info, NSTM, N_op, chunk_size_gb)
      Character (len=*),       intent(in) :: file_info
      Integer,                 intent(in) :: NSTM, N_op
      Real (Kind=Kind(0.d0)),  intent(in) :: chunk_size_gb
      Integer :: n, n1
      Logical :: Toggle, Toggle1

      if (is_root(Group_Comm)) then
         Open (Unit = 50,file=file_info,status="unknown",position="append")
         Write(50,*) 'Sweeps                              : ', get_NSweep()
         If ( abs(get_CPU_MAX()) < Eps_small ) then
            Write(50,*) 'Bins                                : ', get_NBin()
            Write(50,*) 'No CPU-time limitation '
         else
            Write(50,'(" Prog will stop after hours:",2x,F8.4)') get_CPU_MAX()
         endif
         Write(50,*) 'Measure Int.                        : ', get_LOBS_ST(), get_LOBS_EN()
         Write(50,*) 'Stabilization,Wrap                  : ', get_Nwrap()
         Write(50,*) 'Nstm                                : ', NSTM
         Write(50,*) 'Ltau                                : ', get_Ltau()
         Write(50,*) '# of interacting Ops per time slice : ', Size(OP_V,1)
         If ( get_Propose_S0() ) &
              &  Write(50,*) 'Propose Ising moves according to  bare Ising action'
         if ( get_sequential() ) then
             Write(50,*)  '-------  Space-sequential Time-sequential moves -------' 
             Write(50,*) 'Nt_sequential_start: ', get_Nt_sequential_start()
             Write(50,*) 'Nt_sequential_end  : ', get_Nt_sequential_end()
             If (get_Sequential_MALA()) then 
                Write(50,*) 'MALA Update          '
                Write(50,*) 'delta_t            : ', get_Delta_t_MALA_sequential() 
                Write(50,*) 'Max_force          : ', get_Max_Force_MALA_sequential()
             endif
          endif    
          if ( get_Global_tau_moves() .or. get_Global_tau_MALA_moves() ) Then
             Write(50,*)  '-------  Space-global     Time-sequential moves -------' 
             If ( get_Global_tau_moves() ) Then
                Write(50,*) 'Metropolis'
                Write(50,*) 'N_Global_tau       : ', get_N_Global_tau()
             endif
             If ( get_Global_tau_MALA_moves() ) Then
                Write(50,*) 'MALA'
                Write(50,*) 'N_Global_tau_MALA  : ', get_N_Global_tau_MALA()
                Write(50,*) 'delta_t            : ', get_Delta_t_MALA_global_tau()
                Write(50,*) 'Max_Force          : ', get_Max_Force_MALA_global_tau()
             endif
          endif
          if ( get_Global_moves() .or. get_Langevin() .or. get_HMC() .or. get_Global_MALA_moves() ) then 
             Write(50,*)  '-------  Space-global     Time-global     moves -------' 
             if ( get_HMC() ) then
                Write(50,*) ' HMC '
                Write(50,*) 'Leapfrog_Steps    : ', get_Leapfrog_steps()
                Write(50,*) 'delta _t          : ', get_Delta_t_Langevin_HMC()
                Write(50,*) 'HMC_Sweeps        : ', get_N_HMC_sweeps()
             endif
             if ( get_Global_MALA_moves() ) then
                Write(50,*) 'User Defined MALA '
                Write(50,*) '# moves / sweep    : ', get_N_Global_MALA_sweeps()
                Write(50,*) 'delta_t            : ', get_Delta_t_MALA_global()
                Write(50,*) 'Max_Force          : ', get_MAX_Force_MALA_global()
             endif
             If ( get_Global_moves() ) Then
                Write(50,*) 'User Defined Metropolis   '
                Write(50,*) '#  moves / sweep   :', get_N_Global()
             Endif
             if (get_Langevin() ) then
                Write(50,*) 'Langevin '     
                Write(50,*) 'delta_t            : ', get_Delta_t_Langevin_HMC()
                Write(50,*) 'Max_Force          : ', get_Max_Force()
             endif
          endif


         !Write out info  for  amplitude and flip_protocol
         Toggle  = .false.
         Do n = 1,N_op
            if (nsigma%t(n) == 3 .or. nsigma%t(n) == 4)  Toggle = .true.
         Enddo
         if ( Toggle .and. (.not.get_Sequential_MALA()) ) then
            Write(50,*) 'Amplitude  for  t=3,4  vertices is  set to: ', get_Amplitude()
         endif
         Toggle  = .false.
         Do n = 1,N_op
            if (nsigma%t(n) == 4)  Toggle = .true.
         Enddo
         If (Toggle)  then 
            Write(50,"('Flip protocal  for  t=4  vertices is  set to')",advance="no")
            do n1  =  1,4 
               Toggle1 = .false.
               do n =1,N_op
                  if (nsigma%Flip_protocol(n)  ==  n1  .and. nsigma%t(n) == 4  )  Toggle1 = .true.
               enddo
               if  (Toggle1)   Write(50,"(I2,2x)",advance="no")   n1 
            enddo
            Write(50,*)
         endif

#if defined(MPI)
         Write(50,*) 'Number of mpi-processes : ', comm_size(Group_Comm)
         if(use_mpi_shm) Write(50,*) 'Using mpi-shared memory in chunks of ', chunk_size_gb, 'GB.'
#endif
#if defined(GIT)
         Write(50,*) 'This executable represents commit '&
              &      , GIT_COMMIT_HASH , ' of branch ' , GIT_BRANCH , '.'
#endif
#if defined(STAB1)
         Write(50,*) 'STAB1 is defined '
#endif
#if defined(STAB2)
         Write(50,*) 'STAB2 is defined '
#endif
#if defined(STAB3)
         Write(50,*) 'STAB3 is defined '
#endif
#if defined(STABLOG)
         Write(50,*) 'LOG is defined '
#endif
#if defined(QRREF)
         Write(50,*) 'QRREF is defined '
#endif
#if defined(TEMPERING) && !defined(PARALLEL_PARAMS)
         Write(50,*) '# of exchange steps  ', get_N_exchange_steps()
         Write(50,*) 'Tempering frequency  ', get_N_Tempering_frequency()
         Write(50,*) 'Tempering Calc_det   ', get_Tempering_calc_det()
#endif
         close(50)
      endif
    end Subroutine Write_info_parameters


    !--------------------------------------------------------------------
    !> @brief
    !> Start time to the info file and the git status to info_git (group root).
    !--------------------------------------------------------------------
    Subroutine Write_info_start(file_info, file_git_info)
      Character (len=*), intent(in) :: file_info, file_git_info

      if (is_root(Group_Comm)) then
         Open (Unit = 50,file=file_info,status="unknown",position="append")
         write(50,*) "START TIME: " // iso8601_datetime()
         close(50)
#if defined(GIT)
         Open (Unit=50, file=file_git_info, status="unknown", position="append")
         write(50,*) "=================================="
#include "git_status.h"
         write(50,*) "=================================="
         close(50)
#endif
      endif
    end Subroutine Write_info_start


    !--------------------------------------------------------------------
    !> @brief
    !> Effective number of bins (CPU_MAX runs) and end time to the info file.
    !--------------------------------------------------------------------
    Subroutine Write_info_end(file_info, NBin_eff)
      Character (len=*), intent(in) :: file_info
      Integer,           intent(in) :: NBin_eff

      If (is_root(Group_Comm)) then
         Open (Unit=50,file=file_info, status="unknown", position="append")
         if ( abs(get_CPU_MAX()) > Eps_small ) then
            Write(50,*)' Effective number of bins   : ', Nbin_eff
         endif
         write(50,*)'FIN TIME: ' // iso8601_datetime()
         Close(50)
      endif
    end Subroutine Write_info_end


#if defined(HDF5)
    !--------------------------------------------------------------------
    !> @brief
    !> Creates the HDF5 data file if needed and writes the model parameters
    !> to it (group root only).
    !--------------------------------------------------------------------
    Subroutine Init_hdf5_file(file_dat)
      Character (len=*), intent(in) :: file_dat
      INTEGER(HID_T) :: file_id
      Logical :: file_exists
      Integer :: ierr

      if (is_root(Group_Comm)) then
        CALL h5open_f(ierr)
        inquire (file=file_dat, exist=file_exists)
        IF (.not. file_exists) THEN
          ! Create HDF5 file
          CALL h5fcreate_f(file_dat, H5F_ACC_TRUNC_F, file_id, ierr)
          call h5ltset_attribute_string_f(file_id, '/', 'program_name', 'ALF', ierr)
          call h5fclose_f(file_id, ierr)
        endif
        call ham%write_parameters_hdf5(file_dat)

      endif
    end Subroutine Init_hdf5_file
#endif

end Module qmc_setup_mod
