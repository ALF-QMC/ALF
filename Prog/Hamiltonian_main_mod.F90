!  Copyright (C) 2016 - 2023 The ALF project
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
!       http://alf.physik.uni-wuerzburg.de
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module defines the interface between the Hamiltonians (= model and observables definition) and
!> the Monte Carlo core. Hamiltonians are defined as submodules of this module. The Monte Carlo core
!> has only access to public members of this module. For defining a new Hamiltonian named <new_ham_name>,
!> the user has to create the file Hamiltonians/Hamiltonian_<new_ham_name>_smod.F90 and add the line
!> <new_ham_name> to Hamiltonians.list.

!> @details
!> The public variables of this module are the following
!>
!>
!> @param [public] ham
!> \verbatim
!> class(ham_base), allocatable
!> Object that contains all the Hamiltonian-specific procedures needed by the Monte Carlo core.
!> The procedures can be overloaded in the Hamiltonians. \endverbatim
!>
!> @param [public] OP_V
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable
!> List of operators of type=1,2 and 3 describing the sequence of interactions on a time slice.
!> The first index runs over this sequence. The second corresponds to the flavor index.  \endverbatim
!>
!> @param [public] OP_T
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable
!> Sequence of  operators  accounting for the  hopping on a  time slice. This can include  various
!> checkerboard decompositions. The first index runs over this sequence. The second corresponds to
!> the flavor index. \endverbatim
!> *  The progagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n}  \f$.  That is
!> first the hopping and then the potential energy.
!>
!>@param [public] WF_L
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Left trial wave function.  \endverbatim
!>
!> @param [public] WF_R
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Right trial wave function.   For both wave functions the index runs over the flavor index. \endverbatim
!>
!> @param [public]  nsigma
!> \verbatim Type(Fields)
!> Contains all auxiliary fields in the variable f(:,:). The first index runs through the operator
!> sequence. The second through the time slices.   \endverbatim
!
!> @param [public]  Ndim
!> \verbatim Integer
!> Total number of orbitals. e.g. # unit cells * # orbitals per unit cell.  \endverbatim
!
!> @param [public]  N_FL
!> \verbatim Integer
!> # of flavors.  Propagation is block diagonal in flavors.  \endverbatim
!
!> @param [public]  N_SUN
!> \verbatim Integer
!> # of colors.  Propagation is color independent.  \endverbatim
!>
!> @param [public] Ltrot
!> \verbatim Integer
!> Available measurment interval in units of Delta Tau. \endverbatim
!>
!> @param [public] Thtrot
!>  \verbatim Integer
!> Effective projection parameter in units of Delta Tau.  (Only relevant if projective option is turned on) \endverbatim
!>
!> @param [public] Projector
!> \verbatim Logical
!> Flag for projector. If true then the total number of time slices will correspond to Ltrot + 2*Thtrot \endverbatim
!>
!> @param [public] Group_Comm
!> \verbatim Integer
!> Defines MPI communicator  \endverbatim
!
!> @param [public] Symm
!> \verbatim Logical  \endverbatim
!> If set to true then the green functions will be symmetrized
!> before being  sent to the Obser, ObserT subroutines.
!> In particular, the transformation,  \f$ \tilde{G} =  e^{-\Delta \tau T /2 } G e^{\Delta \tau T /2 } \f$
!> will be carried out  and \f$ \tilde{G} \f$  will be sent to the Obser and ObserT subroutines.  Note that
!> if you want to use this  feature, then you have to be sure the hopping and interaction terms are decomposed
!> symmetrically. If Symm is true, the propagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=N_T}^{1}e^{T_n/2} \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n/2}  \f$
!>
!>
!> You still have to add some docu for the other private variables in this module.
!>
!--------------------------------------------------------------------

module Hamiltonian_main
   use runtime_error_mod
   use files_mod
   use Operator_mod, only: operator
   use WaveFunction_mod, only: WaveFunction
   use Observables
   use Fields_mod, only: Fields
   use iso_fortran_env, only: output_unit, error_unit

   implicit none

   private
   public :: Alloc_Ham, ham_base, ham
#ifdef __PGI
   public :: Obs_scal, Obs_eq, Obs_tau
#endif
   type ham_base
   contains
      procedure, nopass :: ham_set => ham_set_base
      procedure, nopass :: Alloc_obs => Alloc_obs_base
      procedure, nopass :: Obser => Obser_base
      procedure, nopass :: ObserT => ObserT_base
      procedure, nopass :: Pr_obs => Pr_obs_base
      procedure, nopass :: Init_obs => Init_obs_base
      procedure, nopass :: Global_move_tau => Global_move_tau_base
      procedure, nopass :: Hamiltonian_set_nsigma => Hamiltonian_set_nsigma_base
      procedure, nopass :: Overide_global_tau_sampling_parameters => Overide_global_tau_sampling_parameters_base
      procedure, nopass :: Global_move => Global_move_base
      procedure, nopass :: Delta_S0_global => Delta_S0_global_base
      procedure, nopass :: S0 => S0_base
      procedure, nopass :: Ham_Langevin_HMC_S0 => Ham_Langevin_HMC_S0_base
      procedure, nopass :: weight_reconstruction => weight_reconstruction_base
      procedure, nopass :: GR_reconstruction => GR_reconstruction_base
      procedure, nopass :: GRT_reconstruction => GRT_reconstruction_base
      procedure, nopass :: Apply_B_HMC => Apply_B_HMC_base
#ifdef HDF5
      procedure, nopass :: write_parameters_hdf5 => write_parameters_hdf5_base
#endif
   end type ham_base

   class(ham_base), allocatable :: ham

   type(operator), dimension(:, :), allocatable, public :: Op_V
   type(operator), dimension(:, :), allocatable, public :: Op_T
   type(WaveFunction), dimension(:), allocatable, public :: WF_L
   type(WaveFunction), dimension(:), allocatable, public :: WF_R
   logical, dimension(:), allocatable, public :: Calc_Fl
   integer, dimension(:), allocatable, public :: Calc_Fl_map
   type(Fields), public        :: nsigma
   integer, public        :: Ndim
   integer, public        :: N_FL, N_FL_eff
   integer, public        :: N_SUN
   integer, public        :: Ltrot
   integer, public        :: Thtrot
   logical, public        :: Projector
   integer, public        :: Group_Comm
   logical, public        :: Symm
   logical, public        :: reconstruction_needed
   logical, public        :: leap_frog_bulk

   !>    Privat Observables
   type(Obser_Vec), dimension(:), allocatable :: Obs_scal
   type(Obser_Latt), dimension(:), allocatable :: Obs_eq
   type(Obser_Latt), dimension(:), allocatable :: Obs_tau

#include "Hamiltonians_interface.h"
!!$      This file will  be dynamically generated and appended
!!$      interface
!!$         module subroutine Ham_Alloc_Kondo()
!!$         end subroutine Ham_Alloc_Kondo
!!$         module subroutine Ham_Alloc_Hubbard()
!!$         end subroutine Ham_Alloc_Hubbard
!!$         module subroutine Ham_Alloc_Hubbard_Plain_Vanilla()
!!$         end subroutine Ham_Alloc_Hubbard_Plain_Vanilla
!!$         module subroutine Ham_Alloc_tV()
!!$         end subroutine Ham_Alloc_tV
!!$         module subroutine Ham_Alloc_LRC()
!!$         end subroutine Ham_Alloc_LRC
!!$         module subroutine Ham_Alloc_Z2_Matter()
!!$         end subroutine Ham_Alloc_Z2_Matter
!!$      end interface

contains

   subroutine Alloc_Ham(ham_name)
      implicit none
      character(len=64), intent(in) :: ham_name

      select case (str_to_upper(ham_name))
#include "Hamiltonians_case.h"
!!$       This file will be dynamically generated and appended
      case default
         write (error_unit, '("A","A","A")') 'Hamiltonian ', ham_name, ' not yet implemented!'
         call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
      end select
   end subroutine Alloc_Ham

   !--------------------------------------------------------------------
   !> @brief
   !> Sets the Hamiltonian.
   !> @details
   !> This has to be overloaded in the Hamiltonian submodule.
   !--------------------------------------------------------------------
   subroutine Ham_Set_base()
      implicit none
      write (error_unit, *) 'Ham_set not defined!'
      call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
   end subroutine Ham_Set_base

   !--------------------------------------------------------------------
   !> @author
   !> ALF Collaboration
   !>
   !> @brief
   !> Single spin flip S0 ratio
   !> @details
   !> S0=exp(-S0(new))/exp(-S0(old)) where the new configuration correpsonds to the old one up to
   !> a spin flip of Operator n on time slice nt
   !> @details
   !--------------------------------------------------------------------
   real(Kind=kind(0.d0)) function S0_base(n, nt, Hs_new)
      implicit none
      !> Operator index
      integer, intent(IN) :: n
      !> Time slice
      integer, intent(IN) :: nt
      !> New local field on time slice nt and operator index n
      complex(Kind=kind(0.d0)), intent(In) :: Hs_new

      S0_base = 1.d0
      if (Op_V(n, 1)%type == 1) then
         write (error_unit, *) 'function S0 not implemented'
         call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
      end if

   end function S0_base

   !--------------------------------------------------------------------
   !> @author
   !> ALF Collaboration
   !>
   !> @brief
   !> Specifiy the equal time and time displaced observables
   !> @details
   !--------------------------------------------------------------------
   subroutine Alloc_obs_base(Ltau)

      implicit none
      !>  Ltau=1 if time displaced correlations are considered.
      integer, intent(In) :: Ltau
      write (error_unit, *) "Warning: Alloc_obs not implemented."
   end subroutine Alloc_obs_base

   !--------------------------------------------------------------------
   !> @author
   !> ALF Collaboration
   !>
   !> @brief
   !> Global moves
   !>
   !> @details
   !>  This routine generates a
   !>  global update  and returns the propability T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)
   !> @param [IN] nsigma_old,  Type(Fields)
   !> \verbatim
   !>  Old configuration. The new configuration is stored in nsigma.
   !> \endverbatim
   !> @param [OUT]  T0_Proposal_ratio Real
   !> \verbatimam
   !>  T0_Proposal_ratio  =  T0( sigma_new -> sigma_old ) /  T0( sigma_old -> sigma_new)
   !> \endverbatim
   !> @param [OUT]  Size_clust Real
   !> \verbatim
   !>  Size of cluster that will be flipped.
   !> \endverbatim
   !-------------------------------------------------------------------
   subroutine Global_move_base(T0_Proposal_ratio, nsigma_old, size_clust)

      implicit none
      real(Kind=kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
      type(Fields), intent(IN)  :: nsigma_old

      write (error_unit, *) 'Global_move not implemented'
      call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)

   end subroutine Global_move_base

   !--------------------------------------------------------------------
   !> @author
   !> ALF Collaboration
   !>
   !> @brief
   !> Computes the ratio exp(S0(new))/exp(S0(old))
   !>
   !> @details
   !> This function computes the ratio \verbatim  e^{-S0(nsigma)}/e^{-S0(nsigma_old)} \endverbatim
   !> @param [IN] nsigma_old,  Type(Fields)
   !> \verbatim
   !>  Old configuration. The new configuration is stored in nsigma.
   !> \endverbatim
   !-------------------------------------------------------------------
   real(Kind=kind(0.d0)) function Delta_S0_global_base(Nsigma_old)

      !  This function computes the ratio:  e^{-S0(nsigma)}/e^{-S0(nsigma_old)}
      implicit none

      ! Arguments
      type(Fields), intent(IN) :: nsigma_old

      logical, save              :: first_call = .true.
      integer                    :: field_id, tau, Nfields, Ntau
      complex(kind=kind(0.0d0)) :: Hs_old

      Delta_S0_global_base = 1.d0
      Nfields = size(nsigma_old%f, 1)
      Ntau = size(nsigma_old%f, 2)
      do tau = 1, Ntau
         do field_id = 1, Nfields
            ! S0 returns S0=exp(-S0(HS_old))/exp(-S0(nsigma))
            ! purposly call ham%S0 instead of S0_base such that S0 may be provided in derived Hamiltonian
            Hs_old = nsigma_old%f(field_id, tau)
            ! note we need exp(-S0(new))/exp(-S0(old)) but nsigma is already the new config and we provide HS_old
            ! in contrast to HS_new. Hence, S0 returns the inverse of whar we need!
            Delta_S0_global_base = Delta_S0_global_base/ham%S0(field_id, tau, Hs_old)
         end do
      end do

      if (first_call) then
         write (output_unit, *)
         write (output_unit, *) "ATTENTION:     The base implementation of Delta_S0_global is used!"
         write (output_unit, *) "This relies on a proper version of S0, but may be inefficient for some models."
         write (output_unit, *) "Consider overwriting this function to benefit from model specific properties."
         write (output_unit, *) "Suppressing further printouts of this message."
         write (output_unit, *)
         first_call = .false.
      end if

   end function Delta_S0_global_base

   !--------------------------------------------------------------------
   !> @author
   !> ALF Collaboration
   !>
   !> @brief
   !> Computes equal time observables
   !> @details
   !> @param [IN] Gr   Complex(:,:,:)
   !> \verbatim
   !>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
   !> \endverbatim
   !> @param [IN] Phase   Complex
   !> \verbatim
   !>  Phase
   !> \endverbatim
   !> @param [IN] Ntau Integer
   !> \verbatim
   !>  Time slice
   !> \endverbatim
   !-------------------------------------------------------------------
   subroutine Obser_base(GR, Phase, Ntau, Mc_step_weight)

      implicit none

      complex(Kind=kind(0.d0)), intent(IN) :: GR(Ndim, Ndim, N_FL)
      complex(Kind=kind(0.d0)), intent(IN) :: PHASE
      integer, intent(IN)          :: Ntau
      real(Kind=kind(0.d0)), intent(IN) :: Mc_step_weight
      logical, save              :: first_call = .true.

      if (first_call) then
         write (error_unit, *) "Warning: Obser not implemented."
         first_call = .false.
      end if
   end subroutine Obser_base

   !--------------------------------------------------------------------
   !> @author
   !> ALF Collaboration
   !>
   !> @brief
   !> Computes time displaced  observables
   !> @details
   !> @param [IN] NT, Integer
   !> \verbatim
   !>  Imaginary time
   !> \endverbatim
   !> @param [IN] GT0, GTT, G00, GTT,  Complex(:,:,:)
   !> \verbatim
   !>  Green functions:
   !>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
   !>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
   !>  G00(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(0  )>
   !>  GTT(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(tau)>
   !> \endverbatim
   !> @param [IN] Phase   Complex
   !> \verbatim
   !>  Phase
   !> \endverbatim
   !-------------------------------------------------------------------
   subroutine ObserT_base(NT, GT0, G0T, G00, GTT, PHASE, Mc_step_weight)
      implicit none

      integer, intent(IN) :: NT
      complex(Kind=kind(0.d0)), intent(IN) :: GT0(Ndim, Ndim, N_FL), G0T(Ndim, Ndim, N_FL)
      complex(Kind=kind(0.d0)), intent(IN) :: G00(Ndim, Ndim, N_FL), GTT(Ndim, Ndim, N_FL)
      complex(Kind=kind(0.d0)), intent(IN) :: Phase
      real(Kind=kind(0.d0)), intent(IN) :: Mc_step_weight
      logical, save              :: first_call = .true.

      if (first_call) then
         write (error_unit, *) "Warning: ObserT not implemented."
         first_call = .false.
      end if

   end subroutine ObserT_base

   !--------------------------------------------------------------------
   !> @author
   !> ALF Collaboration
   !>
   !> @brief
   !> Prints out the bins.  No need to change this routine.
   !-------------------------------------------------------------------
   subroutine Pr_obs_base(LTAU)

      implicit none

      integer, intent(In) ::  Ltau

      !Local
      integer :: I

      if (allocated(Obs_scal)) then
         do I = 1, size(Obs_scal, 1)
            call Print_bin_Vec(Obs_scal(I), Group_Comm)
         end do
      end if
      if (allocated(Obs_eq)) then
         do I = 1, size(Obs_eq, 1)
            call Print_bin_Latt(Obs_eq(I), Group_Comm)
         end do
      end if
      if (allocated(Obs_tau)) then
         do I = 1, size(Obs_tau, 1)
            call Print_bin_Latt(Obs_tau(I), Group_Comm)
         end do
      end if

   end subroutine Pr_obs_base

   !--------------------------------------------------------------------
   !> @author
   !> ALF Collaboration
   !>
   !> @brief
   !> Initializes observables to zero before each bins.  No need to change
   !> this routine.
   !-------------------------------------------------------------------
   subroutine Init_obs_base(Ltau)

      implicit none
      integer, intent(In) :: Ltau

      ! Local
      integer :: I

      if (allocated(Obs_scal)) then
         do I = 1, size(Obs_scal, 1)
            call Obser_vec_Init(Obs_scal(I))
         end do
      end if

      if (allocated(Obs_eq)) then
         do I = 1, size(Obs_eq, 1)
            call Obser_Latt_Init(Obs_eq(I))
         end do
      end if

      if (allocated(Obs_tau)) then
         do I = 1, size(Obs_tau, 1)
            call Obser_Latt_Init(Obs_tau(I))
         end do
      end if

   end subroutine Init_obs_base

   !--------------------------------------------------------------------
   !> @author
   !> ALF Collaboration
   !>
   !> @brief
   !> Specify a global move on a given time slice tau.
   !>
   !> @details
   !> @param[in] ntau Integer
   !> \verbatim
   !>  Time slice
   !> \endverbatim
   !> @param[out] T0_Proposal_ratio, Real
   !> \verbatim
   !>  T0_Proposal_ratio = T0( sigma_new -> sigma ) /  T0( sigma -> sigma_new)
   !> \endverbatim
   !> @param[out] S0_ratio, Real
   !> \verbatim
   !>  S0_ratio = e^( S_0(sigma_new) ) / e^( S_0(sigma) )
   !> \endverbatim
   !> @param[out] Flip_length  Integer
   !> \verbatim
   !>  Number of flips stored in the first  Flip_length entries of the array Flip_values.
   !>  Has to be smaller than NDIM
   !> \endverbatim
   !> @param[out] Flip_list  Integer(Ndim)
   !> \verbatim
   !>  List of spins to be flipped: nsigma%f(Flip_list(1),ntau) ... nsigma%f(Flip_list(Flip_Length),ntau)
   !>  Note that Ndim = size(Op_V,1)
   !> \endverbatim
   !> @param[out] Flip_value  Real(Ndim)
   !> \verbatim
   !>  Flip_value(:)= nsigma%flip(Flip_list(:),ntau)
   !>  Note that Ndim = size(Op_V,1)
   !> \endverbatim
   !--------------------------------------------------------------------
   subroutine Global_move_tau_base(T0_Proposal_ratio, S0_ratio, &
         &                     Flip_list, Flip_length, Flip_value, ntau)

      implicit none
      real(Kind=kind(0.d0)), intent(OUT) :: T0_Proposal_ratio, S0_ratio
      integer, intent(OUT) :: Flip_list(:)
      complex(Kind=kind(0.d0)), intent(OUT) :: Flip_value(:)
      integer, intent(OUT) :: Flip_length
      integer, intent(IN)  :: ntau

      write (error_unit, *) 'Global_move_tau not implemented'
      call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
   end subroutine Global_move_tau_base

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> The user can set the initial field.
!>
!> @details
!> @param[OUT] Initial_field Real(:,:)
!> \verbatim
!>  Upon entry Initial_field is not allocated. If alloacted then it will contain the
!>  the initial field
!> \endverbatim
!--------------------------------------------------------------------
   subroutine Hamiltonian_set_nsigma_base(Initial_field)
      implicit none

      complex(Kind=kind(0.d0)), allocatable, dimension(:, :), intent(INOUT) :: Initial_field

      !  Consider  when we implement  different debugging  levels
!!$             write(output_unit,*)
!!$             write(output_unit,*) "ATTENTION:     Base implementation of Hamiltonian_set_nsigma is getting calling!"
!!$             write(output_unit,*) "This routine does not actually change the initial field configuration."
!!$             write(output_unit,*)

   end subroutine Hamiltonian_set_nsigma_base

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> This routine allows to user to  determine the global_tau sampling parameters at run time
!> It is especially usefull if these parameters are dependent on other parameters.
!>
!> @details
!> \endverbatim
!--------------------------------------------------------------------
   subroutine Overide_global_tau_sampling_parameters_base(Nt_sequential_start, Nt_sequential_end, N_Global_tau)

      implicit none
      integer, intent(INOUT) :: Nt_sequential_start, Nt_sequential_end, N_Global_tau

!!$             write(output_unit,*)
!!$             write(output_unit,*) "ATTENTION:     Base implementation of Overide_global_tau_sampling_parameters is getting calling!"
!!$             write(output_unit,*) "This routine does not actually change the parameters."
!!$             write(output_unit,*)

   end subroutine Overide_global_tau_sampling_parameters_base

   !--------------------------------------------------------------------
   !> @author
   !> ALF Collaboration
   !>
   !> @brief
   !>   Forces_0  = \partial S_0 / \partial s  are calculated and returned to  main program.
   !>
   !-------------------------------------------------------------------
   subroutine Ham_Langevin_HMC_S0_base(Forces_0)

      implicit none

      real(Kind=kind(0.d0)), intent(inout), allocatable :: Forces_0(:, :)

      logical, save                                      :: first_call = .true.

      Forces_0 = 0.d0

      if (first_call) then
         write (output_unit, *)
         write (output_unit, *) "ATTENTION:     Base implementation of Ham_Langevin_HMC_S0 is getting calling!"
         write (output_unit, *) "This assumes trivial S0 action and is likely incorrect!"
         write (output_unit, *) "Consider overwritting this routine according to the model in your Hamiltonian."
         write (output_unit, *) "Suppressing further printouts of this message."
         write (output_unit, *)
         first_call = .false.
      end if
      ! Johannes: I would actually like to terminate the code. I cannot come up with a scenario where Forces_0=0 is correct!

   end subroutine Ham_Langevin_HMC_S0_base

   !--------------------------------------------------------------------
   !> @author
   !> ALF Collaboration
   !>
   !> @brief
   !>   p_HMC are the conjugate momenta to the fields phi in the HMC updating scheme.
   !>   The relevant part of the action reads -1/2 p_HMC^T M^{-1} p_HMC, where M^{-1} = B^T * B is the
   !>   "mass" of the conjugate momenta (tuning parameter in HMC)
   !>   If we define p_tilde_HMC= B*p_HMC then the Forces to update p_tilde_HMC are:
   !>        p_tilde_HMC --> p_tilde_HMC - delta t B del H /del phi .
   !>   and the Forces to update phi are:
   !>        phi --> phi + delta t B^T p_tilde_HMC
   !>   By detault, M=1 => B=1, hence this routine does nothing!
   !>   We  note  that  Forces_HMC(:,:)   has  to be  understood as  a vection  with  superindex x=(n1,n2)
   !>   such  that B =  B(x,y).   Since  we  want  M to be positive definite  the  Kernel of  B has  to consist
   !>   solely of the zero vector.
   !>
   !> @details
   !> @param[inout] Forces_HMC Real(:,:)
   !> \verbatim
   !>  del H /del phi or p_tilde_HMC on input
   !>  Forces_HMS --> B^op Forces_HMC
   !> \endverbatim .
   !> @param[in] ltrans Logical
   !> \verbatim
   !>  ltrans = .true.  : op='T'
   !>  ltrans = .false. : op='N'
   !> \endverbatim
   !>
   !-------------------------------------------------------------------
   subroutine Apply_B_HMC_base(Forces_HMC, ltrans)

      implicit none

      real(Kind=kind(0.d0)), intent(inout), allocatable :: Forces_HMC(:, :)
      logical, intent(in)                 :: ltrans

   end subroutine Apply_B_HMC_base

!--------------------------------------------------------------------
!> @brief
!> Reconstructs dependent flavors of the configuration's weight.
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!--------------------------------------------------------------------
   subroutine weight_reconstruction_base(weight)
      implicit none
      complex(Kind=kind(0.d0)), intent(inout) :: weight(:)

      write (error_unit, *) 'weight_reconstruction not defined!'
      call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)

   end subroutine weight_reconstruction_base

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reconstructs dependent flavors of equal time Greens function
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!> @param [INOUT] Gr   Complex(:,:,:)
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!-------------------------------------------------------------------
   subroutine GR_reconstruction_base(GR)

      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT) :: GR(Ndim, Ndim, N_FL)

      write (error_unit, *) "Warning: GR_reconstruction not implemented."
      call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
   end subroutine GR_reconstruction_base

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Reconstructs dependent flavors of time displaced Greens function G0T and GT0
!> @details
!> This has to be overloaded in the Hamiltonian submodule.
!> @param [INOUT] GT0, G0T,  Complex(:,:,:)
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
!> \endverbatim
!-------------------------------------------------------------------
   subroutine GRT_reconstruction_base(GT0, G0T)
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT) :: GT0(Ndim, Ndim, N_FL), G0T(Ndim, Ndim, N_FL)

      write (error_unit, *) "Warning: GRT_reconstruction not implemented."
      call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
   end subroutine GRT_reconstruction_base

#ifdef HDF5
   subroutine write_parameters_hdf5_base(filename)
      implicit none

      character(len=64), intent(in) :: filename

      write (output_unit, *)
      write (output_unit, *) "ATTENTION:     Base implementation of write_parameters_hdf5 is getting calling!"
      write (output_unit, *) "This routine does not actually store the parameters in the hdf5 file."
      write (output_unit, *) "Usually, a python script is generating a model specific routine, granted that the Hamiltonian"
      write (output_unit, *) "respects the design rule. Without parameters, the HDF5 file's compatibility with pyALF is limited."
      write (output_unit, *)

   end subroutine write_parameters_hdf5_base
#endif

end module Hamiltonian_main
