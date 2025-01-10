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
   use Operator_mod, only: operator
   use WaveFunction_mod, only: WaveFunction
   use Observables
   use Fields_mod, only: Fields
   use iso_fortran_env, only: output_unit, error_unit

   implicit none

   private
   public :: Alloc_Ham, ham_base, ham
#ifdef __PGI
   public :: Obs_scal
#endif
   type ham_base
   contains
      procedure, nopass :: ham_set => ham_set_base
      procedure, nopass :: Alloc_obs => Alloc_obs_base
      procedure, nopass :: Obser => Obser_base
      procedure, nopass :: ObserT => ObserT_base
      procedure, nopass :: E0_local => E0_local_base
      procedure, nopass :: sum_weight => sum_weight_base
      procedure, nopass :: update_fac_norm => update_fac_norm_base
      procedure, nopass :: Pr_obs => Pr_obs_base
      procedure, nopass :: Init_obs => Init_obs_base
      procedure, nopass :: Init_obs_mc => Init_obs_mc_base
      procedure, nopass :: S0 => S0_base
      procedure, nopass :: weight_reconstruction => weight_reconstruction_base
      procedure, nopass :: GR_reconstruction => GR_reconstruction_base
      procedure, nopass :: GRT_reconstruction => GRT_reconstruction_base
      procedure, nopass :: bp_obsert => bp_obsert_base
      procedure, nopass :: obsert_mc => obsert_mc_base
      procedure, nopass :: set_xloc => set_xloc_base

#ifdef HDF5
      procedure, nopass :: write_parameters_hdf5 => write_parameters_hdf5_base
#endif
   end type ham_base

   class(ham_base), allocatable :: ham

   type(operator), dimension(:, :), allocatable, public :: Op_V
   type(operator), dimension(:, :), allocatable, public :: Op_T
   type(WaveFunction), dimension(:, :), allocatable, public :: WF_L
   type(WaveFunction), dimension(:, :), allocatable, public :: WF_R
   logical, dimension(:), allocatable, public :: Calc_Fl
   integer, dimension(:), allocatable, public :: Calc_Fl_map
   type(Fields), dimension(:), allocatable, public :: nsigma_bp
   integer, public        :: Ndim
   integer, public        :: N_FL, N_FL_eff
   integer, public        :: N_SUN
   integer, public        :: N_wlk
   integer, public        :: N_slat
   integer, public        :: N_grc
   integer, public        :: N_wlk_mpi
   integer, public        :: N_grc_mpi
   integer, public        :: ltrot
   integer, public        :: Group_Comm
   logical, public        :: Symm
   logical, public        :: reconstruction_needed

   complex(Kind=kind(0.d0)), public :: fac_norm
   complex(Kind=kind(0.d0)), dimension(:), allocatable, public :: weight_k
   complex(Kind=kind(0.d0)), dimension(:), allocatable, public :: overlap
   complex(Kind=Kind(0.d0)), dimension(:), allocatable, public :: x_local

   !>    Privat Observables
   type(obser_Vec) , dimension(:), allocatable :: obs_scal
   type(obser_Latt), dimension(:), allocatable :: obs_eq
   type(obser_Latt), dimension(:), allocatable :: obs_tau

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

      select case (ham_name)
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
      real(Kind=kind(0.d0)), intent(In) :: Hs_new

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
   subroutine Alloc_obs_base

      implicit none
      !>  Ltau=1 if time displaced correlations are considered.
      write (error_unit, *) "Warning: Alloc_obs not implemented."
   end subroutine Alloc_obs_base

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
   !> @param [IN] Ntau Integer
   !> \verbatim
   !>  Time slice
   !> \endverbatim
   !-------------------------------------------------------------------
   subroutine Obser_base(GR, GR_mix, i_wlk, i_grc, sum_w, sum_o, act_mea)

      implicit none

      complex(Kind=kind(0.d0)), intent(IN) :: GR(Ndim, Ndim, N_FL)
      complex(Kind=kind(0.d0)), intent(IN) :: GR_mix(Ndim, Ndim, N_FL)
      complex(Kind=kind(0.d0)), intent(IN) :: sum_w, sum_o
      integer, intent(IN) :: i_wlk, i_grc, act_mea
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
   !-------------------------------------------------------------------
   subroutine ObserT_base(NT, GT0, G0T, G00, GTT, i_wlk, i_grc, sum_w, sum_o, act_mea)
      implicit none

      integer, intent(IN) :: NT, i_wlk, i_grc, act_mea
      complex(Kind=kind(0.d0)), intent(IN) :: GT0(Ndim, Ndim, N_FL), G0T(Ndim, Ndim, N_FL)
      complex(Kind=kind(0.d0)), intent(IN) :: G00(Ndim, Ndim, N_FL), GTT(Ndim, Ndim, N_FL)
      complex(Kind=kind(0.d0)), intent(IN) :: sum_w, sum_o
      logical, save              :: first_call = .true.

      if (first_call) then
         write (error_unit, *) "Warning: ObserT not implemented."
         first_call = .false.
      end if

   end subroutine ObserT_base

   subroutine bp_obsert_base(i_wlk, sum_w, act_mea)
      implicit none

      integer, intent(in) :: i_wlk, act_mea
      complex(Kind=kind(0.d0)), intent(in) :: sum_w

   end subroutine bp_obsert_base

   subroutine obsert_mc_base(nt, gt0, g0t, g00, gtt, overlap_mc)
      implicit none

      integer, intent(IN) :: nt
      complex(Kind=kind(0.d0)), intent(IN) :: gt0(ndim, ndim, n_fl, n_grc), g0t(ndim, ndim, n_fl, n_grc)
      complex(Kind=kind(0.d0)), intent(IN) :: g00(ndim, ndim, n_fl, n_grc), gtt(ndim, ndim, n_fl, n_grc)
      complex(Kind=kind(0.d0)), intent(IN) :: overlap_mc(n_grc)

   end subroutine obsert_mc_base

   complex(Kind=kind(0.d0)) function E0_local_base(GR)
      implicit none

      complex(Kind=kind(0.d0)), intent(IN) :: GR(Ndim, Ndim, N_FL)

      E0_local_base = cmplx(0.d0, 0.d0, kind(0.d0))

   end function E0_local_base

   subroutine set_xloc_base(GR)
     Implicit none
      
     Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
     
   end subroutine set_xloc_base

   subroutine sum_weight_base(z_sum_weight)
      implicit none

      complex(Kind=kind(0.d0)), intent(out) :: z_sum_weight

      z_sum_weight = cmplx(0.d0, 0.d0, kind(0.d0))

   end subroutine sum_weight_base

   subroutine update_fac_norm_base(GR, ntw)
      implicit none

      complex(Kind=kind(0.d0)), intent(in) :: GR(Ndim, Ndim, N_FL, N_wlk)
      integer, intent(in) :: ntw

   end subroutine update_fac_norm_base

   !--------------------------------------------------------------------
   !> @author
   !> ALF Collaboration
   !>
   !> @brief
   !> Prints out the bins.  No need to change this routine.
   !-------------------------------------------------------------------
   subroutine Pr_obs_base(LTAU)

      implicit none

      integer, intent(In) :: Ltau

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

   subroutine Init_obs_mc_base

      implicit none

   end subroutine Init_obs_mc_base

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
