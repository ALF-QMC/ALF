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
      procedure, nopass :: Count_obs => Count_obs_base
      procedure, nopass :: Init_obs => Init_obs_base
      procedure, nopass :: s0 => s0_base
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
   type(Fields), dimension(:), allocatable, public :: nsigma_bp
   integer, public        :: Ndim
   integer, public        :: N_FL
   integer, public        :: N_SUN
   integer, public        :: N_wlk
   integer, public        :: N_slat
   integer, public        :: N_grc
   integer, public        :: N_wlk_mpi
   integer, public        :: N_grc_mpi
   integer, public        :: ltrot
   integer, public        :: Group_Comm
   logical, public        :: Symm

   complex(Kind=kind(0.d0)), public :: fac_norm
   complex(Kind=kind(0.d0)), dimension(:), allocatable, public :: weight_k
   complex(Kind=kind(0.d0)), dimension(:), allocatable, public :: overlap
   complex(kind=kind(0.d0)), dimension(:), allocatable, public :: x_local

   !>    Privat Observables
   type(Obser_Vec),  dimension(:), allocatable :: obs_scal
   type(Obser_Latt), dimension(:), allocatable :: obs_eq
   type(Obser_Latt), dimension(:), allocatable :: obs_tau

#include "Hamiltonians_interface.h"

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
   subroutine Obser_base(GR, GR_mix, i_grc, re_w, sum_w, sum_o)

      implicit none

      complex(Kind=kind(0.d0)), intent(IN) :: GR(Ndim, Ndim, N_FL)
      complex(Kind=kind(0.d0)), intent(IN) :: GR_mix(Ndim, Ndim, N_FL)
      complex(Kind=kind(0.d0)), intent(IN) :: sum_w, sum_o
      real   (Kind=kind(0.d0)), intent(IN) :: re_w
      integer, intent(IN) :: i_grc
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
   subroutine ObserT_base(NT, GT0, G0T, G00, GTT, i_grc, re_w, sum_w, sum_o)
      implicit none

      integer, intent(IN) :: NT, i_grc
      complex(Kind=kind(0.d0)), intent(IN) :: GT0(Ndim, Ndim, N_FL), G0T(Ndim, Ndim, N_FL)
      complex(Kind=kind(0.d0)), intent(IN) :: G00(Ndim, Ndim, N_FL), GTT(Ndim, Ndim, N_FL)
      complex(Kind=kind(0.d0)), intent(IN) :: sum_w, sum_o
      real   (Kind=kind(0.d0)), intent(IN) :: re_w
      logical, save              :: first_call = .true.

      if (first_call) then
         write (error_unit, *) "Warning: ObserT not implemented."
         first_call = .false.
      end if

   end subroutine ObserT_base

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

      complex(Kind=kind(0.d0)), intent(in) :: GR(Ndim, Ndim, N_FL, N_grc)
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
   
   subroutine Count_obs_base

      implicit none

      !Local
      integer :: I

      if (allocated(obs_scal)) then
         do i = 1, size(obs_scal, 1)
            obs_scal(i)%n = obs_scal(i)%n + 1
            obs_scal(i)%ave_sign = obs_scal(i)%ave_sign + 1.d0
         end do
      end if
      if (allocated(obs_eq)) then
         do i = 1, size(obs_eq, 1)
            obs_eq(i)%n = obs_eq(i)%n + 1
            obs_eq(i)%ave_sign = obs_eq(i)%ave_sign + 1.d0
         end do
      end if
      if (allocated(obs_tau)) then
         do i = 1, size(obs_tau, 1)
            obs_tau(i)%n = obs_tau(i)%n + 1
            obs_tau(i)%ave_sign = obs_tau(i)%ave_sign + 1.d0
         end do
      end if

   end subroutine Count_obs_base

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
