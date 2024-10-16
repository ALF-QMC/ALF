!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module computes and stores the exponential of the hopping matrix.  It also provide routines to carry out multiplications  with
!> \f$  e^{ - \Delta \tau  H_t  }   = \prod_{n=1}^{N} e^{ - \Delta \tau_n  H_t(n) }   \f$,
!> \f$  e^{   \Delta \tau  H_t  }   = \prod_{n=N}^{1} e^{   \Delta \tau_n  H_t(n) }   \f$,
!> \f$  e^{ - \Delta \tau  H_t/2  }   = \prod_{n=1}^{N} e^{ - \Delta \tau_n  H_t(n)/2 }   \f$, and
!> \f$  e^{   \Delta \tau  H_t/2  }   = \prod_{n=N}^{1} e^{   \Delta \tau_n  H_t(n)/2 }   \f$.
!> The last equality are important for the symmetric Trotter option. (See variable: Symm in the Hamiltonian module)
!
!>
!--------------------------------------------------------------------

module Hop_mod

   use Hamiltonian_main
   use Operator_mod
   use Random_wrap
   use DynamicMatrixArray_mod
   use ContainerElementBase_mod
   use OpTTypes_mod
   use iso_fortran_env, only: output_unit, error_unit

   ! Private variables
   type(DynamicMatrixArray), private, allocatable :: ExpOpT_vec(:) ! for now we have for simplicity for each flavour a vector
   integer, private, save ::  Ncheck
   real(Kind=kind(0.d0)), private, save  :: Zero

contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This function serves as a central entry point to collect the
!> processing that occurs in mapping an OpT input matrix to the internal
!> matrix-like data structure.
!
!> @param ExpOpT_vec[inout] a DynamicMatrixArray structure to which we append new elements.
!> @param op[in] an Operator that describes an OpT hopping matrix.
!
!--------------------------------------------------------------------
   subroutine OpT_postprocess(ExpOpT_vec, op)
      use Operator_mod
      implicit none

      type(DynamicMatrixArray), intent(inout) :: ExpOpT_vec
      type(operator), intent(in) :: op

      class(CmplxExpOpT), pointer :: cmplxexp => null()
      class(RealExpOpT), pointer :: realexp => null()

      if (Op_is_real(op)) then
         ! branch for real operators
         allocate (realexp) ! Yep, this is a manifest memory leak. Using the ptr we can allocate onto the same variable
         call realexp%init(op)
         call ExpOpT_vec%pushback(realexp)
      else
         ! branch for complex operators
         allocate (cmplxexp)
         call cmplxexp%init(op)
         call ExpOpT_vec%pushback(cmplxexp)
      end if

   end subroutine

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This functions sets up the exponentiated matrices.
!> We symmetrize the upper part of those matrices.
!
!--------------------------------------------------------------------
   subroutine Hop_mod_init
      use runtime_error_mod
      implicit none

      integer :: nc, nf

      Ncheck = size(Op_T, 1)
      if (size(Op_T, 2) /= N_FL) then
         write (error_unit, *) 'Hop_mod_init: Error in the number of flavors.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      allocate (ExpOpT_vec(N_FL))

      do nf = 1, N_FL
         call ExpOpT_vec(nf)%init()
         do nc = 1, Ncheck
            call OpT_postprocess(ExpOpT_vec(nf), Op_T(nc, nf))
         end do
      end do

      Zero = 1.e-12
   end subroutine Hop_mod_init

!--------------------------------------------------------------------

   subroutine Hop_mod_mmthr(In, nf, t)

      ! InOut:  In = e^{ -dtau T }.IN
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT)  :: IN(:, :)
      integer, intent(IN) :: nf, t

      !Local
      integer :: nc
      class(ContainerElementBase), pointer :: dummy

      do nc = Ncheck, 1, -1
         dummy => ExpOpT_vec(nf)%at(nc)
         call dummy%lmult(In, t)
      end do
   end subroutine Hop_mod_mmthr

   subroutine Hop_mod_mmthr_1D2(In, nf, t)

      ! InOut:  In = e^{ -dtau T }.IN
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT)  :: IN(:, :)
      integer, intent(IN) :: nf, t

      !Local
      integer :: nc
      class(ContainerElementBase), pointer :: dummy

      do nc = Ncheck, 1, -1
         dummy => ExpOpT_vec(nf)%at(nc)
         call dummy%lmult1D2(In, t)
      end do
   end subroutine Hop_mod_mmthr_1D2

   subroutine Hop_mod_mmthr_m1(In, nf, t)

      ! InOut:  In = e^{  dtau T }.IN
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT)  :: IN(:, :)
      integer :: nf
      integer, intent(in) :: t

      !Local
      integer :: nc
      class(ContainerElementBase), pointer :: dummy

      do nc = 1, Ncheck
         dummy => ExpOpT_vec(nf)%at(nc)
         call dummy%lmultinv(In, t)
      end do

   end subroutine Hop_mod_mmthr_m1

   subroutine Hop_mod_mmthr_m1_1D2(In, nf, t)

      ! InOut:  In = e^{  dtau T }.IN
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT)  :: IN(:, :)
      integer :: nf
      integer, intent(in) :: t

      !Local
      integer :: nc
      class(ContainerElementBase), pointer :: dummy

      do nc = 1, Ncheck
         dummy => ExpOpT_vec(nf)%at(nc)
         call dummy%lmultinv1D2(In, t)
      end do

   end subroutine Hop_mod_mmthr_m1_1D2

   subroutine Hop_mod_mmthlc_m1_1D2(In, nf, t)

      ! InOut:  In = e^{  dtau T }.IN
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT)  :: IN(:, :)
      integer :: nf
      integer, intent(in) :: t

      !Local
      integer :: nc
      class(ContainerElementBase), pointer :: dummy

      do nc = Ncheck, 1, -1
         dummy => ExpOpT_vec(nf)%at(nc)
         call dummy%lmultinv1D2(In, t)
      end do

   end subroutine Hop_mod_mmthlc_m1_1D2

!--------------------------------------------------------------------

   subroutine Hop_mod_mmthl(In, nf, t)

      ! InOut:  In = IN * e^{ -dtau T }
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT)  :: IN(:, :)
      integer :: nf
      integer, intent(in) :: t

      !Local
      integer :: nc
      class(ContainerElementBase), pointer :: dummy

      do nc = 1, Ncheck
         dummy => ExpOpT_vec(nf)%at(nc)
         call dummy%rmult(In, t)
      end do

   end subroutine Hop_mod_mmthl

   subroutine Hop_mod_mmthl_1D2(In, nf, t)

      ! InOut:  In = e^{ -dtau T }.IN
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT)  :: IN(:, :)
      integer, intent(IN) :: nf, t

      !Local
      integer :: nc
      class(ContainerElementBase), pointer :: dummy

      do nc = Ncheck, 1, -1
         dummy => ExpOpT_vec(nf)%at(nc)
         call dummy%rmult1D2(In, t)
      end do
   end subroutine Hop_mod_mmthl_1D2

!--------------------------------------------------------------------

   subroutine Hop_mod_mmthlc(In, nf, t)

      ! InOut:  In = IN * e^{ -dtau T }
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT)  :: IN(:, :)
      integer :: nf
      integer, intent(in) :: t

      !Local
      integer :: nc
      class(ContainerElementBase), pointer :: dummy

      do nc = 1, Ncheck
         dummy => ExpOpT_vec(nf)%at(nc)
         call dummy%lmult(In, t)
      end do

   end subroutine Hop_mod_mmthlc

!--------------------------------------------------------------------

   subroutine Hop_mod_mmthlc_1D2(In, nf, t)

      ! InOut:  In = IN * e^{ -dtau T }
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT)  :: IN(:, :)
      integer :: nf
      integer, intent(in) :: t

      !Local
      integer :: nc
      class(ContainerElementBase), pointer :: dummy

      do nc = 1, Ncheck
         dummy => ExpOpT_vec(nf)%at(nc)
         call dummy%lmult1D2(In, t)
      end do

   end subroutine Hop_mod_mmthlc_1D2

!--------------------------------------------------------------------

   subroutine Hop_mod_mmthl_m1(In, nf, t)

      ! InOut:  In = IN * e^{ dtau T }
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT)  :: IN(:, :)
      integer :: nf
      integer, intent(in) :: t

      !Local
      integer :: nc
      class(ContainerElementBase), pointer :: dummy

      do nc = Ncheck, 1, -1
         dummy => ExpOpT_vec(nf)%at(nc)
         call dummy%rmultinv(In, t)
      end do

   end subroutine Hop_mod_mmthl_m1

!--------------------------------------------------------------------

   subroutine Hop_mod_mmthl_m1_1D2(In, nf, t)

      ! InOut:  In = IN * e^{ dtau T }
      implicit none

      complex(Kind=kind(0.d0)), intent(INOUT)  :: IN(:, :)
      integer :: nf
      integer, intent(in) :: t

      !Local
      integer :: nc
      class(ContainerElementBase), pointer :: dummy

      do nc = Ncheck, 1, -1
         dummy => ExpOpT_vec(nf)%at(nc)
         call dummy%rmultinv1D2(In, t)
      end do

   end subroutine Hop_mod_mmthl_m1_1D2

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Out = exp(-Delta tau /2 T) In exp( Delta tau /2 T)
!>
!
!--------------------------------------------------------------------
   subroutine Hop_mod_Symm(Out, In, t1, t2)

      implicit none
      complex(Kind=kind(0.d0)), dimension(:, :, :), intent(Out):: Out
      complex(Kind=kind(0.d0)), dimension(:, :, :), intent(IN):: In
      integer, intent(in) :: t1
      integer, optional, intent(in) :: t2

      integer :: nf, nc, nf_eff
      class(ContainerElementBase), pointer :: dummy

      Out = In
      do nf_eff = 1, N_FL_eff !size(In,3)
         nf = Calc_Fl_map(nf_eff)
         do nc = Ncheck, 1, -1
            dummy => ExpOpT_vec(nf)%at(nc)
            call dummy%adjointaction(Out(:, :, nf), t1, t2)
         end do
      end do

   end subroutine Hop_mod_Symm
end module Hop_mod
