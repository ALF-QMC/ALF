! Declare a common base class for the interface: it multiplies with complex matrices
module ContainerElementBase_mod
   implicit none

   private
   public :: ContainerElementBase

   ! Base for defining the interface
   type, abstract :: ContainerElementBase
   contains
      procedure(rmultinterface), deferred :: rmult
      procedure(lmultinterface), deferred :: lmult
      procedure(rmult1D2interface), deferred :: rmult1D2
      procedure(lmult1D2interface), deferred :: lmult1D2
      procedure(rmultinvinterface), deferred :: rmultinv
      procedure(lmultinvinterface), deferred :: lmultinv
      procedure(rmultinv1D2interface), deferred :: rmultinv1D2
      procedure(lmultinv1D2interface), deferred :: lmultinv1D2
      procedure(adjointactioninterface), deferred :: adjointaction
      procedure(dump), deferred :: dump
      procedure(dealloc), deferred :: dealloc
   end type ContainerElementBase

   abstract interface

      !--------------------------------------------------------------------
      !> @brief
      !> multiplies this with arg from the right.
      !
      !> @param[in] this
      !--------------------------------------------------------------------
      subroutine rmultinterface(this, arg, t)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
         integer, intent(in) :: t
      end subroutine

      subroutine rmult1D2interface(this, arg, t)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
         integer, intent(in) :: t
      end subroutine

      !--------------------------------------------------------------------
      !> @brief
      !> multiplies this with arg from the left.
      !
      !> @param[in] this
      !--------------------------------------------------------------------
      subroutine lmultinterface(this, arg, t)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
         integer, intent(in) :: t
      end subroutine

      subroutine lmult1D2interface(this, arg, t)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
         integer, intent(in) :: t
      end subroutine

      !--------------------------------------------------------------------
      !> @brief
      !> multiplies this^-1 with arg from the right.
      !
      !> @param[in] this
      !--------------------------------------------------------------------
      subroutine rmultinvinterface(this, arg, t)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
         integer, intent(in) :: t
      end subroutine

      subroutine rmultinv1D2interface(this, arg, t)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
         integer, intent(in) :: t
      end subroutine

      !--------------------------------------------------------------------
      !> @brief
      !> multiplies this^-1 with arg from the left.
      !
      !> @param[in] this
      !--------------------------------------------------------------------
      subroutine lmultinvinterface(this, arg, t)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
         integer, intent(in) :: t
      end subroutine

      subroutine lmultinv1D2interface(this, arg, t)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
         integer, intent(in) :: t
      end subroutine

      !--------------------------------------------------------------------
      !> @brief
      !> This dumps the content to the screen.
      !
      !> @param[in] this
      !--------------------------------------------------------------------
      subroutine dump(this)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
      end subroutine

      !--------------------------------------------------------------------
      !> @brief
      !> Free the used memory
      !
      !> @param[in] this
      !--------------------------------------------------------------------
      subroutine dealloc(this)
         import ContainerElementBase
         class(ContainerElementBase), intent(inout) :: this
      end subroutine

      !--------------------------------------------------------------------
      !> @brief
      !> Perform the similarity transform e^{-T/2} arg e^{T/2}
      !
      !> @param[in] this
      !> @param[inout] the matrix that we intend to transform.
      !--------------------------------------------------------------------
      subroutine adjointactioninterface(this, arg, t1, t2)
         import ContainerElementBase
         class(ContainerElementBase), intent(in) :: this
         complex(kind=kind(0.d0)), intent(inout), dimension(:, :) :: arg
         integer, intent(in) :: t1
         integer, optional, intent(in) :: t2
      end subroutine
   end interface
end module ContainerElementBase_mod
