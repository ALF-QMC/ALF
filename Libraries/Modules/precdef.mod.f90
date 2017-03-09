!===============================================================================
   MODULE precdef
!-------------------------------------------------------------------------------
!
   IMPLICIT NONE

   INTEGER, PARAMETER :: &
      byte   = selected_int_kind(2),          & ! -128 ... 127, 1 byte
      long   = selected_int_kind(9),          & ! −2147483648 ... 2147483647, 4 byte
      int64  = selected_int_kind(18),         & ! −9223372036854775808 ... 9223372036854775807 8 byte
      single = selected_real_kind(p=6,r=37),  & ! kind(1.0), 4 byte
      !double = selected_real_kind(p=15,r=307)   ! selected_real_kind(2*precision(1.0_double)), 8 byte
      double = 8   ! selected_real_kind(2*precision(1.0_double)), 8 byte
      
   REAL(kind=Kind(0.d0)), PARAMETER :: &
      rone = 1.0D0, &
      rzero = 0.0D0

   COMPLEX(kind=Kind(0.d0)), PARAMETER :: &
      cone = cmplx(rone,rzero,double), &
      czero = cmplx(rzero,rzero,double)

   END MODULE precdef
