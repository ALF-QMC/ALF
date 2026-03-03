!  Copyright (C) 2018-2026 The ALF project
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

!===============================================================================
   MODULE precdef
!--------------------------------------------------------------------
!> @author ALF-project
!> @brief Precision and numeric type definitions for the ALF project.
!>
!> @details
!> This module centralizes numeric kinds and constants to ensure
!> consistency and portability across platforms.
!>
!> **Integer Kinds:**
!> - byte: 1-byte signed integer (-128 to 127)
!> - long: 4-byte signed integer
!>
!> **Floating Point Precision:**
!> - single: 4-byte (single precision), ~6 digits
!> - double: 8-byte (double precision), ~15 digits (REAL64)
!>
!> **Numeric Constants:**
!> - Real: rone (1.0) and rzero (0.0)
!> - Complex: cone (1.0+0i) and czero (0.0+0i)
!>
!> All definitions use ISO_FORTRAN_ENV for portability.
!--------------------------------------------------------------------
   use, intrinsic :: iso_fortran_env
   IMPLICIT NONE

   INTEGER, PARAMETER :: &
      byte   = selected_int_kind(2),          & ! -128 ... 127, 1 byte
      long   = selected_int_kind(9),          & ! −2147483648 ... 2147483647, 4 byte
!      int64  = selected_int_kind(18),         & ! −9223372036854775808 ... 9223372036854775807 8 byte
      single = selected_real_kind(p=6,r=37),  & ! kind(1.0), 4 byte
      !double = selected_real_kind(p=15,r=307)   ! selected_real_kind(2*precision(1.0_double)), 8 byte
      double = REAL64   ! selected_real_kind(2*precision(1.0_double)), 8 byte
      
   REAL(kind=Kind(0.d0)), PARAMETER :: &
      rone = 1.0D0, &
      rzero = 0.0D0

   COMPLEX(kind=Kind(0.d0)), PARAMETER :: &
      cone = cmplx(rone,rzero,REAL64), &
      czero = cmplx(rzero,rzero,REAL64)

   END MODULE precdef
