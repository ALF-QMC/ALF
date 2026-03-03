!  Copyright (C) 2025-2026 The ALF project
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

module iso8601_datetime_mod
!--------------------------------------------------------------------
!> @author ALF-project (adapted from https://cyber.dabamos.de/programming/modernfortran/date-and-time.html)
!> @brief ISO 8601 date/time formatting utility.
!
!> @details
!> Provides function to retrieve current system date and time formatted
!> according to ISO 8601 standard (YYYY-MM-DDTHH:MM:SS.mmm±HH:MM).
!>
!> Used for timestamping simulation outputs, log files, and data files
!> with unambiguous, internationally recognized datetime format.
!>
!> ISO 8601 advantages:
!> - Unambiguous: avoids MM/DD vs DD/MM confusion
!> - Sortable: lexicographic sort = chronological sort
!> - Machine-readable: widely supported standard
!> - Includes timezone information
!
!--------------------------------------------------------------------
    implicit none
    contains
    
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Get current date and time in ISO 8601 format.
!
!> @details
!> Returns system date/time as ISO 8601 string:
!>   YYYY-MM-DDTHH:MM:SS.mmm±HH:MM
!>
!> Format breakdown:
!> - YYYY: 4-digit year
!> - MM: 2-digit month (01-12)
!> - DD: 2-digit day (01-31)
!> - T: date/time separator
!> - HH: 2-digit hour (00-23)
!> - MM: 2-digit minute (00-59)
!> - SS: 2-digit second (00-59)
!> - mmm: 3-digit millisecond (000-999)
!> - ±HH:MM: timezone offset from UTC
!>
!> Example output: "2026-03-15T14:32:05.123+01:00"
!>
!> Uses Fortran intrinsic date_and_time() to query system clock.
!
!> @return ISO 8601 formatted datetime string (29 characters)
!
!> @note Timezone accuracy depends on system configuration
!
!> @see date_and_time (Fortran intrinsic)
!--------------------------------------------------------------------
    character(len=29) function iso8601_datetime()
        ! ISO 8601 format: YYYY-MM-DDTHH:MM:SS.mmm±HH:MM
        character(len=*), parameter :: ISO_FMT = &
            '(i4, 2("-", i2.2), "T", 2(i0.2, ":"), i0.2, ".", i0.3, a, ":", a)'
        character(len=5)  :: zone      ! Timezone: ±HHMM
        integer           :: dt(8)     ! Date/time values from system

        ! Query system date and time
        ! dt(1)=year, dt(2)=month, dt(3)=day, dt(5)=hour, 
        ! dt(6)=minute, dt(7)=second, dt(8)=millisecond
        call date_and_time(values=dt, zone=zone)

        ! Format as ISO 8601 string
        write (iso8601_datetime, ISO_FMT) dt(1), dt(2), dt(3), dt(5), dt(6), &
                                dt(7), dt(8), zone(1:3), zone(4:5)
    end function iso8601_datetime
end module iso8601_datetime_mod
