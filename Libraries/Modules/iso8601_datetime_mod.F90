module iso8601_datetime_mod
    implicit none
    contains
    character(len=29) function iso8601_datetime()
        ! Return current date and time in ISO 8601 format.
        ! Taken from https://cyber.dabamos.de/programming/modernfortran/date-and-time.html
        character(len=*), parameter :: ISO_FMT = &
            '(i4, 2("-", i2.2), "T", 2(i0.2, ":"), i0.2, ".", i0.3, a, ":", a)'
        character(len=5)  :: zone
        integer           :: dt(8)

        call date_and_time(values=dt, zone=zone)

        write (iso8601_datetime, ISO_FMT) dt(1), dt(2), dt(3), dt(5), dt(6), &
                                dt(7), dt(8), zone(1:3), zone(4:5)
    end function iso8601_datetime
end module iso8601_datetime_mod
