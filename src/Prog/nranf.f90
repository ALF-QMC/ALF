    integer function nranf(N)
      Use Random_wrap
      implicit none
      integer :: N
      
      nranf  = nint(ranf()*dble(N) + 0.5)

      if (nranf .lt. 1 ) nranf = 1
      if (nranf .gt. N ) nranf = N 

    end function nranf

