! beta is inverse temperature
! dtau is imaginary time
! U is hubbard U
! n is the number of time slices
! t_i i time dependent hopping
! h_i is time dpendent interactions
! time lambadi is time dependent interactions

program main
    implicit none
    integer :: n, iter, maxitr,i
    real*8 :: beta, U, dtau, cost, eps, error, gamma
    real*8, dimension(:), allocatable :: hi, ti, lambdai
    ! Input paramers
    beta = 4.0
    U = 2.0
    dtau = 0.1
    n = 5
    ! Call the subroutine to calculate hi, ti, and lambdat
    call time_dependent_int(beta, U, dtau, n, hi, ti, lambdai)
    !print *, 'Optimized hi:'
    do i = 1, n
        print*, i, ti(i), hi(i), lambdai(i), dtau
    end do
contains

    ! Subroutine to calculate hi, ti, and lambdat
    subroutine time_dependent_int(beta, U, dtau, n, hi, ti, lambdai)
        implicit none
        real*8, intent(in) :: beta, U, dtau
        integer, intent(in) :: n
        real*8, dimension(:), allocatable, intent(out) :: hi, ti, lambdai
        real*8 :: gamma,gmin,gmax,fun_val
        ! Local variables
        integer :: i, itr, maxitr
        real*8 :: cost, eps, error
        !real, external:: fun

        allocate(hi(n), ti(n+1), lambdai(n))
        ! Calculate threshhold
        cost = beta / dtau / 2.0
        eps = 1.0e-14
        maxitr = 100
        error = 1.0

        if (cost > real(n)) then
            gmin = 0.0
            gmax = 1.0
            itr = 0 
            do while (error > eps .and. itr < maxitr)
                gamma = (gmin + gmax) / 2.0
                fun_val= fun(gamma, n)
                if (fun_val > cost) then
                    gmin = gamma
                else
                    gmax = gamma
                end if
                error = abs(fun_val - cost)
                itr = itr + 1
                print *, ' itr,  error =', itr, error, gmin, gmax
            end do
        
            ! Calculate hi
            hi(n) = dtau
            do i = n - 1, 1, -1
                hi(i) = hi(i + 1) / gamma
            end do
        else
            gamma = 1.0
            hi = beta / real(n) / 2.0
        endif

        do i=1,n
            !print*,'i,hi(i)=',i,hi(i)
            lambdai(i) = sqrt(U * hi(i))
        enddo

        ti(1) = hi(1) / 2.0
        do i = 1, n-1
            ti(i+1) = (hi(i+1) + hi(i)) / 2.0
        end do
        ti(n + 1) = hi(n)
        !do i = 1, n
        !   print*, i, ti(i), hi(i), lambdai(i), dtau
        ! end do
        print*, ('t=',ti(i),i=1,n+1)
        print*, ('h=',hi(i),i=1,n)
        print*, ('l=',lambdai(i),i=1,n)
    end subroutine time_dependent_int


    ! Function to fix constraint on time-dependent coupling
    real*8 function fun(p, n)
        integer :: n
        real*8 :: p
        cost=p**(n-1)
        fun=(1.d0-p*cost)/(1.d0-p)/cost
    end function fun
end program main
