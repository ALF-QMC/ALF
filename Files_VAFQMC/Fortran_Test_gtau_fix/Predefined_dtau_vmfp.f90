module Predefined_dtau_vmfp
    implicit none
    
contains

    ! Subroutine to calculate hi, ti, and lambdat
    subroutine time_dependent_int(beta, U, dtau, n, hi, ti, lambdai)
        implicit none
        integer :: i, itr, maxitr
        integer, intent(in) :: n
        Real (Kind=Kind(0.d0)), intent(in) :: beta, U, dtau
        Real (Kind=Kind(0.d0)), allocatable :: hi(:), ti(:), lambdai(:)
        Real (Kind=Kind(0.d0)) :: gamma, gmin, gmax, fun_val
        Real (Kind=Kind(0.d0)) :: cost, eps, error
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
        !print*, ('t=',ti(i),i=1,n+1)
        !print*, ('h=',hi(i),i=1,n)
        !print*, ('l=',lambdai(i),i=1,n)
    end subroutine time_dependent_int

    ! Function to fix constraint on time-dependent coupling
    Real (Kind=Kind(0.d0)) function fun(p, n)
        integer :: n
        Real (Kind=Kind(0.d0)) :: p
        Real (Kind=Kind(0.d0)) :: cost
        cost = p**(n-1)
        fun = (1.d0 - p * cost) / (1.d0 - p) / cost
    end function fun

end module Predefined_dtau_vmfp
