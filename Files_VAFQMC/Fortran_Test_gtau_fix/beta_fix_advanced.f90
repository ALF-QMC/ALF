program main
    implicit none
    integer :: n,nt,i
    real*8 :: beta, Ur, dtau, gamma
    real*8 :: betan, dtaun, gamman
    real*8, dimension(:), allocatable :: hi, ti, lambdai
    real*8, dimension(:), allocatable :: hin, tin, lambdain
    logical :: hirsch

    beta = 4.0
    Ur = 2.0
    dtau = 0.1
    n = 5
    nt = 2 * n + 1
    hirsch =.false.

    call  time_dependent_int(n, beta, Ur, dtau, hi, gamma, lambdai, ti, hirsch)

    print*, 't=',(ti(i),i=1, nt)
    print*, 'h=',(hi(i),i=1,n)
    print*, 'l=',(lambdai(i),i=1, nt)

    tin=  ti
    hin =  hi
    lambdain =  lambdai
    call time_dependent_int_n(Ur, n, beta, betan, dtau, dtaun, hi, hin, gamma,lambdai, lambdain, ti, tin, hirsch)
    print*, 'beta, betan, gamma =',beta, betan, gamma
    !print*, 'tn=',(tin(i),i=1, nt)
    !print*, 'hn=',(hin(i),i=1,n)
    !print*, 'ln=',(lambdain(i),i=1, nt)

    call fixfirst(nt, beta, Ur, lambdai, hirsch)
    print*, 'lc=',(lambdai(i),i=1, nt)
contains

subroutine time_dependent_int(n, beta, Ur, dtau, hi, gamma, lambdai, ti, hirsch)
    implicit none
    real*8, intent(in) :: beta, Ur, dtau
    real*8, intent(out) :: gamma
    integer, intent(in) :: n
    integer :: i,  nt, maxitr, itr
    real*8 ::  U, error, eps, cost, gmin, gmax, fun_val, costu
    logical, intent(in) :: hirsch
    real*8, dimension(:), allocatable, intent(out) :: hi, ti, lambdai
    !real, external:: fun
    allocate(hi(n))
    allocate(ti(2*n+1), lambdai(2*n+1))
    
    U = abs(Ur)
    nt = 2 * n + 1
    cost = beta / dtau / 2.d0
    error = 1.d0
    eps = 1d-14
    maxitr = 100
    
    if (cost > real(n)) then
        gmin = 0.d0
        gmax = 1.d0
    elseif (cost > 1.d0) then
        gmin = 1.d0
        gmax = cost / (cost - 1.d0)
    else
        gamma = 1.d0
        error = 0.d0
    endif
    itr = 0
    do while (error > eps .and. itr < maxitr)
        gamma = (gmin + gmax) / 2.d0
        fun_val = fun(gamma, n)
        if (fun_val > cost) then
            gmin = gamma
        else
            gmax = gamma
        endif
        error = abs(fun_val - cost)
        itr = itr + 1
!        write(*, *) ' itr  error =', itr, error, pmin, pmax
    end do

    hi(n) = dtau
    do i = n - 1, 1, -1
        hi(i) = hi(i + 1) / gamma
    end do

    lambdai(nt) = 0.d0
    if (hirsch) then
        do i = 1, n
            costu = exp(Ur * hi(i) / 2.d0)
            lambdai(i) = log(costu + sqrt(costu**2 - 1.d0))
            lambdai(nt - i) = lambdai(i)
        end do
    else
        do i = 1, n
            lambdai(i) = sqrt(U * hi(i))
            lambdai(nt - i) = lambdai(i)
        end do
    endif
    ti(1) = hi(1) / 2.d0
    do i = 1, n - 1
        ti(i + 1) = (hi(i + 1) + hi(i)) / 2.d0
    end do
    ti(n + 1) = hi(n)
    do i = 1, n
        ti(nt - i + 1) = ti(i)
    end do
end subroutine




subroutine time_dependent_int_n(Ur, n, beta, betan, dtau, dtaun, hi, hin, gamma, &
    lambdai, lambdain, ti, tin, hirsch)
    implicit none
    real*8, intent(in) :: beta, Ur, dtau, gamma
    integer, intent(in) :: n
    real*8, intent(out) :: betan
    real*8, dimension(n), intent(in) :: hi
    real*8, dimension(2*n+1), intent(in) :: ti, lambdai
    real*8, dimension(n), intent(out) :: hin
    real*8, dimension(2*n+1), intent(out) :: tin, lambdain
    integer :: i, nt
    real*8 :: dtaun, gamman, U
    real*8 :: cost, costn, der, costu, costun, cost_der
    logical, intent(in) :: hirsch

    U = abs(Ur)
    cost = beta / dtau / 2.d0
    nt = 2 * n + 1
    hin = 0.d0    
    gamman = 0.d0 

    do i = n, 1, -1
        tin(i) = tin(i) + tin(nt - i + 1)
        tin(nt - i + 1) = 0.d0
    end do

    hin(n) = hin(n) + tin(n + 1)
    tin(n + 1) = 0.d0

    do i = n - 1, 1, -1
        hin(i) = hin(i) + tin(i + 1) / 2.d0
        hin(i + 1) = hin(i + 1) + tin(i + 1) / 2.d0
        tin(i + 1) = 0.d0
    end do

    hin(1) = hin(1) + tin(1) / 2.d0
    tin(1) = 0.d0

    if (hirsch) then
        do i = n, 1, -1
            lambdain(i) = lambdain(i) + lambdain(nt - i)
            lambdain(nt - i) = 0.d0
            costu = exp(U * hi(i) / 2.d0)
            costun = 1.d0 / sqrt(costu**2 - 1.d0) * lambdain(i)
            hin(i) = hin(i) + U / 2.d0 * costu * costun
            lambdain(i) = 0.d0
        end do
    else
        do i = n, 1, -1
            lambdain(i) = lambdain(i) + lambdain(nt - i)
            lambdain(nt - i) = 0.d0
            hin(i) = hin(i) + 0.5d0 / lambdai(i) * U * lambdain(i)
            lambdain(i) = 0.d0
        end do
    endif 

    lambdain(nt) = 0.d0
    do i = 1, n - 1
        hin(i + 1) = hin(i + 1) + hin(i) / gamma
        gamman = gamman - hi(i + 1) / gamma**2 * hin(i)
        hin(i) = 0.d0
    end do
    dtaun = dtaun + hin(n)
    hin(n) = 0.d0

    ! Reverse of input cost output gamma der = dcost/dp|gamma
    cost_der = (n - 1.d0 - n * gamma + gamma**n)
    if (gamma .ne. 1.d0 .and. cost_der .ne. 0.d0) then
        der = -cost_der / (gamma - 1.d0)**2 / gamma**n
    else
        der = -((n - 1) * n) / 2.d0
    endif
    costn = gamman / der
    gamman = 0.d0
    ! Reverse of cost = beta / dtmin / 2.d0
    dtaun = dtaun - cost / dtau * costn
    betan = betan + costn / dtau / 2.d0
    costn = 0.d0
    print*, 'dtaun, betan, costn,cost, gamman, der =', dtaun, betan, costn, cost, gamman, der
    print*, 'tb1=', (tin(i), i=1, nt)
    print*, 'hb1=', (hin(i), i=1, n)
    print*, 'lb1=', (lambdain(i), i=1, nt)
end subroutine



subroutine fixfirst(nt, beta, Ur, lambdai, hirsch)
    implicit none
    integer :: nt, i
    real*8, intent(in) :: beta, Ur
    logical, intent(in) :: hirsch
    real*8 ::  dtt, cost, costi, lambdai(*)

    if (hirsch) then
        cost = beta
        do i = 2, (nt - 1) / 2
            dtt = exp(lambdai(i))
            costi = 2.d0 * log((1.d0 + dtt**2) / 2.d0 / dtt)
            cost = cost - costi / Ur
        end do
        if (cost > 0.d0) then
            cost = exp(cost * Ur / 2.d0)
            dtt = cost + sqrt(cost**2 - 1.d0)
            lambdai(1) = log(dtt)
        else
            write(*, *) ' ERROR beta inconsistency !!! '
            lambdai(1) = 0.d0
        endif
    else
        cost = beta
        do i = 2, (nt - 1) / 2
            cost = cost - lambdai(i)**2 / Ur
        end do
        if (cost > 0.d0) then
            lambdai(1) = sqrt(cost * Ur)
        else
            write(*, *) ' ERROR beta inconsistency !!! '
            lambdai(1) = 0.d0
        endif
    endif
    do i = 1, (nt - 1) / 2
        lambdai(nt - i) = lambdai(i)
    end do
end subroutine

function fun(p, n)
    implicit none
    real*8 :: fun, p, cost
    integer :: n 
    cost = p**(n - 1)
    fun = (1.d0 - p * cost) / (1.d0 - p) / cost
    return
    end function
end program main