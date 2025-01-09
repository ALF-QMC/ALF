module test_qr
    implicit none
    contains
    subroutine zfill_matrix(m, n, A)
        integer :: m, n
        complex(kind=kind(0.d0)) :: A(:,:)
        
        integer :: i, j
        real*8 :: re, im
        
        do j = 1, n
            do i = 1, m
                call random_number( re )
                call random_number( im )
                A(i,j) = cmplx(re, im, kind=kind(0.d0))
            end do
        end do
    end subroutine
end module test_qr
program main
    use test_qr
    use QDRP_mod
    implicit none

    integer :: Ndim, N_part, nargs
    ! integer :: N_reps=10, i_rep
    integer :: Lwork
    complex(kind=kind(0.d0)), allocatable :: Mat(:,:), tau(:), work(:), D(:)
    integer, allocatable :: IPVT(:)
    CHARACTER(len=32) :: arg

    nargs = COMMAND_ARGUMENT_COUNT()
    if (nargs > 0) then
        CALL GET_COMMAND_ARGUMENT(1, arg)
        read(arg,*) Ndim
    else
        Ndim = 2
    endif
    if (nargs > 1) then
        CALL GET_COMMAND_ARGUMENT(2, arg)
        read(arg,*) N_part
    else
        N_part = Ndim
    endif

    print*, Ndim, N_part


    allocate(mat(Ndim, N_part), D(N_part), IPVT(N_part), tau(N_part))

    OPEN (UNIT = 10, FILE='mat', STATUS='UNKNOWN', ACTION='READ')
    read(10, *) Mat
    close(10)
    OPEN (UNIT = 10, FILE='IPVT', STATUS='UNKNOWN', ACTION='READ')
    read(10, *) IPVT
    close(10)

    ! do i_rep=1, N_reps
    !    print*, "fill matrix", i_rep
    !    call zfill_matrix(Ndim, N_part, mat)
    !    IPVT = 0

        print*, 'decompose'
        call QDRP_decompose(Ndim, N_part, Mat, D, IPVT, TAU, WORK, LWORK)
        deallocate(work)
    ! enddo
End program
