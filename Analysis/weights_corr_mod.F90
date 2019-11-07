Module Weights_corr
  use :: iso_fortran_env

  implicit none

contains

  integer function get_norb(filename)
    implicit none
    character(len=*), intent(in) :: filename

    real(kind=kind(0.d0)) :: X
    integer :: unit, ierr

    open(newunit=unit, file=filename, status='old', iostat=ierr)
    if (ierr /= 0) then
       write(error_unit,'(3A)') 'unable to open ', filename, ' for reading'
       stop
    end if
    read(unit,*) X, get_norb
    close(unit)
  end function get_norb
  
  function get_weights_corr(Norb) result(w)
    implicit none
    integer, intent(in) :: norb
    real (kind=kind(0.d0)), allocatable :: w(:,:)

    real, parameter :: SENTINEL=HUGE(0.)
    real (kind=kind(0.d0)), allocatable :: weights(:), weights_corr(:,:)
    integer :: unit, ierr, nread, nread_corr, i, j
    
    namelist /VAR_weights_corr/ weights, weights_corr

    allocate(weights(Norb), weights_corr(Norb, Norb), w(Norb, Norb))
    weights = SENTINEL
    weights_corr = SENTINEL
    open(newunit=unit, file='parameters', status='old', iostat=ierr)
    if (ierr /= 0) then
       write(error_unit, '(A)') 'unable to open <parameters>'
       stop
    end if
    read(unit, nml=VAR_weights_corr)
    close(unit)

    nread = 0
    do i = 1, Norb
       if (weights(i) /= SENTINEL) then
          nread = nread + 1
       end if
    end do
    if ((nread > 0).and.(nread < Norb)) then
       write(error_unit, '(A,I3,A,I3)') 'Error: read ', nread, ' values for weights, expected ', Norb
       stop
    end if

    nread_corr = 0
    do i = 1, Norb
       do j = 1, Norb
          if (weights_corr(j, i) /= SENTINEL) then
             nread_corr = nread_corr + 1
          end if
       end do
    end do
    if ((nread_corr > 0).and.(nread_corr < (Norb**2))) then
       write(error_unit, '(A,I3,A,I3)') 'Error: read ', nread_corr, ' values for weights_corr, expected ', Norb**2
       stop
    end if

    if (nread > 0) then
       ! weights_corr cannot be given if weights is present
       if (nread_corr > 0) then
          write(error_unit, '(A)') 'Error: weights and weights_corr cannot be provided at the same time'
          stop
       else
          ! Here weights is given and we construct weights_corr from weights
          do i = 1, Norb
             do j = 1, Norb
                w(j, i) = weights(i) * weights(j)
             end do
          end do
       end if ! if (nread_corr > 0)
    else if (nread_corr > 0) then
       ! Here weights_corr is given
       w = weights_corr
    else
       ! Here neither weights nor weights_corr is given
       ! By default, the result is the identity matrix
       w = 0.d0
       do i = 1, Norb
          w(i, i) = 1.d0
       end do
    end if ! if (nread > 0)
    
    deallocate(weights, weights_corr)
  end function get_weights_corr
end Module Weights_corr
