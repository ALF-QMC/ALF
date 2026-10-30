!  Copyright (C) 2026 The ALF project
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

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Benchmarks the BLAS/LAPACK wrappers of alf_linalg_mod on the CPU against
!> offloading to the GPU and determines the offload thresholds.
!>
!> Usage: gpu_linalg_bench.out [NMAX [ROUTINE ...]]
!>
!> For every routine and size N = 16, 23, 32, 45, ..., NMAX (default 4096)
!> a square problem of size N is timed on the CPU (threshold -1) and on the
!> GPU (threshold 0), including the host-device transfers of the call.
!> The GPU result is compared to the CPU result.
!> The recommended threshold is the smallest N above which the GPU is faster
!> for all tested sizes.
!--------------------------------------------------------------------
program gpu_linalg_bench

  use alf_linalg_mod, only: alf_linalg_init, alf_linalg_set_threshold, alf_linalg_offload_active
  use runtime_error_mod
  use iso_fortran_env, only: int64, output_unit
  implicit none

  integer, parameter :: dp = kind(0.d0)
  ! The wrappers are external procedures like BLAS/LAPACK
  real(dp), external :: ALF_DZNRM2
  integer, parameter :: N_ALL = 25
  character(len=6), parameter :: all_names(N_ALL) = [character(len=6) :: &
     & 'ZGEMM',  'DGEMM',  'ZHEMM',  'ZTRSM',  'ZTRMM',  'ZGEMV',  'DGEMV',  'ZGERU', &
     & 'ZSCAL',  'ZCOPY',  'ZSWAP',  'DZNRM2', 'ZGETRF', 'DGETRF', 'ZGETRS', 'ZGETRI', &
     & 'DGETRI', 'ZGEQRF', 'DGEQRF', 'ZUNGQR', 'DORGQR', 'ZUNMQR', 'ZHEEV',  'DSYEV', &
     & 'ZGESVD']
  ! Minimal total time per measurement, minimal and maximal number of repetitions
  real(dp), parameter :: t_min = 0.2_dp
  integer,  parameter :: rep_min = 3, rep_max = 1000
  ! Relative deviation of GPU from CPU result that is reported as failure
  real(dp), parameter :: tol = 1.d-8

  ! Input data, copied to the work arrays before every call
  complex(dp), allocatable :: za0(:,:), zb0(:,:), zc0(:,:), zh0(:,:), zlu0(:,:), zqr0(:,:), ztau0(:)
  real(dp),    allocatable :: da0(:,:), db0(:,:), dc0(:,:), ds0(:,:), dqr0(:,:), dtau0(:)
  integer,     allocatable :: ipiv0(:)
  ! Work arrays
  complex(dp), allocatable :: za(:,:), zb(:,:), zc(:,:), ztau(:), zwork(:), zu(:,:), zvt(:,:)
  real(dp),    allocatable :: da(:,:), db(:,:), dc(:,:), dtau(:), dwork(:), rwork(:), w(:)
  integer,     allocatable :: ipiv(:)
  ! Result of the last call, for comparing CPU and GPU
  complex(dp), allocatable :: res_cpu(:), res_gpu(:)

  character(len=6), allocatable :: routines(:)
  character(len=32) :: arg
  integer :: nmax, nsizes, ir, is, n, crossover, nfail
  integer, allocatable :: sizes(:)
  real(dp), allocatable :: t_cpu(:), t_gpu(:)
  real(dp) :: err

  nmax = 4096
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    read(arg, *) nmax
  endif
  if (command_argument_count() >= 2) then
    allocate(routines(command_argument_count() - 1))
    do ir = 1, size(routines)
      call get_command_argument(ir + 1, arg)
      routines(ir) = upcase(trim(arg))
    enddo
  else
    allocate(routines(N_ALL))
    routines = all_names
  endif

  ! Sizes 16 * sqrt(2)**k up to nmax
  nsizes = 0
  do while (nint(16*sqrt(2._dp)**nsizes) <= nmax)
    nsizes = nsizes + 1
  enddo
  allocate(sizes(nsizes), t_cpu(nsizes), t_gpu(nsizes))
  do is = 1, nsizes
    sizes(is) = nint(16*sqrt(2._dp)**(is - 1))
  enddo

  call alf_linalg_init(0)
  if (.not. alf_linalg_offload_active()) then
    write(output_unit, '(A)') ' GPU offloading is not active, compile with configure.sh option GPU.'
    write(output_unit, '(A)') ' Only CPU timings are measured.'
  endif

  nfail = 0
  do ir = 1, size(routines)
    write(output_unit, '(/,A)') ' ' // trim(routines(ir))
    write(output_unit, '(A8,2A14,A10,A12)') 'N', 't_cpu [s]', 't_gpu [s]', 'speedup', 'rel. dev.'
    do is = 1, nsizes
      n = sizes(is)
      call setup(n)
      t_cpu(is) = measure(routines(ir), n, -1)
      res_cpu = result(routines(ir))
      if (alf_linalg_offload_active()) then
        t_gpu(is) = measure(routines(ir), n, 0)
        res_gpu = result(routines(ir))
        err = maxval(abs(res_gpu - res_cpu)) / max(maxval(abs(res_cpu)), tiny(1._dp))
      else
        t_gpu(is) = huge(1._dp)
        err = 0
      endif
      if (err > tol) nfail = nfail + 1
      write(output_unit, '(I8,2ES14.4,F10.3,ES12.3,A)') n, t_cpu(is), t_gpu(is), t_cpu(is)/t_gpu(is), &
         & err, merge('  FAILED', '        ', err > tol)
      flush(output_unit)
    enddo
    ! Smallest size above which the GPU is always faster
    crossover = -1
    do is = nsizes, 1, -1
      if (t_gpu(is) >= t_cpu(is)) exit
      crossover = sizes(is)
    enddo
    write(output_unit, '(A,A,A,I0)') ' Recommended: export ALF_GPU_THRESH_', trim(routines(ir)), '=', crossover
  enddo

  if (nfail > 0) then
    write(output_unit, '(/,A,I0,A)') ' ', nfail, ' GPU results deviated from the CPU results!'
    call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  endif

contains

  function upcase(s) result(u)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: u
    integer :: i
    u = s
    do i = 1, len(s)
      if (s(i:i) >= 'a' .and. s(i:i) <= 'z') u(i:i) = achar(iachar(s(i:i)) - 32)
    enddo
  end function upcase

  real(dp) function wtime()
    integer(int64) :: count, rate
    call system_clock(count, rate)
    wtime = real(count, dp) / real(rate, dp)
  end function wtime

  !> Allocates and fills the input data for size n.
  subroutine setup(n)
    integer, intent(in) :: n
    real(dp), allocatable :: x(:,:), y(:,:)
    integer :: i, info

    if (allocated(za0)) deallocate(za0, zb0, zc0, zh0, zlu0, zqr0, ztau0, da0, db0, dc0, ds0, dqr0, &
       & dtau0, ipiv0, za, zb, zc, ztau, zwork, zu, zvt, da, db, dc, dtau, dwork, rwork, w, ipiv)
    allocate(x(n,n), y(n,n))
    call random_number(x); call random_number(y)
    za0 = cmplx(x - 0.5_dp, y - 0.5_dp, dp)
    call random_number(x); call random_number(y)
    zb0 = cmplx(x - 0.5_dp, y - 0.5_dp, dp)
    call random_number(x); call random_number(y)
    zc0 = cmplx(x - 0.5_dp, y - 0.5_dp, dp)
    ! Well conditioned matrices for factorizations and triangular solves
    do i = 1, n
      za0(i,i) = za0(i,i) + n
    enddo
    zh0 = za0 + conjg(transpose(za0))
    call random_number(x); da0 = x - 0.5_dp
    call random_number(x); db0 = x - 0.5_dp
    call random_number(x); dc0 = x - 0.5_dp
    do i = 1, n
      da0(i,i) = da0(i,i) + n
    enddo
    ds0 = da0 + transpose(da0)

    allocate(za(n,n), zb(n,n), zc(n,n), ztau(n), zu(n,n), zvt(n,n), da(n,n), db(n,n), dc(n,n), &
       & dtau(n), ipiv(n), w(n), rwork(5*n), ipiv0(n), ztau0(n), dtau0(n))
    ! Workspace large enough for all routines on the CPU
    allocate(zwork(64*n), dwork(64*n))
    ! LU and QR factors as input for GETRS, GETRI, UNGQR, ORGQR and UNMQR, computed on the CPU
    call alf_linalg_set_threshold('ZGETRF', -1)
    call alf_linalg_set_threshold('ZGEQRF', -1)
    call alf_linalg_set_threshold('DGEQRF', -1)
    zlu0 = za0
    call ALF_ZGETRF(n, n, zlu0, n, ipiv0, info)
    zqr0 = za0
    call ALF_ZGEQRF(n, n, zqr0, n, ztau0, zwork, size(zwork), info)
    dqr0 = da0
    call ALF_DGEQRF(n, n, dqr0, n, dtau0, dwork, size(dwork), info)
  end subroutine setup

  !> Median time of a call of routine name with threshold thr (-1: CPU, 0: GPU).
  real(dp) function measure(name, n, thr)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n, thr
    real(dp) :: times(rep_max), t_total
    integer :: nrep, i
    call alf_linalg_set_threshold(name, thr)
    ! Warm-up
    call reset()
    times(1) = call_routine(name, n)
    nrep = 0
    t_total = 0
    do while (nrep < rep_max .and. (nrep < rep_min .or. t_total < t_min))
      nrep = nrep + 1
      call reset()
      times(nrep) = call_routine(name, n)
      t_total = t_total + times(nrep)
    enddo
    ! Median by sorting
    do i = 2, nrep
      call insert(times, i)
    enddo
    measure = times((nrep + 1)/2)
    call alf_linalg_set_threshold(name, -1)
  end function measure

  subroutine insert(t, i)
    real(dp), intent(inout) :: t(:)
    integer, intent(in) :: i
    real(dp) :: x
    integer :: j
    x = t(i)
    j = i - 1
    do while (j >= 1)
      if (t(j) <= x) exit
      t(j+1) = t(j)
      j = j - 1
    enddo
    t(j+1) = x
  end subroutine insert

  !> Copies the input data to the work arrays.
  subroutine reset()
    za = za0; zb = zb0; zc = zc0
    da = da0; db = db0; dc = dc0
    ipiv = ipiv0
  end subroutine reset

  !> Calls the wrapper of routine name and returns the elapsed time.
  real(dp) function call_routine(name, n) result(t)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    complex(dp), parameter :: z1 = (1._dp, 0._dp), zh = (0.5_dp, 0.25_dp)
    real(dp) :: t0, nrm
    integer :: info

    select case (name)
    case ('ZHEEV')
      za = zh0
    case ('DSYEV')
      da = ds0
    case ('ZGETRS', 'ZGETRI')
      za = zlu0
    case ('DGETRI')
      call ALF_DGETRF(n, n, da, n, ipiv, info)
    case ('ZUNGQR', 'ZUNMQR')
      za = zqr0
      ztau = ztau0
    case ('DORGQR')
      da = dqr0
      dtau = dtau0
    end select

    t0 = wtime()
    select case (name)
    case ('ZGEMM');  call ALF_ZGEMM('N', 'N', n, n, n, z1, za, n, zb, n, zh, zc, n)
    case ('DGEMM');  call ALF_DGEMM('N', 'N', n, n, n, 1._dp, da, n, db, n, 0.5_dp, dc, n)
    case ('ZHEMM');  call ALF_ZHEMM('L', 'U', n, n, z1, za, n, zb, n, zh, zc, n)
    case ('ZTRSM');  call ALF_ZTRSM('L', 'U', 'N', 'N', n, n, z1, za, n, zb, n)
    case ('ZTRMM');  call ALF_ZTRMM('L', 'U', 'N', 'N', n, n, z1, za, n, zb, n)
    case ('ZGEMV');  call ALF_ZGEMV('N', n, n, z1, za, n, zb, 1, zh, zc, 1)
    case ('DGEMV');  call ALF_DGEMV('N', n, n, 1._dp, da, n, db, 1, 0.5_dp, dc, 1)
    case ('ZGERU');  call ALF_ZGERU(n, n, z1, zb, 1, zc, 1, za, n)
    case ('ZSCAL');  call ALF_ZSCAL(n, zh, za, 1)
    case ('ZCOPY');  call ALF_ZCOPY(n, za, 1, zb, 1)
    case ('ZSWAP');  call ALF_ZSWAP(n, za, 1, zb, 1)
    case ('DZNRM2'); nrm = ALF_DZNRM2(n, za, 1); za(1,1) = nrm
    case ('ZGETRF'); call ALF_ZGETRF(n, n, za, n, ipiv, info)
    case ('DGETRF'); call ALF_DGETRF(n, n, da, n, ipiv, info)
    case ('ZGETRS'); call ALF_ZGETRS('N', n, n, za, n, ipiv, zb, n, info)
    case ('ZGETRI'); call ALF_ZGETRI(n, za, n, ipiv, zwork, size(zwork), info)
    case ('DGETRI'); call ALF_DGETRI(n, da, n, ipiv, dwork, size(dwork), info)
    case ('ZGEQRF'); call ALF_ZGEQRF(n, n, za, n, ztau, zwork, size(zwork), info)
    case ('DGEQRF'); call ALF_DGEQRF(n, n, da, n, dtau, dwork, size(dwork), info)
    case ('ZUNGQR'); call ALF_ZUNGQR(n, n, n, za, n, ztau, zwork, size(zwork), info)
    case ('DORGQR'); call ALF_DORGQR(n, n, n, da, n, dtau, dwork, size(dwork), info)
    case ('ZUNMQR'); call ALF_ZUNMQR('L', 'C', n, n, n, za, n, ztau, zb, n, zwork, size(zwork), info)
    case ('ZHEEV');  call ALF_ZHEEV('V', 'U', n, za, n, w, zwork, size(zwork), rwork, info)
    case ('DSYEV');  call ALF_DSYEV('V', 'U', n, da, n, w, dwork, size(dwork), info)
    case ('ZGESVD'); call ALF_ZGESVD('A', 'A', n, n, za, n, w, zu, n, zvt, n, zwork, size(zwork), rwork, info)
    case default
      write(output_unit, *) 'Unknown routine ', name
      call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
    end select
    t = wtime() - t0
  end function call_routine

  !> Output of the last call of routine name, independent of sign conventions.
  !> Copies via the helpers below, since large temporaries overflow the stack.
  function result(name) result(r)
    character(len=*), intent(in) :: name
    complex(dp), allocatable :: r(:)
    select case (name)
    case ('ZGEMM', 'ZHEMM', 'ZGEMV');                   r = flat_z(zc)
    case ('DGEMM', 'DGEMV');                            r = flat_d(dc)
    case ('ZTRSM', 'ZTRMM', 'ZGETRS', 'ZUNMQR', 'ZCOPY'); r = flat_z(zb)
    case ('ZSWAP')
      allocate(r(size(za) + size(zb)))
      r(1:size(za)) = flat_z(za)
      r(size(za)+1:) = flat_z(zb)
    case ('ZGERU', 'ZSCAL', 'DZNRM2', 'ZGETRF', 'ZGETRI', 'ZUNGQR')
                                                        r = flat_z(za)
    case ('DGETRF', 'DGETRI', 'DORGQR');                r = flat_d(da)
    ! Householder vectors may differ in sign, compare the absolute values of R
    case ('ZGEQRF');                                    r = upper_abs_z(za)
    case ('DGEQRF');                                    r = upper_abs_d(da)
    ! Eigen- and singular vectors are not unique, compare the values only
    case ('ZHEEV', 'DSYEV', 'ZGESVD');                  r = cmplx(w, 0._dp, dp)
    case default;                                       r = [(0._dp, 0._dp)]
    end select
  end function result

  function flat_z(a) result(r)
    complex(dp), intent(in) :: a(:,:)
    complex(dp), allocatable :: r(:)
    integer :: j
    allocate(r(size(a)))
    do j = 1, size(a,2)
      r((j-1)*size(a,1)+1 : j*size(a,1)) = a(:,j)
    enddo
  end function flat_z

  function flat_d(a) result(r)
    real(dp), intent(in) :: a(:,:)
    complex(dp), allocatable :: r(:)
    integer :: i, j
    allocate(r(size(a)))
    do j = 1, size(a,2)
      do i = 1, size(a,1)
        r(i + (j-1)*size(a,1)) = cmplx(a(i,j), 0._dp, dp)
      enddo
    enddo
  end function flat_d

  !> Absolute values of the upper triangle of a, flattened.
  function upper_abs_z(a) result(r)
    complex(dp), intent(in) :: a(:,:)
    complex(dp), allocatable :: r(:)
    integer :: i, j
    allocate(r(size(a)))
    r = 0
    do j = 1, size(a,2)
      do i = 1, min(j, size(a,1))
        r(i + (j-1)*size(a,1)) = abs(a(i,j))
      enddo
    enddo
  end function upper_abs_z

  function upper_abs_d(a) result(r)
    real(dp), intent(in) :: a(:,:)
    complex(dp), allocatable :: r(:)
    integer :: i, j
    allocate(r(size(a)))
    r = 0
    do j = 1, size(a,2)
      do i = 1, min(j, size(a,1))
        r(i + (j-1)*size(a,1)) = abs(a(i,j))
      enddo
    enddo
  end function upper_abs_d

end program gpu_linalg_bench
