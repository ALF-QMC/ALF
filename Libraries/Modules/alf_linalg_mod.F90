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

! Compile time defaults of the offload thresholds, see module documentation.
! A negative value means that the routine is never offloaded.
#ifndef ALF_GPU_THRESH_ZGEMM
#define ALF_GPU_THRESH_ZGEMM -1
#endif
#ifndef ALF_GPU_THRESH_DGEMM
#define ALF_GPU_THRESH_DGEMM -1
#endif
#ifndef ALF_GPU_THRESH_ZHEMM
#define ALF_GPU_THRESH_ZHEMM -1
#endif
#ifndef ALF_GPU_THRESH_ZTRSM
#define ALF_GPU_THRESH_ZTRSM -1
#endif
#ifndef ALF_GPU_THRESH_ZTRMM
#define ALF_GPU_THRESH_ZTRMM -1
#endif
#ifndef ALF_GPU_THRESH_ZGEMV
#define ALF_GPU_THRESH_ZGEMV -1
#endif
#ifndef ALF_GPU_THRESH_DGEMV
#define ALF_GPU_THRESH_DGEMV -1
#endif
#ifndef ALF_GPU_THRESH_ZGERU
#define ALF_GPU_THRESH_ZGERU -1
#endif
#ifndef ALF_GPU_THRESH_ZSCAL
#define ALF_GPU_THRESH_ZSCAL -1
#endif
#ifndef ALF_GPU_THRESH_ZCOPY
#define ALF_GPU_THRESH_ZCOPY -1
#endif
#ifndef ALF_GPU_THRESH_ZSWAP
#define ALF_GPU_THRESH_ZSWAP -1
#endif
#ifndef ALF_GPU_THRESH_DZNRM2
#define ALF_GPU_THRESH_DZNRM2 -1
#endif
#ifndef ALF_GPU_THRESH_ZGETRF
#define ALF_GPU_THRESH_ZGETRF -1
#endif
#ifndef ALF_GPU_THRESH_DGETRF
#define ALF_GPU_THRESH_DGETRF -1
#endif
#ifndef ALF_GPU_THRESH_ZGETRS
#define ALF_GPU_THRESH_ZGETRS -1
#endif
#ifndef ALF_GPU_THRESH_ZGETRI
#define ALF_GPU_THRESH_ZGETRI -1
#endif
#ifndef ALF_GPU_THRESH_DGETRI
#define ALF_GPU_THRESH_DGETRI -1
#endif
#ifndef ALF_GPU_THRESH_ZGEQRF
#define ALF_GPU_THRESH_ZGEQRF -1
#endif
#ifndef ALF_GPU_THRESH_DGEQRF
#define ALF_GPU_THRESH_DGEQRF -1
#endif
#ifndef ALF_GPU_THRESH_ZUNGQR
#define ALF_GPU_THRESH_ZUNGQR -1
#endif
#ifndef ALF_GPU_THRESH_DORGQR
#define ALF_GPU_THRESH_DORGQR -1
#endif
#ifndef ALF_GPU_THRESH_ZUNMQR
#define ALF_GPU_THRESH_ZUNMQR -1
#endif
#ifndef ALF_GPU_THRESH_ZHEEV
#define ALF_GPU_THRESH_ZHEEV -1
#endif
#ifndef ALF_GPU_THRESH_DSYEV
#define ALF_GPU_THRESH_DSYEV -1
#endif
#ifndef ALF_GPU_THRESH_ZGESVD
#define ALF_GPU_THRESH_ZGESVD -1
#endif

module alf_linalg_mod

  !--------------------------------------------------------------------
  !> @author
  !> ALF-project
  !
  !> @brief
  !> Configuration of the wrappers around all BLAS and LAPACK routines called by ALF.
  !>
  !> The wrappers ALF_<ROUTINE> in alf_linalg_wrappers.F90 have the same argument
  !> list as the BLAS/LAPACK routine <ROUTINE>, so that a call can be converted by
  !> renaming it. Like BLAS/LAPACK they are external procedures with implicit
  !> interfaces, since ALF passes e.g. elements of assumed-shape arrays as start of
  !> a matrix, which is not allowed for procedures with explicit interfaces.
  !> The function ALF_DZNRM2 has to be declared as REAL(Kind=Kind(0.d0)).
  !>
  !> If ALF is compiled with -DALF_GPU (configure.sh option GPU), calls are
  !> offloaded to the GPU via OpenMP and oneMKL, if the size of the call is
  !> at least the threshold N of the routine. N is the size of an equivalent
  !> square problem: A call is offloaded if its work is at least N**3 for
  !> BLAS level 3 and LAPACK routines, N**2 for BLAS level 2 routines and
  !> N for BLAS level 1 routines. A negative N means "never offload".
  !>
  !> The threshold of e.g. ZGEMM is set at compile time with the macro
  !> ALF_GPU_THRESH_ZGEMM (e.g. ALF_FLAGS_EXT="-DALF_GPU_THRESH_ZGEMM=512")
  !> and can be overwritten at runtime with the environment variable of the
  !> same name. The environment variable ALF_GPU_OFFLOAD=0 disables offloading.
  !>
  !> Every offloaded call copies its arguments to the device and back, so the
  !> thresholds should be determined with the benchmark in Benchmarks/.
  !>
  !> Routines without OpenMP offload support in oneMKL are pure pass-throughs.
  !--------------------------------------------------------------------

#ifdef ALF_GPU
  use omp_lib
#endif
  use iso_fortran_env, only: int64, output_unit
  use runtime_error_mod
  implicit none
  private

  public :: alf_linalg_init, alf_linalg_set_threshold, alf_linalg_get_threshold, &
     &      alf_linalg_print_thresholds, alf_linalg_offload_active
  ! Used by the wrappers in alf_linalg_wrappers.F90
  public :: dp
#ifdef ALF_GPU
  public :: use_gpu, mspan, vspan, is_char
  public :: ID_ZGEMM, ID_DGEMM, ID_ZHEMM, ID_ZTRSM, ID_ZTRMM, ID_ZGEMV, ID_DGEMV, ID_ZGERU, &
     &      ID_ZSCAL, ID_ZCOPY, ID_ZSWAP, ID_DZNRM2, ID_ZGETRF, ID_DGETRF, ID_ZGETRS, ID_ZGETRI, &
     &      ID_DGETRI, ID_ZGEQRF, ID_DGEQRF, ID_ZUNGQR, ID_DORGQR, ID_ZUNMQR, ID_ZHEEV, ID_DSYEV, &
     &      ID_ZGESVD
#endif

  integer, parameter :: dp = kind(0.d0)

  integer, parameter :: ID_ZGEMM  =  1, ID_DGEMM  =  2, ID_ZHEMM  =  3, ID_ZTRSM  =  4, &
     &                  ID_ZTRMM  =  5, ID_ZGEMV  =  6, ID_DGEMV  =  7, ID_ZGERU  =  8, &
     &                  ID_ZSCAL  =  9, ID_ZCOPY  = 10, ID_ZSWAP  = 11, ID_DZNRM2 = 12, &
     &                  ID_ZGETRF = 13, ID_DGETRF = 14, ID_ZGETRS = 15, ID_ZGETRI = 16, &
     &                  ID_DGETRI = 17, ID_ZGEQRF = 18, ID_DGEQRF = 19, ID_ZUNGQR = 20, &
     &                  ID_DORGQR = 21, ID_ZUNMQR = 22, ID_ZHEEV  = 23, ID_DSYEV  = 24, &
     &                  ID_ZGESVD = 25
  integer, parameter :: N_ROUTINES = 25

  character(len=6), parameter :: names(N_ROUTINES) = [character(len=6) :: &
     & 'ZGEMM',  'DGEMM',  'ZHEMM',  'ZTRSM',  'ZTRMM',  'ZGEMV',  'DGEMV',  'ZGERU', &
     & 'ZSCAL',  'ZCOPY',  'ZSWAP',  'DZNRM2', 'ZGETRF', 'DGETRF', 'ZGETRS', 'ZGETRI', &
     & 'DGETRI', 'ZGEQRF', 'DGEQRF', 'ZUNGQR', 'DORGQR', 'ZUNMQR', 'ZHEEV',  'DSYEV', &
     & 'ZGESVD']
  ! Exponent relating the threshold N to the work of a call
  integer, parameter :: powers(N_ROUTINES) = [3, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1, 1, &
     &                                        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]

  integer, save :: thresh(N_ROUTINES) = [ &
     & ALF_GPU_THRESH_ZGEMM,  ALF_GPU_THRESH_DGEMM,  ALF_GPU_THRESH_ZHEMM,  ALF_GPU_THRESH_ZTRSM,  &
     & ALF_GPU_THRESH_ZTRMM,  ALF_GPU_THRESH_ZGEMV,  ALF_GPU_THRESH_DGEMV,  ALF_GPU_THRESH_ZGERU,  &
     & ALF_GPU_THRESH_ZSCAL,  ALF_GPU_THRESH_ZCOPY,  ALF_GPU_THRESH_ZSWAP,  ALF_GPU_THRESH_DZNRM2, &
     & ALF_GPU_THRESH_ZGETRF, ALF_GPU_THRESH_DGETRF, ALF_GPU_THRESH_ZGETRS, ALF_GPU_THRESH_ZGETRI, &
     & ALF_GPU_THRESH_DGETRI, ALF_GPU_THRESH_ZGEQRF, ALF_GPU_THRESH_DGEQRF, ALF_GPU_THRESH_ZUNGQR, &
     & ALF_GPU_THRESH_DORGQR, ALF_GPU_THRESH_ZUNMQR, ALF_GPU_THRESH_ZHEEV,  ALF_GPU_THRESH_DSYEV,  &
     & ALF_GPU_THRESH_ZGESVD]
  ! thresh**powers, compared against the work of a call
  integer(int64), save :: thresh_work(N_ROUTINES)
  logical, save :: initialized = .false.
  logical, save :: offload_active = .false.

contains

  !--------------------------------------------------------------------
  !> @brief
  !> Reads the thresholds from the environment and checks for GPU devices.
  !> Is called automatically by the first wrapper call, call explicitly to
  !> print the active configuration.
  !>
  !> @param[in] irank  Optional MPI rank, the configuration is printed on rank 0.
  !--------------------------------------------------------------------
  subroutine alf_linalg_init(irank)
    integer, intent(in), optional :: irank

    character(len=64) :: val
    integer :: id, length, stat, ios, n

    initialized = .true.
    do id = 1, N_ROUTINES
      call get_environment_variable('ALF_GPU_THRESH_'//trim(names(id)), val, length, stat)
      if (stat == 0 .and. length > 0) then
        read(val, *, iostat=ios) n
        if (ios /= 0) then
          write(output_unit, *) 'alf_linalg_init: Cannot read ALF_GPU_THRESH_'//trim(names(id))//'=', trim(val)
          call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
        endif
        thresh(id) = n
      endif
    enddo
    call update_thresh_work()

#ifdef ALF_GPU
    offload_active = omp_get_num_devices() > 0
    call get_environment_variable('ALF_GPU_OFFLOAD', val, length, stat)
    if (stat == 0 .and. length > 0) then
      if (trim(val) == '0') offload_active = .false.
    endif
#endif

    if (present(irank)) then
      if (irank == 0) call alf_linalg_print_thresholds()
    endif
  end subroutine alf_linalg_init

  subroutine update_thresh_work()
    integer :: id
    do id = 1, N_ROUTINES
      if (thresh(id) < 0) then
        thresh_work(id) = -1
      else
        thresh_work(id) = int(thresh(id), int64)**powers(id)
      endif
    enddo
  end subroutine update_thresh_work

  integer function name_to_id(name) result(id)
    character(len=*), intent(in) :: name
    character(len=len(name)) :: upper
    integer :: i, c
    upper = name
    do i = 1, len(name)
      c = iachar(name(i:i))
      if (c >= iachar('a') .and. c <= iachar('z')) upper(i:i) = achar(c - 32)
    enddo
    do id = 1, N_ROUTINES
      if (trim(names(id)) == trim(upper)) return
    enddo
    write(output_unit, *) 'alf_linalg: No offloadable routine named ', trim(name)
    call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
  end function name_to_id

  !> @brief Sets the offload threshold N of the routine called name, e.g. 'ZGEMM'.
  subroutine alf_linalg_set_threshold(name, n)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    if (.not. initialized) call alf_linalg_init()
    thresh(name_to_id(name)) = n
    call update_thresh_work()
  end subroutine alf_linalg_set_threshold

  !> @brief Returns the offload threshold N of the routine called name, e.g. 'ZGEMM'.
  integer function alf_linalg_get_threshold(name)
    character(len=*), intent(in) :: name
    if (.not. initialized) call alf_linalg_init()
    alf_linalg_get_threshold = thresh(name_to_id(name))
  end function alf_linalg_get_threshold

  !> @brief Returns true if calls above the thresholds are offloaded to a GPU.
  logical function alf_linalg_offload_active()
    if (.not. initialized) call alf_linalg_init()
    alf_linalg_offload_active = offload_active
  end function alf_linalg_offload_active

  subroutine alf_linalg_print_thresholds()
    integer :: id
#ifdef ALF_GPU
    if (offload_active) then
      write(output_unit, '(A,I0,A,I0)') ' GPU offloading active, devices: ', omp_get_num_devices(), &
         & ', default device: ', omp_get_default_device()
      write(output_unit, '(A)') ' GPU offload thresholds (equivalent square size, -1: never):'
      do id = 1, N_ROUTINES
        write(output_unit, '(3X,A6,1X,I0)') names(id), thresh(id)
      enddo
    else
      write(output_unit, '(A)') ' GPU offloading compiled in, but inactive (no device or ALF_GPU_OFFLOAD=0)'
    endif
#else
    id = 0
#endif
  end subroutine alf_linalg_print_thresholds

#ifdef ALF_GPU
  !> True if a call of routine id with the given work is offloaded.
  logical function use_gpu(id, work)
    integer, intent(in) :: id
    integer(int64), intent(in) :: work
    if (.not. initialized) call alf_linalg_init()
    use_gpu = offload_active .and. thresh_work(id) >= 0 .and. work > 0 .and. work >= thresh_work(id)
  end function use_gpu

  !> Number of elements touched by a rows x cols matrix with leading dimension ld.
  !> At least one, so that the base address of an empty matrix is mapped as well.
  pure integer(int64) function mspan(ld, rows, cols)
    integer, intent(in) :: ld, rows, cols
    if (rows <= 0 .or. cols <= 0) then
      mspan = 1
    else
      mspan = int(ld, int64)*(cols - 1) + rows
    endif
  end function mspan

  !> Number of elements touched by a vector of length n with increment inc.
  pure integer(int64) function vspan(n, inc)
    integer, intent(in) :: n, inc
    if (n <= 0) then
      vspan = 1
    else
      vspan = int(abs(inc), int64)*(n - 1) + 1
    endif
  end function vspan

  pure logical function is_char(c, upper)
    character, intent(in) :: c, upper
    is_char = c == upper .or. iachar(c) == iachar(upper) + 32
  end function is_char
#endif

end module alf_linalg_mod
