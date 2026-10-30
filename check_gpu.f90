! Small check if BLAS calls can be offloaded to a GPU with OpenMP and oneMKL.
include "mkl_blas_omp_offload_lp64.f90"

program check_gpu

  use iso_fortran_env
  use omp_lib
  use onemkl_blas_omp_offload_lp64
  implicit none
  integer, parameter :: N=64
  complex (kind=kind(0.d0)) :: A(N,N), B(N,N), C(N,N), C_ref(N,N)
  real    (kind=kind(0.d0)) :: X(N,N)

  if (omp_get_num_devices() < 1) then
    write(error_unit, *) "No GPU device found"
    error stop 1
  endif
  write(output_unit, '(A,I0,A,I0)') " GPU devices found: ", omp_get_num_devices(), &
     & ", default device: ", omp_get_default_device()

  call random_number(X)
  A = cmplx(X, 1.d0 - X, kind(0.d0))
  call random_number(X)
  B = cmplx(1.d0 - X, X, kind(0.d0))
  C = (0.d0, 0.d0)
  C_ref = matmul(A, B)

  !$omp target data map(to: A, B) map(tofrom: C)
  !$omp dispatch
  call ZGEMM('N', 'N', N, N, N, (1.d0, 0.d0), A, N, B, N, (0.d0, 0.d0), C, N)
  !$omp end target data

  if (.not. maxval(abs(C - C_ref)) < 1d-10) then
    write(error_unit, *) "Problem with offloaded ZGEMM"
    error stop 1
  endif

end program check_gpu
