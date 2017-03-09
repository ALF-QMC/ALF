!compile with gfortran 12-UDV-R.f90 -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ ../../Libraries/Modules/modules_90.a -llapack -lblas ../../Libraries/myold/libnag.a

Program UDVR
Use MyMats
implicit none
REAL (KIND=KIND(0.D0)), dimension(:, :), allocatable :: a
REAL (KIND=KIND(0.D0)), dimension(:, :), allocatable :: u, v, test
REAL (KIND=KIND(0.D0)), dimension(:), allocatable :: d
INTEGER :: ncon, i, j, ndim, nr
REAL (KIND=KIND(0.D0)) :: Z

ndim = 30
ncon = 1
allocate(a(ndim, ndim), u(ndim, ndim), v(ndim, ndim), d(ndim), test(ndim, ndim))
! initialize A to the Hilbert matrix
do i = 1, ndim
do j = 1, ndim
a(i, j) = 1.D0/(i+j)
enddo
enddo
call UDV1_R(a, u, d, v, ncon)
! check decomposition
           DO J = 1,ndim
              DO I = 1,ndim
                 Z = 0.D0
                 DO NR = 1,ndim 
                    Z = Z + U(I,NR)*D(NR)*V(NR,J)
                 ENDDO
                 TEST(I,J) = Z
              ENDDO
           ENDDO
do i = 1, ndim
do j = 1, ndim
if (Abs(a(i,j) - test(i,j)) > max(abs(a(i,j)), abs(test(i,j)))*1D-14) then
write (*, *) i, j, a(i, j),  " != ", test(i, j), "diff: ", Abs(a(i,j) - test(i,j))
STOP 2
end if
enddo
enddo
deallocate(a, u, d, v, test)
end program
