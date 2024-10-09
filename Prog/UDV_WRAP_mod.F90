!  Copyright (C) 2016 - 2018 The ALF project
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
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
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
!> This module contains two version of the stabilization.  To switch between the two schemes
!> you should   define  STAB1  in the set_env.sh file.    The defaut scheme is quick  and
!> gernerically works better.
!
!--------------------------------------------------------------------

module UDV_Wrap_mod

   use MyMats
   use Files_mod

contains

!***************************************************************

#if defined(STAB1)
   subroutine UDV_Wrap_Pivot(A, U, D, V, NCON, N1, N2)
      use QDRP_mod

      implicit none
      complex(Kind=kind(0.d0)), intent(IN), dimension(:, :) :: A
      complex(Kind=kind(0.d0)), intent(INOUT), dimension(:, :) :: U, V
      complex(Kind=kind(0.d0)), intent(INOUT), dimension(:) :: D
      integer, intent(IN) :: NCON
      integer, intent(IN) :: N1, N2

      ! Locals
      real(Kind=kind(0.d0)) :: VHELP(N2), XNORM(N2), XMAX, XMEAN
      integer :: IVPT(N2), IVPTM1(N2), I, J, K, IMAX
      complex(Kind=kind(0.d0))  :: A1(N1, N2), A2(N1, N2), V1(N2, N2), phase, beta

      do I = 1, N2
         XNORM(I) = 0.d0
         do J = 1, N1
            XNORM(I) = XNORM(I) + dble(A(J, I)*conjg(A(J, I)))
         end do
      end do
      VHELP = XNORM

      do I = 1, N2
         XMAX = VHELP(1)
         IMAX = 1
         do J = 2, N2
            if (VHELP(J) .gt. XMAX) then
               IMAX = J
               XMAX = VHELP(J)
            end if
         end do
         VHELP(IMAX) = -1.d0
         IVPTM1(IMAX) = I
         IVPT(I) = IMAX
      end do
      do I = 1, N2
         K = IVPT(I)
         A1(:, I) = A(:, K)
      end do

      call UDV_Wrap(A1, U, D, V, NCON)
      V1 = V
      Phase = Det_C(V1, N2)
      call Pivot_Phase(phase, IVPT, size(D, 1))
      beta = 1/Phase
      !scale first row of R with 1/phase to set Det(R)=1 [=Det(V)]
      call ZSCAL(size(V, 2), beta, V(1, 1), size(V, 1))
      ! scale first column of U to correct the scaling in V such that UDV is not changed
      call ZSCAL(size(U, 1), phase, U(1, 1), 1)

      V1 = V
      do I = 1, N2
         K = IVPTM1(I)
         V(:, I) = V1(:, K)
      end do

      if (NCON == 1) then
         !Check the result  A = U D V
         do J = 1, N2
            A1(:, J) = D*V(:, J)
         end do
         write (6, *) 'Here'
         call MMULT(A2, U, A1)
         call COMPARE(A, A2, XMAX, XMEAN)
         write (6, *) 'Check afer  Pivoting', XMAX
      end if

   end subroutine UDV_Wrap_Pivot
#else
   subroutine UDV_Wrap_Pivot(A, U, D, V, NCON, N1, N2)
      use QDRP_mod

      implicit none
      complex(Kind=kind(0.d0)), intent(IN), dimension(:, :) :: A
      complex(Kind=kind(0.d0)), intent(INOUT), dimension(:, :) :: U, V
      complex(Kind=kind(0.d0)), intent(INOUT), dimension(:) :: D
      integer, intent(IN) :: NCON
      integer, intent(IN) :: N1, N2

      ! Locals
      real(Kind=kind(0.d0)) :: VHELP(N2), XNORM(N2), XMAX, XMEAN
      integer :: IVPT(N2), IVPTM1(N2), I, J, K, IMAX
      complex(Kind=kind(0.d0))  :: A1(N1, N2), A2(N1, N2), V1(N2, N2), Z, phase, beta

      do I = 1, N2
         XNORM(I) = 0.d0
         do J = 1, N1
            XNORM(I) = XNORM(I) + dble(A(J, I)*conjg(A(J, I)))
         end do
      end do
      VHELP = XNORM

      do I = 1, N2
         XMAX = VHELP(1)
         IMAX = 1
         do J = 2, N2
            if (VHELP(J) .gt. XMAX) then
               IMAX = J
               XMAX = VHELP(J)
            end if
         end do
         VHELP(IMAX) = -1.d0
         IVPTM1(IMAX) = I
         IVPT(I) = IMAX
      end do
      do I = 1, N2
         K = IVPT(I)
         A1(:, I) = A(:, K)/cmplx(XNORM(K), 0.d0, kind(0.d0))
      end do

      call UDV(A1, U, D, V1, NCON)
      Phase = cmplx(1.d0, 0.d0, kind(0.d0))
      do i = 1, size(D, 1)
         Phase = Phase*V1(i, i)
      end do
      call Pivot_Phase(phase, IVPT, size(D, 1))
      beta = 1/Phase
      !scale first row of R with 1/phase to set Det(R)=1 [=Det(V)]
      call ZSCAL(size(V1, 2), beta, V1(1, 1), size(V1, 1))
      ! scale first column of U to correct the scaling in V such that UDV is not changed
      call ZSCAL(size(U, 1), phase, U(1, 1), 1)

      ! Finish the pivotting.
      do I = 1, N2
         D(I) = D(I)*XNORM(IVPT(I))
      end do
      do I = 1, N2 - 1
         Z = 1.d0/XNORM(IVPT(I))
         do J = I + 1, N2
            V1(I, J) = V1(I, J)*XNORM(IVPT(J))*Z
         end do
      end do

      do J = 1, N2
         do I = 1, N2
            V(I, J) = V1(I, IVPTM1(J))
         end do
      end do

      if (NCON == 1) then
         !Check the result  A = U D V
         do J = 1, N2
            A1(:, J) = D*V(:, J)
         end do
         !Write(6,*) 'Here'
         call MMULT(A2, U, A1)
         !Call MMULT (A2,U,V)
         call COMPARE(A, A2, XMAX, XMEAN)
         write (6, *) 'Check afer  Pivoting', XMAX
      end if

   end subroutine UDV_Wrap_Pivot
#endif

!***************************************************************
   subroutine UDV_Wrap(A, U, D, V, NCON)
#ifdef MPI
      use mpi
#endif

      implicit none

      complex(Kind=kind(0.d0)), intent(IN), dimension(:, :) :: A
      complex(Kind=kind(0.d0)), intent(INOUT), dimension(:, :) :: U, V
      complex(Kind=kind(0.d0)), intent(INOUT), dimension(:) :: D
      integer, intent(IN) :: NCON

      !Local
      complex(Kind=kind(0.d0)), allocatable ::  A1(:, :), U1(:, :)
      integer :: N
      character(len=64) :: file_sr, File
#ifdef MPI
      integer :: Isize, Irank, Ierr

      call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
#endif

      File_sr = "SDV"
#ifdef MPI
      File = File_i(File_sr, Irank)
#else
      File = File_sr
#endif
      call QR(A, U, V, NCON)
      N = size(A, 1)
      allocate (A1(N, N), U1(N, N))
      !DO I = 1,N
      !   DO J = 1,N
      !      WRITE(6,*) I,J,V(I,J)
      !   ENDDO
      !ENDDO
      A1 = V
      call SVD(A1, U1, D, V, NCON)
      call MMULT(A1, U, U1)
      U = A1

   end subroutine UDV_Wrap

end module UDV_Wrap_mod
