!  Copyright (C) 2016 - 2022 The ALF project
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

module cgr1_mod
   implicit none
contains

   subroutine CGR(PHASE, NVAR, GRUP, udvr, udvl)

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!>    Computes  GRUP = (1 + UR*DR*VR*VL*DL*UL)^-1
!>    and      PHASE = det(1 + UR*DR*VR*VL*DL*UL) / abs(det(1 + UR*DR*VR*VL*DL*UL))
!>    NVAR = 1 Big scales are in DL
!>    NVAR = 2 Big scales are in DR
!> Implementation note: we calculate the Phase as:
!> NVAR = 1 : Phase = det(URUP * ULUP)/ |det(URUP * ULUP)| * det(P) * det(R) *det(Q)/ |det(R) det(Q)|
!> NVAR = 2 : Phase = det(URUP * ULUP)/ |det(URUP * ULUP)| * det(P) * det^*(R) *det^*(Q)/ |det(R) det(Q)|
!> If STAB3 is selected the following tweak is applied
!> We seperate D as D^+ * D^- where D^+ (D^-) contains the scales larger (smaller) then 1.0
!> Also, we use (DR^+^-1 UR^* UL^* DL^+^-1 + DR^- VR VL DL^- )^-1 = UL^* DL^+^-1 GRUP Dr^+^-1 UR^* .
!
!--------------------------------------------------------------------

      use UDV_State_mod

#if (defined(STAB2) || defined(STAB1)) && !defined(STABLOG)
      use UDV_Wrap_mod

      implicit none

      !Arguments.
      class(UDV_State), intent(IN) :: udvl, udvr
      complex(Kind=kind(0.d0)), dimension(:, :), intent(INOUT) :: GRUP
      complex(Kind=kind(0.d0)) :: PHASE
      integer         :: NVAR

      !Local
      logical, save :: Scale_warning_message = .true.
      type(UDV_State) :: udvlocal
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable :: TPUP, TPUP1, TPUPM1
      integer, dimension(:), allocatable :: IPVT
      complex(Kind=kind(0.d0)) ::  ZDUP1, ZDDO1, ZDUP2, ZDDO2, Z1, ZUP, ZDO, alpha, beta
      integer :: I, J, N_size, NCON, info
      real(Kind=kind(0.d0)) :: X, Xmax, sv

      if (.not. allocated(UDVL%V)) then
         !call projector cgr
         call cgrp(phase, grup, udvr, udvl)
         return
      end if

      if (udvl%side .ne. "L" .and. udvl%side .ne. "l") then
         write (*, *) "calling wrong decompose"
      end if
      if (udvr%side .ne. "R" .and. udvr%side .ne. "r") then
         write (*, *) "calling wrong decompose"
      end if

      N_size = udvl%Ndim
      NCON = 0
      alpha = 1.d0
      beta = 0.d0
      allocate (TPUP(N_size, N_size), TPUP1(N_size, N_size), TPUPM1(N_size, N_size), IPVT(N_size))
      call udvlocal%alloc(N_size)
      !Write(6,*) 'In CGR', N_size
      call MMULT(udvlocal%V, udvr%V, udvl%V)
      if (dble(udvr%D(1)*udvlocal%V(1, 1)*udvl%D(1)) > 0.1*huge(1.0d0) .and. Scale_warning_message) then
         write (error_unit, *)
         write (error_unit, *) "Warning: Large number encountered; Generation of NaN's is imminent"
         write (error_unit, *) "         Switching to LOG stablilization scheme is likely to help"
         Scale_warning_message = .false.
      end if
      do J = 1, N_size
         TPUP(:, J) = udvr%D(:)*udvlocal%V(:, J)*udvl%D(J)
      end do
      call ZGEMM('C', 'N', N_size, N_size, N_size, alpha, udvr%U(1, 1), N_size, udvl%U(1, 1), N_size, alpha, TPUP, N_size)
      !>  Syntax
      !>  ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      !>  C := alpha*op( A )*op( B ) + beta*C
      !>  TPUP =  (URUP)^(dagger) ULUP^(dagger) + TPUP
      if (NVAR .eq. 1) then
         !WRITE(6,*) 'UDV of U + DR * V * DL'
         call UDV_WRAP_Pivot(TPUP, udvlocal%U, udvlocal%D, udvlocal%V, NCON, N_size, N_Size)
         !CALL UDV(TPUP,UUP,DUP,VUP,NCON)
!            CALL MMULT(TPUP,udvlocal%V, udvl%U)
         call ZGEMM('N', 'C', N_size, N_size, N_size, alpha, udvlocal%V(1, 1), N_size, udvl%U(1, 1), N_size, beta, TPUP, N_size)
         !Do I = 1,N_size
         !   Write(6,*) DLUP(I)
         !enddo
         call INV(TPUP, TPUPM1, ZDUP1)
         call MMULT(TPUP1, udvr%U, udvlocal%U)
         call ZGETRF(N_size, N_size, TPUP1, N_size, IPVT, info)
         !>  TPUP1 = P * L * U   LU-decomposition
         Z1 = ZDUP1
         do i = 1, N_size
            if (IPVT(i) .ne. i) then
               Z1 = -Z1
            end if
            Z1 = Z1*TPUP1(I, I)

         end do
      else
         !WRITE(6,*) 'UDV of (U + DR * V * DL)^{*}'
         TPUP1 = CT(TPUP)
         call UDV_WRAP_Pivot(TPUP1, udvlocal%U, udvlocal%D, udvlocal%V, NCON, N_size, N_size)
         !CALL UDV(TPUP1,UUP,DUP,VUP,NCON)
         call ZGEMM('N', 'N', N_size, N_size, N_size, alpha, udvl%U(1, 1), N_size, udvlocal%U, N_size, beta, TPUPM1, N_size)
         call ZGEMM('N', 'C', N_size, N_size, N_size, alpha, udvr%U(1, 1), N_size, udvlocal%V, N_size, beta, TPUP1, N_size)
         call ZGETRF(N_size, N_size, TPUP1, N_size, IPVT, info)
         !>  TPUP1 = P * L * U   LU-decomposition
         ZDUP2 = 1.d0
         do i = 1, N_size
            ZDUP2 = ZDUP2*TPUP1(I, I)
            if (IPVT(i) .ne. i) then
               ZDUP2 = -ZDUP2
            end if
         end do
         TPUP = TPUPM1
         ZDUP1 = DET_C(TPUP, N_size)! Det destroys its argument
         Z1 = ZDUP2/ZDUP1
      end if
      do J = 1, N_size
         sv = dble(udvlocal%D(J))
         X = abs(sv)
         if (J == 1) Xmax = X
         if (X < Xmax) Xmax = X
         sv = 1.d0/sv
         do I = 1, N_size
            udvlocal%U(J, I) = TPUPM1(I, J)*sv
         end do
      end do
      call ZGETRS('T', N_size, N_size, TPUP1, N_size, IPVT, udvlocal%U, N_size, info)
      !> Syntax
      !> ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      !> Op(A) * X = B
      !> On output  B = X
      GRUP = transpose(udvlocal%U)
      PHASE = Z1/abs(Z1)
      call udvlocal%dealloc
      deallocate (TPUP, TPUP1, TPUPM1, IPVT)

#else

      use MyMats
      use QDRP_mod

      implicit none
      !Arguments.
!         COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(IN)   ::  URUP, VRUP, ULUP, VLUP
!         COMPLEX(Kind=Kind(0.d0)), Dimension(:),   Intent(IN)   ::  DLUP, DRUP
      class(UDV_State), intent(IN) :: udvl, udvr
      complex(Kind=kind(0.d0)), dimension(:, :), intent(INOUT) :: GRUP
      complex(Kind=kind(0.d0)), intent(INOUT) :: PHASE
      integer         :: NVAR

      !Local
      logical, save :: Scale_warning_message = .true.
      complex(Kind=kind(0.d0)), dimension(:, :), allocatable ::  TPUP, RHS
      complex(Kind=kind(0.d0)), dimension(:), allocatable ::  DUP
      integer, dimension(:), allocatable :: IPVT, VISITED
      complex(Kind=kind(0.d0)) ::  alpha, beta, Z, DLJ
      integer :: I, J, N_size, info, LWORK, next, L
      real(Kind=kind(0.d0)) :: X, Xmax, sv

      complex(Kind=kind(0.d0)), allocatable, dimension(:) :: TAU, WORK
      logical :: FORWRD

      if (udvl%side .ne. "L" .and. udvl%side .ne. "l") then
         write (*, *) "calling wrong decompose"
      end if
      if (udvr%side .ne. "R" .and. udvr%side .ne. "r") then
         write (*, *) "calling wrong decompose"
      end if

      if (.not. allocated(UDVL%V)) then
         !call projector cgr
         call cgrp(phase, grup, udvr, udvl)
         return
      end if

      N_size = udvl%ndim
      alpha = 1.d0
      beta = 0.d0
      allocate (TPUP(N_size, N_size), RHS(N_size, N_size), IPVT(N_size), TAU(N_size), DUP(N_size))
      !Write(6,*) 'In CGR', N_size
      ! can be inserted again once we are sure that we may assume that UR and UL stem from householder reflectors
!        CALL ZGEMM('C', 'C', N_size, N_size, N_size, alpha, URUP, N_size, ULUP, N_size, alpha, TPUP, N_size)
      call ZGEMM('C', 'N', N_size, N_size, N_size, alpha, udvr%U, N_size, udvl%U, N_size, beta, RHS(1, 1), N_size)

      call MMULT(TPUP, udvr%V, udvl%V)
#if !(defined(STAB3) || defined(STABLOG))
      if (dble(udvr%D(1)*TPUP(1, 1)*udvl%D(1) + RHS(1, 1)) > 0.1*huge(1.0d0) .and. Scale_warning_message) then
         write (error_unit, *)
         write (error_unit, *) "Warning: Large number encountered; Generation of NaN's is imminent"
         write (error_unit, *) "         Switching to LOG stablilization scheme is likely to help"
         Scale_warning_message = .false.
      end if
      do J = 1, N_size
         TPUP(:, J) = udvr%D(:)*TPUP(:, J)*udvl%D(J)
      end do
      TPUP = TPUP + RHS
#else
#if ! defined(STABLOG)
      !missuse DUP(I) as DR(I) for temporary storage
      !scales in D are assumed to be real and positive
      do I = 1, N_size
         if (dble(udvr%D(I)) <= 1.d0) then
            DUP(I) = udvr%D(I)
         else
            DUP(I) = 1.d0/udvr%D(I)
         end if
      end do
      do J = 1, N_size
         if (dble(udvl%D(J)) <= 1.d0) then
            DLJ = udvl%D(J)
            do I = 1, N_size
               if (dble(udvr%D(I)) <= 1.d0) then
                  TPUP(I, J) = RHS(I, J) + udvr%D(I)*udvl%D(J)*TPUP(I, J)
               else
                  TPUP(I, J) = DUP(I)*RHS(I, J) + DLJ*TPUP(I, J)
               end if
            end do
         else
            DLJ = 1.d0/udvl%D(J)
            do I = 1, N_size
               if (dble(udvr%D(I)) <= 1.d0) then
                  TPUP(I, J) = DLJ*RHS(I, J) + DUP(I)*TPUP(I, J)
               else
                  TPUP(I, J) = RHS(I, J)/udvr%D(I)/udvl%D(J) + TPUP(I, J)
               end if
            end do
         end if
      end do
#else
      !missuse DUP(I) as DR(I) for temporary storage
      do I = 1, N_size
         if (udvr%L(I) <= 0.d0) then
            DUP(I) = cmplx(exp(udvr%L(I)), 0.d0, kind(0.d0))
         else
            DUP(I) = cmplx(exp(-udvr%L(I)), 0.d0, kind(0.d0))
         end if
      end do
      do J = 1, N_size
         if (udvl%L(J) <= 0.d0) then
            DLJ = cmplx(exp(udvl%L(J)), 0.d0, kind(0.d0))
            do I = 1, N_size
               if (udvr%L(I) <= 0.d0) then
                  TPUP(I, J) = RHS(I, J) + cmplx(exp(udvr%L(I) + udvl%L(J)), 0.d0, kind(0.d0))*TPUP(I, J)
               else
                  TPUP(I, J) = DUP(I)*RHS(I, J) + DLJ*TPUP(I, J)
               end if
            end do
         else
            DLJ = cmplx(exp(-udvl%L(J)), 0.d0, kind(0.d0))
            do I = 1, N_size
               if (udvr%L(I) <= 0.d0) then
                  TPUP(I, J) = DLJ*RHS(I, J) + DUP(I)*TPUP(I, J)
               else
                  TPUP(I, J) = cmplx(exp(-udvr%L(I) - udvl%L(J)), 0.d0, kind(0.d0))*RHS(I, J) + TPUP(I, J)
               end if
            end do
         end if
      end do
#endif
#endif
      ! calculate determinant of UR*UL
      ! as the D's are real and positive, they do not contribute the the phase of det so they can be ignored
      PHASE = conjg(DET_C(RHS, N_size))
      PHASE = PHASE/abs(PHASE)
      IPVT = 0
      if (NVAR .ne. 1) then
         TPUP = conjg(transpose(TPUP))
      end if
      call QDRP_decompose(N_size, udvl%N_part, TPUP, DUP, IPVT, TAU, WORK, LWORK)
      allocate (VISITED(N_size))
      ! Calculate the sign of the permutation from the pivoting. Somehow the format used by the QR decomposition of lapack
      ! is different from that of the LU decomposition of lapack
      VISITED = 0
      do i = 1, N_size
         if (VISITED(i) .eq. 0) then
            next = i
            L = 0
            do while (VISITED(next) .eq. 0)
               L = L + 1
               VISITED(next) = 1
               next = IPVT(next)
            end do
            if (mod(L, 2) .eq. 0) then
               PHASE = -PHASE
            end if
         end if
      end do
      !calculate the determinant of the unitary matrix Q and the upper triangular matrix R
      do i = 1, N_size
         Z = TAU(i)
         if (NVAR .eq. 1) then
            PHASE = PHASE*TPUP(i, i)/abs(TPUP(i, i))
         else
            ! the diagonal should be real, but let's be safe....
            PHASE = PHASE*conjg(TPUP(i, i))/abs(TPUP(i, i))
            Z = conjg(Z) ! conjugate the elementary reflector
         end if
         if (Z .ne. cmplx(0.d0, 0.d0, Kind=kind(0.d0))) then
            ! here we calculate the determinant of a single householder reflector: det(1 - tau * v v* ) = 1 - tau * v^* v
            ! In lapack the scalar tau and the vector v are scaled such that |tau|^2 |v|^2 = 2 Re(tau)
            ! The complete determinant det(Q) is the product of all reflectors. See http://www.netlib.org/lapack/lug/node128.html
            X = abs(Z)
            Z = 1.d0 - 2.d0*(Z/X)*(dble(Z)/X)
            PHASE = PHASE*Z/abs(Z)
         end if
      end do
      if (NVAR .eq. 1) then
         ! This is supposed to solve the system
         ! URUP U D V P^dagger ULUP G = 1
         ! initialize the rhs with CT(URUP)
         RHS = CT(udvr%U)
#if (defined(STAB3) || defined(STABLOG))
         !scale RHS=R_+^-1*RHS
         do J = 1, N_size
#if !defined(STABLOG)
            if (dble(UDVR%D(J)) > 1.d0) call ZSCAL(N_size, 1.d0/UDVR%D(J), RHS(J, 1), N_size)
#else
            if (UDVR%L(J) > 0.d0) call ZSCAL(N_size, cmplx(exp(-UDVR%L(J)), 0.d0, kind(0.d0)), RHS(J, 1), N_size)
#endif
         end do
#endif
         ! RHS = U^dagger * RHS
         call ZUNMQR('L', 'C', N_size, N_size, N_size, TPUP(1, 1), N_size, TAU(1), RHS(1, 1), N_size, WORK(1), LWORK, INFO)
         deallocate (TAU, WORK)
         !apply inverse of D to RHS from the left
         do J = 1, N_size
            sv = dble(DUP(J))
            X = abs(sv)
            if (J == 1) Xmax = X
            if (X < Xmax) Xmax = X
            do I = 1, N_size
               RHS(I, J) = RHS(I, J)/DUP(I)
            end do
         end do
         ! We solve the equation
         !  A * G = RHS for G with A = R * P^dagger * ULUP
         ! first we solve R *y = RHS. The solution is afterwards in RHS
         call ZTRSM('L', 'U', 'N', 'N', N_size, N_size, alpha, TPUP(1, 1), N_size, RHS(1, 1), N_size)
         ! apply permutation matrix
         FORWRD = .false.
         call ZLAPMR(FORWRD, N_size, N_size, RHS(1, 1), N_size, IPVT(1))
#if (defined(STAB3) || defined(STABLOG))
         !scale RHS=L_+^-1*RHS
         do J = 1, N_size
#if !defined(STABLOG)
            if (dble(UDVL%D(J)) > 1.d0) call ZSCAL(N_size, 1.d0/UDVL%D(J), RHS(J, 1), N_size)
#else
            if (UDVL%L(J) > 0.d0) call ZSCAL(N_size, cmplx(exp(-UDVL%L(J)), 0.d0, kind(0.d0)), RHS(J, 1), N_size)
#endif
         end do
#endif
         ! perform multiplication with ULUP and store in GRUP
         call ZGEMM('N', 'N', N_size, N_size, N_size, alpha, udvl%U(1, 1), N_size, RHS(1, 1), N_size, beta, GRUP(1, 1), N_size)
      else
         ! This solves the system G * URUP * P * R^dagger * D * U^dagger * ULUP = 1

         ! RHS = ULUP * UUP
         RHS = udvl%U !CT(udvl%U)
#if (defined(STAB3) || defined(STABLOG))
         !scale RHS=RHS*L_+^-1
         do J = 1, N_size
#if !defined(STABLOG)
            if (dble(UDVL%D(J)) > 1.d0) call ZSCAL(N_size, 1.d0/UDVL%D(J), RHS(1, J), 1)
#else
            if (UDVL%L(J) > 0.d0) call ZSCAL(N_size, cmplx(exp(-UDVL%L(J)), 0.d0, kind(0.d0)), RHS(1, J), 1)
#endif
         end do
#endif
         call ZUNMQR('R', 'N', N_size, N_size, N_size, TPUP(1, 1), N_size, TAU(1), RHS(1, 1), N_size, WORK(1), LWORK, INFO)
         deallocate (TAU, WORK)
         ! apply D^-1 to RHS from the right
         do J = 1, N_size
            sv = dble(DUP(J))
            X = abs(sv)
            if (J == 1) Xmax = X
            if (X < Xmax) Xmax = X
            sv = 1.d0/sv
            do I = 1, N_size
               RHS(I, J) = RHS(I, J)*sv
            end do
         end do

         ! We solve the equation
         ! G * A = RHS for G with A = URUP * P * R^dagger
         ! first we solve y * R^dagger = RHS
         call ZTRSM('R', 'U', 'C', 'N', N_size, N_size, alpha, TPUP(1, 1), N_size, RHS(1, 1), N_size)
         ! apply inverse permutation matrix
         FORWRD = .false.
         call ZLAPMT(FORWRD, N_size, N_size, RHS(1, 1), N_size, IPVT(1))
#if (defined(STAB3) || defined(STABLOG))
         ! first scale RHS=RHS*R_+^-1
         do J = 1, N_size
#if !defined(STABLOG)
            if (dble(UDVR%D(J)) > 1.d0) call ZSCAL(N_size, 1.d0/UDVR%D(J), RHS(1, J), 1)
#else
            if (UDVR%L(J) > 0.d0) call ZSCAL(N_size, cmplx(exp(-UDVR%L(J)), 0.d0, kind(0.d0)), RHS(1, J), 1)
#endif
         end do
#endif
         ! perform multiplication with URUP
         call ZGEMM('N', 'C', N_size, N_size, N_size, alpha, RHS(1, 1), N_size, udvr%U(1, 1), N_size, beta, GRUP(1, 1), N_size)
      end if
      deallocate (TPUP, DUP, IPVT, VISITED, RHS)
#endif

   end subroutine CGR

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Computes the Green's function in the projective implementation.
!
!> @param[out] PHASE
!> @param[out] GRUP
!> @param[in] udvr
!> @param[in] udvl
!
!--------------------------------------------------------------------
   subroutine CGRP(phase, GRUP, udvr, udvl)
      use UDV_State_mod
      class(UDV_State), intent(IN) :: udvl, udvr
      complex(Kind=kind(0.d0)), dimension(:, :), intent(OUT) :: GRUP
      complex(Kind=kind(0.d0)), intent(OUT) :: phase

      complex(Kind=kind(0.d0)), allocatable, dimension(:, :) :: sMat, rMat
      integer, allocatable :: ipiv(:)
      complex(Kind=kind(0.d0)) :: alpha, beta
      integer :: Ndim, N_part, info, n

      if ((udvl%side .ne. "L") .and. (udvl%side .ne. "l")) then
         write (*, *) "cgrp: udvl is not of type left"
         write (*, *) "cgrp: actual side is ", udvl%side
      end if
      if ((udvr%side .ne. "R") .and. (udvr%side .ne. "r")) then
         write (*, *) "cgrp: udvr is not of type right"
         write (*, *) "cgrp: actual side is ", udvr%side
      end if

      Ndim = udvl%ndim
      N_part = udvl%n_part
      allocate (sMat(N_part, N_part), ipiv(N_part), rMat(N_part, Ndim))

      ! Gr = Ur (Ul Ur)^-1 Ul
      ! Phase = 1 + Ur (Ul Ur)^-1 Ul
      ! Ul = udvl%U ^dag
      alpha = 1.d0
      beta = 0.d0
      call ZGEMM('C', 'N', N_part, N_part, Ndim, alpha, udvl%U(1, 1), Ndim, udvr%U(1, 1), Ndim, beta, sMat(1, 1), N_part)

      ! ZGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call ZGETRF(N_part, N_part, sMat, N_part, ipiv, info)
      phase = 1.d0
      do n = 1, N_part
         if (ipiv(n) .ne. n) then
            phase = -phase*sMat(n, n)/abs(sMat(n, n))
         else
            phase = phase*sMat(n, n)/abs(sMat(n, n))
         end if
      end do
      rMat = conjg(transpose(udvl%U))
      call zgetrs('N', N_part, Ndim, sMat(1, 1), N_part, ipiv, rMat(1, 1), N_part, info)
      alpha = -1.d0
      call ZGEMM('N', 'N', Ndim, Ndim, N_part, alpha, udvr%U(1, 1), Ndim, rMat(1, 1), N_part, beta, GRUP(1, 1), Ndim)
      do n = 1, Ndim
         GRUP(n, n) = GRUP(n, n) + cmplx(1.d0, 0.d0, kind(0.d0))
      end do
      deallocate (sMat, rMat, ipiv)

   end subroutine CGRP
end module cgr1_mod
