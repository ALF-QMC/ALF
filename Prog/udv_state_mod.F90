!  Copyright (C) 2017, 2018 The ALF project
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
!       https://alf.physik.uni-wuerzburg.de .
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage https://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Handles UDV decompositions
!>
!>
!--------------------------------------------------------------------

module UDV_State_mod
   use runtime_error_mod
   use iso_fortran_env, only: output_unit, error_unit

   implicit none
   private
   public :: UDV_State

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> Handles UDV decompositions
!>
!> @details
!> The UDV instance contains the following variables
!> @param Ndim Integer
!> \verbatim
!> The total number of orbitals
!> \endverbatim
!> @param N_part Integer
!> \verbatim
!> For the finite temperature code N_part=Ndim. For the projective code, N_part corresponds to the number of particles per flavor.
!> \endverbatim
!> @param U(Ndim,N\_part) Complex
!> @param D(N\_part) Complex
!> @param V(N\_part,N\_part) Complex
!> @param side Char
!> \verbatim
!>  side = R  for right propagation :   U(tau,0) * P_R = U_R d_r v_r
!>  side = L  for left  propagation :   (P_L * U(Theta,tau) )^{dag} =  U_L d_l {v_l}^{dag}
!>  For the finite temperature code, P_L and _P_R are unit matrices
!>  For the projective code the matrices v_l and v_r are not allocated
!> \endverbatim
!> @param L(N\_part) Real
!> \verbatim
!>  Space for the logscale option D=e^{L}
!> \endverbatim
!--------------------------------------------------------------------

   type UDV_State
      complex(Kind=kind(0.d0)), allocatable :: U(:, :), V(:, :)
#if !defined(STABLOG)
      complex(Kind=kind(0.d0)), allocatable :: D(:)
#else
      real(Kind=kind(0.d0)), allocatable :: L(:)
#endif
      integer   :: ndim, n_part  ! ndim: number of orbitals per flavor. n_part: number of particles per flavor
      character :: side          ! side = R  for right propagation :   B       * P_R = U d v
      ! side = L  for lesft propagation :   (P_L B)^{dag} = U d v^{dag}
   contains
      procedure :: alloc => alloc_UDV_state
      procedure :: init => init_UDV_state
      procedure :: dealloc => dealloc_UDV_state
      procedure :: reset => reset_UDV_state
      procedure :: assign => assign_UDV_state
      procedure :: decompose => decompose_UDV_state
      procedure :: print => print_UDV_state
      procedure :: setscale => setscale_UDV_state
      procedure :: getscale => getscale_UDV_state
      procedure :: testscale => testscale_UDV_state
#if defined(MPI)
      procedure :: MPI_Sendrecv => MPI_Sendrecv_UDV_state
#endif
      generic :: assignment(=) => assign
   end type UDV_State

   logical, save :: trigger_scale_warning = .true.

contains
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This function initializes the memory of an object.
!>
!> @param [inout] this class(UDV_state)
!> \verbatim The object to be  allocated \endverbatim
!> @param [in] t Integer
!> \verbatim Number of orbitals \endverbatim
!> @param [in]  t_part optional   Integer
!> \verbatim Number of particles (projective code)  or number of orbitals (finite temperature code).
!>If not present then t_part is set to t \endverbatim
!>
!-------------------------------------------------------------------
   subroutine alloc_UDV_state(this, t, t_part)
      implicit none
      class(UDV_State), intent(INOUT) :: this
      integer, intent(IN) :: t
      integer, intent(IN), optional :: t_part

      this%ndim = t
      if (present(t_part)) then
         this%N_part = t_part
         allocate (this%U(this%ndim, this%N_part))
      else
         this%N_part = t
         allocate (this%U(this%ndim, this%N_part), this%V(this%N_part, this%N_part))
      end if
#if !defined(STABLOG)
      allocate (this%D(this%N_part))
#else
      allocate (this%L(this%N_part))
#endif
   end subroutine alloc_UDV_state

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This function allocates  the memory of an object and initializes
!> the U=P and V=1 and D=1.
!>
!> @param [inout] this  class(UDV_state)
!> @param [in]  t Integer
!> \verbatim Number of orbitals \endverbatim
!> @param [in] side Character
!> @param [in] P(:,:) ,optional,  Complex
!> \verbatim If present sets this%U = P  \endverbatim
!>
!-------------------------------------------------------------------
   subroutine init_UDV_state(this, t, side, P)
      implicit none
      class(UDV_State), intent(INOUT) :: this
      integer, intent(IN) :: t
      character, intent(IN) :: side
      complex(kind=kind(0.d0)), intent(IN), optional :: P(:, :)

      this%side = side
      if (present(P)) then
         if (t .ne. size(P, 1)) then
            write (error_unit, *) "Mismatching Ndim between explicitly provided argument and implicitly provided size(P,1)"
            call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
         end if
         if (t < size(P, 2) .or. size(P, 2) < 0) then
            write (error_unit, *) "Illegal number of particles provided as size(P,2) (0 <= N_part <= Ndim)"
            call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
         end if
         call this%alloc(t, size(P, 2))
         call this%reset(side, P)
      else
         call this%alloc(t)
         call this%reset(side)
      end if
   end subroutine init_UDV_state

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This function initializes the scales of an object: this%%D(scale_idx)=scale_val
!>
!> @param [inout] this Class(UDV_State)
!> @param [in]  scale_val, Complex
!> @param [in]  scale_idx, Integer
!>
!-------------------------------------------------------------------
   subroutine setscale_UDV_state(this, scale_val, scale_idx)
      implicit none
      class(UDV_State), intent(INOUT) :: this
      complex(Kind=kind(0.d0)), intent(IN) :: scale_val
      integer, intent(IN) :: scale_idx

#if !defined(STABLOG)
      this%D(scale_idx) = scale_val
#else
      this%L(scale_idx) = log(dble(scale_val))
#endif
   end subroutine setscale_UDV_state

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This function returns the scales of an object. scale_val=this%D(scale_idx)
!>
!> @param [in] this Class(UDV_state)
!> @param [out]  scale_val, Complex
!> @param [in]  scale_idx, Integer
!-------------------------------------------------------------------
   subroutine getscale_UDV_state(this, scale_val, scale_idx)
      implicit none
      class(UDV_State), intent(IN) :: this
      complex(Kind=kind(0.d0)), intent(out) :: scale_val
      integer, intent(IN) :: scale_idx

#if !defined(STABLOG)
      scale_val = this%D(scale_idx)
#else
      scale_val = cmplx(exp(this%L(scale_idx)), 0.d0, kind(0.d0))
#endif
   end subroutine getscale_UDV_state

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This function prints a warning message if the scales in D grow too large or too small.
!> A switch to STABLOG is recommended for this scenario which usually solves it and prevents NaN's from getting generated.
!> The warning is only printed once and the test does not apply to the STABLOG option.
!>
!> @param [in] this Class(UDV_state)
!-------------------------------------------------------------------
   subroutine testscale_UDV_state(this)
      use runtime_error_mod
      implicit none
      class(UDV_State), intent(IN) :: this

#if !defined(STABLOG)
      real(Kind=kind(this%D(1))) :: dummy_dp

      ! Check if any scale is NaN
      if (any(this%D /= this%D)) then
         write (error_unit, *)
         write (error_unit, *) "Error: At least one scale is NaN."
         write (error_unit, *) "       Switch to LOG is required."
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      ! ATTENTION, the test assumes a (mostly) sorted array D [real and positive numbers]!
      ! Check if largest scale is approaching the largest representable value
      if (dble(this%D(1)) > 0.1*huge(dummy_dp) .and. trigger_scale_warning) then
         write (error_unit, *)
         write (error_unit, *) "Warning: Largest scale is approaching the largest representable value."
         write (error_unit, *) "         Consider switching to LOG."
         trigger_scale_warning = .false.
      end if
      ! Check if myVariable is approaching the smallest representable value
      if (dble(this%D(this%n_part)) < 10.0*tiny(dummy_dp) .and. trigger_scale_warning) then
         write (error_unit, *)
         write (error_unit, *) "Warning: Smallest scale is approaching the smalles representable value."
         write (error_unit, *) "         Consider switching to LOG."
         trigger_scale_warning = .false.
      end if
#endif
   end subroutine testscale_UDV_state

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This function deallocates the occupied memory.
!>
!> @param [inout] this Class(UDV_State)
!-------------------------------------------------------------------
   subroutine dealloc_UDV_state(this)
      implicit none
      class(UDV_State), intent(INOUT) :: this

      !V is only allocated in finite temperature version
      if (allocated(this%V)) deallocate (this%V)
      deallocate (this%U)
#if !defined(STABLOG)
      deallocate (this%D)
#else
      deallocate (this%L)
#endif
   end subroutine dealloc_UDV_state

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This function reinitializes the opject to
!> U=P, V=1 and D=1.  If P is not present then U=1.
!>
!> @param [inout] this Class(UDV_state)
!> @param [IN] side Character
!> @param [IN] P(:,:), optional   Complex
!-------------------------------------------------------------------
   subroutine reset_UDV_state(this, side, P)
      implicit none
      class(UDV_State), intent(INOUT) :: this
      character, intent(IN) ::side
      complex(Kind=kind(0.d0)), optional :: P(:, :)
      complex(Kind=kind(0.d0)) :: alpha, beta

      alpha = 0.d0
      beta = 1.d0
      this%side = side
      if (present(P)) then
         if (size(P, 1) .ne. this%ndim .or. size(P, 2) .ne. this%N_part) then
            call this%dealloc
            call this%alloc(size(P, 1), size(P, 2))
         end if
         call ZLACPY('A', this%ndim, this%N_part, P(1, 1), this%ndim, this%U(1, 1), this%ndim)
      else
         call ZLASET('A', this%ndim, this%ndim, alpha, beta, this%U(1, 1), this%ndim)
         call ZLASET('A', this%ndim, this%ndim, alpha, beta, this%V(1, 1), this%ndim)
      end if
#if !defined(STABLOG)
      this%D = beta
#else
      this%L = 0.d0
#endif
   end subroutine reset_UDV_state

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> A helper function to print the state.
!>
!> @param [inout] this Class(UVD_state)
!-------------------------------------------------------------------
   subroutine print_UDV_state(this)
      implicit none
      class(UDV_State), intent(IN) :: this
      integer :: i

      write (*, *) "Side = ", this%side
      write (*, *) "NDim = ", this%ndim
      write (*, *) "N_part = ", this%N_part
      do i = 1, this%ndim
         write (*, *) this%U(i, :)
      end do
      write (*, *) "======================"
      if (allocated(this%V)) then
         do i = 1, this%n_part
            write (*, *) this%V(i, :)
         end do
      else
         write (*, *) "V is only stored in finite temperature version"
      end if
      write (*, *) "======================"
#if !defined(STABLOG)
      write (*, *) this%D(:)
#else
      write (*, *) this%L(:)
#endif
   end subroutine print_UDV_state

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Assign this=src
!>
!> @param [inout] this  Class(UDV_state)
!> @param [in] src Class(UDV_state)
!-------------------------------------------------------------------
#if __INTEL_COMPILER_BUILD_DATE == 20190206 || __INTEL_COMPILER_BUILD_DATE == 20190416 || __INTEL_COMPILER_BUILD_DATE == 20190815
   ! Handle bug in ifort 19.3, 19.4 and 19.5, that breaks ASSIGNMENT(=), IMPURE is an Intel keyword.
   IMPURE elemental subroutine assign_UDV_state(this, src)
#else
      subroutine assign_UDV_state(this, src)
#endif
         implicit none
         class(UDV_State), intent(INOUT) :: this
         class(UDV_State), intent(IN) :: src

         if (this%ndim .ne. src%ndim .or. this%n_part .ne. src%n_part) call this%dealloc
         this%ndim = src%ndim
         this%n_part = src%n_part
         this%side = src%side

         if (.not. allocated(this%U)) allocate (this%U(this%ndim, this%n_part))
         if (.not. allocated(this%V) .and. allocated(src%V)) allocate (this%V(this%n_part, this%n_part))
         associate (ndim => src%ndim)
            call ZLACPY('A', ndim, this%n_part, src%U(1, 1), ndim, this%U(1, 1), ndim)
            if (allocated(src%V)) call ZLACPY('A', this%n_part, this%n_part, src%V(1, 1), &
                 & this%n_part, this%V(1, 1), this%n_part)
         end associate
#if !defined(STABLOG)
         if (.not. allocated(this%D)) allocate (this%D(this%n_part))
         this%D = src%D
#else
         if (.not. allocated(this%L)) allocate (this%L(this%n_part))
         this%L = src%L
#endif
      end subroutine assign_UDV_state

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> UDV  decomposition based on QR
!>
!> @details
!> @param [inout] UDV Class(UDV_State)
      !> \verbatim
      !>  if UDV%side = r then
      !>  On input  IN  = A * UDV%D * UDV%V  and A is  an arbitrary matrix stored in UDV%U.
      !>                  UDV%D  and  UDV%V stem from previous calls to this routine.
      !>  On outut  IN  = UDV%U * UDV%D * UDV%V
      !>                  Here Det(V) = 1,  D is a real diagonal matrix, and U column orthornormal
      !>
      !>  if UDV%side = l then
      !>  On input  IN  = A * UDV%D * (UDV%V)^{dag}  and A is  an arbitrary matrix stored in UDV%U.
      !>                  UDV%D and  UDV%V stem from previous calls to this routine.
      !>  On outut  IN  = UDV%U * UDV%D * (UDV%V)^{dag}
      !>                  Here Det(V) = 1,  D is a real diagonal matrix, and U column orthornormal
      !> \endverbatim
!>
!-------------------------------------------------------------------
      subroutine decompose_UDV_state(UDVR)
         use QDRP_mod
         use MyMats
         implicit none
         class(UDV_State), intent(inout) :: UDVR
         complex(Kind=kind(0.d0)), allocatable, dimension(:) :: TAU, WORK
         complex(Kind=kind(0.d0)) ::  Z_ONE, beta, phase
         integer :: INFO, i, LWORK, Ndim, N_part
         integer, allocatable, dimension(:) :: IPVT
#if defined(STABLOG)
         real(Kind=kind(0.d0)), allocatable, dimension(:) :: tmpnorm
         real(Kind=kind(0.d0)) :: tmpL, DZNRM2
         integer :: J, PVT
         complex(Kind=kind(0.d0)), allocatable, dimension(:) :: D
#else
         logical :: FORWRD
#endif

         ! QR(TMP * U * D) * V
         Z_ONE = cmplx(1.d0, 0.d0, kind(0.d0))
         Ndim = UDVR%ndim
         N_part = UDVR%n_part
         allocate (TAU(N_part), IPVT(N_part))
#if !defined(STABLOG)
         ! TMP1 = TMP1 * D
         if (allocated(UDVR%V)) then
            do i = 1, N_part
               UDVR%U(:, i) = UDVR%U(:, i)*UDVR%D(i)
            end do
         end if
         !use lapack internal pivoting
         IPVT = 0
         call QDRP_decompose(Ndim, N_part, UDVR%U, UDVR%D, IPVT, TAU, WORK, LWORK)
         Phase = cmplx(1.d0, 0.d0, kind(0.d0))
         do i = 1, N_part
            Phase = Phase*UDVR%U(i, i)
         end do
         call Pivot_Phase(phase, IPVT, N_part)
         if (udvr%side == "L" .or. udvr%side == "l") then
            Phase = conjg(Phase)
         end if
         beta = 1/Phase
         if (allocated(UDVR%V)) then
            !scale first row of R with 1/phase to set Det(R)=1 [=Det(V)]
            call ZSCAL(N_part, beta, UDVR%U(1, 1), Ndim)
            ! Permute V. Since we multiply with V from the right we have to permute the rows.
            ! A V = A P P^-1 V = Q R P^-1 V
            FORWRD = .true.
            if (udvr%side == "R" .or. udvr%side == "r") then
               call ZLAPMR(FORWRD, N_part, N_part, UDVR%V, N_part, IPVT(1)) ! lapack 3.3
            else
               call ZLAPMT(FORWRD, N_part, N_part, UDVR%V, N_part, IPVT(1))
            end if
         end if
#else
         !manually perform pivoting (using the logscale if LOG is defined)
         allocate (tmpnorm(N_part), D(N_part))
         if (allocated(UDVR%V)) then
            do i = 1, N_part
               tmpnorm(i) = log(DZNRM2(Ndim, UDVR%U(1, I), 1)) + UDVR%L(I)
            end do
         else
            do i = 1, N_part
               tmpnorm(i) = log(DZNRM2(Ndim, UDVR%U(1, I), 1))
            end do
         end if
         ! TmpMat=UDVr%V
         ! phase=det_c(tmpmat,ndim)
         ! write(*,*) "Phase in:",phase
         Phase = cmplx(1.d0, 0.d0, kind(0.d0))
         do i = 1, N_part
            PVT = I
            do j = I + 1, N_part
               if (tmpnorm(J) > tmpnorm(PVT)) PVT = J
            end do
            IPVT(I) = PVT
            if (PVT .ne. I) then
               call ZSWAP(ndim, UDVR%U(1, PVT), 1, UDVR%U(1, I), 1)
               if (allocated(UDVR%V)) then
                  if (udvr%side == "R" .or. udvr%side == "r") then
                     call ZSWAP(N_part, UDVR%V(PVT, 1), N_part, UDVR%V(I, 1), N_part)
                  else
                     call ZSWAP(N_part, UDVR%V(1, PVT), 1, UDVR%V(1, I), 1)
                  end if
               end if
               tmpL = UDVR%L(I)
               UDVR%L(I) = UDVR%L(PVT)
               UDVR%L(PVT) = tmpL
               tmpnorm(PVT) = tmpnorm(I)
               phase = -phase
            end if
         end do
         !disable lapack internal pivoting
         IPVT = 1
         call QDRP_decompose(Ndim, N_part, UDVR%U, D, IPVT, TAU, WORK, LWORK)
         do i = 1, N_part
            Phase = Phase*UDVR%U(i, i)
         end do
         if (udvr%side == "L" .or. udvr%side == "l") then
            Phase = conjg(Phase)
         end if
         beta = 1/Phase
         if (allocated(UDVR%V)) then
            !scale first row of R with 1/phase to set Det(R)=1 [=Det(V)]
            call ZSCAL(N_part, beta, UDVR%U(1, 1), Ndim)
            do i = 1, N_part
               do j = i + 1, N_part
                  UDVR%U(i, j) = UDVR%U(i, j)*cmplx(exp(UDVR%L(j) - UDVR%L(I)), 0.d0, kind(0.d0))
               end do
               UDVR%L(I) = log(dble(D(I))) + UDVR%L(I)
            end do
         else
            do i = 1, N_part
               UDVR%L(I) = log(dble(D(I)))
            end do
         end if
         deallocate (D, tmpnorm)
#endif
         if (allocated(UDVR%V)) then
            if (UDVR%side == "R" .or. UDVR%side == "r") then
               ! V = R * V
               call ZTRMM('L', 'U', 'N', 'N', N_part, N_part, Z_ONE, UDVR%U, Ndim, UDVR%V, N_part)
            else
               ! V = V * R^dagger
               call ZTRMM('R', 'U', 'C', 'N', N_part, N_part, Z_ONE, UDVR%U, Ndim, UDVR%V, N_part)
            end if
         end if
         ! Generate explicitly U in the previously abused storage of U
         call ZUNGQR(Ndim, N_part, N_part, UDVR%U, Ndim, TAU, WORK, LWORK, INFO)
         ! scale first column of U to correct the scaling in V such that UDV is not changed
         call ZSCAL(Ndim, phase, UDVR%U(1, 1), 1)

         deallocate (TAU, WORK, IPVT)

      end subroutine decompose_UDV_state

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This function sends the UDV-State to MPI-process of rank dest
!> and replaces the state with the one it receives from MPI-process of rank source
!
!> @param [inout] this The Class(UDV_state)
!> @param [in] dest MPI-rank of process where this will be sent
!> @param [in] sendtag
!> @param [in] source MPI-rank of process from which the new state will be received
!> @param [in] recvtag
!> @param [out] STATUS
!> @param [out] IERR
      !-------------------------------------------------------------------

#if defined(MPI)
      subroutine MPI_Sendrecv_UDV_state(this, dest, sendtag, source, recvtag, STATUS, IERR)
         use mpi
         implicit none

         class(UDV_State), intent(INOUT) :: this
         integer, intent(in)  :: dest, sendtag, source, recvtag
         integer, intent(out) :: STATUS(MPI_STATUS_SIZE), IERR
         integer :: n

         n = this%ndim*this%ndim
         call MPI_Sendrecv_replace(this%U, n, MPI_COMPLEX16, dest, sendtag, &
              &                source, recvtag, MPI_COMM_WORLD, STATUS, IERR)
         call MPI_Sendrecv_replace(this%V, n, MPI_COMPLEX16, dest, sendtag, &
              &                source, recvtag, MPI_COMM_WORLD, STATUS, IERR)
#if !defined(STABLOG)
         call MPI_Sendrecv_replace(this%D, this%ndim, MPI_COMPLEX16, dest, sendtag, &
              &                source, recvtag, MPI_COMM_WORLD, STATUS, IERR)
#else
         call MPI_Sendrecv_replace(this%L, this%ndim, MPI_REAL8, dest, sendtag, &
              &                source, recvtag, MPI_COMM_WORLD, STATUS, IERR)
#endif
      end subroutine MPI_Sendrecv_UDV_state
#endif

      end module UDV_State_mod
