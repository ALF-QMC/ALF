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

MODULE UDV_State_mod
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: UDV_State
    TYPE UDV_State
        COMPLEX (Kind=Kind(0.d0)), allocatable :: U(:, :), V(:, :)
#if !defined(LOG)
        COMPLEX (Kind=Kind(0.d0)), allocatable :: D(:)
#else
        REAL    (Kind=Kind(0.d0)), allocatable :: L(:)
#endif
        INTEGER :: ndim, n_part
        CHARACTER :: side

        CONTAINS
            PROCEDURE :: alloc => alloc_UDV_state
            PROCEDURE :: init => init_UDV_state
            PROCEDURE :: dealloc => dealloc_UDV_state
            PROCEDURE :: reset => reset_UDV_state
            PROCEDURE :: assign => assign_UDV_state
!             PROCEDURE :: left_decompose => left_decompose_UDV_state
            PROCEDURE :: decompose => decompose_UDV_state
            PROCEDURE :: print => print_UDV_state
            PROCEDURE :: setscale => setscale_UDV_state
            PROCEDURE :: getscale => getscale_UDV_state
#if defined(MPI)
            PROCEDURE :: MPI_Sendrecv => MPI_Sendrecv_UDV_state
#endif
            GENERIC :: ASSIGNMENT(=) => assign
    END TYPE UDV_State

CONTAINS
!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function initializes the memory of an object.
!>
!> @param [inout] this The object to be modified.
!> @param [in] t the size of the involved matrices.
!-------------------------------------------------------------------
SUBROUTINE alloc_UDV_state(this, t, t_part)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: t
    INTEGER, INTENT(IN), OPTIONAL :: t_part

    this%ndim = t
    if( present(t_part) ) then
      this%N_part=t_part
      ALLOCATE(this%U(this%ndim, this%N_part))
    else
      this%N_part=t
      ALLOCATE(this%U(this%ndim, this%N_part), this%V(this%N_part, this%N_part))
    endif
#if !defined(LOG)
    ALLOCATE(this%D(this%N_part))
#else
    ALLOCATE(this%L(this%N_part))
#endif
END SUBROUTINE alloc_UDV_state

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function initializes the memory of an object and initializes
!> the U and V matrices to unit matrices and sets the D vector to one.
!>
!> @param [inout] this The object to be modified.
!> @param [in] t the size of the involved matrices.
!-------------------------------------------------------------------
SUBROUTINE init_UDV_state(this, t, side, P)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: t
    CHARACTER, INTENT(IN) :: side
    COMPLEX(kind=kind(0.d0)), INTENT(IN), OPTIONAL :: P(:,:)
    
    this%side=side
    if( present(P)) then
      if ( t .ne. size(P,1) ) then
        write(*,*) "Mismatching Ndim between explicitly provided argument and implicitly provided size(P,1)"
        stop 1
      endif
      if ( t < size(P,2) .or. size(P,2) < 0 ) then
        write(*,*) "Illegal number of particles provided as size(P,2) (0 <= N_part <= Ndim)"
        stop 1
      endif
      CALL this%alloc(t,size(P,2))
      CALL this%reset(side,P)
    else
      CALL this%alloc(t)
      CALL this%reset(side)
    endif
END SUBROUTINE init_UDV_state

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function initializes the scales of an object.
!>
!> @param [inout] this The object to be modified.
!> @param [in] t the size of the involved matrices.
!-------------------------------------------------------------------
SUBROUTINE setscale_UDV_state(this, scale_val, scale_idx)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    COMPLEX (Kind=Kind(0.d0)), INTENT(IN) :: scale_val
    INTEGER, INTENT(IN) :: scale_idx

#if !defined(LOG)
    this%D(scale_idx)=scale_val
#else
    this%L(scale_idx)=log(dble(scale_val))
#endif
END SUBROUTINE setscale_UDV_state

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function returns the scales of an object.
!>
!> @param [inout] this The object to be modified.
!> @param [in] t the size of the involved matrices.
!-------------------------------------------------------------------
SUBROUTINE getscale_UDV_state(this, scale_val, scale_idx)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    COMPLEX (Kind=Kind(0.d0)), INTENT(out) :: scale_val
    INTEGER, INTENT(IN) :: scale_idx

#if !defined(LOG)
    scale_val=this%D(scale_idx)
#else
    scale_val=cmplx(exp(this%L(scale_idx)),0.d0,kind(0.d0))
#endif
END SUBROUTINE getscale_UDV_state

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function deallocates the occupied memory.
!>
!> @param [inout] this The object to be modified.
!-------------------------------------------------------------------
SUBROUTINE dealloc_UDV_state(this)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    
    !V is only allocated in finite temperature version
    IF(ALLOCATED(this%V)) DEALLOCATE(this%V)
    DEALLOCATE(this%U)
#if !defined(LOG)
    DEALLOCATE(this%D)
#else
    DEALLOCATE(this%L)
#endif
END SUBROUTINE dealloc_UDV_state

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function reinitializes the
!> U and V matrices to unit matrices and sets the D vector to one.
!>
!> @param [inout] this The object to be reset.
!-------------------------------------------------------------------
SUBROUTINE reset_UDV_state(this, side, P)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    CHARACTER, INTENT(IN) ::side
    COMPLEX (Kind=Kind(0.d0)), OPTIONAL :: P(:,:)
    COMPLEX (Kind=Kind(0.d0)) :: alpha, beta

    alpha = 0.D0
    beta = 1.D0
    this%side=side
    if( present(P) ) then
      if(size(P,1) .ne. this%ndim .or. size(P,2) .ne. this%N_part) then
        CALL this%dealloc
        CALL this%alloc(size(P,1),size(P,2))
      endif
      CALL ZLACPY('A', this%ndim, this%N_part, P(1, 1), this%ndim, this%U(1, 1), this%ndim)
    else
      CALL ZLASET('A', this%ndim, this%ndim, alpha, beta, this%U(1, 1), this%ndim)
      CALL ZLASET('A', this%ndim, this%ndim, alpha, beta, this%V(1, 1), this%ndim)
    endif
#if !defined(LOG)
    this%D = beta
#else
    this%L = 0.d0
#endif
END SUBROUTINE reset_UDV_state

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> A helper function to print the state.
!>
!> @param [inout] this The object to be modified.
!-------------------------------------------------------------------
SUBROUTINE print_UDV_state(this)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(IN) :: this
    INTEGER :: i

    WRITE(*,*) "Side = ", this%side
    WRITE(*,*) "NDim = ", this%ndim
    WRITE(*,*) "N_part = ", this%N_part
    DO i = 1, this%ndim
        WRITE(*,*) this%U(i, :)
    ENDDO
    WRITE(*,*) "======================"
    if( ALLOCATED(this%V)) then
      DO i = 1, this%n_part
          WRITE(*,*) this%V(i, :)
      ENDDO
    else
      WRITE(*,*) "V is only stored in finite temperature version"
    endif
    WRITE(*,*) "======================"
#if !defined(LOG)
    WRITE(*,*) this%D(:)
#else
    WRITE(*,*) this%L(:)
#endif
END SUBROUTINE print_UDV_state

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function sets this to the values stored in src,
!> thereby implementing assignment.
!>
!> @param [inout] this It gets overwritten with values in src.
!> @param [inout] src The matrices are taken from here
!-------------------------------------------------------------------
SUBROUTINE assign_UDV_state(this, src)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    CLASS(UDV_State), INTENT(IN) :: src
    
    IF(this%ndim .ne. src%ndim .or. this%n_part .ne. src%n_part) call this%dealloc
    this%ndim = src%ndim
    this%n_part = src%n_part
    this%side = src%side
    
    IF(.not. ALLOCATED(this%U)) ALLOCATE(this%U(this%ndim, this%n_part))
    IF(.not. ALLOCATED(this%V) .and. ALLOCATED(src%V)) ALLOCATE(this%V(this%n_part, this%n_part))
    ASSOCIATE(ndim => src%ndim)
        CALL ZLACPY('A', ndim, this%n_part, src%U(1, 1), ndim, this%U(1, 1), ndim)
        if (ALLOCATED(src%V)) CALL ZLACPY('A', this%n_part, this%n_part, src%V(1, 1), &
            & this%n_part, this%V(1, 1), this%n_part)
    END ASSOCIATE
#if !defined(LOG)
    IF(.not. ALLOCATED(this%D)) ALLOCATE(this%D(this%n_part))
    this%D = src%D
#else
    IF(.not. ALLOCATED(this%L)) ALLOCATE(this%L(this%n_part))
    this%L = src%L
#endif
END SUBROUTINE assign_UDV_state

! !--------------------------------------------------------------------
! !> @author 
! !> ALF-project
! !
! !> @brief 
! !> This function updates the UDV matrices with the new matrix stored in TMP:
! !
! !> @param [inout] UDV The UDV object which we update
! !> @param [in] TMP A full matrix
! !> @param [in] TMP1 temporary storage
! !> @param [in] NCON wether we check.
! !-------------------------------------------------------------------
!  SUBROUTINE right_decompose_UDV_state(UDV)!, TMP, TMP1, NCON)
!         Use QDRP_mod
!         Use MyMats
!         Implicit None
! !         INTEGER, intent(in) :: NCON
! !         COMPLEX (Kind=Kind(0.d0)), intent(in), allocatable, Dimension(: ,:) :: TMP
! !         COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, Dimension(:, :) :: TMP1
!         CLASS(UDV_State), intent(inout) :: UDV
!         COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:) :: TAU, WORK, D
!         REAL (Kind=Kind(0.d0)), allocatable, Dimension(:) :: tmpnorm
!         REAL (Kind=Kind(0.d0)) :: tmpL, DZNRM2
!         COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta, tmpD, phase, TmpMat(udv%ndim,udv%ndim)
!         INTEGER, allocatable, Dimension(:) :: IPVT
!         INTEGER :: INFO, i, j, LWORK, Ndim, PVT,n_part
!         LOGICAL :: FORWRD
!         
!         if(udv%side .ne. "L" .and. udv%side .ne. "l" ) then
!           write(*,*) "calling wrong decompose"
!         endif
! 
!         Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
!         beta = 0.D0
!         Ndim = UDV%ndim
!         N_part = UDV%n_part
!         ALLOCATE(TAU(N_part), IPVT(N_part))
! #if !defined(LOG)
!         ! TMP1 = TMP1 * D
!         DO i = 1,NDim
!             UDV%U(:, i) = UDV%U(:, i) * UDV%D(i)
!         ENDDO
!         IPVT = 0
!         call QDRP_decompose(Ndim, N_part, UDV%U, UDV%D, IPVT, TAU, WORK, LWORK)
!         Phase=cmplx(1.d0,0.d0,kind(0.d0))
!         do i=1,N_part
!           Phase=Phase*UDV%U(i,i)
!         enddo
!         Call Pivot_Phase(phase,IPVT,N_part)
!         beta=1/Phase
!         !scale first row of R with 1/phase to set Det(R)=1 [=Det(V)]
!         call ZSCAL(N_part,beta,UDV%U(1,1),size(UDV%U,1))
!         ! Permute V, since we multiply with V from the left we have to permute its columns
!         FORWRD = .true.
!         CALL ZLAPMT(FORWRD, N_part, N_part, UDV%V, N_part, IPVT)
! #else
!         ALLOCATE(tmpnorm(Ndim),D(Ndim))
!         Do i=1,Ndim
!             tmpnorm(i) = log(DZNRM2( Ndim, UDV%U( 1, I ), 1 ))+UDV%L(I)
!         enddo
! !         TmpMat=UDV%V
! !         phase=det_c(tmpmat,ndim)
! !         write(*,*) "Phase in:",phase
!         Phase=cmplx(1.d0,0.d0,kind(0.d0))
!         do i=1,Ndim
!             PVT = I
!             do j=I+1,Ndim
!               if( tmpnorm(J)>tmpnorm(PVT) ) PVT=J
!             enddo
!             IPVT(I)=PVT
!             IF( PVT.NE.I ) THEN
!                 CALL ZSWAP( ndim, UDV%U( 1, PVT ), 1, UDV%U( 1, I ), 1 )
!                 CALL ZSWAP( ndim, UDV%V( 1, PVT ), 1, UDV%V( 1, I ), 1 )
!                 tmpL=UDV%L(I)
!                 UDV%L(I)=UDV%L(PVT)
!                 UDV%L(PVT)=tmpL
!                 tmpnorm( PVT ) = tmpnorm( I )
!                 phase=-phase
!             END IF
!         enddo
! !         TmpMat=UDV%V
! !         write(*,*) "Phase after pivot:",det_c(tmpmat,ndim)
! !         write(*,*) "Phase pivot:",phase
!         IPVT = 1
!         call QDRP_decompose(Ndim, N_part, UDV%U, D, IPVT, TAU, WORK, LWORK)
!         do i=1,N_part
!           Phase=Phase*UDV%U(i,i)
!         enddo
!         Phase=CONJG(Phase)
!         beta=1/Phase
!         !scale first row of R with 1/phase to set Det(R)=1 [=Det(V)]
!         call ZSCAL(N_part,beta,UDV%U(1,1),size(UDV%U,1))
!         ! Permute V, since we multiply with V from the left we have to permute its columns
!         do i=1,Ndim
!           do j=i+1,Ndim
!             UDV%U(i,j)=UDV%U(i,j)*cmplx(exp(UDV%L(j)-UDV%L(I)),0.d0,kind(0.d0))
!           enddo
!           !D contains absolute values, hence imag. part is zero
!           UDV%L(I)=log(DBLE(D(I))) + UDV%L(I)
!         enddo
!         DEALLOCATE(D, tmpnorm)
! #endif
!         ! V = V * R^dagger
!         CALL ZTRMM('R', 'U', 'C', 'N', N_part, N_part, Z_ONE, UDV%U, Ndim, UDV%V, N_part)
!         ! create explicitly U in the storage already present for it
!         CALL ZUNGQR(Ndim, N_part, N_part, UDV%U, Ndim, TAU, WORK, LWORK, INFO)
!         call ZSCAL(size(UDV%U,1),phase,UDV%U(1,1),1)
! !         UDV%U = CONJG(TRANSPOSE(UDV%U ))
!         DEALLOCATE(TAU, WORK, IPVT)
! END SUBROUTINE right_decompose_UDV_state

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function updates the UDV matrices with the new matrix stored in TMP.
!> Essentially we calculate the product TMP * U * D * V
!> For the result we generate a new decomposition in the form U, D, V
!
!> @param [inout] UDVR The UDV object which we update
!> @param [in] TMP A full matrix
!> @param [in] TMP1 temporary storage
!> @param [in] NCON wether we check.(TODO: currently not used)
!-------------------------------------------------------------------
 SUBROUTINE decompose_UDV_state(UDVR)!, TMP, TMP1, NCON)
        Use QDRP_mod
        Use MyMats
        Implicit None
!         INTEGER, intent(in) :: NCON
!         COMPLEX (Kind=Kind(0.d0)), intent(in), allocatable, dimension(:, :) :: TMP
!         COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:, :) :: TMP1
        CLASS(UDV_State), intent(inout) :: UDVR
        COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:) :: TAU, WORK, D
        REAL (Kind=Kind(0.d0)), allocatable, Dimension(:) :: tmpnorm
        REAL (Kind=Kind(0.d0)) :: tmpL, DZNRM2
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta, tmpD, phase, TmpMat(udvr%ndim,udvr%ndim)
        INTEGER :: INFO, i, j, LWORK, Ndim, PVT, N_part
        INTEGER, allocatable, Dimension(:) :: IPVT
        LOGICAL :: FORWRD
        
!         if(udvr%side .ne. "R" .and. udvr%side .ne. "r" ) then
!           write(*,*) "calling wrong decompose"
!         endif
        
        ! QR(TMP * U * D) * V
        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        Ndim = UDVR%ndim
        N_part = UDVR%n_part
        ALLOCATE(TAU(N_part), IPVT(N_part))
#if !defined(LOG)
        ! TMP1 = TMP1 * D
        If( ALLOCATED(UDVR%V) ) then
          DO i = 1,N_part
              UDVR%U(:, i) = UDVR%U(:, i)*UDVR%D(i)
          ENDDO
        endif
        !use lapack internal pivoting
        IPVT = 0
        call QDRP_decompose(Ndim, N_part, UDVR%U, UDVR%D, IPVT, TAU, WORK, LWORK)
        Phase=cmplx(1.d0,0.d0,kind(0.d0))
        do i=1,N_part
          Phase=Phase*UDVR%U(i,i)
        enddo
        Call Pivot_Phase(phase,IPVT,N_part)
        if(udvr%side == "L" .or. udvr%side == "l" ) then
          Phase=CONJG(Phase)
        endif
        beta=1/Phase
        If( ALLOCATED(UDVR%V) ) then
          !scale first row of R with 1/phase to set Det(R)=1 [=Det(V)]
          call ZSCAL(N_part,beta,UDVR%U(1,1),Ndim)
          ! Permute V. Since we multiply with V from the right we have to permute the rows.
          ! A V = A P P^-1 V = Q R P^-1 V
          FORWRD = .true.
          if(udvr%side == "R" .or. udvr%side == "r" ) then
            CALL ZLAPMR(FORWRD, N_part, N_part, UDVR%V, N_part, IPVT(1)) ! lapack 3.3
          else
            CALL ZLAPMT(FORWRD, N_part, N_part, UDVR%V, N_part, IPVT(1))
          endif
        endif
#else
        !manually perform pivoting (using the logscale if LOG is defined)
        ALLOCATE(tmpnorm(N_part),D(N_part))
        if ( ALLOCATED(UDVR%V) ) then
          Do i=1,N_part
              tmpnorm(i) = log(DZNRM2( Ndim, UDVR%U( 1, I ), 1 ))+UDVR%L(I)
          enddo
        else
          Do i=1,N_part
              tmpnorm(i) = log(DZNRM2( Ndim, UDVR%U( 1, I ), 1 ))
          enddo
        endif
!         TmpMat=UDVr%V
!         phase=det_c(tmpmat,ndim)
!         write(*,*) "Phase in:",phase
        Phase=cmplx(1.d0,0.d0,kind(0.d0))
        do i=1,N_part
            PVT = I
            do j=I+1,N_part
              if( tmpnorm(J)>tmpnorm(PVT) ) PVT=J
            enddo
            IPVT(I)=PVT
            IF( PVT.NE.I ) THEN
                CALL ZSWAP( ndim, UDVR%U( 1, PVT ), 1, UDVR%U( 1, I ), 1 )
                If( ALLOCATED(UDVR%V) ) then
                  if(udvr%side == "R" .or. udvr%side == "r" ) then
                    CALL ZSWAP( N_part, UDVR%V( PVT, 1 ), N_part, UDVR%V( I, 1 ), N_part )
                  else
                    CALL ZSWAP( N_part, UDVR%V( 1, PVT ), 1, UDVR%V( 1, I ), 1 )
                  endif
                endif
                tmpL=UDVR%L(I)
                UDVR%L(I)=UDVR%L(PVT)
                UDVR%L(PVT)=tmpL
                tmpnorm( PVT ) = tmpnorm( I )
                phase=-phase
            END IF
        enddo
        !disable lapack internal pivoting
        IPVT = 1
        call QDRP_decompose(Ndim, N_part, UDVR%U, D, IPVT, TAU, WORK, LWORK)
        do i=1,N_part
          Phase=Phase*UDVR%U(i,i)
        enddo
        if(udvr%side == "L" .or. udvr%side == "l" ) then
          Phase=CONJG(Phase)
        endif
        beta=1/Phase
        If( ALLOCATED(UDVR%V) ) then
          !scale first row of R with 1/phase to set Det(R)=1 [=Det(V)]
          call ZSCAL(N_part,beta,UDVR%U(1,1),Ndim)
          do i=1,N_part
            do j=i+1,N_part
              UDVR%U(i,j)=UDVR%U(i,j)*cmplx(exp(UDVR%L(j)-UDVR%L(I)),0.d0,kind(0.d0))
            enddo
            UDVR%L(I)=log(dble(D(I))) + UDVR%L(I)
          enddo
        else
          do i=1,N_part
            UDVR%L(I)=log(dble(D(I)))
          enddo
        endif
        DEALLOCATE(D, tmpnorm)
#endif
        If( ALLOCATED(UDVR%V) ) then
          if(UDVR%side == "R" .or. UDVR%side == "r" ) then
            ! V = R * V
            CALL ZTRMM('L', 'U', 'N', 'N', N_part, N_part, Z_ONE, UDVR%U, Ndim, UDVR%V, N_part)
          else
            ! V = V * R^dagger
            CALL ZTRMM('R', 'U', 'C', 'N', N_part, N_part, Z_ONE, UDVR%U, Ndim, UDVR%V, N_part)
          endif
        endif
        ! Generate explicitly U in the previously abused storage of U
        CALL ZUNGQR(Ndim, N_part, N_part, UDVR%U, Ndim, TAU, WORK, LWORK, INFO)
        ! scale first column of U to correct the scaling in V such that UDV is not changed
        call ZSCAL(Ndim,phase,UDVR%U(1,1),1)
        DEALLOCATE(TAU, WORK, IPVT)
END SUBROUTINE decompose_UDV_state

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function sends the UDV-State to MPI-process of rank dest
!> and replaces the state with the one it receives from MPI-process of rank source 
!
!> @param [inout] this The object to be modified.
!> @param [in] dest MPI-rank of process where this will be sent
!> @param [in] sendtag
!> @param [in] source MPI-rank of process from which the new state will be received
!> @param [in] recvtag
!> @param [out] STATUS
!> @param [out] IERR
!-------------------------------------------------------------------
#if defined(MPI) 
 SUBROUTINE MPI_Sendrecv_UDV_state(this, dest, sendtag, source, recvtag, STATUS, IERR)
        Use mpi
        Implicit None

        CLASS(UDV_State), INTENT(INOUT) :: this
        INTEGER, intent(in)  :: dest, sendtag, source, recvtag
        Integer, intent(out) :: STATUS(MPI_STATUS_SIZE), IERR
        INTEGER :: n

        n = this%ndim * this%ndim
        CALL MPI_Sendrecv_replace(this%U, n, MPI_COMPLEX16, dest, sendtag, &
                 &                source, recvtag, MPI_COMM_WORLD, STATUS, IERR)
        CALL MPI_Sendrecv_replace(this%V, n, MPI_COMPLEX16, dest, sendtag, &
                 &                source, recvtag, MPI_COMM_WORLD, STATUS, IERR)
#if !defined(LOG)
        CALL MPI_Sendrecv_replace(this%D, this%ndim, MPI_COMPLEX16, dest, sendtag, &
                 &                source, recvtag, MPI_COMM_WORLD, STATUS, IERR)
#else
        CALL MPI_Sendrecv_replace(this%L, this%ndim, MPI_REAL8, dest, sendtag, &
                 &                source, recvtag, MPI_COMM_WORLD, STATUS, IERR)
#endif
END SUBROUTINE MPI_Sendrecv_UDV_state
#endif

END MODULE UDV_State_mod
