!  Copyright (C) 2017 The ALF project
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
        INTEGER :: ndim

        CONTAINS
            PROCEDURE :: alloc => alloc_UDV_state
            PROCEDURE :: init => init_UDV_state
            PROCEDURE :: dealloc => dealloc_UDV_state
            PROCEDURE :: reset => reset_UDV_state
            PROCEDURE :: assign => assign_UDV_state
            PROCEDURE :: left_decompose => left_decompose_UDV_state
            PROCEDURE :: right_decompose => right_decompose_UDV_state
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
SUBROUTINE alloc_UDV_state(this, t)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: t

    this%ndim = t
    ALLOCATE(this%U(this%ndim, this%ndim), this%V(this%ndim, this%ndim))
#if !defined(LOG)
    ALLOCATE(this%D(this%ndim))
#else
    ALLOCATE(this%L(this%ndim))
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
SUBROUTINE init_UDV_state(this, t)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    INTEGER :: t

    CALL this%alloc(t)
    CALL this%reset
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

    DEALLOCATE(this%U, this%V)
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
SUBROUTINE reset_UDV_state(this)
    IMPLICIT NONE
    CLASS(UDV_State), INTENT(INOUT) :: this
    COMPLEX (Kind=Kind(0.d0)) :: alpha, beta

    alpha = 0.D0
    beta = 1.D0
    CALL ZLASET('A', this%ndim, this%ndim, alpha, beta, this%U(1, 1), this%ndim)
    CALL ZLASET('A', this%ndim, this%ndim, alpha, beta, this%V(1, 1), this%ndim)
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

    WRITE(*,*) "NDim = ", this%ndim
    DO i = 1, this%ndim
        WRITE(*,*) this%U(i, :)
    ENDDO
    WRITE(*,*) "======================"
    DO i = 1, this%ndim
        WRITE(*,*) this%V(i, :)
    ENDDO
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
    
    this%ndim = src%ndim
    IF(.not. ALLOCATED(this%U)) ALLOCATE(this%U(this%ndim, this%ndim))
    IF(.not. ALLOCATED(this%V)) ALLOCATE(this%V(this%ndim, this%ndim))
    ASSOCIATE(ndim => src%ndim)
        CALL ZLACPY('A', ndim, ndim, src%U(1, 1), ndim, this%U(1, 1), ndim)
        CALL ZLACPY('A', ndim, ndim, src%V(1, 1), ndim, this%V(1, 1), ndim)
    END ASSOCIATE
#if !defined(LOG)
    this%D = src%D
#else
    this%L = src%L
#endif
END SUBROUTINE assign_UDV_state

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> This function updates the UDV matrices with the new matrix stored in TMP:
!
!> @param [inout] UDVL The UDV object which we update
!> @param [in] TMP A full matrix
!> @param [in] TMP1 temporary storage
!> @param [in] NCON wether we check.
!-------------------------------------------------------------------
 SUBROUTINE right_decompose_UDV_state(UDVL)!, TMP, TMP1, NCON)
        Use QDRP_mod
        Implicit None
!         INTEGER, intent(in) :: NCON
!         COMPLEX (Kind=Kind(0.d0)), intent(in), allocatable, Dimension(: ,:) :: TMP
!         COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, Dimension(:, :) :: TMP1
        CLASS(UDV_State), intent(inout) :: UDVL
        COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:) :: TAU, WORK, D
        REAL (Kind=Kind(0.d0)), allocatable, Dimension(:) :: tmpnorm
        REAL (Kind=Kind(0.d0)) :: tmpL, DZNRM2
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta, tmpD, phase
        INTEGER, allocatable, Dimension(:) :: IPVT
        INTEGER :: INFO, i, j, LWORK, Ndim, PVT
        LOGICAL :: FORWRD

        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        Ndim = UDVL%ndim
        ! TMP1 = TMP^dagger * U^dagger
!         CALL ZGEMM('C', 'C', Ndim, Ndim, Ndim, Z_ONE, TMP(1, 1), Ndim, UDVL%U, Ndim, beta, TMP1(1, 1), Ndim)
        ALLOCATE(TAU(Ndim), IPVT(Ndim))
#if !defined(LOG)
        ! TMP1 = TMP1 * D
        DO i = 1,NDim
            UDVL%U(:, i) = UDVL%U(:, i) * UDVL%D(i)
        ENDDO
        IPVT = 0
        call QDRP_decompose(Ndim, UDVL%U, UDVL%D, IPVT, TAU, WORK, LWORK)
        Phase=cmplx(1.d0,0.d0,kind(0.d0))
        do i=1,size(UDVL%D,1)
          Phase=Phase*UDVL%U(i,i)
        enddo
        Call Pivot_Phase(phase,IPVT,size(UDVL%D,1))
        beta=1/Phase
        !scale first row of R with 1/phase to set Det(R)=1 [=Det(V)]
        call ZSCAL(size(UDVL%D,1),beta,UDVL%U(1,1),size(UDVL%U,1))
        ! Permute V, since we multiply with V from the left we have to permute its columns
        FORWRD = .true.
        CALL ZLAPMT(FORWRD, Ndim, Ndim, UDVL%V, Ndim, IPVT)
#else
        ALLOCATE(tmpnorm(Ndim),D(Ndim))
        Do i=1,Ndim
            tmpnorm(i) = log(DZNRM2( Ndim, UDVL%U( 1, I ), 1 ))+UDVL%L(I)
        enddo
        do i=1,Ndim
            PVT = I
            do j=I+1,Ndim
              if( tmpnorm(J)>tmpnorm(PVT) ) PVT=J
            enddo
            IF( PVT.NE.I ) THEN
                CALL ZSWAP( ndim, UDVL%U( 1, PVT ), 1, UDVL%U( 1, I ), 1 )
!                 CALL ZCOPY( ndim, UDVL%U( 1, I   ), 1, UDVL%U( 1, PVT ), 1 )
                CALL ZSWAP( ndim, UDVL%V( 1, PVT ), 1, UDVL%V( 1, I ), 1 )
                tmpL=UDVL%L(I)
                UDVL%L(I)=UDVL%L(PVT)
                UDVL%L(PVT)=tmpL
                tmpnorm( PVT ) = tmpnorm( I )
!             ELSE
!                 CALL ZCOPY( ndim, UDVL%U( 1, I ), 1, UDVL%U( 1, I ), 1 )
            END IF
        enddo
        IPVT = 1
        call QDRP_decompose(Ndim, UDVL%U, D, IPVT, TAU, WORK, LWORK)
        Phase=cmplx(1.d0,0.d0,kind(0.d0))
        do i=1,size(D,1)
          Phase=Phase*UDVL%U(i,i)
        enddo
        beta=1/Phase
        !scale first row of R with 1/phase to set Det(R)=1 [=Det(V)]
        call ZSCAL(size(D,1),beta,UDVL%U(1,1),size(UDVL%U,1))
        ! Permute V, since we multiply with V from the left we have to permute its columns
        do i=1,Ndim
          do j=i+1,Ndim
            UDVL%U(i,j)=UDVL%U(i,j)*cmplx(exp(UDVL%L(j)-UDVL%L(I)),0.d0,kind(0.d0))
          enddo
          !D contains absolute values, hence imag. part is zero
          UDVL%L(I)=log(DBLE(D(I))) + UDVL%L(I)
        enddo
        DEALLOCATE(D, tmpnorm)
#endif
        ! V = V * R^dagger
        CALL ZTRMM('R', 'U', 'C', 'N', Ndim, Ndim, Z_ONE, UDVL%U, Ndim, UDVL%V, Ndim)
        ! create explicitly U in the storage already present for it
        CALL ZUNGQR(Ndim, Ndim, Ndim, UDVL%U, Ndim, TAU, WORK, LWORK, INFO)
        call ZSCAL(size(UDVL%U,1),phase,UDVL%U(1,1),1)
        UDVL%U = CONJG(TRANSPOSE(UDVL%U ))
        DEALLOCATE(TAU, WORK, IPVT)
END SUBROUTINE right_decompose_UDV_state

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
 SUBROUTINE left_decompose_UDV_state(UDVR)!, TMP, TMP1, NCON)
        Use QDRP_mod
        Implicit None
!         INTEGER, intent(in) :: NCON
!         COMPLEX (Kind=Kind(0.d0)), intent(in), allocatable, dimension(:, :) :: TMP
!         COMPLEX (Kind=Kind(0.d0)), intent(inout), allocatable, dimension(:, :) :: TMP1
        CLASS(UDV_State), intent(inout) :: UDVR
        COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:) :: TAU, WORK, D
        REAL (Kind=Kind(0.d0)), allocatable, Dimension(:) :: tmpnorm
        REAL (Kind=Kind(0.d0)) :: tmpL, DZNRM2
        COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta, tmpD, phase
        INTEGER :: INFO, i, j, LWORK, Ndim, PVT
        INTEGER, allocatable, Dimension(:) :: IPVT
        LOGICAL :: FORWRD
        
        ! QR(TMP * U * D) * V
        Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
        beta = 0.D0
        Ndim = UDVR%ndim
        ALLOCATE(TAU(Ndim), IPVT(Ndim))
#if !defined(LOG)
        ! TMP1 = TMP1 * D
        DO i = 1,NDim
            UDVR%U(:, i) = UDVR%U(:, i)*UDVR%D(i)
        ENDDO
        !use lapack internal pivoting
        IPVT = 0
        call QDRP_decompose(Ndim, UDVR%U, UDVR%D, IPVT, TAU, WORK, LWORK)
        Phase=cmplx(1.d0,0.d0,kind(0.d0))
        do i=1,size(UDVR%D,1)
          Phase=Phase*UDVR%U(i,i)
        enddo
        Call Pivot_Phase(phase,IPVT,size(UDVR%D,1))
        beta=1/Phase
        !scale first row of R with 1/phase to set Det(R)=1 [=Det(V)]
        call ZSCAL(size(UDVR%D,1),beta,UDVR%U(1,1),size(UDVR%U,1))
        ! Permute V. Since we multiply with V from the right we have to permute the rows.
        ! A V = A P P^-1 V = Q R P^-1 V
        FORWRD = .true.
        CALL ZLAPMR(FORWRD, Ndim, Ndim, UDVR%V, Ndim, IPVT(1)) ! lapack 3.3
#else
        !manually perform pivoting (using the logscale if LOG is defined)
        ALLOCATE(tmpnorm(Ndim),D(Ndim))
        Do i=1,Ndim
            tmpnorm(i) = log(DZNRM2( Ndim, UDVR%U( 1, I ), 1 ))+UDVR%L(I)
        enddo
        do i=1,Ndim
            PVT = I
            do j=I+1,Ndim
              if( tmpnorm(J)>tmpnorm(PVT) ) PVT=J
            enddo
            IF( PVT.NE.I ) THEN
                CALL ZSWAP( ndim, UDVR%U( 1, PVT ), 1, UDVR%U( 1, I ), 1 )
!                 CALL ZCOPY( ndim, UDVR%U( 1, I   ), 1, UDVR%U( 1, PVT ), 1 )
                CALL ZSWAP( ndim, UDVR%V( PVT, 1 ), Ndim, UDVR%V( I, 1 ), Ndim )
                tmpL=UDVR%L(I)
                UDVR%L(I)=UDVR%L(PVT)
                UDVR%L(PVT)=tmpL
                tmpnorm( PVT ) = tmpnorm( I )
!             ELSE
!                 CALL ZCOPY( ndim, UDVR%U( 1, I ), 1, UDVR%U( 1, I ), 1 )
            END IF
        enddo
        !disable lapack internal pivoting
        IPVT = 1
        call QDRP_decompose(Ndim, UDVR%U, D, IPVT, TAU, WORK, LWORK)
        Phase=cmplx(1.d0,0.d0,kind(0.d0))
        do i=1,size(D,1)
          Phase=Phase*UDVR%U(i,i)
        enddo
        beta=1/Phase
        !scale first row of R with 1/phase to set Det(R)=1 [=Det(V)]
        call ZSCAL(size(D,1),beta,UDVR%U(1,1),size(UDVR%U,1))
        do i=1,Ndim
          do j=i+1,Ndim
            UDVR%U(i,j)=UDVR%U(i,j)*cmplx(exp(UDVR%L(j)-UDVR%L(I)),0.d0,kind(0.d0))
          enddo
          UDVR%L(I)=log(dble(D(I))) + UDVR%L(I)
        enddo
        DEALLOCATE(D, tmpnorm)
#endif
        ! V = R * V
        CALL ZTRMM('L', 'U', 'N', 'N', Ndim, Ndim, Z_ONE, UDVR%U, Ndim, UDVR%V, Ndim)
        ! Generate explicitly U in the previously abused storage of U
        CALL ZUNGQR(Ndim, Ndim, Ndim, UDVR%U, Ndim, TAU, WORK, LWORK, INFO)
        ! scale first column of U to correct the scaling in V such that UDV is not changed
        call ZSCAL(size(UDVR%U,1),phase,UDVR%U(1,1),1)
        DEALLOCATE(TAU, WORK, IPVT)
END SUBROUTINE left_decompose_UDV_state

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
