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


MODULE UDV_State_mod
    use runtime_error_mod
    use iso_fortran_env, only: output_unit, error_unit

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: UDV_State

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

    TYPE UDV_State
        COMPLEX (Kind=Kind(0.d0)), allocatable :: U(:, :), V(:, :)
#if !defined(STABLOG)
        COMPLEX (Kind=Kind(0.d0)), allocatable :: D(:)
#else
        REAL    (Kind=Kind(0.d0)), allocatable :: L(:)
#endif
        INTEGER   :: ndim, n_part  ! ndim: number of orbitals per flavor. n_part: number of particles per flavor
        CHARACTER :: side          ! side = R  for right propagation :   B       * P_R = U d v
                                   ! side = L  for lesft propagation :   (P_L B)^{dag} = U d v^{dag}
        CONTAINS
            PROCEDURE :: alloc => alloc_UDV_state
            PROCEDURE :: init => init_UDV_state
            PROCEDURE :: dealloc => dealloc_UDV_state
            PROCEDURE :: reset => reset_UDV_state
            PROCEDURE :: assign => assign_UDV_state
            PROCEDURE :: decompose => decompose_UDV_state
            PROCEDURE :: print => print_UDV_state
            PROCEDURE :: setscale => setscale_UDV_state
            PROCEDURE :: getscale => getscale_UDV_state
            PROCEDURE :: testscale => testscale_UDV_state
#if defined(MPI)
            PROCEDURE :: MPI_Sendrecv => MPI_Sendrecv_UDV_state
#endif
            GENERIC :: ASSIGNMENT(=) => assign
    END TYPE UDV_State

    logical, save :: trigger_scale_warning = .true.

CONTAINS
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
#if !defined(STABLOG)
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
     SUBROUTINE init_UDV_state(this, t, side, P)
       IMPLICIT NONE
       CLASS(UDV_State), INTENT(INOUT) :: this
       INTEGER,   INTENT(IN) :: t
       CHARACTER, INTENT(IN) :: side
       COMPLEX(kind=kind(0.d0)), INTENT(IN), OPTIONAL :: P(:,:)

       this%side=side
       if( present(P)) then
          if ( t .ne. size(P,1) ) then
             write(error_unit,*) "Mismatching Ndim between explicitly provided argument and implicitly provided size(P,1)"
             CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          endif
          if ( t < size(P,2) .or. size(P,2) < 0 ) then
             write(error_unit,*) "Illegal number of particles provided as size(P,2) (0 <= N_part <= Ndim)"
             CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
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
!> This function initializes the scales of an object: this%%D(scale_idx)=scale_val
!>
!> @param [inout] this Class(UDV_State)
!> @param [in]  scale_val, Complex
!> @param [in]  scale_idx, Integer
!>
!-------------------------------------------------------------------
     SUBROUTINE setscale_UDV_state(this, scale_val, scale_idx)
       IMPLICIT NONE
       CLASS(UDV_State), INTENT(INOUT) :: this
       COMPLEX (Kind=Kind(0.d0)), INTENT(IN) :: scale_val
       INTEGER, INTENT(IN) :: scale_idx

#if !defined(STABLOG)
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
!> This function returns the scales of an object. scale_val=this%D(scale_idx)
!>
!> @param [in] this Class(UDV_state)
!> @param [out]  scale_val, Complex
!> @param [in]  scale_idx, Integer
!-------------------------------------------------------------------
     SUBROUTINE getscale_UDV_state(this, scale_val, scale_idx)
       IMPLICIT NONE
       CLASS(UDV_State), INTENT(IN) :: this
       COMPLEX (Kind=Kind(0.d0)), INTENT(out) :: scale_val
       INTEGER, INTENT(IN) :: scale_idx

#if !defined(STABLOG)
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
!> This function prints a warning message if the scales in D grow too large or too small.
!> A switch to STABLOG is recommended for this scenario which usually solves it and prevents NaN's from getting generated.
!> The warning is only printed once and the test does not apply to the STABLOG option.
!>
!> @param [in] this Class(UDV_state)
!-------------------------------------------------------------------
     SUBROUTINE testscale_UDV_state(this)
       use runtime_error_mod
       IMPLICIT NONE
       CLASS(UDV_State), INTENT(IN) :: this

#if !defined(STABLOG)
       real (Kind=Kind(this%D(1))) :: dummy_dp

       ! Check if any scale is NaN
       if ( any(this%D /= this%D) )  then
          write(error_unit,*) 
          write(error_unit,*) "Error: At least one scale is NaN."
          write(error_unit,*) "       Switch to LOG is required."
          call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
       end if 

       ! ATTENTION, the test assumes a (mostly) sorted array D [real and positive numbers]!
       ! Check if largest scale is approaching the largest representable value
       if ( dble(this%D(1)) > 0.1*huge(dummy_dp) .and. trigger_scale_warning) then
          write(error_unit,*) 
          write(error_unit,*) "Warning: Largest scale is approaching the largest representable value."
          write(error_unit,*) "         Consider switching to LOG."
          trigger_scale_warning = .false.
       end if  
       ! Check if myVariable is approaching the smallest representable value
       if ( dble(this%D(this%n_part)) < 10.0*tiny(dummy_dp) .and. trigger_scale_warning) then
          write(error_unit,*) 
          write(error_unit,*) "Warning: Smallest scale is approaching the smalles representable value."
          write(error_unit,*) "         Consider switching to LOG."
          trigger_scale_warning = .false.
       end if
#endif
     END SUBROUTINE testscale_UDV_state

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This function deallocates the occupied memory.
!>
!> @param [inout] this Class(UDV_State)
!-------------------------------------------------------------------
     SUBROUTINE dealloc_UDV_state(this)
       IMPLICIT NONE
       CLASS(UDV_State), INTENT(INOUT) :: this

       !V is only allocated in finite temperature version
       IF(ALLOCATED(this%V)) DEALLOCATE(this%V)
       DEALLOCATE(this%U)
#if !defined(STABLOG)
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
!> This function reinitializes the opject to
!> U=P, V=1 and D=1.  If P is not present then U=1.
!>
!> @param [inout] this Class(UDV_state)
!> @param [IN] side Character
!> @param [IN] P(:,:), optional   Complex
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
#if !defined(STABLOG)
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
!> @param [inout] this Class(UVD_state)
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
#if !defined(STABLOG)
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
!> Assign this=src
!>
!> @param [inout] this  Class(UDV_state)
!> @param [in] src Class(UDV_state)
!-------------------------------------------------------------------
#if __INTEL_COMPILER_BUILD_DATE == 20190206 || __INTEL_COMPILER_BUILD_DATE == 20190416 || __INTEL_COMPILER_BUILD_DATE == 20190815
     ! Handle bug in ifort 19.3, 19.4 and 19.5, that breaks ASSIGNMENT(=), IMPURE is an Intel keyword.
     IMPURE ELEMENTAL SUBROUTINE assign_UDV_state(this, src)
#else
     SUBROUTINE assign_UDV_state(this, src)
#endif
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
#if !defined(STABLOG)
       IF(.not. ALLOCATED(this%D)) ALLOCATE(this%D(this%n_part))
       this%D = src%D
#else
       IF(.not. ALLOCATED(this%L)) ALLOCATE(this%L(this%n_part))
       this%L = src%L
#endif
     END SUBROUTINE assign_UDV_state

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
     SUBROUTINE decompose_UDV_state(UDVR)
       Use QDRP_mod
       Use MyMats
       Implicit None
       CLASS(UDV_State), intent(inout) :: UDVR
       COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:) :: TAU, WORK
       COMPLEX (Kind=Kind(0.d0)) ::  Z_ONE, beta, phase
       INTEGER :: INFO, i, LWORK, Ndim, N_part
       INTEGER, allocatable, Dimension(:) :: IPVT
#if defined(STABLOG)
       REAL (Kind=Kind(0.d0)), allocatable, Dimension(:) :: tmpnorm
       REAL (Kind=Kind(0.d0)) :: tmpL, DZNRM2
       INTEGER :: J, PVT
       COMPLEX (Kind=Kind(0.d0)), allocatable, Dimension(:) :: D
#else
       LOGICAL :: FORWRD
#endif

       ! QR(TMP * U * D) * V
       Z_ONE = cmplx(1.d0, 0.d0, kind(0.D0))
       Ndim = UDVR%ndim
       N_part = UDVR%n_part
       ALLOCATE(TAU(N_part), IPVT(N_part))
#if !defined(STABLOG)
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
       ! TmpMat=UDVr%V
       ! phase=det_c(tmpmat,ndim)
       ! write(*,*) "Phase in:",phase
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
!> @param [inout] this The Class(UDV_state)
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
#if !defined(STABLOG)
       CALL MPI_Sendrecv_replace(this%D, this%ndim, MPI_COMPLEX16, dest, sendtag, &
            &                source, recvtag, MPI_COMM_WORLD, STATUS, IERR)
#else
       CALL MPI_Sendrecv_replace(this%L, this%ndim, MPI_REAL8, dest, sendtag, &
            &                source, recvtag, MPI_COMM_WORLD, STATUS, IERR)
#endif
     END SUBROUTINE MPI_Sendrecv_UDV_state
#endif

   END MODULE UDV_State_mod
