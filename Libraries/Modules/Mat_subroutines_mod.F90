!  Copyright (C) 2016 - 2022 The ALF project
!
!  This file is part of the ALF project.
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


module mat_subroutines
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Matrix operations for operator type.
!> type.
!
!--------------------------------------------------------------------
use runtime_error_mod
implicit none

contains

subroutine ZSLGEMM(side, op, N, M1, M2, A, P, Mat)
! Small Large general matrix multiplication

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> !!!!!  Side = L
!>   M = op( P^T A P) * M
!>   Side = R
!>   M =  M * op( P^T A P)
!>   On input: P =  Op%P and A = Op%O
!> !!!!!   type
!>   op = N  -->  None
!>   op = T  -->  Transposed
!>   op = C  -->  Transposed + Complex conjugation
!> !!!!! Mat has dimensions M1,M2, N is dimension of A
!>
!--------------------------------------------------------------------
        use iso_fortran_env, only: output_unit, error_unit

        IMPLICIT NONE
        CHARACTER (1)            , INTENT(IN) :: side, op
        INTEGER                  , INTENT(IN) :: N, M1, M2
        COMPLEX (KIND=KIND(0.D0)), INTENT(IN)   , DIMENSION(N,N) :: A
        COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(M1,M2) :: Mat
        INTEGER                  , INTENT(IN)   , DIMENSION(N)   :: P

        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: WORK, WORK2
        Complex (Kind = Kind(0.D0)) :: alpha, beta, Z(8)
        INTEGER :: I, L, IDX, NUMBLOCKS, op_id
        INTEGER, DIMENSION(:), ALLOCATABLE :: IDXLIST, DIMLIST
        LOGICAL :: COMPACT, LEFT

        IF ( side == 'L' .or. side == 'l' ) THEN
          LEFT=.true.
        ELSEIF ( side == 'R' .or. side == 'r' ) THEN
          LEFT=.false.
        ELSE
          write(error_unit,*) 'ZSLGEMM: Illegal argument for side=',side,': It is not one of [R,r,L,l] !'
          Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        ENDIF

        IF ( op == 'N' .or. op == 'n' ) THEN
          op_id=0
        ELSEIF ( op == 'T' .or. op == 't' ) THEN
          op_id=1
        ELSEIF ( op == 'C' .or. op == 'c' ) THEN
          op_id=2
        ELSE
          write(error_unit,*) 'ZSLGEMM: Illegal argument for op=',op,': It is not one of [N,n,T,t,C,c] !'
          Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        ENDIF

        alpha = 1.D0
        beta = 0.D0

        !identify possible block structure
        !only used in default case for n>4
        IF(N > 8) THEN
          COMPACT = .TRUE.
          L = 1
          IDX = 1
          ALLOCATE(IDXLIST(N),DIMLIST(N))
          NUMBLOCKS=0
          DO I=1,N-1
            IF ( P(I)+1 .ne. P(I+1) ) THEN
              COMPACT = .FALSE.
              NUMBLOCKS=NUMBLOCKS+1
              IDXLIST(NUMBLOCKS)=IDX
              DIMLIST(NUMBLOCKS)=L
              IDX=IDX+L
              L=1
            ELSE
              L=L+1
            ENDIF
          ENDDO
          IF(IDX<N+1) THEN
            ! last block if need be
              NUMBLOCKS=NUMBLOCKS+1
              IDXLIST(NUMBLOCKS)=IDX
              DIMLIST(NUMBLOCKS)=L
          ENDIF
        ELSEIF(N>1) THEN
          ALLOCATE(WORK(N,N))
          IF( op_id==0) THEN
            CALL ZLACPY('A', N, N, A(1,1), N, WORK(1,1), N)
          ELSEIF( op_id==1 ) then
            DO I=1,N
            DO L=1,N
              WORK(I,L)=A(L,I)
            ENDDO
            ENDDO
          ELSE
            DO I=1,N
            DO L=1,N
              WORK(I,L)=conjg(A(L,I))
            ENDDO
            ENDDO
          ENDIF
        ENDIF

        IF ( LEFT ) THEN
          ! multiply op(A) from the left  [ Mat = op(A)*Mat ]

          SELECT CASE(N)
          CASE (1)
            ! Here only one row is rescaled
            IF(op_id == 0 .or. op_id == 1) then
              CALL ZSCAL(M2,A(1,1),Mat(P(1),1),M1)
            else
              CALL ZSCAL(M2,conjg(A(1,1)),Mat(P(1),1),M1)
            endif
          CASE (2)
            ! perform inplace matmult
!             L=M2/4
!             DO I=1,4*L-1,4
!               Z(1)=Mat(P(1),I)
!               Z(2)=Mat(P(2),I)
!               Z(3)=Mat(P(1),I+1)
!               Z(4)=Mat(P(2),I+1)
!               Z(5)=Mat(P(1),I+2)
!               Z(6)=Mat(P(2),I+2)
!               Z(7)=Mat(P(1),I+3)
!               Z(8)=Mat(P(2),I+3)
!               Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)
!               Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)
!               Mat(P(1),I+1)=WORK(1,1)*Z(3)+WORK(1,2)*Z(4)
!               Mat(P(2),I+1)=WORK(2,1)*Z(3)+WORK(2,2)*Z(4)
!               Mat(P(1),I+2)=WORK(1,1)*Z(5)+WORK(1,2)*Z(6)
!               Mat(P(2),I+2)=WORK(2,1)*Z(5)+WORK(2,2)*Z(6)
!               Mat(P(1),I+3)=WORK(1,1)*Z(7)+WORK(1,2)*Z(8)
!               Mat(P(2),I+3)=WORK(2,1)*Z(7)+WORK(2,2)*Z(8)
!             ENDDO
            DO I=1,M2!4*L,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)
            ENDDO
            DEALLOCATE(WORK)
          CASE (3)
            ! perform inplace matmult
!             L=M2/2
!             DO I=1,2*L-1,2
!               Z(1)=Mat(P(1),I)
!               Z(2)=Mat(P(2),I)
!               Z(3)=Mat(P(3),I)
!               Z(4)=Mat(P(1),I+1)
!               Z(5)=Mat(P(2),I+1)
!               Z(6)=Mat(P(3),I+1)
!               Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)
!               Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)
!               Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)
!               Mat(P(1),I+1)=WORK(1,1)*Z(4)+WORK(1,2)*Z(5)+WORK(1,3)*Z(6)
!               Mat(P(2),I+1)=WORK(2,1)*Z(4)+WORK(2,2)*Z(5)+WORK(2,3)*Z(6)
!               Mat(P(3),I+1)=WORK(3,1)*Z(4)+WORK(3,2)*Z(5)+WORK(3,3)*Z(6)
!             ENDDO
            DO I=1,M2!2*L,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)
            ENDDO
            DEALLOCATE(WORK)
          CASE (4)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)+WORK(1,4)*Z(4)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)+WORK(2,4)*Z(4)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)+WORK(3,4)*Z(4)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)+WORK(4,4)*Z(4)
            ENDDO
            DEALLOCATE(WORK)
          CASE (5)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Z(5)=Mat(P(5),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)+WORK(1,4)*Z(4)+WORK(1,5)*Z(5)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)+WORK(2,4)*Z(4)+WORK(2,5)*Z(5)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)+WORK(3,4)*Z(4)+WORK(3,5)*Z(5)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)+WORK(4,4)*Z(4)+WORK(4,5)*Z(5)
              Mat(P(5),I)=WORK(5,1)*Z(1)+WORK(5,2)*Z(2)+WORK(5,3)*Z(3)+WORK(5,4)*Z(4)+WORK(5,5)*Z(5)
            ENDDO
            DEALLOCATE(WORK)
          CASE (6)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Z(5)=Mat(P(5),I)
              Z(6)=Mat(P(6),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)&
                &+WORK(1,4)*Z(4)+WORK(1,5)*Z(5)+WORK(1,6)*Z(6)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)&
                &+WORK(2,4)*Z(4)+WORK(2,5)*Z(5)+WORK(2,6)*Z(6)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)&
                &+WORK(3,4)*Z(4)+WORK(3,5)*Z(5)+WORK(3,6)*Z(6)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)&
              &+WORK(4,4)*Z(4)+WORK(4,5)*Z(5)+WORK(4,6)*Z(6)
              Mat(P(5),I)=WORK(5,1)*Z(1)+WORK(5,2)*Z(2)+WORK(5,3)*Z(3)&
                &+WORK(5,4)*Z(4)+WORK(5,5)*Z(5)+WORK(5,6)*Z(6)
              Mat(P(6),I)=WORK(6,1)*Z(1)+WORK(6,2)*Z(2)+WORK(6,3)*Z(3)&
                &+WORK(6,4)*Z(4)+WORK(6,5)*Z(5)+WORK(6,6)*Z(6)
            ENDDO
            DEALLOCATE(WORK)
          CASE (7)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Z(5)=Mat(P(5),I)
              Z(6)=Mat(P(6),I)
              Z(7)=Mat(P(7),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)&
                &+WORK(1,4)*Z(4)+WORK(1,5)*Z(5)+WORK(1,6)*Z(6)+WORK(1,7)*Z(7)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)&
                &+WORK(2,4)*Z(4)+WORK(2,5)*Z(5)+WORK(2,6)*Z(6)+WORK(2,7)*Z(7)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)&
                &+WORK(3,4)*Z(4)+WORK(3,5)*Z(5)+WORK(3,6)*Z(6)+WORK(3,7)*Z(7)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)&
                &+WORK(4,4)*Z(4)+WORK(4,5)*Z(5)+WORK(4,6)*Z(6)+WORK(4,7)*Z(7)
              Mat(P(5),I)=WORK(5,1)*Z(1)+WORK(5,2)*Z(2)+WORK(5,3)*Z(3)&
                &+WORK(5,4)*Z(4)+WORK(5,5)*Z(5)+WORK(5,6)*Z(6)+WORK(5,7)*Z(7)
              Mat(P(6),I)=WORK(6,1)*Z(1)+WORK(6,2)*Z(2)+WORK(6,3)*Z(3)&
                &+WORK(6,4)*Z(4)+WORK(6,5)*Z(5)+WORK(6,6)*Z(6)+WORK(6,7)*Z(7)
              Mat(P(7),I)=WORK(7,1)*Z(1)+WORK(7,2)*Z(2)+WORK(7,3)*Z(3)&
                &+WORK(7,4)*Z(4)+WORK(7,5)*Z(5)+WORK(7,6)*Z(6)+WORK(7,7)*Z(7)
            ENDDO
            DEALLOCATE(WORK)
          CASE (8)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Z(5)=Mat(P(5),I)
              Z(6)=Mat(P(6),I)
              Z(7)=Mat(P(7),I)
              Z(8)=Mat(P(8),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)+WORK(1,4)*Z(4)&
                &+WORK(1,5)*Z(5)+WORK(1,6)*Z(6)+WORK(1,7)*Z(7)+WORK(1,8)*Z(8)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)+WORK(2,4)*Z(4)&
                &+WORK(2,5)*Z(5)+WORK(2,6)*Z(6)+WORK(2,7)*Z(7)+WORK(2,8)*Z(8)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)+WORK(3,4)*Z(4)&
                &+WORK(3,5)*Z(5)+WORK(3,6)*Z(6)+WORK(3,7)*Z(7)+WORK(3,8)*Z(8)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)+WORK(4,4)*Z(4)&
                &+WORK(4,5)*Z(5)+WORK(4,6)*Z(6)+WORK(4,7)*Z(7)+WORK(4,8)*Z(8)
              Mat(P(5),I)=WORK(5,1)*Z(1)+WORK(5,2)*Z(2)+WORK(5,3)*Z(3)+WORK(5,4)*Z(4)&
                &+WORK(5,5)*Z(5)+WORK(5,6)*Z(6)+WORK(5,7)*Z(7)+WORK(5,8)*Z(8)
              Mat(P(6),I)=WORK(6,1)*Z(1)+WORK(6,2)*Z(2)+WORK(6,3)*Z(3)+WORK(6,4)*Z(4)&
                &+WORK(6,5)*Z(5)+WORK(6,6)*Z(6)+WORK(6,7)*Z(7)+WORK(6,8)*Z(8)
              Mat(P(7),I)=WORK(7,1)*Z(1)+WORK(7,2)*Z(2)+WORK(7,3)*Z(3)+WORK(7,4)*Z(4)&
                &+WORK(7,5)*Z(5)+WORK(7,6)*Z(6)+WORK(7,7)*Z(7)+WORK(7,8)*Z(8)
              Mat(P(8),I)=WORK(8,1)*Z(1)+WORK(8,2)*Z(2)+WORK(8,3)*Z(3)+WORK(8,4)*Z(4)&
                &+WORK(8,5)*Z(5)+WORK(8,6)*Z(6)+WORK(8,7)*Z(7)+WORK(8,8)*Z(8)
            ENDDO
            DEALLOCATE(WORK)
          CASE DEFAULT
            ! allocate memory and copy blocks of Mat to work
            ALLOCATE(WORK(N,M2))
            DO I=1,NUMBLOCKS
              CALL ZLACPY('A', DIMLIST(I), M2, Mat(P(IDXLIST(I)),1), M1, WORK(IDXLIST(I),1), N)
            ENDDO

            ! Perform Mat multiplication
            IF(COMPACT) THEN
              !write result directly into mat
              CALL ZGEMM(op,'N', N, M2, N, alpha, A(1, 1), N, WORK(1, 1), N, beta, Mat(P(1), 1), M1)
            ELSE
              !additional space for result
              ALLOCATE(WORK2(N,M2))
              CALL ZGEMM(op,'N', N, M2, N, alpha, A(1, 1), N, WORK(1, 1), N, beta, WORK2(1, 1), N)
              !distribute result back into mat using blocks
              DO I=1,NUMBLOCKS
                CALL ZLACPY('A', DIMLIST(I), M2, WORK2(IDXLIST(I),1), N, Mat(P(IDXLIST(I)),1), M1)
              ENDDO
              !free result memory
              DEALLOCATE(WORK2)
            ENDIF
            !free memory of first mat copy
            DEALLOCATE(WORK,IDXLIST,DIMLIST)
          END SELECT

       ELSE
          ! multiply op(A) from the right [ Mat = Mat*op(A) ]

          SELECT CASE(N)
          CASE (1)
            ! Here only one column is rescaled
            IF(op_id == 0 .or. op_id == 1) then
              CALL ZSCAL(M1,A(1,1),Mat(1,P(1)),1)
            ELSE
              CALL ZSCAL(M1,conjg(A(1,1)),Mat(1,P(1)),1)
            ENDIF
          CASE (2)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)
            ENDDO
            DEALLOCATE(WORK)
          CASE (3)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)
            ENDDO
            DEALLOCATE(WORK)
          CASE (4)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)+WORK(4,1)*Z(4)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)+WORK(4,2)*Z(4)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)+WORK(4,3)*Z(4)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)+WORK(4,4)*Z(4)
            ENDDO
            DEALLOCATE(WORK)
          CASE (5)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Z(5)=Mat(I,P(5))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)+WORK(4,1)*Z(4)+WORK(5,1)*Z(5)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)+WORK(4,2)*Z(4)+WORK(5,2)*Z(5)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)+WORK(4,3)*Z(4)+WORK(5,3)*Z(5)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)+WORK(4,4)*Z(4)+WORK(5,4)*Z(5)
              Mat(I,P(5))=WORK(1,5)*Z(1)+WORK(2,5)*Z(2)+WORK(3,5)*Z(3)+WORK(4,5)*Z(4)+WORK(5,5)*Z(5)
            ENDDO
            DEALLOCATE(WORK)
          CASE (6)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Z(5)=Mat(I,P(5))
              Z(6)=Mat(I,P(6))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)&
                &+WORK(4,1)*Z(4)+WORK(5,1)*Z(5)+WORK(6,1)*Z(6)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)&
                &+WORK(4,2)*Z(4)+WORK(5,2)*Z(5)+WORK(6,2)*Z(6)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)&
                &+WORK(4,3)*Z(4)+WORK(5,3)*Z(5)+WORK(6,3)*Z(6)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)&
                &+WORK(4,4)*Z(4)+WORK(5,4)*Z(5)+WORK(6,4)*Z(6)
              Mat(I,P(5))=WORK(1,5)*Z(1)+WORK(2,5)*Z(2)+WORK(3,5)*Z(3)&
                &+WORK(4,5)*Z(4)+WORK(5,5)*Z(5)+WORK(6,5)*Z(6)
              Mat(I,P(6))=WORK(1,6)*Z(1)+WORK(2,6)*Z(2)+WORK(3,6)*Z(3)&
                &+WORK(4,6)*Z(4)+WORK(5,6)*Z(5)+WORK(6,6)*Z(6)
            ENDDO
            DEALLOCATE(WORK)
          CASE (7)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Z(5)=Mat(I,P(5))
              Z(6)=Mat(I,P(6))
              Z(7)=Mat(I,P(7))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)&
                &+WORK(4,1)*Z(4)+WORK(5,1)*Z(5)+WORK(6,1)*Z(6)+WORK(7,1)*Z(7)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)&
                &+WORK(4,2)*Z(4)+WORK(5,2)*Z(5)+WORK(6,2)*Z(6)+WORK(7,2)*Z(7)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)&
                &+WORK(4,3)*Z(4)+WORK(5,3)*Z(5)+WORK(6,3)*Z(6)+WORK(7,3)*Z(7)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)&
                &+WORK(4,4)*Z(4)+WORK(5,4)*Z(5)+WORK(6,4)*Z(6)+WORK(7,4)*Z(7)
              Mat(I,P(5))=WORK(1,5)*Z(1)+WORK(2,5)*Z(2)+WORK(3,5)*Z(3)&
                &+WORK(4,5)*Z(4)+WORK(5,5)*Z(5)+WORK(6,5)*Z(6)+WORK(7,5)*Z(7)
              Mat(I,P(6))=WORK(1,6)*Z(1)+WORK(2,6)*Z(2)+WORK(3,6)*Z(3)&
                &+WORK(4,6)*Z(4)+WORK(5,6)*Z(5)+WORK(6,6)*Z(6)+WORK(7,6)*Z(7)
              Mat(I,P(7))=WORK(1,7)*Z(1)+WORK(2,7)*Z(2)+WORK(3,7)*Z(3)&
                &+WORK(4,7)*Z(4)+WORK(5,7)*Z(5)+WORK(6,7)*Z(6)+WORK(7,7)*Z(7)
            ENDDO
            DEALLOCATE(WORK)
          CASE (8)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Z(5)=Mat(I,P(5))
              Z(6)=Mat(I,P(6))
              Z(7)=Mat(I,P(7))
              Z(8)=Mat(I,P(8))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)+WORK(4,1)*Z(4)&
                &+WORK(5,1)*Z(5)+WORK(6,1)*Z(6)+WORK(7,1)*Z(7)+WORK(8,1)*Z(8)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)+WORK(4,2)*Z(4)&
                &+WORK(5,2)*Z(5)+WORK(6,2)*Z(6)+WORK(7,2)*Z(7)+WORK(8,2)*Z(8)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)+WORK(4,3)*Z(4)&
                &+WORK(5,3)*Z(5)+WORK(6,3)*Z(6)+WORK(7,3)*Z(7)+WORK(8,3)*Z(8)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)+WORK(4,4)*Z(4)&
                &+WORK(5,4)*Z(5)+WORK(6,4)*Z(6)+WORK(7,4)*Z(7)+WORK(8,4)*Z(8)
              Mat(I,P(5))=WORK(1,5)*Z(1)+WORK(2,5)*Z(2)+WORK(3,5)*Z(3)+WORK(4,5)*Z(4)&
                &+WORK(5,5)*Z(5)+WORK(6,5)*Z(6)+WORK(7,5)*Z(7)+WORK(8,5)*Z(8)
              Mat(I,P(6))=WORK(1,6)*Z(1)+WORK(2,6)*Z(2)+WORK(3,6)*Z(3)+WORK(4,6)*Z(4)&
                &+WORK(5,6)*Z(5)+WORK(6,6)*Z(6)+WORK(7,6)*Z(7)+WORK(8,6)*Z(8)
              Mat(I,P(7))=WORK(1,7)*Z(1)+WORK(2,7)*Z(2)+WORK(3,7)*Z(3)+WORK(4,7)*Z(4)&
                &+WORK(5,7)*Z(5)+WORK(6,7)*Z(6)+WORK(7,7)*Z(7)+WORK(8,7)*Z(8)
              Mat(I,P(8))=WORK(1,8)*Z(1)+WORK(2,8)*Z(2)+WORK(3,8)*Z(3)+WORK(4,8)*Z(4)&
                &+WORK(5,8)*Z(5)+WORK(6,8)*Z(6)+WORK(7,8)*Z(7)+WORK(8,8)*Z(8)
            ENDDO
            DEALLOCATE(WORK)
          CASE DEFAULT
            ! allocate memory and copy blocks of Mat to work
            ALLOCATE(WORK(M1,N))
            DO I=1,NUMBLOCKS
              CALL ZLACPY('A', M1, DIMLIST(I), Mat(1,P(IDXLIST(I))), M1, WORK(1,IDXLIST(I)), M1)
            ENDDO

            ! Perform Mat multiplication
            IF(COMPACT) THEN
              !write result directly into mat
              CALL ZGEMM('N',op, M1, N, N, alpha, WORK(1, 1), M1, A(1, 1), N, beta, Mat(1, P(1)), M1)
            ELSE
              !additional space for result
              ALLOCATE(WORK2(M1,N))
              CALL ZGEMM('N',op, M1, N, N, alpha, WORK(1, 1), M1, A(1, 1), N, beta, WORK2(1, 1), M1)
              !distribute result back into mat using blocks
              DO I=1,NUMBLOCKS
                CALL ZLACPY('A', M1, DIMLIST(I), WORK2(1,IDXLIST(I)), M1, Mat(1,P(IDXLIST(I))), M1)
              ENDDO
              !free result memory
              DEALLOCATE(WORK2)
            ENDIF
            !free memory of first mat copy
            DEALLOCATE(WORK,IDXLIST,DIMLIST)
          END SELECT

        ENDIF

end subroutine ZSLGEMM


subroutine ZSLHEMM(side, uplo, N, M1, M2, A, P, Mat)
! Small Large  hermitian matrix multiplication

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> P^T A P is hermitian
!> !!!!!  Side = L
!>   M = op( P^T A P) * M
!>   Side = R
!>   M =  M * op( P^T A P)
!>   On input: P =  Op%P and A = Op%O
!> !!!!!   type
!>   op = N  -->  None
!>   op = T  -->  Transposed
!>   op = C  -->  Transposed + Complex conjugation. Same as N.
!> !!!!! Mat has dimensions M1,M2
!>
!--------------------------------------------------------------------
        use iso_fortran_env, only: output_unit, error_unit



        IMPLICIT NONE
        CHARACTER (1)            , INTENT(IN) :: side, uplo
        INTEGER                  , INTENT(IN) :: N, M1, M2
        COMPLEX (KIND=KIND(0.D0)), INTENT(IN)   , DIMENSION(N,N) :: A
        COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(M1,M2) :: Mat
        INTEGER                  , INTENT(IN)   , DIMENSION(N)   :: P

        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: WORK, WORK2
        Complex (Kind = Kind(0.D0)) :: alpha, beta, Z(8)
        INTEGER :: I,L,IDX, NUMBLOCKS
        INTEGER, DIMENSION(:), ALLOCATABLE :: IDXLIST, DIMLIST
        LOGICAL :: COMPACT, ALLOC

        alpha = 1.D0
        beta = 0.D0

        !identify possible block structure
        !only used in default case for n>4
        IF(N > 8) THEN
          COMPACT = .TRUE.
          L = 1
          IDX = 1
          ALLOCATE(IDXLIST(N),DIMLIST(N))
          ALLOC = .TRUE.
          NUMBLOCKS=0
          DO I=1,N-1
            IF ( P(I)+1 .ne. P(I+1) ) THEN
              COMPACT = .FALSE.
              NUMBLOCKS=NUMBLOCKS+1
              IDXLIST(NUMBLOCKS)=IDX
              DIMLIST(NUMBLOCKS)=L
              IDX=IDX+L
              L=1
            ELSE
              L=L+1
            ENDIF
          ENDDO
          IF(IDX<N+1) THEN
            ! last block if need be
              NUMBLOCKS=NUMBLOCKS+1
              IDXLIST(NUMBLOCKS)=IDX
              DIMLIST(NUMBLOCKS)=L
          ENDIF
        ELSEIF(N>1) THEN
          ALLOCATE(WORK(N,N))
          CALL ZLACPY(uplo, N, N, A(1,1), N, WORK(1,1), N)
          ! Fill the rest of WORK (thereby ignoring uplo)
          IF(uplo=='U' .or. uplo=='u') THEN
            DO L=1,N
            DO I=1,L-1
              WORK(L,I)=conjg(WORK(I,L))
            ENDDO
            ENDDO
          ELSE
            DO L=1,N
            DO I=L+1,N
              WORK(L,I)=conjg(WORK(I,L))
            ENDDO
            ENDDO
          ENDIF
        ENDIF

        IF ( side == 'L' .or. side == 'l' ) THEN
          ! multiply op(A) from the left  [ Mat = op(A)*Mat ]

          SELECT CASE(N)
          CASE (1)
            ! Here only one row is rescaled
            ! uplo and transpositions as well as conjugation has no effect for 1x1 herm. Matrices
            CALL ZSCAL(M2,A(1,1),Mat(P(1),1),M1)
          CASE (2)
            ! perform inplace matmult
!             L=M2/4
!             DO I=1,4*L-1,4
!               Z(1)=Mat(P(1),I)
!               Z(2)=Mat(P(2),I)
!               Z(3)=Mat(P(1),I+1)
!               Z(4)=Mat(P(2),I+1)
!               Z(5)=Mat(P(1),I+2)
!               Z(6)=Mat(P(2),I+2)
!               Z(7)=Mat(P(1),I+3)
!               Z(8)=Mat(P(2),I+3)
!               Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)
!               Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)
!               Mat(P(1),I+1)=WORK(1,1)*Z(3)+WORK(1,2)*Z(4)
!               Mat(P(2),I+1)=WORK(2,1)*Z(3)+WORK(2,2)*Z(4)
!               Mat(P(1),I+2)=WORK(1,1)*Z(5)+WORK(1,2)*Z(6)
!               Mat(P(2),I+2)=WORK(2,1)*Z(5)+WORK(2,2)*Z(6)
!               Mat(P(1),I+3)=WORK(1,1)*Z(7)+WORK(1,2)*Z(8)
!               Mat(P(2),I+3)=WORK(2,1)*Z(7)+WORK(2,2)*Z(8)
!             ENDDO
            DO I=1,M2!4*L,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)
            ENDDO
            DEALLOCATE(WORK)
          CASE (3)
            ! perform inplace matmult
!             L=M2/2
!             DO I=1,2*L-1,2
!               Z(1)=Mat(P(1),I)
!               Z(2)=Mat(P(2),I)
!               Z(3)=Mat(P(3),I)
!               Z(4)=Mat(P(1),I+1)
!               Z(5)=Mat(P(2),I+1)
!               Z(6)=Mat(P(3),I+1)
!               Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)
!               Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)
!               Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)
!               Mat(P(1),I+1)=WORK(1,1)*Z(4)+WORK(1,2)*Z(5)+WORK(1,3)*Z(6)
!               Mat(P(2),I+1)=WORK(2,1)*Z(4)+WORK(2,2)*Z(5)+WORK(2,3)*Z(6)
!               Mat(P(3),I+1)=WORK(3,1)*Z(4)+WORK(3,2)*Z(5)+WORK(3,3)*Z(6)
!             ENDDO
            DO I=1,M2!2*L,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)
            ENDDO
            DEALLOCATE(WORK)
          CASE (4)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)+WORK(1,4)*Z(4)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)+WORK(2,4)*Z(4)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)+WORK(3,4)*Z(4)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)+WORK(4,4)*Z(4)
            ENDDO
            DEALLOCATE(WORK)
          CASE (5)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Z(5)=Mat(P(5),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)+WORK(1,4)*Z(4)+WORK(1,5)*Z(5)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)+WORK(2,4)*Z(4)+WORK(2,5)*Z(5)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)+WORK(3,4)*Z(4)+WORK(3,5)*Z(5)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)+WORK(4,4)*Z(4)+WORK(4,5)*Z(5)
              Mat(P(5),I)=WORK(5,1)*Z(1)+WORK(5,2)*Z(2)+WORK(5,3)*Z(3)+WORK(5,4)*Z(4)+WORK(5,5)*Z(5)
            ENDDO
            DEALLOCATE(WORK)
          CASE (6)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Z(5)=Mat(P(5),I)
              Z(6)=Mat(P(6),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)&
                &+WORK(1,4)*Z(4)+WORK(1,5)*Z(5)+WORK(1,6)*Z(6)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)&
                &+WORK(2,4)*Z(4)+WORK(2,5)*Z(5)+WORK(2,6)*Z(6)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)&
                &+WORK(3,4)*Z(4)+WORK(3,5)*Z(5)+WORK(3,6)*Z(6)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)&
                &+WORK(4,4)*Z(4)+WORK(4,5)*Z(5)+WORK(4,6)*Z(6)
              Mat(P(5),I)=WORK(5,1)*Z(1)+WORK(5,2)*Z(2)+WORK(5,3)*Z(3)&
                &+WORK(5,4)*Z(4)+WORK(5,5)*Z(5)+WORK(5,6)*Z(6)
              Mat(P(6),I)=WORK(6,1)*Z(1)+WORK(6,2)*Z(2)+WORK(6,3)*Z(3)&
                &+WORK(6,4)*Z(4)+WORK(6,5)*Z(5)+WORK(6,6)*Z(6)
            ENDDO
            DEALLOCATE(WORK)
          CASE (7)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Z(5)=Mat(P(5),I)
              Z(6)=Mat(P(6),I)
              Z(7)=Mat(P(7),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)&
                &+WORK(1,4)*Z(4)+WORK(1,5)*Z(5)+WORK(1,6)*Z(6)+WORK(1,7)*Z(7)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)&
                &+WORK(2,4)*Z(4)+WORK(2,5)*Z(5)+WORK(2,6)*Z(6)+WORK(2,7)*Z(7)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)&
                &+WORK(3,4)*Z(4)+WORK(3,5)*Z(5)+WORK(3,6)*Z(6)+WORK(3,7)*Z(7)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)&
                &+WORK(4,4)*Z(4)+WORK(4,5)*Z(5)+WORK(4,6)*Z(6)+WORK(4,7)*Z(7)
              Mat(P(5),I)=WORK(5,1)*Z(1)+WORK(5,2)*Z(2)+WORK(5,3)*Z(3)&
                &+WORK(5,4)*Z(4)+WORK(5,5)*Z(5)+WORK(5,6)*Z(6)+WORK(5,7)*Z(7)
              Mat(P(6),I)=WORK(6,1)*Z(1)+WORK(6,2)*Z(2)+WORK(6,3)*Z(3)&
                &+WORK(6,4)*Z(4)+WORK(6,5)*Z(5)+WORK(6,6)*Z(6)+WORK(6,7)*Z(7)
              Mat(P(7),I)=WORK(7,1)*Z(1)+WORK(7,2)*Z(2)+WORK(7,3)*Z(3)&
                &+WORK(7,4)*Z(4)+WORK(7,5)*Z(5)+WORK(7,6)*Z(6)+WORK(7,7)*Z(7)
            ENDDO
            DEALLOCATE(WORK)
          CASE (8)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Z(5)=Mat(P(5),I)
              Z(6)=Mat(P(6),I)
              Z(7)=Mat(P(7),I)
              Z(8)=Mat(P(8),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)+WORK(1,4)*Z(4)&
                &+WORK(1,5)*Z(5)+WORK(1,6)*Z(6)+WORK(1,7)*Z(7)+WORK(1,8)*Z(8)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)+WORK(2,4)*Z(4)&
                &+WORK(2,5)*Z(5)+WORK(2,6)*Z(6)+WORK(2,7)*Z(7)+WORK(2,8)*Z(8)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)+WORK(3,4)*Z(4)&
                &+WORK(3,5)*Z(5)+WORK(3,6)*Z(6)+WORK(3,7)*Z(7)+WORK(3,8)*Z(8)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)+WORK(4,4)*Z(4)&
                &+WORK(4,5)*Z(5)+WORK(4,6)*Z(6)+WORK(4,7)*Z(7)+WORK(4,8)*Z(8)
              Mat(P(5),I)=WORK(5,1)*Z(1)+WORK(5,2)*Z(2)+WORK(5,3)*Z(3)+WORK(5,4)*Z(4)&
                &+WORK(5,5)*Z(5)+WORK(5,6)*Z(6)+WORK(5,7)*Z(7)+WORK(5,8)*Z(8)
              Mat(P(6),I)=WORK(6,1)*Z(1)+WORK(6,2)*Z(2)+WORK(6,3)*Z(3)+WORK(6,4)*Z(4)&
                &+WORK(6,5)*Z(5)+WORK(6,6)*Z(6)+WORK(6,7)*Z(7)+WORK(6,8)*Z(8)
              Mat(P(7),I)=WORK(7,1)*Z(1)+WORK(7,2)*Z(2)+WORK(7,3)*Z(3)+WORK(7,4)*Z(4)&
                &+WORK(7,5)*Z(5)+WORK(7,6)*Z(6)+WORK(7,7)*Z(7)+WORK(7,8)*Z(8)
              Mat(P(8),I)=WORK(8,1)*Z(1)+WORK(8,2)*Z(2)+WORK(8,3)*Z(3)+WORK(8,4)*Z(4)&
                &+WORK(8,5)*Z(5)+WORK(8,6)*Z(6)+WORK(8,7)*Z(7)+WORK(8,8)*Z(8)
            ENDDO
            DEALLOCATE(WORK)
          CASE DEFAULT
            ! allocate memory and copy blocks of Mat to work
            ALLOCATE(WORK(N,M2))
            DO I=1,NUMBLOCKS
              CALL ZLACPY('A', DIMLIST(I), M2, Mat(P(IDXLIST(I)),1), M1, WORK(IDXLIST(I),1), N)
            ENDDO

            ! Perform Mat multiplication
            IF(COMPACT) THEN
              !write result directly into mat
              CALL ZHEMM(side, uplo, N, M2, alpha, A(1, 1), N, WORK(1, 1), N, beta, Mat(P(1), 1), M1)
            ELSE
              !additional space for result
              ALLOCATE(WORK2(N,M2))
              CALL ZHEMM(side, uplo, N, M2, alpha, A(1, 1), N, WORK(1, 1), N, beta, WORK2(1, 1), N)
              !distribute result back into mat using blocks
              DO I=1,NUMBLOCKS
                CALL ZLACPY('A', DIMLIST(I), M2, WORK2(IDXLIST(I),1), N, Mat(P(IDXLIST(I)),1), M1)
              ENDDO
              !free result memory
              DEALLOCATE(WORK2)
            ENDIF
            !free memory of first mat copy
            DEALLOCATE(WORK,IDXLIST,DIMLIST)
          END SELECT

        ELSEIF ( side == 'R' .or. side == 'r' ) THEN
          ! multiply op(A) from the right [ Mat = Mat*op(A) ]

          SELECT CASE(N)
          CASE (1)
            ! Here only one column is rescaled
            ! uplo and transposition as well as conjugation has no effect for 1x1 herm. matrices
            CALL ZSCAL(M1,A(1,1),Mat(1,P(1)),1)
          CASE (2)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)
            ENDDO
            DEALLOCATE(WORK)
          CASE (3)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)
            ENDDO
            DEALLOCATE(WORK)
          CASE (4)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)+WORK(4,1)*Z(4)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)+WORK(4,2)*Z(4)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)+WORK(4,3)*Z(4)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)+WORK(4,4)*Z(4)
            ENDDO
            DEALLOCATE(WORK)
          CASE (5)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Z(5)=Mat(I,P(5))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)+WORK(4,1)*Z(4)+WORK(5,1)*Z(5)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)+WORK(4,2)*Z(4)+WORK(5,2)*Z(5)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)+WORK(4,3)*Z(4)+WORK(5,3)*Z(5)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)+WORK(4,4)*Z(4)+WORK(5,4)*Z(5)
              Mat(I,P(5))=WORK(1,5)*Z(1)+WORK(2,5)*Z(2)+WORK(3,5)*Z(3)+WORK(4,5)*Z(4)+WORK(5,5)*Z(5)
            ENDDO
            DEALLOCATE(WORK)
          CASE (6)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Z(5)=Mat(I,P(5))
              Z(6)=Mat(I,P(6))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)&
                &+WORK(4,1)*Z(4)+WORK(5,1)*Z(5)+WORK(6,1)*Z(6)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)&
                &+WORK(4,2)*Z(4)+WORK(5,2)*Z(5)+WORK(6,2)*Z(6)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)&
                &+WORK(4,3)*Z(4)+WORK(5,3)*Z(5)+WORK(6,3)*Z(6)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)&
                &+WORK(4,4)*Z(4)+WORK(5,4)*Z(5)+WORK(6,4)*Z(6)
              Mat(I,P(5))=WORK(1,5)*Z(1)+WORK(2,5)*Z(2)+WORK(3,5)*Z(3)&
                &+WORK(4,5)*Z(4)+WORK(5,5)*Z(5)+WORK(6,5)*Z(6)
              Mat(I,P(6))=WORK(1,6)*Z(1)+WORK(2,6)*Z(2)+WORK(3,6)*Z(3)&
                &+WORK(4,6)*Z(4)+WORK(5,6)*Z(5)+WORK(6,6)*Z(6)
            ENDDO
            DEALLOCATE(WORK)
          CASE (7)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Z(5)=Mat(I,P(5))
              Z(6)=Mat(I,P(6))
              Z(7)=Mat(I,P(7))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)&
                &+WORK(4,1)*Z(4)+WORK(5,1)*Z(5)+WORK(6,1)*Z(6)+WORK(7,1)*Z(7)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)&
                &+WORK(4,2)*Z(4)+WORK(5,2)*Z(5)+WORK(6,2)*Z(6)+WORK(7,2)*Z(7)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)&
                &+WORK(4,3)*Z(4)+WORK(5,3)*Z(5)+WORK(6,3)*Z(6)+WORK(7,3)*Z(7)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)&
                &+WORK(4,4)*Z(4)+WORK(5,4)*Z(5)+WORK(6,4)*Z(6)+WORK(7,4)*Z(7)
              Mat(I,P(5))=WORK(1,5)*Z(1)+WORK(2,5)*Z(2)+WORK(3,5)*Z(3)&
                &+WORK(4,5)*Z(4)+WORK(5,5)*Z(5)+WORK(6,5)*Z(6)+WORK(7,5)*Z(7)
              Mat(I,P(6))=WORK(1,6)*Z(1)+WORK(2,6)*Z(2)+WORK(3,6)*Z(3)&
                &+WORK(4,6)*Z(4)+WORK(5,6)*Z(5)+WORK(6,6)*Z(6)+WORK(7,6)*Z(7)
              Mat(I,P(7))=WORK(1,7)*Z(1)+WORK(2,7)*Z(2)+WORK(3,7)*Z(3)&
                &+WORK(4,7)*Z(4)+WORK(5,7)*Z(5)+WORK(6,7)*Z(6)+WORK(7,7)*Z(7)
            ENDDO
            DEALLOCATE(WORK)
          CASE (8)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Z(5)=Mat(I,P(5))
              Z(6)=Mat(I,P(6))
              Z(7)=Mat(I,P(7))
              Z(8)=Mat(I,P(8))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)+WORK(4,1)*Z(4)&
                &+WORK(5,1)*Z(5)+WORK(6,1)*Z(6)+WORK(7,1)*Z(7)+WORK(8,1)*Z(8)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)+WORK(4,2)*Z(4)&
                &+WORK(5,2)*Z(5)+WORK(6,2)*Z(6)+WORK(7,2)*Z(7)+WORK(8,2)*Z(8)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)+WORK(4,3)*Z(4)&
                &+WORK(5,3)*Z(5)+WORK(6,3)*Z(6)+WORK(7,3)*Z(7)+WORK(8,3)*Z(8)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)+WORK(4,4)*Z(4)&
                &+WORK(5,4)*Z(5)+WORK(6,4)*Z(6)+WORK(7,4)*Z(7)+WORK(8,4)*Z(8)
              Mat(I,P(5))=WORK(1,5)*Z(1)+WORK(2,5)*Z(2)+WORK(3,5)*Z(3)+WORK(4,5)*Z(4)&
                &+WORK(5,5)*Z(5)+WORK(6,5)*Z(6)+WORK(7,5)*Z(7)+WORK(8,5)*Z(8)
              Mat(I,P(6))=WORK(1,6)*Z(1)+WORK(2,6)*Z(2)+WORK(3,6)*Z(3)+WORK(4,6)*Z(4)&
                &+WORK(5,6)*Z(5)+WORK(6,6)*Z(6)+WORK(7,6)*Z(7)+WORK(8,6)*Z(8)
              Mat(I,P(7))=WORK(1,7)*Z(1)+WORK(2,7)*Z(2)+WORK(3,7)*Z(3)+WORK(4,7)*Z(4)&
                &+WORK(5,7)*Z(5)+WORK(6,7)*Z(6)+WORK(7,7)*Z(7)+WORK(8,7)*Z(8)
              Mat(I,P(8))=WORK(1,8)*Z(1)+WORK(2,8)*Z(2)+WORK(3,8)*Z(3)+WORK(4,8)*Z(4)&
                &+WORK(5,8)*Z(5)+WORK(6,8)*Z(6)+WORK(7,8)*Z(7)+WORK(8,8)*Z(8)
            ENDDO
            DEALLOCATE(WORK)
          CASE DEFAULT
            ! allocate memory and copy blocks of Mat to work
            ALLOCATE(WORK(M1,N))
            DO I=1,NUMBLOCKS
              CALL ZLACPY('A', M1, DIMLIST(I), Mat(1,P(IDXLIST(I))), M1, WORK(1,IDXLIST(I)), M1)
            ENDDO

            ! Perform Mat multiplication
            IF(COMPACT) THEN
              !write result directly into mat
              CALL ZHEMM(side, uplo, M1, N, alpha, A(1, 1), N, WORK(1, 1), M1, beta, Mat(1, P(1)), M1)
            ELSE
              !additional space for result
              ALLOCATE(WORK2(M1,N))
              CALL ZHEMM(side, uplo, M1, N, alpha, A(1, 1), N, WORK(1, 1), M1, beta, WORK2(1, 1), M1)
              !distribute result back into mat using blocks
              DO I=1,NUMBLOCKS
                CALL ZLACPY('A', M1, DIMLIST(I), WORK2(1,IDXLIST(I)), M1, Mat(1,P(IDXLIST(I))), M1)
              ENDDO
              !free result memory
              DEALLOCATE(WORK2)
            ENDIF
            !free memory of first mat copy
            DEALLOCATE(WORK,IDXLIST,DIMLIST)
          END SELECT

        ELSE
          write(error_unit,*) 'ZSLHEMM: Illegal argument for side: It is not one of [R,r,L,l] !'
          Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        ENDIF

end subroutine ZSLHEMM


subroutine ZDSLSYMM(side, uplo, N, M1, M2, A, P, Mat)
! Mixed Complex-Real Small Large symmetric matrix multiplication

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> P^T A P is symmetric
!> !!!!!  Side = L
!>   M = op( P^T A P) * M
!>   Side = R
!>   M =  M * op( P^T A P)
!>   On input: P =  Op%P and A = Op%O
!> !!!!!   type
!>   op = N  -->  None
!>   op = T  -->  Transposed
!>   op = C  -->  Transposed + Complex conjugation. Same as N.
!> !!!!! Mat has dimensions M1,M2
!>
!--------------------------------------------------------------------
        use iso_fortran_env, only: output_unit, error_unit

        IMPLICIT NONE
        CHARACTER (1)            , INTENT(IN) :: side, uplo
        INTEGER                  , INTENT(IN) :: N, M1, M2
        REAL (KIND=KIND(0.D0)), INTENT(IN)   , DIMENSION(N,N) :: A
        COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(M1,M2) :: Mat
        INTEGER                  , INTENT(IN)   , DIMENSION(N)   :: P

        REAL (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: WORK
        REAL (KIND=KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: RWORK ! required for the ZLACRM calls
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: WORK2, WORKCMPLX
        Complex (Kind = Kind(0.D0)) :: alpha, beta, Z(8)
        INTEGER :: I,L,IDX, NUMBLOCKS
        INTEGER, DIMENSION(:), ALLOCATABLE :: IDXLIST, DIMLIST
        LOGICAL :: COMPACT, ALLOC

        alpha = 1.D0
        beta = 0.D0

        !identify possible block structure
        !only used in default case for n>4
        IF(N > 8) THEN
          IF(uplo=='L' .or. uplo=='l') THEN ! uplo == l unimplemented and never used in this case
            Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          ENDIF
          COMPACT = .TRUE.
          L = 1
          IDX = 1
          ALLOCATE(IDXLIST(N),DIMLIST(N))
          ALLOC = .TRUE.
          NUMBLOCKS=0
          DO I=1,N-1
            IF ( P(I)+1 .ne. P(I+1) ) THEN
              COMPACT = .FALSE.
              NUMBLOCKS=NUMBLOCKS+1
              IDXLIST(NUMBLOCKS)=IDX
              DIMLIST(NUMBLOCKS)=L
              IDX=IDX+L
              L=1
            ELSE
              L=L+1
            ENDIF
          ENDDO
          IF(IDX<N+1) THEN
            ! last block if need be
              NUMBLOCKS=NUMBLOCKS+1
              IDXLIST(NUMBLOCKS)=IDX
              DIMLIST(NUMBLOCKS)=L
          ENDIF
        ELSEIF(N>1) THEN
          ALLOCATE(WORK(N,N))
          CALL DLACPY(uplo, N, N, A(1,1), N, WORK(1,1), N)
          ! Fill the rest of WORK (thereby ignoring uplo)
          IF(uplo=='U' .or. uplo=='u') THEN
            DO L=1,N
            DO I=1,L-1
              WORK(L,I)=WORK(I,L)
            ENDDO
            ENDDO
          ELSE
            DO L=1,N
            DO I=L+1,N
              WORK(L,I)=WORK(I,L)
            ENDDO
            ENDDO
          ENDIF
        ENDIF

        IF ( side == 'L' .or. side == 'l' ) THEN
          ! multiply op(A) from the left  [ Mat = op(A)*Mat ]

          SELECT CASE(N)
          CASE (1)
            ! Here only one row is rescaled
            ! uplo and transpositions have no effect for 1x1 herm. Matrices
            CALL ZDSCAL(M2,A(1,1),Mat(P(1),1),M1)
          CASE (2)
            ! perform inplace matmult
!             L=M2/4
!             DO I=1,4*L-1,4
!               Z(1)=Mat(P(1),I)
!               Z(2)=Mat(P(2),I)
!               Z(3)=Mat(P(1),I+1)
!               Z(4)=Mat(P(2),I+1)
!               Z(5)=Mat(P(1),I+2)
!               Z(6)=Mat(P(2),I+2)
!               Z(7)=Mat(P(1),I+3)
!               Z(8)=Mat(P(2),I+3)
!               Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)
!               Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)
!               Mat(P(1),I+1)=WORK(1,1)*Z(3)+WORK(1,2)*Z(4)
!               Mat(P(2),I+1)=WORK(2,1)*Z(3)+WORK(2,2)*Z(4)
!               Mat(P(1),I+2)=WORK(1,1)*Z(5)+WORK(1,2)*Z(6)
!               Mat(P(2),I+2)=WORK(2,1)*Z(5)+WORK(2,2)*Z(6)
!               Mat(P(1),I+3)=WORK(1,1)*Z(7)+WORK(1,2)*Z(8)
!               Mat(P(2),I+3)=WORK(2,1)*Z(7)+WORK(2,2)*Z(8)
!             ENDDO
            DO I=1,M2!4*L,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)
            ENDDO
            DEALLOCATE(WORK)
          CASE (3)
            ! perform inplace matmult
!             L=M2/2
!             DO I=1,2*L-1,2
!               Z(1)=Mat(P(1),I)
!               Z(2)=Mat(P(2),I)
!               Z(3)=Mat(P(3),I)
!               Z(4)=Mat(P(1),I+1)
!               Z(5)=Mat(P(2),I+1)
!               Z(6)=Mat(P(3),I+1)
!               Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)
!               Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)
!               Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)
!               Mat(P(1),I+1)=WORK(1,1)*Z(4)+WORK(1,2)*Z(5)+WORK(1,3)*Z(6)
!               Mat(P(2),I+1)=WORK(2,1)*Z(4)+WORK(2,2)*Z(5)+WORK(2,3)*Z(6)
!               Mat(P(3),I+1)=WORK(3,1)*Z(4)+WORK(3,2)*Z(5)+WORK(3,3)*Z(6)
!             ENDDO
            DO I=1,M2!2*L,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)
            ENDDO
            DEALLOCATE(WORK)
          CASE (4)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)+WORK(1,4)*Z(4)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)+WORK(2,4)*Z(4)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)+WORK(3,4)*Z(4)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)+WORK(4,4)*Z(4)
            ENDDO
            DEALLOCATE(WORK)
          CASE (5)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Z(5)=Mat(P(5),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)+WORK(1,4)*Z(4)+WORK(1,5)*Z(5)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)+WORK(2,4)*Z(4)+WORK(2,5)*Z(5)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)+WORK(3,4)*Z(4)+WORK(3,5)*Z(5)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)+WORK(4,4)*Z(4)+WORK(4,5)*Z(5)
              Mat(P(5),I)=WORK(5,1)*Z(1)+WORK(5,2)*Z(2)+WORK(5,3)*Z(3)+WORK(5,4)*Z(4)+WORK(5,5)*Z(5)
            ENDDO
            DEALLOCATE(WORK)
          CASE (6)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Z(5)=Mat(P(5),I)
              Z(6)=Mat(P(6),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)&
                &+WORK(1,4)*Z(4)+WORK(1,5)*Z(5)+WORK(1,6)*Z(6)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)&
                &+WORK(2,4)*Z(4)+WORK(2,5)*Z(5)+WORK(2,6)*Z(6)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)&
                &+WORK(3,4)*Z(4)+WORK(3,5)*Z(5)+WORK(3,6)*Z(6)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)&
                &+WORK(4,4)*Z(4)+WORK(4,5)*Z(5)+WORK(4,6)*Z(6)
              Mat(P(5),I)=WORK(5,1)*Z(1)+WORK(5,2)*Z(2)+WORK(5,3)*Z(3)&
                &+WORK(5,4)*Z(4)+WORK(5,5)*Z(5)+WORK(5,6)*Z(6)
              Mat(P(6),I)=WORK(6,1)*Z(1)+WORK(6,2)*Z(2)+WORK(6,3)*Z(3)&
                &+WORK(6,4)*Z(4)+WORK(6,5)*Z(5)+WORK(6,6)*Z(6)
            ENDDO
            DEALLOCATE(WORK)
          CASE (7)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Z(5)=Mat(P(5),I)
              Z(6)=Mat(P(6),I)
              Z(7)=Mat(P(7),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)&
                &+WORK(1,4)*Z(4)+WORK(1,5)*Z(5)+WORK(1,6)*Z(6)+WORK(1,7)*Z(7)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)&
                &+WORK(2,4)*Z(4)+WORK(2,5)*Z(5)+WORK(2,6)*Z(6)+WORK(2,7)*Z(7)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)&
                &+WORK(3,4)*Z(4)+WORK(3,5)*Z(5)+WORK(3,6)*Z(6)+WORK(3,7)*Z(7)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)&
                &+WORK(4,4)*Z(4)+WORK(4,5)*Z(5)+WORK(4,6)*Z(6)+WORK(4,7)*Z(7)
              Mat(P(5),I)=WORK(5,1)*Z(1)+WORK(5,2)*Z(2)+WORK(5,3)*Z(3)&
                &+WORK(5,4)*Z(4)+WORK(5,5)*Z(5)+WORK(5,6)*Z(6)+WORK(5,7)*Z(7)
              Mat(P(6),I)=WORK(6,1)*Z(1)+WORK(6,2)*Z(2)+WORK(6,3)*Z(3)&
                &+WORK(6,4)*Z(4)+WORK(6,5)*Z(5)+WORK(6,6)*Z(6)+WORK(6,7)*Z(7)
              Mat(P(7),I)=WORK(7,1)*Z(1)+WORK(7,2)*Z(2)+WORK(7,3)*Z(3)&
                &+WORK(7,4)*Z(4)+WORK(7,5)*Z(5)+WORK(7,6)*Z(6)+WORK(7,7)*Z(7)
            ENDDO
            DEALLOCATE(WORK)
          CASE (8)
            ! perform inplace matmult
            DO I=1,M2
              Z(1)=Mat(P(1),I)
              Z(2)=Mat(P(2),I)
              Z(3)=Mat(P(3),I)
              Z(4)=Mat(P(4),I)
              Z(5)=Mat(P(5),I)
              Z(6)=Mat(P(6),I)
              Z(7)=Mat(P(7),I)
              Z(8)=Mat(P(8),I)
              Mat(P(1),I)=WORK(1,1)*Z(1)+WORK(1,2)*Z(2)+WORK(1,3)*Z(3)+WORK(1,4)*Z(4)&
                &+WORK(1,5)*Z(5)+WORK(1,6)*Z(6)+WORK(1,7)*Z(7)+WORK(1,8)*Z(8)
              Mat(P(2),I)=WORK(2,1)*Z(1)+WORK(2,2)*Z(2)+WORK(2,3)*Z(3)+WORK(2,4)*Z(4)&
                &+WORK(2,5)*Z(5)+WORK(2,6)*Z(6)+WORK(2,7)*Z(7)+WORK(2,8)*Z(8)
              Mat(P(3),I)=WORK(3,1)*Z(1)+WORK(3,2)*Z(2)+WORK(3,3)*Z(3)+WORK(3,4)*Z(4)&
                &+WORK(3,5)*Z(5)+WORK(3,6)*Z(6)+WORK(3,7)*Z(7)+WORK(3,8)*Z(8)
              Mat(P(4),I)=WORK(4,1)*Z(1)+WORK(4,2)*Z(2)+WORK(4,3)*Z(3)+WORK(4,4)*Z(4)&
                &+WORK(4,5)*Z(5)+WORK(4,6)*Z(6)+WORK(4,7)*Z(7)+WORK(4,8)*Z(8)
              Mat(P(5),I)=WORK(5,1)*Z(1)+WORK(5,2)*Z(2)+WORK(5,3)*Z(3)+WORK(5,4)*Z(4)&
                &+WORK(5,5)*Z(5)+WORK(5,6)*Z(6)+WORK(5,7)*Z(7)+WORK(5,8)*Z(8)
              Mat(P(6),I)=WORK(6,1)*Z(1)+WORK(6,2)*Z(2)+WORK(6,3)*Z(3)+WORK(6,4)*Z(4)&
                &+WORK(6,5)*Z(5)+WORK(6,6)*Z(6)+WORK(6,7)*Z(7)+WORK(6,8)*Z(8)
              Mat(P(7),I)=WORK(7,1)*Z(1)+WORK(7,2)*Z(2)+WORK(7,3)*Z(3)+WORK(7,4)*Z(4)&
                &+WORK(7,5)*Z(5)+WORK(7,6)*Z(6)+WORK(7,7)*Z(7)+WORK(7,8)*Z(8)
              Mat(P(8),I)=WORK(8,1)*Z(1)+WORK(8,2)*Z(2)+WORK(8,3)*Z(3)+WORK(8,4)*Z(4)&
                &+WORK(8,5)*Z(5)+WORK(8,6)*Z(6)+WORK(8,7)*Z(7)+WORK(8,8)*Z(8)
            ENDDO
            DEALLOCATE(WORK)
          CASE DEFAULT
            ! allocate memory and copy blocks of Mat to work
            ALLOCATE(WORKCMPLX(N,M2))
            DO I=1,NUMBLOCKS
              CALL ZLACPY('A', DIMLIST(I), M2, Mat(P(IDXLIST(I)),1), M1, WORKCMPLX(IDXLIST(I),1), N)
            ENDDO

            ! Perform Mat multiplication
            IF(COMPACT) THEN
              !write result directly into mat
              allocate (RWORK(2*N*M2))
              call ZLARCM(N, M2, A(1, 1), N, WORKCMPLX(1,1), N, Mat(P(1), 1), M1, RWORK)
              
!              CALL ZHEMM(side, uplo, N, M2, alpha, A(1, 1), N, WORKCMPLX(1, 1), N, beta, Mat(P(1), 1), M1)
              DEALLOCATE(RWORK)
            ELSE
              !additional space for result
              ALLOCATE(WORK2(N,M2), RWORK(2*N*M2))
              call ZLARCM(N, M2, A(1, 1), N, WORKCMPLX(1,1), N, WORK2(1,1), N, RWORK)
!               CALL ZHEMM(side, uplo, N, M2, alpha, A(1, 1), N, WORKCMPLX(1, 1), N, beta, WORK2(1, 1), N)
              !distribute result back into mat using blocks
              DO I=1,NUMBLOCKS
                CALL ZLACPY('A', DIMLIST(I), M2, WORK2(IDXLIST(I),1), N, Mat(P(IDXLIST(I)),1), M1)
              ENDDO
              !free result memory
              DEALLOCATE(WORK2, RWORK)
            ENDIF
            !free memory of first mat copy
            DEALLOCATE(WORKCMPLX,IDXLIST,DIMLIST)
          END SELECT

        ELSEIF ( side == 'R' .or. side == 'r' ) THEN
          ! multiply op(A) from the right [ Mat = Mat*op(A) ]

          SELECT CASE(N)
          CASE (1)
            ! Here only one column is rescaled
            ! uplo and transpositions have no effect for 1x1 matrices
            CALL ZDSCAL(M1,A(1,1),Mat(1,P(1)),1)
          CASE (2)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)
            ENDDO
            DEALLOCATE(WORK)
          CASE (3)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)
            ENDDO
            DEALLOCATE(WORK)
          CASE (4)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)+WORK(4,1)*Z(4)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)+WORK(4,2)*Z(4)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)+WORK(4,3)*Z(4)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)+WORK(4,4)*Z(4)
            ENDDO
            DEALLOCATE(WORK)
          CASE (5)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Z(5)=Mat(I,P(5))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)+WORK(4,1)*Z(4)+WORK(5,1)*Z(5)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)+WORK(4,2)*Z(4)+WORK(5,2)*Z(5)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)+WORK(4,3)*Z(4)+WORK(5,3)*Z(5)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)+WORK(4,4)*Z(4)+WORK(5,4)*Z(5)
              Mat(I,P(5))=WORK(1,5)*Z(1)+WORK(2,5)*Z(2)+WORK(3,5)*Z(3)+WORK(4,5)*Z(4)+WORK(5,5)*Z(5)
            ENDDO
            DEALLOCATE(WORK)
          CASE (6)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Z(5)=Mat(I,P(5))
              Z(6)=Mat(I,P(6))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)&
                &+WORK(4,1)*Z(4)+WORK(5,1)*Z(5)+WORK(6,1)*Z(6)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)&
                &+WORK(4,2)*Z(4)+WORK(5,2)*Z(5)+WORK(6,2)*Z(6)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)&
                &+WORK(4,3)*Z(4)+WORK(5,3)*Z(5)+WORK(6,3)*Z(6)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)&
                &+WORK(4,4)*Z(4)+WORK(5,4)*Z(5)+WORK(6,4)*Z(6)
              Mat(I,P(5))=WORK(1,5)*Z(1)+WORK(2,5)*Z(2)+WORK(3,5)*Z(3)&
                &+WORK(4,5)*Z(4)+WORK(5,5)*Z(5)+WORK(6,5)*Z(6)
              Mat(I,P(6))=WORK(1,6)*Z(1)+WORK(2,6)*Z(2)+WORK(3,6)*Z(3)&
                &+WORK(4,6)*Z(4)+WORK(5,6)*Z(5)+WORK(6,6)*Z(6)
            ENDDO
            DEALLOCATE(WORK)
          CASE (7)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Z(5)=Mat(I,P(5))
              Z(6)=Mat(I,P(6))
              Z(7)=Mat(I,P(7))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)&
                &+WORK(4,1)*Z(4)+WORK(5,1)*Z(5)+WORK(6,1)*Z(6)+WORK(7,1)*Z(7)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)&
                &+WORK(4,2)*Z(4)+WORK(5,2)*Z(5)+WORK(6,2)*Z(6)+WORK(7,2)*Z(7)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)&
                &+WORK(4,3)*Z(4)+WORK(5,3)*Z(5)+WORK(6,3)*Z(6)+WORK(7,3)*Z(7)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)&
                &+WORK(4,4)*Z(4)+WORK(5,4)*Z(5)+WORK(6,4)*Z(6)+WORK(7,4)*Z(7)
              Mat(I,P(5))=WORK(1,5)*Z(1)+WORK(2,5)*Z(2)+WORK(3,5)*Z(3)&
                &+WORK(4,5)*Z(4)+WORK(5,5)*Z(5)+WORK(6,5)*Z(6)+WORK(7,5)*Z(7)
              Mat(I,P(6))=WORK(1,6)*Z(1)+WORK(2,6)*Z(2)+WORK(3,6)*Z(3)&
                &+WORK(4,6)*Z(4)+WORK(5,6)*Z(5)+WORK(6,6)*Z(6)+WORK(7,6)*Z(7)
              Mat(I,P(7))=WORK(1,7)*Z(1)+WORK(2,7)*Z(2)+WORK(3,7)*Z(3)&
                &+WORK(4,7)*Z(4)+WORK(5,7)*Z(5)+WORK(6,7)*Z(6)+WORK(7,7)*Z(7)
            ENDDO
            DEALLOCATE(WORK)
          CASE (8)
            ! perform inplace matmult
            DO I=1,M1
              Z(1)=Mat(I,P(1))
              Z(2)=Mat(I,P(2))
              Z(3)=Mat(I,P(3))
              Z(4)=Mat(I,P(4))
              Z(5)=Mat(I,P(5))
              Z(6)=Mat(I,P(6))
              Z(7)=Mat(I,P(7))
              Z(8)=Mat(I,P(8))
              Mat(I,P(1))=WORK(1,1)*Z(1)+WORK(2,1)*Z(2)+WORK(3,1)*Z(3)+WORK(4,1)*Z(4)&
                &+WORK(5,1)*Z(5)+WORK(6,1)*Z(6)+WORK(7,1)*Z(7)+WORK(8,1)*Z(8)
              Mat(I,P(2))=WORK(1,2)*Z(1)+WORK(2,2)*Z(2)+WORK(3,2)*Z(3)+WORK(4,2)*Z(4)&
                &+WORK(5,2)*Z(5)+WORK(6,2)*Z(6)+WORK(7,2)*Z(7)+WORK(8,2)*Z(8)
              Mat(I,P(3))=WORK(1,3)*Z(1)+WORK(2,3)*Z(2)+WORK(3,3)*Z(3)+WORK(4,3)*Z(4)&
                &+WORK(5,3)*Z(5)+WORK(6,3)*Z(6)+WORK(7,3)*Z(7)+WORK(8,3)*Z(8)
              Mat(I,P(4))=WORK(1,4)*Z(1)+WORK(2,4)*Z(2)+WORK(3,4)*Z(3)+WORK(4,4)*Z(4)&
                &+WORK(5,4)*Z(5)+WORK(6,4)*Z(6)+WORK(7,4)*Z(7)+WORK(8,4)*Z(8)
              Mat(I,P(5))=WORK(1,5)*Z(1)+WORK(2,5)*Z(2)+WORK(3,5)*Z(3)+WORK(4,5)*Z(4)&
                &+WORK(5,5)*Z(5)+WORK(6,5)*Z(6)+WORK(7,5)*Z(7)+WORK(8,5)*Z(8)
              Mat(I,P(6))=WORK(1,6)*Z(1)+WORK(2,6)*Z(2)+WORK(3,6)*Z(3)+WORK(4,6)*Z(4)&
                &+WORK(5,6)*Z(5)+WORK(6,6)*Z(6)+WORK(7,6)*Z(7)+WORK(8,6)*Z(8)
              Mat(I,P(7))=WORK(1,7)*Z(1)+WORK(2,7)*Z(2)+WORK(3,7)*Z(3)+WORK(4,7)*Z(4)&
                &+WORK(5,7)*Z(5)+WORK(6,7)*Z(6)+WORK(7,7)*Z(7)+WORK(8,7)*Z(8)
              Mat(I,P(8))=WORK(1,8)*Z(1)+WORK(2,8)*Z(2)+WORK(3,8)*Z(3)+WORK(4,8)*Z(4)&
                &+WORK(5,8)*Z(5)+WORK(6,8)*Z(6)+WORK(7,8)*Z(7)+WORK(8,8)*Z(8)
            ENDDO
            DEALLOCATE(WORK)
          CASE DEFAULT
            ! allocate memory and copy blocks of Mat to work
            ALLOCATE(WORKCMPLX(M1,N))
            DO I=1,NUMBLOCKS
              CALL ZLACPY('A', M1, DIMLIST(I), Mat(1,P(IDXLIST(I))), M1, WORKCMPLX(1,IDXLIST(I)), M1)
            ENDDO

            ! Perform Mat multiplication
            IF(COMPACT) THEN ! B*A + C
              ALLOCATE(RWORK(2*M1*N))
              !write result directly into mat
              call ZLACRM(M1, N, WORKCMPLX(1,1), M1, A(1, 1), N, Mat(1, P(1)), M1, RWORK)
!              CALL ZHEMM(side, uplo, M1, N, alpha, A(1, 1), N, WORKCMPLX(1, 1), M1, beta, Mat(1, P(1)), M1)
              DEALLOCATE(RWORK)
            ELSE
              !additional space for result
              ALLOCATE(WORK2(M1,N), RWORK(2*M1*N))
              call ZLACRM(M1, N, WORKCMPLX(1,1), M1, A(1, 1), N, WORK2(1, 1), M1, RWORK)
!              CALL ZHEMM(side, uplo, M1, N, alpha, A(1, 1), N, WORKCMPLX(1, 1), M1, beta, WORK2(1, 1), M1)
              !distribute result back into mat using blocks
              DO I=1,NUMBLOCKS
                CALL ZLACPY('A', M1, DIMLIST(I), WORK2(1,IDXLIST(I)), M1, Mat(1,P(IDXLIST(I))), M1)
              ENDDO
              !free result memory
              DEALLOCATE(WORK2, RWORK)
            ENDIF
            !free memory of first mat copy
            DEALLOCATE(WORKCMPLX,IDXLIST,DIMLIST)
          END SELECT

        ELSE
          write(error_unit,*) 'ZDSLSYMM: Illegal argument for side: It is not one of [R,r,L,l] !'
          Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        ENDIF

end subroutine ZDSLSYMM

end module
