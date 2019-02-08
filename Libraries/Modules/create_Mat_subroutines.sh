#!/bin/bash
sub_ZSLGEMM_LEFT()
{
    N=$1
    cat <<ABC
subroutine ZSLGEMM_LEFT_${N}(M1, M2, P, WORK, Mat)
    IMPLICIT NONE
    INTEGER                  , INTENT(IN)    :: M1, M2, P(${N})
    COMPLEX (KIND=KIND(0.D0)), INTENT(IN)    :: WORK(${N},${N})
    COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT) :: Mat(M1,M2)
    
    INTEGER :: I, J
    COMPLEX (KIND=KIND(0.D0)) :: Z(${N})
    ! perform inplace matmult
    DO I=1,M2
ABC
    for i in $(seq $N); do
        echo "      Z(${i})=Mat(P(${i}),I)"
    done

    for i in $(seq $N); do
        s="      Mat(P(${i}),I)=WORK(${i},1)*Z(1)"
        for j in $(seq 2 $N); do
            if (( $j % 4 == 1 )); then
                s="${s}&\n                &"
            fi
            s="${s}+WORK(${i},${j})*Z(${j})"
        done
        echo -e "${s}"
    done
    echo "    ENDDO"
    echo "end subroutine ZSLGEMM_LEFT_${N}"
    echo ""
}

sub_ZSLGEMM_RIGHT()
{
    N=$1
    cat <<ABC
subroutine ZSLGEMM_RIGHT_${N}(M1, M2, P, WORK, Mat)
    IMPLICIT NONE
    INTEGER                  , INTENT(IN)    :: M1, M2, P(${N})
    COMPLEX (KIND=KIND(0.D0)), INTENT(IN)    :: WORK(${N},${N})
    COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT) :: Mat(M1,M2)
    
    INTEGER :: I, J
    COMPLEX (KIND=KIND(0.D0)) :: Z(${N})
    ! perform inplace matmult
    DO I=1,M1
ABC
    for i in $(seq $N); do
        echo "      Z(${i})=Mat(I,P(${i}))"
    done

    for i in $(seq $N); do
        s="      Mat(I,P(${i}))=WORK(1,${i})*Z(1)"
        for j in $(seq 2 $N); do
            if (( $j % 4 == 1 )); then
                s="${s}&\n                &"
            fi
            s="${s}+WORK(${j},${i})*Z(${j})"
        done
        echo -e "${s}"
    done
    echo "    ENDDO"
    echo "end subroutine ZSLGEMM_RIGHT_${N}"
    echo ""
}

N_max="12"
filename="Mat_subroutines.f90"
rm -f $filename

cat <<ABC >> $filename
!  Copyright (C) 2016 - 2018 The ALF project
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


!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Matrix operations for operator type.   
!> type. 
!
!--------------------------------------------------------------------

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
!> !!!!! Mat has dimensions M1,M2
!>
!--------------------------------------------------------------------
  
        IMPLICIT NONE
        CHARACTER (1)            , INTENT(IN) :: side, op
        INTEGER                  , INTENT(IN) :: N, M1, M2
        COMPLEX (KIND=KIND(0.D0)), INTENT(IN)   , DIMENSION(N,N) :: A
        COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(M1,M2) :: Mat
        INTEGER                  , INTENT(IN)   , DIMENSION(N)   :: P
        
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: WORK, WORK2
        Complex (Kind = Kind(0.D0)) :: alpha, beta
        INTEGER :: I, L, IDX, NUMBLOCKS, op_id
        INTEGER, DIMENSION(:), ALLOCATABLE :: IDXLIST, DIMLIST
        LOGICAL :: COMPACT, LEFT
        
        IF ( side == 'L' .or. side == 'l' ) THEN
          LEFT=.true.
        ELSEIF ( side == 'R' .or. side == 'r' ) THEN
          LEFT=.false.
        ELSE
          write(*,*) 'Illegal argument for side=',side,': It is not one of [R,r,L,l] !'
          stop 2
        ENDIF
        
        IF ( op == 'N' .or. op == 'n' ) THEN
          op_id=0
        ELSEIF ( op == 'T' .or. op == 't' ) THEN
          op_id=1
        ELSEIF ( op == 'C' .or. op == 'c' ) THEN
          op_id=2
        ELSE
          write(*,*) 'Illegal argument for op=',op,': It is not one of [N,n,T,t,C,c] !'
          stop 2
        ENDIF
        
        alpha = 1.D0
        beta = 0.D0
        
        !identify possible block structure
        !only used in default case for n>4
        IF(N > ${N_max}) THEN
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
ABC

for N in $(seq 2 ${N_max}); do
    cat <<ABC >> $filename
          CASE (${N})
            call ZSLGEMM_LEFT_${N}(M1, M2, P, WORK, Mat)
            DEALLOCATE(WORK)
ABC
done

cat <<ABC >> $filename
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
ABC

for N in $(seq 2 ${N_max}); do
    cat <<ABC >> $filename
          CASE (${N})
            call ZSLGEMM_RIGHT_${N}(M1, M2, P, WORK, Mat)
            DEALLOCATE(WORK)
ABC
done

cat <<ABC >> $filename
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
  

  
        IMPLICIT NONE
        CHARACTER (1)            , INTENT(IN) :: side, uplo
        INTEGER                  , INTENT(IN) :: N, M1, M2
        COMPLEX (KIND=KIND(0.D0)), INTENT(IN)   , DIMENSION(N,N) :: A
        COMPLEX (KIND=KIND(0.D0)), INTENT(INOUT), DIMENSION(M1,M2) :: Mat
        INTEGER                  , INTENT(IN)   , DIMENSION(N)   :: P
        
        COMPLEX (KIND=KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: WORK, WORK2
        Complex (Kind = Kind(0.D0)) :: alpha, beta, Z(${N_max})
        INTEGER :: I,L,IDX, NUMBLOCKS
        INTEGER, DIMENSION(:), ALLOCATABLE :: IDXLIST, DIMLIST
        LOGICAL :: COMPACT, ALLOC
        
        alpha = 1.D0
        beta = 0.D0
        
        !identify possible block structure
        !only used in default case for n>4
        IF(N > ${N_max}) THEN
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
            ! uplo and transpositions as well as conjugatio has no effect for 1x1 herm. Matrices
            CALL ZSCAL(M2,A(1,1),Mat(P(1),1),M1)
ABC

for N in $(seq 2 ${N_max}); do
    cat <<ABC >> $filename
          CASE (${N})
            call ZSLGEMM_LEFT_${N}(M1, M2, P, WORK, Mat)
            DEALLOCATE(WORK)
ABC
done

cat <<ABC >> $filename
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
            ! uplo and transpositions as well as conjugatio has no effect for 1x1 herm. Matrices
            CALL ZSCAL(M1,A(1,1),Mat(1,P(1)),1)
ABC

for N in $(seq 2 ${N_max}); do
    cat <<ABC >> $filename
          CASE (${N})
            call ZSLGEMM_RIGHT_${N}(M1, M2, P, WORK, Mat)
            DEALLOCATE(WORK)
ABC
done


cat <<ABC >> $filename
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
          write(*,*) 'Illegal argument for side: It is not one of [R,r,L,l] !'
          stop 1
        ENDIF

end subroutine ZSLHEMM

ABC


for N in $(seq 2 ${N_max}); do
    sub_ZSLGEMM_LEFT $N >> $filename
done
for N in $(seq 2 ${N_max}); do
    sub_ZSLGEMM_RIGHT $N >> $filename
done
