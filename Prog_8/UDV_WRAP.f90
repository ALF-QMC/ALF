!  Copyright (C) 2016 The ALF project
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
!       http://alf.physik.uni-wuerzburg.de .
!       
!     - We require the preservation of the above copyright notice and this license in all original files.
!     
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain 
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
! 
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.
   Module UDV_Wrap_mod

     Use MyMats
     Use Files_mod

   Contains

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

!***************************************************************

#if defined(STAB1) 
     Subroutine UDV_Wrap_Pivot(A,U,D,V,NCON,N1,N2)
       
       Implicit NONE
       COMPLEX (Kind=Kind(0.d0)), INTENT(IN),    DIMENSION(:,:) :: A
       COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:,:) :: U,V
       COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:) :: D
       INTEGER, INTENT(IN) :: NCON
       INTEGER, INTENT(IN) :: N1,N2
       
       ! Locals
       REAL (Kind=Kind(0.d0)) :: VHELP(N2), XNORM(N2), XMAX, XMEAN
       INTEGER :: IVPT(N2), IVPTM1(N2), I, J, K, IMAX
       COMPLEX (Kind=Kind(0.d0))  :: A1(N1,N2), A2(N1,N2), V1(N2,N2)
       
       DO I = 1,N2
          XNORM(I) = 0.D0
          DO J = 1,N1
             XNORM(I) = XNORM(I) + DBLE( A(J,I) * CONJG( A(J,I) ) )
          ENDDO
       ENDDO
       VHELP = XNORM
      
       DO I = 1,N2
          XMAX = VHELP(1)
          IMAX = 1
          DO J = 2,N2
             IF (VHELP(J).GT.XMAX)  THEN
                IMAX = J
                XMAX = VHELP(J)
             ENDIF
          ENDDO
          VHELP(IMAX) = -1.D0
          IVPTM1(IMAX)=  I
          IVPT(I) = IMAX
       ENDDO
       DO I = 1,N2
          K = IVPT(I)
          A1(:, I) = A(:, K)
       ENDDO
       
       CALL UDV_Wrap(A1,U,D,V,NCON)
       
       V1 = V
       DO I = 1,N2
          K = IVPTM1(I)
          V(:, I) = V1(:,K)
       ENDDO
       

       IF (NCON == 1) THEN
          !Check the result  A = U D V
          DO J = 1,N2
             A1(:, J) = D * V(:, J)
          ENDDO
          Write(6,*) 'Here'
          Call MMULT (A2,U,A1)
          CALL COMPARE(A,A2,XMAX,XMEAN)
          Write (6,*) 'Check afer  Pivoting', XMAX
       ENDIF
       
     End Subroutine UDV_Wrap_Pivot
#else
     Subroutine UDV_Wrap_Pivot(A,U,D,V,NCON,N1,N2)
       
       Implicit NONE
       COMPLEX (Kind=Kind(0.d0)), INTENT(IN),    DIMENSION(:,:) :: A
       COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:,:) :: U,V
       COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:) :: D
       INTEGER, INTENT(IN) :: NCON
       INTEGER, INTENT(IN) :: N1,N2
       
       ! Locals
       REAL (Kind=Kind(0.d0)) :: VHELP(N2), XNORM(N2), XMAX, XMEAN
       INTEGER :: IVPT(N2), IVPTM1(N2), I, J, K, IMAX
       COMPLEX (Kind=Kind(0.d0))  :: A1(N1,N2), A2(N1,N2), V1(N2,N2), U1(N2,N2), Z
       
       DO I = 1,N2
          XNORM(I) = 0.D0
          DO J = 1,N1
             XNORM(I) = XNORM(I) + DBLE( A(J,I) * CONJG( A(J,I) ) )
          ENDDO
       ENDDO
       VHELP = XNORM
      
       DO I = 1,N2
          XMAX = VHELP(1)
          IMAX = 1
          DO J = 2,N2
             IF (VHELP(J).GT.XMAX)  THEN
                IMAX = J
                XMAX = VHELP(J)
             ENDIF
          ENDDO
          VHELP(IMAX) = -1.D0
          IVPTM1(IMAX)=  I
          IVPT(I) = IMAX
       ENDDO
       DO I = 1,N2
          K = IVPT(I)
          A1(:, I) = A(:, K)/CMPLX(XNORM(K),0.d0,Kind(0.d0))
       ENDDO
       
       CALL UDV(A1,U,D,V1,NCON)
       

       ! Finish the pivotting.
       DO I = 1,N2
          D(I) = D(I)*XNORM(IVPT(I))
       ENDDO
       DO I = 1,N2-1
          Z = 1.D0/XNORM(IVPT(I))
          DO J = I+1,N2
             V1(I,J) = V1(I,J)*XNORM(IVPT(J))*Z
          ENDDO
       ENDDO

       DO J = 1,N2
          DO I = 1,N2
             V(I,J) = V1(I,IVPTM1(J))
          ENDDO
       ENDDO
       

       IF (NCON == 1) THEN
          !Check the result  A = U D V
          DO J = 1,N2
             A1(:, J) = D * V(:, J)
          ENDDO
          !Write(6,*) 'Here'
          Call MMULT (A2,U,A1)
          !Call MMULT (A2,U,V)
          CALL COMPARE(A,A2,XMAX,XMEAN)
          Write (6,*) 'Check afer  Pivoting', XMAX
       ENDIF

       
       
     End Subroutine UDV_Wrap_Pivot
#endif

!***************************************************************
     Subroutine UDV_Wrap(A,U,D,V,NCON)


       Implicit None
#ifdef MPI            
       INCLUDE 'mpif.h'
#endif
       COMPLEX (Kind=Kind(0.d0)), INTENT(IN),    DIMENSION(:,:) :: A
       COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:,:) :: U,V
       COMPLEX (Kind=Kind(0.d0)), INTENT(INOUT), DIMENSION(:) :: D
       INTEGER, INTENT(IN) :: NCON
       
       !Local 
       Complex (Kind=Kind(0.d0)), Allocatable ::  A1(:,:),U1(:,:)
       Integer :: N,I,J
       character (len=64) :: file_sr, File
#ifdef MPI  
       INTEGER :: STATUS(MPI_STATUS_SIZE)
       INTEGER :: Isize, Irank,Ierr
            
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
#endif
   
       File_sr = "SDV"
#ifdef MPI 
       File = File_i(File_sr, Irank)
#else
       File = File_sr
#endif 
       CALL QR(A,U,V,NCON)
       N = Size(A,1)
       Allocate (A1(N,N),U1(N,N))
       !DO I = 1,N
       !   DO J = 1,N
       !      WRITE(6,*) I,J,V(I,J)
       !   ENDDO
       !ENDDO
       A1 = V
       CALL SVD(A1,U1,D,V,NCON)
       Call MMULT(A1,U,U1)
       U = A1
       
     End Subroutine UDV_Wrap
     
   End Module UDV_Wrap_mod
   
