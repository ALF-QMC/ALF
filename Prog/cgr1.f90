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

      SUBROUTINE CGR(PHASE,NVAR, GRUP, URUP,DRUP,VRUP, ULUP,DLUP,VLUP)

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!>    Cmputes  GRUP = (1 + UR*DR*VR*VL*DL*UL)^-1
!>    NVAR = 1 Big scales are in DL
!>    NVAR = 2 Big scales are in DR
!
!--------------------------------------------------------------------

        Use UDV_Wrap_mod

        Implicit None


	!Arguments.
        COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(IN)   ::  URUP, VRUP, ULUP, VLUP
        COMPLEX(Kind=Kind(0.d0)), Dimension(:),   Intent(In)   ::  DLUP, DRUP
        COMPLEX(Kind=Kind(0.d0)), Dimension(:,:), Intent(INOUT) :: GRUP
        COMPLEX(Kind=Kind(0.d0)) :: PHASE
        INTEGER         :: NVAR
 
        !Local
        COMPLEX (Kind=Kind(0.d0)), Dimension(:,:), Allocatable ::  UUP, VUP, TPUP, TPUP1, TPUPM1
        COMPLEX (Kind=Kind(0.d0)), Dimension(:) , Allocatable ::  DUP
        INTEGER, Dimension(:), Allocatable :: IPVT
        COMPLEX (Kind=Kind(0.d0)) ::  ZDUP1, ZDDO1, ZDUP2, ZDDO2, Z1, ZUP, ZDO, alpha, beta
        Integer :: I,J, N_size, NCON, info
        Real (Kind=Kind(0.D0)) :: X, Xmax, sv
        
        N_size = SIZE(DLUP,1)
        NCON = 0
        alpha = 1.D0
        beta = 0.D0
        Allocate( UUP(N_size,N_size),  VUP(N_size,N_size), TPUP(N_size,N_size), TPUP1(N_size,N_size), &
             & TPUPM1(N_size,N_size), DUP(N_size), IPVT(N_size) )
        !Write(6,*) 'In CGR', N_size
        CALL MMULT(VUP,VRUP,VLUP)
        DO J = 1,N_size
            TPUP(:,J) = DRUP(:)*VUP(:,J)*DLUP(J)
        ENDDO
        CALL ZGEMM('C', 'C', N_size, N_size, N_size, alpha, URUP, N_size, ULUP, N_size, alpha, TPUP, N_size)
        IF (NVAR.EQ.1) THEN
           !WRITE(6,*) 'UDV of U + DR * V * DL'
           CALL UDV_WRAP_Pivot(TPUP,UUP,DUP,VUP,NCON,N_size,N_Size)
           !CALL UDV(TPUP,UUP,DUP,VUP,NCON)
           CALL MMULT(TPUP,VUP,ULUP)
           !Do I = 1,N_size
           !   Write(6,*) DLUP(I)
           !enddo
           CALL INV(TPUP,TPUPM1,ZDUP1)
           CALL MMULT(TPUP1,URUP,UUP)
           CALL ZGETRF(N_size, N_size, TPUP1, N_size, IPVT, info)
           Z1 = ZDUP1
           Do i = 1, N_size
           IF (IPVT(i) .ne. i) THEN
            Z1 = -Z1
           endif
           Z1 = Z1 * TPUP1(I, I)
           enddo
        ELSE
           !WRITE(6,*) 'UDV of (U + DR * V * DL)^{*}'
           TPUP1 = CT(TPUP)
           CALL UDV_WRAP_Pivot(TPUP1,UUP,DUP,VUP,NCON,N_size,N_size)
           !CALL UDV(TPUP1,UUP,DUP,VUP,NCON)
           CALL ZGEMM('C', 'N', N_size, N_size, N_size, alpha, ULUP, N_size, UUP, N_size, beta, TPUPM1, N_size)
           CALL ZGEMM('N', 'C', N_size, N_size, N_size, alpha, URUP, N_size, VUP, N_size, beta, TPUP1, N_size)
           CALL ZGETRF(N_size, N_size, TPUP1, N_size, IPVT, info)
           ZDUP2 = 1.D0
           do i = 1, N_size
           ZDUP2 = ZDUP2 * TPUP1(I,I)
           IF (IPVT(i) .ne. i) THEN
            ZDUP2 = -ZDUP2
           endif
           enddo
           TPUP = TPUPM1
           ZDUP1 = DET_C(TPUP, N_size)! Det destroys its argument
           Z1 = ZDUP2/ZDUP1
        ENDIF
        DO J = 1, N_size
           sv = DBLE(DUP(J))
           X = ABS(sv)
           if (J == 1)  Xmax = X
           if ( X  < Xmax ) Xmax = X
           sv = 1.D0/sv
           DO I = 1, N_size
              UUP(J, I) = TPUPM1(I, J) * sv
           ENDDO
        ENDDO
        call ZGETRS('T', N_size, N_size, TPUP1, N_size, IPVT, UUP, N_size, info)
        GRUP = TRANSPOSE(UUP)
        PHASE = Z1/ABS(Z1)
        Deallocate(UUP, VUP, TPUP,TPUP1,TPUPM1, DUP, IPVT )

      END SUBROUTINE CGR
