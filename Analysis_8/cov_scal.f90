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

      Program Cov_vec

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Analysis program for scalar observables
!
!--------------------------------------------------------------------

        
        Use ERRORS
        Implicit none

        REAL    (Kind=Kind(0.d0)), DIMENSION(:,:), ALLOCATABLE :: OBS
        REAL    (Kind=Kind(0.d0)), DIMENSION(:),   ALLOCATABLE :: EN, SIGN
        REAL    (Kind=Kind(0.d0)) :: XM, XERR, X

        Complex (Kind=Kind(0.d0)) Z1,Z2,Z3,Z4,Z5
        Complex (Kind=Kind(0.d0)), Allocatable  :: Tmp(:)
        Integer :: Nobs 
        Integer :: NST, NS, NS1, NS2, NSTEP, NC, NP,  Nbins, NP_EFF, ISEED, I, IOBS, N_Back
        Integer :: N,N1, NBIN

        Integer :: n_skip, N_rebin, N_Cov, ierr
        NAMELIST /VAR_errors/   n_skip, N_rebin, N_Cov, N_Back


         OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'unable to open <parameters>',ierr
            STOP
         END IF
         READ(5,NML=VAR_errors)
         CLOSE(5)


        ! Count the number of bins
        Open (Unit=10, File="Var_scal", status="unknown")
        Read(10,*) NOBS
        allocate (Tmp(NOBS) )
        rewind(10)
        nbins = 0
        do
            read(10,*,End=10) N,  (Tmp(I), I=1,size(Tmp,1)-1), X
            Tmp(NOBS) = cmplx(X,0.d0,kind(0.d0))
            nbins = nbins + 1
         enddo
10       continue
         Write(6,*) "# of bins: ", Nbins
         Close(10) 
         
         NP = NBINS
         
         ALLOCATE(OBS(NP,NOBS))

         NST = N_skip; NS1 = N_Rebin; NS2 = N_Rebin; NSTEP = 1
         OPEN (UNIT=20, FILE='Var_scal', STATUS='old')
         NC = 0
         DO N = 1,NP
            IF (N.GE.NST) THEN
               NC = NC + 1
               READ(20,*) N1, (Tmp(I), I=1,size(Tmp,1)-1),X 
               Tmp(NOBS) = cmplx(X,0.d0,kind(0.d0))
               OBS(NC,:) = dble(Tmp(:)) 
            ELSE
               READ(20,*) N1, (Tmp(I), I=1,size(Tmp,1)-1),X 
            ENDIF
         ENDDO
         CLOSE(20)
2100     FORMAT(I6,2X,F16.8)
         
         OPEN (UNIT=21, FILE='Var_scalJ', STATUS='unknown')
         WRITE(21,*) 'Effective number of bins, and bins: ', NC, NP
         NP_EFF = NC
         ALLOCATE (EN(NP_EFF), SIGN(NP_EFF))
         DO IOBS = 1,NOBS
            WRITE(21,*)
            DO I = 1,NP_EFF
               EN  (I) = Real(OBS(I,IOBS), kind(0.d0))
               SIGN(I) = Real(OBS(I,NOBS), kind(0.d0))
            ENDDO
            DO NBIN = NS1, NS2, NSTEP
               if (NBIN.gt.0) then
                  IF (IOBS.EQ.NOBS  ) then 
                     CALL ERRCALCJ(EN,XM,XERR,NBIN)
                  else
                     CALL ERRCALCJ(EN,SIGN,XM,XERR,NBIN)
                  endif
                  WRITE(21,2001) IOBS, XM,  XERR
                  ! NBOOT = 40
                  ! CALL BOOTSTRAP( EN,XM_BS,XERR_BS,NBOOT,ISEED)
                  ! WRITE(21,2001) IOBS, XM_BS,  XERR_BS
                  ! IF (IOBS == 4) Write(22,"(F14.7,2x,F14.7)")  XM/dble(L*L), XERR/dble(L*L)
               endif
            ENDDO
         ENDDO
         CLOSE(21)
2001     FORMAT('OBS : ', I4,4x,F12.6,2X, F12.6)
         
         DEALLOCATE (EN,SIGN,OBS)
         
       END Program Cov_vec
       
