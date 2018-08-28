      
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

       Program Cov_eq

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Analysis program for equal time observables.
!> 
!
!--------------------------------------------------------------------
         Use Errors
         Use MyMats
         Use Matrix
         Use Lattices_v3 

         Implicit none


         Integer      :: Nunit, Norb, ierr, qmin(6), qinst, nmin, avmin, avmax, avdiff, avmax2, avdiff2
         Integer      :: no, no1, n, n1,m,  nbins, nbins1, n_skip, nb, N_rebin, N_cov, N_Back
         real (Kind=Kind(0.d0)):: X, Y , qmin_norm
         real (Kind=Kind(0.d0)), allocatable :: Ratio1r(:), Ratio2r(:)
         Complex (Kind=Kind(0.d0)), allocatable :: Phase(:), Ratio1(:), Ratio2(:)
         Complex (Kind=Kind(0.d0)), allocatable :: U(:,:)
         Type  (Mat_C), allocatable :: Bins (:,:), Bins_R(:,:)
         Type  (Mat_C), allocatable :: Bins1 (:,:), Bins_R1(:,:)
         Type  (Mat_C), allocatable :: Bins1_tmp (:), Bins_R1_tmp(:)
         Real (Kind=Kind(0.d0)), allocatable :: EW(:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins0(:,:)
         Complex (Kind=Kind(0.d0)) :: Z, Z1, Z2, Xmean,Xerr, Xmean_r, Xerr_r
         Real (Kind=Kind(0.d0)) :: Xm,Xe
         Real    (Kind=Kind(0.d0)) :: Xk_p(2), XR_p(2) , XR1_p(2), XK_p_norm
         Complex (Kind=Kind(0.d0)), allocatable :: V_help(:), V_help_TR(:)
         Real (Kind=Kind(0.d0)) :: Pi, a1_p(2), a2_p(2), L1_p(2), L2_p(2), del_p(2)
         Real (Kind=Kind(0.d0)), allocatable :: AutoCorr(:),En(:)

         Integer             :: L1, L2, I, N_auto, N_SUN, nav=2
         Character (len=64)  :: Model, Lattice_type
         Type (Lattice)      :: Latt
         Character (len=64)  :: File_out
         Logical             :: Checkerboard	 

         NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model, N_SUN, Checkerboard
         NAMELIST /VAR_errors/   n_skip, N_rebin, N_Cov, N_Back, N_auto



         Checkerboard = .false.
         N_SUN  = 1
         N_Back = 1
         N_auto = 0
         OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'unable to open <parameters>',ierr
            STOP
         END IF
         READ(5,NML=VAR_lattice)
         READ(5,NML=VAR_errors)
         CLOSE(5)

         If ( Lattice_type =="BipartiteSquare" ) then
            a1_p(1) =  1.D0/sqrt(2.D0)  ; a1_p(2) =  1.D0/sqrt(2.D0)
            a2_p(1) =  1.D0/sqrt(2.D0)  ; a2_p(2) = -1.D0/sqrt(2.D0)
            L1_p    =  dble(L1)*a1_p
            L2_p    =  dble(L2)*a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
         elseif ( Lattice_type =="Square" ) then
            a1_p(1) =  1.0  ; a1_p(2) =  0.d0
            a2_p(1) =  0.0  ; a2_p(2) =  1.d0
            L1_p    =  dble(L1)*a1_p
            L2_p    =  dble(L2)*a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
         elseif ( Lattice_type=="Honeycomb" .or. Lattice_type=="Kagome" ) then
            a1_p(1) =  1.d0   ; a1_p(2) =  0.d0
            a2_p(1) =  0.5d0  ; a2_p(2) =  sqrt(3.d0)/2.d0
            del_p   =  (a2_p - 0.5*a1_p ) * 2.0/3.0
            L1_p    =  dble(L1) * a1_p
            L2_p    =  dble(L2) * a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
         elseif ( Lattice_type == "Pi_Flux" ) then 
             a1_p(1) =  1.D0   ; a1_p(2) =   1.d0
             a2_p(1) =  1.D0   ; a2_p(2) =  -1.d0
             !del_p   =  (a2_p - 0.5*a1_p ) * 2.0/3.0
             L1_p    =  dble(L1) * (a1_p - a2_p)/2.d0
             L2_p    =  dble(L2) * (a1_p + a2_p)/2.d0
             Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
         else
            Write(6,*) "Lattice not yet implemented!"
            Stop
         endif
         
         ! Determine the number of bins. 
         Pi = acos(-1.d0)
         Open ( Unit=10, File="ineq", status="unknown" ) 
         nbins = 0
         do
            Read(10,*,End=10) X,Norb,Nunit
            do n = 1,Norb
               Read(10,*) Z
            enddo
            do n = 1,Nunit
               Read(10,*) X,Y
               do no = 1,Norb
                  do no1 = 1,Norb
                     read(10,*) Z
                  enddo
               enddo
            enddo
            nbins = nbins + 1
         enddo
10       continue
         Close(10) 
         Write(6,*) "# of bins: ", Nbins
         nbins  = Nbins - n_skip
         Write(6,*) "Effective # of bins: ", Nbins
         N_auto=min(N_auto,Nbins-3)
         if(Nbins <= 1) then
           write (*,*) "Effective # of bins smaller then 2. Analysis impossible!"
           stop 1
         endif

         ! Allocate  space
         nbins1=nbins/N_rebin
         Allocate ( bins(Nunit,Nbins), bins_r(Nunit,Nbins), Phase(Nbins), Ratio1(Nbins), Ratio2(Nbins), V_help(Nbins), V_help_TR(Nbins), Bins0(Nbins,Norb))
         Do n = 1,Nunit
            do nb = 1,nbins
               Call Make_Mat(bins  (n,nb),Norb)
               Call Make_Mat(bins_r(n,nb),Norb)
               bins_r(n,nb)%el = cmplx(0.d0,0.d0,kind(0.d0))
               bins  (n,nb)%el = cmplx(0.d0,0.d0,kind(0.d0))
            Enddo
         Enddo
         Bins0 = cmplx(0.d0,0.d0,kind(0.d0))
         Open ( Unit=10, File="ineq", status="unknown" ) 
         do nb = 1, nbins + n_skip
            if (nb > n_skip ) then
               Z1=cmplx(0.d0,0.d0,kind(0.d0))
               Z2=cmplx(0.d0,0.d0,kind(0.d0))
               Read(10,*,End=10) X,no,no1
               Phase(nb-n_skip) = cmplx(X,0.d0,kind(0.d0))
               Do no = 1,Norb
                  Read(10,*) Z
                  if (N_Back == 1 ) Bins0(nb-n_skip,no) = Z
               enddo
               do n = 1,Nunit
                  Read(10,*) Xk_p(1), Xk_p(2)
                  m = Inv_K(Xk_p,Latt)
                  do no = 1,norb
                     do no1 = 1,Norb
                        read(10,*) bins(m,nb-n_skip)%el(no,no1) 
                     enddo
                  enddo
                  if ( sqrt(Xk_p(1)**2 + Xk_p(2)**2) < 1.D-6 ) then
                     do no = 1,norb
                        do no1 = 1,Norb
                           bins(m,nb-n_skip)%el(no,no1)  =  bins(m,nb-n_skip)%el(no,no1) -  &
                                &        cmplx(dble(Latt%N),0.d0,kind(0.d0))*Bins0(nb-n_skip,no)*Bins0(nb-n_skip,no1) &
                                &        /Phase(nb-n_skip)
                        enddo
                     enddo
                  endif
                  if (m==qinst) then
                    do no = 1,norb
                      do no1 = 1,Norb
                          Z1 = Z1 + bins(m,nb-n_skip)%el(no,no1) 
                      enddo
                    enddo
                  endif
                  do n1=1,nmin
                    if (m==qmin(n1)) then
                      do no = 1,norb
                        do no1 = 1,Norb
                            Z2 = Z2 + bins(m,nb-n_skip)%el(no,no1) 
                        enddo
                      enddo
                    endif
                  enddo
                  Ratio1(nb-n_skip) = Z2/dble(nmin)
                  Ratio2(nb-n_skip) = Z1 
               enddo
            else
               Read(10,*,End=10) X,no,no1
               Do no = 1,Norb
                  Read(10,*) Z
               enddo
               do n = 1,Nunit
                  Read(10,*) X,Y
                  do no = 1,Norb
                     do no1 = 1,Norb
                        read(10,*) Z
                     enddo
                  enddo
               enddo
            endif
         enddo
         close(10)
        N_auto=min(N_auto,Nbins-3)
        
        
         Allocate ( bins1(Nunit,Nbins1), bins_r1(Nunit,Nbins1), EW(Nunit,Nbins1,Norb), U(Norb,Norb))
         Allocate ( bins1_tmp(Nunit), Ratio1r(Nbins1), Ratio2r(Nbins1))
         Do n = 1,Nunit
            Call Make_Mat(bins1_tmp  (n),Norb)
            do nb = 1,nbins1
               Call Make_Mat(bins1  (n,nb),Norb)
               Call Make_Mat(bins_r1(n,nb),Norb)
               bins_r1(n,nb)%el = cmplx(0.d0,0.d0,kind(0.d0))
               bins1  (n,nb)%el = cmplx(0.d0,0.d0,kind(0.d0))
            Enddo
         Enddo
         Do n = 1,Nunit
            do nb = 1,nbins1*N_rebin
              bins1_tmp(n)%el = bins1_tmp(n)%el + bins(n,nb)%el
            Enddo
         Enddo
         Do n = 1,Nunit
            do nb = 1,nbins1
              bins1(n,nb)%el = bins1_tmp(n)%el
              do n1=(nb-1)*N_rebin+1,nb*N_rebin
                 bins1(n,nb)%el = bins1(n,nb)%el -bins(n,n1)%el
              enddo
              bins1(n,nb)%el = bins1(n,nb)%el/((Nbins1-1)*N_rebin)
              call Diag(bins1(n,nb)%el,U,EW(n,nb,:))
            Enddo
         Enddo

         Open (Unit=33,File="equalJ_Ratio"        ,status="unknown")
         Do n = 1,Nunit
            Xk_p = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p 
            Ratio1r=EW(n,:,norb)
            Ratio2r=0.d0
            n1=0
            if (sqrt(Latt%b1_p(1)**2.d0+Latt%b1_p(2)**2.d0) <= sqrt(Latt%b2_p(1)**2.d0+Latt%b2_p(2)**2.d0)*1.2d0 ) then
              n1=n1+1
              I=Latt%nnlistk(n,1,0)
              Ratio2r=Ratio2r+EW(I,:,norb)
              n1=n1+1
              I=Latt%nnlistk(n,-1,0)
              Ratio2r=Ratio2r+EW(I,:,norb)
            endif
            if (sqrt(Latt%b2_p(1)**2.d0+Latt%b2_p(2)**2.d0) <= sqrt(Latt%b1_p(1)**2.d0+Latt%b1_p(2)**2.d0)*1.2d0 ) then
              n1=n1+1
              I=Latt%nnlistk(n,0,1)
              Ratio2r=Ratio2r+EW(I,:,norb)
              n1=n1+1
              I=Latt%nnlistk(n,0,-1)
              Ratio2r=Ratio2r+EW(I,:,norb)
            endif
            if (L1 == L2 ) then
              n1=n1+1
              I=Latt%nnlistk(n,1,1)
              Ratio2r=Ratio2r+EW(I,:,norb)
              n1=n1+1
              I=Latt%nnlistk(n,-1,-1)
              Ratio2r=Ratio2r+EW(I,:,norb)
            endif
            Ratio1r=1.d0-Ratio2r/Ratio1r/n1
            call ERRCALC(Ratio1r,x,y)
            y = y*DBLE(size(Ratio1r))
            Write(33,"(F12.6,2x,F12.6,2x,F16.8,2x,F16.8)")  Xk_p(1), Xk_p(2), x, y
!             Write(*,*)  n, x, y
         Enddo
         Close(33)

      
       end Program Cov_eq
