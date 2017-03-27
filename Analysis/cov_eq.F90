      
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


         Integer      :: Nunit, Norb, ierr
         Integer      :: no, no1, n, n1,m,  nbins, n_skip, nb, N_rebin, N_cov, N_Back
         real (Kind=Kind(0.d0)):: X, Y 
         Complex (Kind=Kind(0.d0)), allocatable :: Phase(:)
         Type  (Mat_C), allocatable :: Bins (:,:), Bins_R(:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins0(:,:)
         Complex (Kind=Kind(0.d0)) :: Z, Xmean,Xerr, Xmean_r, Xerr_r
         Real    (Kind=Kind(0.d0)) :: Xk_p(2), XR_p(2) , XR1_p(2)
         Complex (Kind=Kind(0.d0)), allocatable :: V_help(:), V_help_R(:)
         Real (Kind=Kind(0.d0)) :: Pi, a1_p(2), a2_p(2), L1_p(2), L2_p(2), del_p(2)

         Integer             :: L1, L2, I
         Character (len=64)  :: Model, Lattice_type
         Type (Lattice)      :: Latt
         

         NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model
         NAMELIST /VAR_errors/   n_skip, N_rebin, N_Cov, N_Back



         
         N_Back = 1
         OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'unable to open <parameters>',ierr
            STOP
         END IF
         READ(5,NML=VAR_lattice)
         READ(5,NML=VAR_errors)
         CLOSE(5)

         If ( Lattice_type =="Square" ) then
            a1_p(1) =  1.0  ; a1_p(2) =  0.d0
            a2_p(1) =  0.0  ; a2_p(2) =  1.d0
            L1_p    =  dble(L1)*a1_p
            L2_p    =  dble(L2)*a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
         elseif ( Lattice_type=="Honeycomb" ) then
            a1_p(1) =  1.d0   ; a1_p(2) =  0.d0
            a2_p(1) =  0.5d0  ; a2_p(2) =  sqrt(3.d0)/2.d0
            del_p   =  (a2_p - 0.5*a1_p ) * 2.0/3.0
            L1_p    =  dble(L1) * a1_p
            L2_p    =  dble(L2) * a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
            !  This will print the  honeycomb lattice. 
            Open (Unit=10,File="Lattice", status="unknown")
            do I = 1,Latt%N
               Xr_p = dble(Latt%list (I,1))*Latt%a1_p + dble(Latt%list (I,2))*Latt%a2_p
               Do n = 1,3
                  if (n==1) Xr1_p = Xr_p - del_p
                  if (n==2) Xr1_p = Xr_p - del_p - a1_p + a2_p 
                  if (n==3) Xr1_p = Xr_p + a2_p  - del_p
                  Write(10,"(F14.7,2x,F14.7)") Xr_p (1), Xr_p (2)
                  Write(10,"(F14.7,2x,F14.7)") Xr1_p(1), Xr1_p(2)
                  Write(10,*)
               enddo
            enddo
            close(10)
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


         ! Allocate  space
         Allocate ( bins(Nunit,Nbins), bins_r(Nunit,Nbins), Phase(Nbins),  V_help(Nbins), V_help_R(Nbins), Bins0(Nbins,Norb))
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
                                &        cmplx(dble(Latt%N),0.d0,kind(0.d0))*Bins0(nb-n_skip,no)*Bins0(nb-n_skip,no1)
                        enddo
                     enddo
                  endif
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
         

         Call Fourier_K_to_R(bins,bins_r,Latt)

         ! Setup symmetries for square lattice. 
#ifdef test
         do n = 1,Nunit
            n1 = n
            Write(6,*) Xk_p(1,n1), Xk_p(2,n1)
            do m = 1,4
               n1 = Rot90(n1, Xk_p, Nunit)
               Write(6,*) n1, Xk_p(1,n1), Xk_p(2,n1)
            enddo
            Write(6,*)
         enddo
#endif
         Open (Unit=33,File="equalJ"        ,status="unknown")
         Open (Unit=34,File="equalJR"       ,status="unknown")
         Do n = 1,Nunit
            Xk_p = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p 
            Xr_p = dble(Latt%list (n,1))*Latt%a1_p + dble(Latt%list (n,2))*Latt%a2_p 
            Write(33,"(F12.6,2x,F12.6)")  Xk_p(1), Xk_p(2)
            Write(34,"(F12.6,2x,F12.6)")  Xr_p(1), Xr_p(2)
            Do no = 1,Norb
               do no1 = 1,Norb
                  do nb = 1,Nbins
                     V_help(nb) = bins  (n,nb)%el(no,no1)
                  enddo
                  call ERRCALCJ( V_help, Phase,XMean, XERR, N_rebin ) 
                  Write(33,"(I3,2x,I3,2x,F16.8,2x,F16.8,2x,F16.8,2x,F16.8)") &
                       &  no,no1, dble(XMean), dble(XERR), aimag(XMean), aimag(XERR)
                  do nb = 1,Nbins
                     V_help(nb) = bins_r(n,nb)%el(no,no1)
                  enddo
                  call ERRCALCJ( V_help,Phase, XMean_r, XERR_r, N_rebin ) 
                  Write(34,"(I3,2x,I3,2x,F16.8,2x,F16.8,2x,F16.8,2x,F16.8)") &
                       &  no,no1, dble(XMean_r), dble(XERR_r), aimag(XMean_r), aimag(XERR_r)
               enddo
            enddo
         enddo
!!$         If (Norb > 1 ) then 
!!$            !Compute susecptibility 
!!$            Xk_p = 0.d0
!!$            n = Inv_K(Xk_p,Latt)
!!$            V_help   = 0.d0
!!$            do nb = 1,Nbins 
!!$               do no = 1,Norb
!!$                  Do no1 = 1,Norb
!!$                     V_help  (nb) = V_help  (nb) + bins(n,nb)%el(no,no1)
!!$                  enddo
!!$               enddo
!!$            enddo
!!$            call ERRCALCJ(V_help,   XMean, XERR, N_rebin ) 
!!$            Write(33,"('# Suscpetibility: ', F12.6,2x,F12.6)")  dble(Xmean  ), dble(Xerr  )
!!$         endif

         Close(33)
         Close(34)


      
       end Program Cov_eq


!!$         Interface
!!$            Integer function Rot90(n, Xk_p, Nunit)
!!$              Implicit none
!!$              Integer, INTENT(IN)       :: Nunit,n
!!$              Real (Kind=Kind(0.d0)), INTENT(IN) :: Xk_p(2,Nunit)
!!$            end function Rot90
!!$         end Interface
!!$       Integer function Rot90(n, Xk_p, Nunit)
!!$
!!$         Implicit none
!!$         Integer, INTENT(IN)       :: Nunit,n
!!$         Real (Kind=Kind(0.d0)), INTENT(IN) :: Xk_p(2,Nunit)
!!$
!!$         !Local
!!$         real (Kind=Kind(0.d0)) :: X1_p(2), Zero, pi, X
!!$         Integer :: m
!!$
!!$         Zero = 1.D-4
!!$         pi = acos(-1.d0)
!!$         X1_p(1)  =  Xk_p(2,n)   
!!$         X1_p(2)  = -Xk_p(1,n)   
!!$         if (X1_p(1) < -pi + Zero )  X1_p(1) = X1_p(1) + 2.0*pi
!!$         if (X1_p(2) < -pi + Zero )  X1_p(2) = X1_p(2) + 2.0*pi
!!$         
!!$         Rot90 = 0
!!$         Do m = 1,Nunit
!!$            X = sqrt( (X1_p(1) -Xk_p(1,m))**2 +  (X1_p(2) -Xk_p(2,m))**2 )
!!$            If ( X < Zero) then
!!$               Rot90 = m
!!$               exit
!!$            endif
!!$         Enddo
!!$
!!$       end function Rot90
