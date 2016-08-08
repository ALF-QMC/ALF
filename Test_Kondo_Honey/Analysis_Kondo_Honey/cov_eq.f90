       Program Cov_eq

         Use Errors
         Use MyMats
         Use Matrix
         Use Lattices_v3 
         !   This version of the analysis program requires the information of the lattice, for fourier transforms 
         !   and for rotations. 

         Implicit none



         Interface
            Integer function Rot90(n, Xk_p, Ndim)
              Implicit none
              Integer, INTENT(IN)       :: Ndim,n
              Real (Kind=8), INTENT(IN) :: Xk_p(2,Ndim)
            end function Rot90
         end Interface

         Integer :: Ndim, Norb, nr, nx, ny,nk, ierr
         Integer :: no, no1, n, n1,m,  nbins, n_skip, nb, N_rebin
         real (Kind=8):: X, Y 
         Real (Kind=8), allocatable :: Phase(:)
         Type  (Mat_C), allocatable :: Bins(:,:), Bins_R(:,:)
         Complex (Kind=8), allocatable :: Bins0(:,:)
         Complex (Kind=8) :: Z, Xmean,Xerr, Xmean_r,Xerr_r
         Real    (Kind=8) :: Xk_p(2), XR_p(2) , XR1_p(2)
         Complex (Kind=8), allocatable :: V_help(:), V_help_R(:)
         Real (Kind=8) :: Pi, a1_p(2), a2_p(2), L1_p(2), L2_p(2), del_p(2)

         Integer             :: L1, L2, I
         Character (len=64)  :: Model, Lattice_type
         Type (Lattice)      :: Latt
         
         NAMELIST /VAR_lattice/  L1, L2, Lattice_type, Model


         OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read',IOSTAT=ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'unable to open <parameters>',ierr
            STOP
         END IF
         READ(5,NML=VAR_lattice)
         CLOSE(5)

         If ( Lattice_type =="Square" ) then
            a1_p(1) =  1.0  ; a1_p(2) =  0.d0
            a2_p(1) =  0.0  ; a2_p(2) =  1.d0
            L1_p    =  dble(L1)*a1_p
            L2_p    =  dble(L2)*a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
         elseif ( Lattice_type=="Honeycomb" ) then
            a1_p(1) =  1.0  ; a1_p(2) =  0.d0
            a2_p(1) =  0.5  ; a2_p(2) =  sqrt(3.0)/2.0
            del_p   =  (a2_p - 0.5*a1_p ) * 2.0/3.0
            L1_p    =  dble(L1) * a1_p
            L2_p    =  dble(L2) * a2_p
            Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
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
            Read(10,*,End=10) X,Norb,Ndim
            do n = 1,Norb
               Read(10,*) Z
            enddo
            do n = 1,Ndim
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
         n_skip = 1
         nbins  = Nbins - n_skip
         Write(6,*) "Effective # of bins: ", Nbins


         ! Allocate  space
         Allocate ( bins(Ndim,Nbins), bins_r(Ndim,Nbins), Phase(Nbins),  V_help(Nbins), V_help_R(Nbins), Bins0(Nbins,Norb))
         Do n = 1,Ndim
            do nb = 1,nbins
               Call Make_Mat(bins  (n,nb),Norb)
               Call Make_Mat(bins_r(n,nb),Norb)
               bins_r(n,nb)%el = 0.d0
            Enddo
         Enddo
         Open ( Unit=10, File="ineq", status="unknown" ) 
         do nb = 1, nbins + n_skip
            if (nb > n_skip ) then
               Read(10,*,End=10) Phase(nb-n_skip),no,no1
               Do no = 1,Norb
                  Read(10,*) Bins0(nb-n_skip,no)
               enddo
               do n = 1,Ndim
                  Read(10,*) Xk_p(1), Xk_p(2)
                  m = Inv_K(Xk_p,Latt)
                  !Write(6,*) m
                  do no = 1,norb
                     do no1 = 1,Norb
                        read(10,*) bins(m,nb-n_skip)%el(no,no1) 
                     enddo
                  enddo
                  if ( sqrt(Xk_p(1)**2 + Xk_p(2)**2) < 1.D-6 ) then
                     do no = 1,norb
                        do no1 = 1,Norb
                           bins(m,nb-n_skip)%el(no,no1)  =  bins(m,nb-n_skip)%el(no,no1) -  &
                                &        cmplx(dble(Latt%N),0.d0)*Bins0(nb-n_skip,no)*Bins0(nb-n_skip,no1)
                        enddo
                     enddo
                  endif
               enddo
            else
               Read(10,*,End=10) X,no,no1
               Do no = 1,Norb
                  Read(10,*) Z
               enddo
               do n = 1,Ndim
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

         ! Setup symmetries for C4v lattice
#ifdef test
         do n = 1,Ndim
            n1 = n
            Write(6,*) Xk_p(1,n1), Xk_p(2,n1)
            do m = 1,4
               n1 = Rot90(n1, Xk_p, Ndim)
               Write(6,*) n1, Xk_p(1,n1), Xk_p(2,n1)
            enddo
            Write(6,*)
         enddo
#endif
         Open (Unit=33,File="equalJ"       ,status="unknown")
         N_rebin = 1
         Do n1 = 1,Ndim
            n = n1
            do m = 1,1
               V_help   = 0.d0
               !n = Rot90(n, Xk_p, Ndim)
               do nb = 1,Nbins 
                  do no = 1,Norb
                     V_help  (nb) = V_help  (nb) + bins(n,nb)%el(no,no)
                  enddo
               enddo
               !V_help = V_help/dble(Norb)
               call ERRCALCJ(V_help,   XMean, XERR, N_rebin ) 
               Xk_p = dble(Latt%listk(n1,1))*Latt%b1_p + dble(Latt%listk(n1,2))*Latt%b2_p 
               Write(33,"(F12.6,2x,F12.6,2x,F12.6,2x,F12.6)") &
                    & Xk_p(1), Xk_p(2), dble(Xmean  ), dble(Xerr  )
            enddo
         enddo
         Close(33)


         Open (Unit=33,File="equalJ_Rc" , status="unknown")
         N_rebin = 1
         Do n1 = 1,Ndim
            n = n1
            do m = 1,1
               do no = 1,2
                  V_help_R = 0.d0
                  !n = Rot90(n, Xk_p, Ndim)
                  do nb = 1,Nbins 
                     V_help_r(nb) = V_help_r(nb) + bins_r(n,nb)%el(1,no)
                  enddo
                  !V_help = V_help/dble(Norb)
                  call ERRCALCJ(V_help_r, XMean_r, XERR_r, N_rebin ) 
                  Xr_p = dble(Latt%list (n1,1))*Latt%a1_p + dble(Latt%list (n1,2))*Latt%a2_p
                  if (no == 2) Xr_p = Xr_p - del_p 
                  Write(33,"(F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,F12.6)") &
                       &  Xr_p(1), Xr_p(2), dble(XMean_r), dble(Xerr_r)
               enddo
            enddo
         enddo
         Close(33)
         
         If (Norb > 2 ) then 
            Open (Unit=33,File="equalJ_Rf" , status="unknown")
            N_rebin = 1
            Do n1 = 1,Ndim
               n = n1
               do m = 1,1
                  do no = 3,4
                     V_help_R = 0.d0
                     !n = Rot90(n, Xk_p, Ndim)
                     do nb = 1,Nbins 
                        V_help_r(nb) = V_help_r(nb) + bins_r(n,nb)%el(3,no)
                     enddo
                     !V_help = V_help/dble(Norb)
                     call ERRCALCJ(V_help_r, XMean_r, XERR_r, N_rebin ) 
                     Xr_p = dble(Latt%list (n1,1))*Latt%a1_p + dble(Latt%list (n1,2))*Latt%a2_p
                     if (no == 4) Xr_p = Xr_p - del_p 
                     Write(33,"(F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,F12.6)") &
                          &  Xr_p(1), Xr_p(2), dble(XMean_r), dble(Xerr_r)
                  enddo
               enddo
            enddo
            Close(33)
         endif
         
#ifdef test
#endif
       end Program Cov_eq

       Integer function Rot90(n, Xk_p, Ndim)

         Implicit none
         Integer, INTENT(IN)       :: Ndim,n
         Real (Kind=8), INTENT(IN) :: Xk_p(2,Ndim)

         !Local
         real (Kind=8) :: X1_p(2), Zero, pi, X
         Integer :: m

         Zero = 1.D-4
         pi = acos(-1.d0)
         X1_p(1)  =  Xk_p(2,n)   
         X1_p(2)  = -Xk_p(1,n)   
         if (X1_p(1) < -pi + Zero )  X1_p(1) = X1_p(1) + 2.0*pi
         if (X1_p(2) < -pi + Zero )  X1_p(2) = X1_p(2) + 2.0*pi
         
         Rot90 = 0
         Do m = 1,Ndim
            X = sqrt( (X1_p(1) -Xk_p(1,m))**2 +  (X1_p(2) -Xk_p(2,m))**2 )
            If ( X < Zero) then
               Rot90 = m
               exit
            endif
         Enddo

       end function Rot90
