       Program Cov_tau

         Use Errors
         Use MyMats
         Use Matrix
         Use Precdef

         Implicit none

         Interface
            Integer function Rot90(n, Xk_p, Ndim)
              Implicit none
              Integer, INTENT(IN)       :: Ndim,n
              Real (Kind=8), INTENT(IN) :: Xk_p(2,Ndim)
            end function Rot90
         end Interface

         Integer :: Ndim, Norb
         Integer :: no, no1, n, nbins, n_skip, nb, nT, Lt
         real (Kind=8):: X, Y,  dtau, X_diag
         real (Kind=8), allocatable :: Xmean(:), Xcov(:,:)
         Complex (Kind=8) :: Z
         Real (Kind=8) :: Zero=1.D-8
         Real (Kind=8), allocatable :: Phase(:)
         Real    (Kind=8), allocatable :: Bins(:,:,:)
         Complex (Kind=8), allocatable :: Bins0(:,:)
         Real (Kind=8), allocatable :: Xk_p(:,:)
         Real (Kind=8), allocatable :: V_help(:,:)
         Character (len=64) :: File_out


         ! Determine the number of bins. 
         Open ( Unit=10, File="intau", status="unknown" ) 
         nbins = 0
         do
            Read(10,*,End=10) X,Norb,Ndim, LT, dtau
            Do no = 1,Norb
               Read(10,*) Z
            enddo
            do n = 1,Ndim
               Read(10,*) X,Y
               do nt = 1,LT
                  do no = 1,Norb
                     do no1 = 1,Norb
                        read(10,*) Z
                     enddo
                  enddo
               enddo
            enddo
            Write(6,*) nbins
            nbins = nbins + 1
         enddo
10       continue
         Close(10) 
         Write(6,*) "# of bins: ", Nbins
         n_skip = 1
         nbins  = Nbins - n_skip
         Write(6,*) "Effective # of bins: ", Nbins

         


         ! Allocate  space
         Allocate ( bins(Ndim,Lt,Nbins), Phase(Nbins), Xk_p(2,ndim), V_help(lt,Nbins), bins0(Nbins,Norb))
         Allocate (Xmean(Lt), Xcov(Lt,Lt))
         bins  = 0.d0
         bins0 = cmplx(0.d0,0.d0,Kind(0.d0))
         Open ( Unit=10, File="intau", status="unknown" ) 
         do nb = 1, nbins + n_skip
            if (nb > n_skip ) then
               Read(10,*,End=10) Phase(nb-n_skip),no,no1,n, X
               X_diag = 0.d0
               Do no = 1,Norb
                  Read(10,*)  bins0(nb-n_skip,no)
                  X_diag =  X_diag + real(bins0(nb-n_skip,no)*bins0(nb-n_skip,no),kind=8)
               Enddo
               do n = 1,Ndim
                  Read(10,*) Xk_p(1,n), Xk_p(2,n)
                  do nt = 1,Lt
                     do no = 1,norb
                        do no1 = 1,Norb
                           read(10,*) Z
                           if (no == no1) bins(n,nt,nb-n_skip) = bins(n,nt,nb-n_skip) +  real(Z,Kind=8) 
                        enddo
                     enddo
                  enddo
                  if ( sqrt(Xk_p(1,n)**2 + Xk_p(2,n)**2) < 1.D-6 ) then
                     Do nt = 1,Lt
                        bins(n,nt,nb-n_skip) =   bins(n,nt,nb-n_skip) - X_diag
                     enddo
                  endif
               enddo
            else
               Read(10,*,End=10) X,no,no1,n,Y
               do no = 1,Norb
                  Read(10,*) Z
               enddo
               do n = 1,Ndim
                  Read(10,*) X,Y
                  do nt = 1,LT
                     do no = 1,Norb
                        do no1 = 1,Norb
                           read(10,*) Z
                        enddo
                     enddo
                  enddo
               enddo
            endif
         enddo
         close(10)

         do n = 1,Nbins
            Write(6,*) Phase(n)
         Enddo
         do n = 1,Ndim
            V_help = 0.d0
            !n1 = n
            if (  Xk_p(1,n) >= -zero .and. XK_p(2,n) >= -zero ) then
               !do m = 1,4
               !   n1 = Rot90(n1, Xk_p, Ndim)
               !   do nt = 1,LT
               !      do nb = 1,nbins
               !         V_help(nt,nb) = V_help(nt,nb) + bins (n1,nt,nb)
               !      enddo
               !   enddo
               !enddo
               !V_help = V_help/4.d0
               call COV(bins(n,:,:), phase, Xcov, Xmean )
               write(File_out,'("g_",F4.2,"_",F4.2)')  Xk_p(1,n), Xk_p(2,n)
               Open (Unit=10,File=File_out,status="unknown")
               do nt = 1, LT
                  Write(10,"(F14.7,2x,F16.8,2x,F16.8)") &
                       & dble(nt-1)*dtau,  Xmean(nt), sqrt(abs(dble(Xcov(nt,nt))))
               enddo
               close(10)
            endif
         enddo

         V_help = 0.d0
         do n = 1,Ndim
            do nb = 1,nbins
               do nt = 1,LT
                  V_help(nt,nb) = V_help(nt,nb) +  bins(n,nt,nb)
               enddo
            enddo
         enddo
         V_help = V_help/dble(Ndim)
         call COV(V_help, phase, Xcov, Xmean )
         write(File_out,'("g_R0")') 
         Open (Unit=10,File=File_out,status="unknown")
         do nt = 1, LT
            Write(10,"(F14.7,2x,F16.8,2x,F16.8)") &
                 & dble(nt-1)*dtau,  Xmean(nt), sqrt(abs(dble(Xcov(nt,nt))))
         enddo
         close(10)
         

       end Program Cov_tau

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
