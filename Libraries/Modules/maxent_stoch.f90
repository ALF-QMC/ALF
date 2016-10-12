
Module MaxEnt_stoch_mod

       Use MyMats
       Use Random_Wrap
       Use Files_mod
       
       Integer, private :: NTAU, nt, Ngamma, ng, Ndis, nd,  L_seed
       Integer, private, allocatable:: Iseed_vec(:)
       Real (Kind=8), private :: Delta, Delta2, OM_st_1, Om_en_1, DeltaXMAX, Beta, Pi, Dom
       Real (Kind=8), allocatable, private :: XQMC1(:)
       
       ! You can still optimize a bit for by redefining the Kernel table to:
       ! xker_table(nt,nw) -> xker_table(nt,nw) / sigma(nt)
       ! This will save quite a lot of divisions in the
       ! MC routine. And this is where all the time goes now.
       
       CONTAINS
         
         Subroutine MaxEnt_stoch(XQMC, Xtau, COV,Xmom1, XKER, Back_Trans_Aom, Beta_1, Alpha_tot,&
              & Ngamma_1, OM_ST, OM_EN, Ndis_1, Nsweeps, NBins, NWarm, L_cov)

           Implicit None

           Real (Kind=8), Dimension(:) :: XQMC, Xtau, Alpha_tot
           Real (Kind=8), Dimension(:,:) :: COV
           Real (Kind=8), External :: XKER, Back_trans_Aom
           Real (Kind=8) :: CHISQ, OM_ST, OM_EN, Beta_1, Xmom1, Err
           Integer :: Nsweeps, NBins, Ngamma_1, Ndis_1, nw, nt1
           Integer, optional :: L_cov

           ! Local
           Integer NSims, ns, nb, nc, Nwarm, nalp1, nalp2, Nex, p_star, Ndis_table, &
                & io_error, io_error1, i, n, nc1
           Real (Kind=8), Allocatable :: Xn_M_tot(:,:), En_M_tot(:), Xn_E_tot(:,:), En_E_tot(:), &
                & Xn_tot(:,:,:), En_tot(:)
           Real (Kind=8), Allocatable :: G_Mean(:), Xn_m(:), Xn_e(:), Xn(:,:), Vhelp(:)
           Real (Kind=8) :: En_M, X, Alpha, Acc_1, Acc_2, En, DeltaE, Ratio, D
           Real (Kind=8) :: Aom, om, XMAX, tau
           Real (Kind=8) :: CPUT, CPUTM
           Integer :: ICPU_1, ICPU_2, N_P_SEC
           Character (64) :: File_root, File1, File_conf, File_Aom
           Real (Kind=8), allocatable :: Xker_table(:,:), U(:,:), sigma(:)
           ! Space for moments.
           Real (Kind=8), allocatable:: Mom_M_tot(:,:), Mom_E_tot(:,:)
           
           Pi = acos(-1.d0)
           NDis = Ndis_1
           DeltaXMAX = 0.01
           delta = 0.001
           delta2 = delta*delta
           Ngamma = Ngamma_1
           Beta = Beta_1 ! Physical temperature for calculation of the kernel.
           Ntau = Size(xqmc,1)
           NSims = Size(Alpha_tot,1)
           Allocate (Xn_tot(Ngamma,2,NSims))
           Allocate (En_m_tot(NSims), En_e_tot(NSims), En_tot(NSims) )
           Allocate (Xn_m_tot(NDis,NSims), Xn_e_tot(NDis,NSims) )
           Allocate (Mom_M_tot(2,Nsims), Mom_E_tot(2,Nsims) )
           Allocate (Xn(Ngamma,2))
           Allocate (Xn_m(NDis), Xn_e(NDis) )
           Om_st_1 = OM_st; Om_en_1 = OM_en
           ! Setup table for the Kernel
           Ndis_table = 50000
           Dom = (OM_EN_1 - OM_ST_1)/dble(Ndis_table-1)
           Allocate ( Xker_table(Ntau, Ndis_table) )
           do nt = 1,Ntau
              do nw = 1,Ndis_table
                 tau = xtau(nt)
                 Om = OM_st + dble(nw-1)*dom
                 Xker_table(nt,nw) = Xker(tau,om,beta)
              enddo
           enddo
           ! Normalize data to have zeroth moment of unity.
           xqmc = xqmc / XMOM1
           cov  = cov / ((XMOM1)**2)
           ! Diagonalize the covariance
           Allocate( U(ntau,ntau), Sigma(ntau), xqmc1(Ntau) )
           If ( Present(L_cov) ) then
              Call Diag(cov,U,sigma)
              ! Write(6,*) " Cov Used"
           else
              Write(6,*) "No Cov Used"
              U = 0.d0
              do nt = 1,ntau
                 U(nt,nt) = 1.d0
                 sigma(nt) = cov(nt,nt)
              enddo
           endif
           do nt = 1,ntau
              sigma(nt) = sqrt(sigma(nt))
           enddo
           xqmc1 = 0.d0
           do nt1 = 1,ntau
              do nt = 1,ntau
                 xqmc1(nt1) = xqmc1(nt1) + xqmc(nt)*U(nt,nt1)
              enddo
              xqmc1(nt1) = xqmc1(nt1)/sigma(nt1)
           enddo
           ! Transform the Kernel
           allocate ( Vhelp(Ntau) )
           do nw = 1,Ndis_table
              Vhelp = 0.d0
              do nt1 = 1,Ntau
                 Vhelp(nt1) = Vhelp(nt1) + dot_product(Xker_table(:, nw), U(:,nt1))
              enddo
              Xker_table(:, nw) = Vhelp/sigma
           enddo
           deallocate( U, Sigma )
           Allocate ( G_Mean(Ntau) )
           G_mean = 0.d0
           ! write(6,*) ' There are ', Ngamma,' delta-functions for a spectrum'
           ! Write(6,*) ' Initializing'
           ! D(om) = 1/(Om_en_1 - Om_st_1)
           D = 1.d0 / (Om_en_1 - Om_st_1)
           
           Call Get_seed_Len(L_seed)
           Allocate(Iseed_vec(L_seed))
           Iseed_vec = 0
           Call Ranset(Iseed_vec)
           File_conf = "dump_conf"
           File_Aom = "dump_Aom"
           Open(unit=41,file=File_conf,status='old',action='read', iostat=io_error)
           Open(unit=42,file=File_Aom, status='old',action='read', iostat=io_error1)
           If (io_error == 0 .and. io_error1 == 0 ) then
              Nwarm = 0
              read(41,*) (Iseed_vec(I), I = 1,size(Iseed_vec,1))
              do ns = 1,Nsims
                 do ng = 1,Ngamma
                    read(41,*) Xn_tot(ng,1,ns), Xn_tot(ng,2,ns)
                 enddo
                 read(41,*) En_m_tot(ns), En_e_tot(ns)
              enddo
              read(42,*) nc
              do ns = 1,Nsims
                 do nd = 1,Ndis
                    read(42,*) Xn_m_tot(nd,ns), Xn_e_tot(nd,ns)
                 enddo
              enddo
              Open (Unit=44,File='Max_stoch_log', Status="unknown", position="append")
              Write(44,*) 'Read from dump: nc = ', nc
              close(44)
           else
              !Iseed is already set.
              Do Ns = 1,NSims
                 do ng = 1,NGamma 
                    Xn_tot(ng,1,ns) = ranf_wrap()
                    Xn_tot(ng,2,ns) = 1.d0/dble(Ngamma)
                 enddo
              enddo
              Xn_m_tot = 0.d0
              En_m_tot = 0.d0
              Xn_e_tot = 0.d0
              En_e_tot = 0.d0
              nc = 0
              Open (Unit=44,File='Max_stoch_log', Status="unknown", position="append")
              Write(44,*) ' No dump data '
              close(44)
           endif
           close(41)
           close(42)
           nc1 = 0
           Mom_M_tot = 0.d0
           Mom_E_tot = 0.d0
           CALL SYSTEM_CLOCK(COUNT_RATE=N_P_SEC)
           CALL SYSTEM_CLOCK(COUNT=ICPU_1)
           ! Start Simulations.
           do Nb = 1,Nbins
              do ns = 1,NSims
                 do ng = 1,Ngamma
                    Xn(ng,1) = Xn_tot(ng,1,ns)
                    Xn(ng,2) = Xn_tot(ng,2,ns)
                 enddo
                 Alpha = Alpha_tot(ns)
                 Call MC(Xtau, Xker_table, Xn, Alpha, NSweeps, Xn_m, En, En_m, &
                      & Acc_1, Acc_2 ) ! Just one bin
                 do ng = 1,Ngamma
                    Xn_tot(ng,1,ns) = Xn(ng,1)
                    Xn_tot(ng,2,ns) = Xn(ng,2)
                 enddo
                 En_tot(ns) = En ! this is the energy of the configuration Xn_tot for simulation ns
                 Open (Unit=44,File='Max_stoch_log', Status="unknown", position="append")
                 Write(44,2003) 1.d0/Alpha, En_m, Acc_1, Acc_2
                 close(44)
                 if (nb.gt.nwarm) then
                    if (ns.eq.1) nc = nc + 1
                    do nd = 1,NDis
                       Xn_m(nd) = Xn_m(nd) * D * dble(Ndis)
                       Xn_m_tot(nd,ns) = Xn_m_tot(nd,ns) + Xn_m(nd)
                       Xn_e_tot(nd,ns) = Xn_e_tot(nd,ns) + Xn_m(nd)*Xn_m(nd)
                    enddo
                    En_m_tot(ns) = En_m_tot(ns) + En_m
                    En_e_tot(ns) = En_e_tot(ns) + En_m*En_m
                    ! Compute moments
                    if (ns.eq.1) nc1 = nc1 + 1
                    do n = 1,Size(Mom_M_tot,1)
                       x = 0.d0
                       do ng = 1,Ngamma
                          X = X + ( Phim1(Xn_tot(ng,1,ns))**(n-1) ) * Xn_tot(ng,2,Ns)
                       enddo
                       Mom_M_tot(n,ns) = Mom_M_tot(n,ns) + X
                       Mom_E_tot(n,ns) = Mom_E_tot(n,ns) + X*X
                    enddo
                 endif
              enddo
              ! Exchange
              Acc_1 = 0.d0
              Do Nex = 1, 2*NSims
                 nalp1= nint( ranf_wrap()*dble(NSims-1) + 0.5 ) ! 1..(NSims-1)
                 nalp2 = nalp1 + 1
                 DeltaE = (Alpha_tot(nalp1)*En_tot(nalp2) + Alpha_tot(nalp2)*En_tot(nalp1))&
                      & -(Alpha_tot(nalp1)*En_tot(nalp1) + Alpha_tot(nalp2)*En_tot(nalp2))
                 Ratio = exp(-DeltaE)
                 if (Ratio.gt.ranf_wrap()) Then
                    Acc_1 = Acc_1 + 1.0
                    !Switch confs an Energies.
                    do ng = 1,Ngamma
                       Xn(ng,1) = Xn_tot(ng,1,nalp1)
                       Xn(ng,2) = Xn_tot(ng,2,nalp1)
                    enddo
                    do ng = 1,Ngamma
                       Xn_tot(ng,1,nalp1) = Xn_tot(ng,1,nalp2)
                       Xn_tot(ng,2,nalp1) = Xn_tot(ng,2,nalp2)
                       Xn_tot(ng,1,nalp2) = Xn(ng,1)
                       Xn_tot(ng,2,nalp2) = Xn(ng,2)
                    enddo
                    En_m = En_tot(nalp1)
                    En_tot(nalp1) = En_tot(nalp2)
                    En_tot(nalp2) = En_m
                 endif
              enddo
              Acc_1 = Acc_1/dble(Nex)
              Open (Unit=44,File='Max_stoch_log', Status="unknown", position="append")
              Write(44,*) 'Acc Exchange: ', Acc_1
              close(44)
           enddo
           CALL SYSTEM_CLOCK(COUNT=ICPU_2)
           CPUT = 0.D0
           CPUT = DBLE(ICPU_2 - ICPU_1)/DBLE(N_P_SEC)
           Open (Unit=44,File='Max_stoch_log', Status="unknown", position="append")
           Write(44,*) 'Total time: ', CPUT
           close(44)
           ! dump so as to restart.
           Open(unit=41,file=File_conf,status='unknown')
           Open(unit=42,file=File_Aom, status='unknown')
           Call Ranget(Iseed_vec)
           write(41,*) (Iseed_vec(I), I = 1,size(Iseed_vec,1))
           do ns = 1,Nsims
              do ng = 1,Ngamma
                 write(41,*) Xn_tot(ng,1,ns), Xn_tot(ng,2,ns)
              enddo
              write(41,*) En_m_tot(ns), En_e_tot(ns)
           enddo
           write(42,*) nc
           do ns = 1,Nsims
              do nd = 1,Ndis
                 write(42,*) Xn_m_tot(nd,ns), Xn_e_tot(nd,ns)
              enddo
           enddo
           close(41)
           close(42)
           ! Stop dump
           Open(Unit=66,File="energies",status="unknown")
           do ns = 1,Nsims
              En_m_tot(ns) = En_m_tot(ns) / dble(nc)
              En_e_tot(ns) = En_e_tot(ns) / dble(nc)
              En_e_tot(ns) = ( En_e_tot(ns) - En_m_tot(ns)**2)/dble(nc)
              if ( En_e_tot(ns) .gt. 0.d0) then
                 En_e_tot(ns) = sqrt(En_e_tot(ns))
              else
                 En_e_tot(ns) = 0.d0
              endif
              write(66,*) Alpha_tot(ns), En_m_tot(ns), En_e_tot(ns)
           enddo
           close(66)
           Open(Unit=66,File="moments",status="unknown")
           do ns = 1,Nsims
              do n = 1,Size(Mom_m_tot,1)
                 Mom_m_tot(n,ns) = Mom_m_tot(n,ns) / dble(nc1)
                 Mom_e_tot(n,ns) = Mom_e_tot(n,ns) / dble(nc1)
                 Mom_e_tot(n,ns) = ( Mom_e_tot(n,ns) - Mom_m_tot(n,ns)**2)/dble(nc1)
                 if ( Mom_e_tot(n,ns) .gt. 0.d0) then
                    Mom_e_tot(n,ns) = sqrt(Mom_e_tot(n,ns))
                 else
                    Mom_e_tot(n,ns) = 0.d0
                 endif
              enddo
              write(66,"(F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,F12.6)") &
                   & Alpha_tot(ns), Mom_m_tot(1,ns), Mom_e_tot(1,ns), &
                   & Mom_m_tot(2,ns), Mom_e_tot(2,ns)
           enddo
           close(66)
           File_root = "Aom"
           do ns = 1,Nsims
              File1 = File_i(File_root,ns)
              Open(Unit=66,File=File1,status="unknown")
              do nd = 1,Ndis
                 Xn_m_tot(nd,ns) = Xn_m_tot(nd,ns) / dble(nc) ! * delta /(dble(nc)*pi)
                 Xn_e_tot(nd,ns) = Xn_e_tot(nd,ns) / dble(nc) ! * delta /(dble(nc)*pi)
                 Xn_e_tot(nd,ns) = (Xn_e_tot(nd,ns) - Xn_m_tot(nd,ns)* Xn_m_tot(nd,ns))/dble(nc)
                 if (Xn_e_tot(nd,ns).gt.0.d0) then
                    Xn_e_tot(nd,ns) = sqrt(Xn_e_tot(nd,ns))
                 else
                    Xn_e_tot(nd,ns) = 0.d0
                 endif
                 om = PhiM1(dble(nd)/dble(NDis))
                 Aom = Xn_m_tot(nd,ns) * Xmom1
                 Err = Xn_e_tot(nd,ns) * Xmom1
                 write(66,2001) om, Back_Trans_Aom(Aom,Beta,om), Back_Trans_Aom(Err,Beta,om)
                 ! PhiM1(dble(nd)/dble(NDis)), Xn_m_tot(nd,ns)
              enddo
              Close(66)
           enddo
           ! Now do the averaging.
           File_root ="Aom_ps"
           do p_star = 1,NSims - 10
              Xn_m = 0.0
              Xn_e = 0.0
              do ns = p_star, NSims-1
                 do nd = 1, NDis
                    Xn_m(nd) = Xn_m(nd) + (En_m_tot(ns) - En_m_tot(ns+1))*Xn_m_tot(nd,ns)
                    Xn_e(nd) = Xn_e(nd) + (En_m_tot(ns) - En_m_tot(ns+1))*Xn_e_tot(nd,ns)
                 enddo
              enddo
              do nd = 1,NDis
                 Xn_m(nd) = Xn_m(nd) / (En_m_tot(p_star) - En_m_tot(NSims))
                 Xn_e(nd) = Xn_e(nd) / (En_m_tot(p_star) - En_m_tot(NSims))
              enddo
              File1 = File_i(File_root,p_star)
              Open(Unit=66,File=File1,status="unknown")
              XMAX = 0.d0
              Do nd = 1,Ndis
                 om = PhiM1(dble(nd)/dble(NDis))
                 Aom = Xn_m(nd) * Xmom1
                 Err = Xn_e(nd) * Xmom1
                 Xn_m(nd) = Back_Trans_Aom(Aom,Beta,om)
                 Xn_e(nd) = Back_Trans_Aom(Err,Beta,om)
                 IF (Xn_m(nd) .gt. XMAX ) XMAX = Xn_m(nd)
              enddo
              do nd = 1,Ndis
                 om = PhiM1(dble(nd)/dble(NDis))
                 write(66,2005) om, Xn_m(nd), Xn_e(nd), Xn_m(nd)/XMAX, Xn_e(nd)/XMAX
                 ! PhiM1(dble(nd)/dble(NDis)), Xn_m(nd)
              enddo
              close(66)
           enddo
           Open (Unit=41,File='Best_fit', Status="unknown")
           do ng = 1,Ngamma
              Write(41,*) Phim1(Xn_tot(ng,1,Nsims)) , Xn_tot(ng,2,Nsims)
           enddo
           close(41)
           Deallocate (Iseed_vec)
           DeAllocate (Xn_tot)
           DeAllocate (En_m_tot, En_e_tot, En_tot )
           DeAllocate (Xn_m_tot, Xn_e_tot )
           DeAllocate (Xn)
           DeAllocate (Xn_m, Xn_e)
           DeAllocate( G_Mean )
           DeAllocate( xqmc1 )
           Deallocate( Xker_table )
2001       format(F14.7,2x,F14.7,2x,F14.7)
2004       format(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)
2005       format(F14.7,2x,F14.7,2x,F14.7,2x,F14.7,2x,F14.7)
2003       format('Alpha, En_m, Acc ', F14.7,2x,F24.12,2x,F14.7,2x,F14.7,2x,F14.7)
           
         end Subroutine MaxEnt_stoch
!-------------------
            
         Subroutine MaxEnt_stoch_fit(XQMC, Xtau, COV, Lcov, XKER, Xmom1, Beta_1, Alpha_tot,&
                 & Ngamma_1, OM_ST, OM_EN, Nsweeps, NBins, NWarm, Aom_res,&
                 & xom_res, Chisq )
              
           Implicit None
           
           Real (Kind=8), Dimension(:) :: XQMC, Xtau, Alpha_tot, Aom_res, Xom_res
           Real (Kind=8), Dimension(:,:) :: COV
           Real (Kind=8), external :: XKER
           Real (Kind=8) :: CHISQ, OM_ST, OM_EN, Beta_1, Xmom1, Err
           Integer :: Nsweeps, NBins, Ngamma_1, nw, nt1, Lcov
           ! Local
           Integer NSims, ns, nb, nc, Nwarm, nalp1, nalp2, Nex
           Real (Kind=8), Allocatable :: Xn_M_tot(:,:), En_M_tot(:), Xn_E_tot(:,:), En_E_tot(:), &
                & Xn_tot(:,:,:), En_tot(:), Xker_table(:,:)
           Real (Kind=8), Allocatable :: G_Mean(:), Xn_m(:), Xn_e(:), Xn(:,:), Vhelp(:)
           Real (Kind=8) :: En_M, Alpha, Acc_1, Acc_2, En, DeltaE, Ratio, D
           Real (Kind=8) :: Aom, om, tau
           Real (Kind=8), allocatable :: U(:,:), sigma(:)
           Pi = acos(-1.d0)
           Call Get_seed_Len(L_seed)
           Allocate(Iseed_vec(L_seed))
           Iseed_vec = 0
           Call Ranset(Iseed_vec)
           NDis = Size(Aom_res,1)
           DeltaXMAX = 0.01
           delta = 0.001
           delta2 = delta*delta
           Ngamma = Ngamma_1
           Beta = Beta_1 ! Physical temperature for calculation of the kernel.
           Ntau = Size(xqmc,1)
           NSims = Size(Alpha_tot,1)
           Allocate (Xn_tot(Ngamma,2,NSims))
           Allocate (En_m_tot(NSims), En_e_tot(NSims), En_tot(NSims) )
           Allocate (Xn_m_tot(NDis,NSims), Xn_e_tot(NDis,NSims) )
           Allocate (Xn(Ngamma,2))
           Allocate (Xn_m(NDis), Xn_e(NDis) )
           Om_st_1 = OM_st; Om_en_1 = OM_en
           ! Setup table for the Kernel
           Ndis = Size(Aom_res)
           Dom = (OM_EN_1 - OM_ST_1)/dble(Ndis-1)
           Allocate ( Xker_table(Ntau, Ndis) )
           Dom = (OM_EN - OM_ST)/dble(Ndis-1)
           do nt = 1,Ntau
              do nw = 1,Ndis
                 tau = xtau(nt)
                 Om = OM_st + dble(nw-1)*dom
                 Xker_table(nt,nw) = Xker(tau,om,beta)
              enddo
           enddo
           ! Normalize data to have zeroth moment of unity.
           xqmc = xqmc / XMOM1
           cov = cov / ((XMOM1)**2)
           ! Diagonalize the covariance
           If (Lcov.eq.1) then
              Allocate( U(ntau,ntau), Sigma(ntau), xqmc1(Ntau) )
              Call Diag(cov,U,sigma)
              do nt = 1,ntau
                 sigma(nt) = sqrt(sigma(nt))
              enddo
              xqmc1 = 0.d0
              do nt1 = 1,ntau
                 do nt = 1,ntau
                    xqmc1(nt1) = xqmc1(nt1) + xqmc(nt)*U(nt,nt1)
                 enddo
                 xqmc1(nt1) = xqmc1(nt1)/sigma(nt1)
              enddo
              ! Transform the Kernel
              allocate ( Vhelp(Ntau) )
              do nw = 1,Ndis
                 Vhelp = 0.d0
                 do nt1 = 1,Ntau
                    do nt = 1,Ntau
                       Vhelp(nt1) = Vhelp(nt1) + Xker_table(nt,nw)*U(nt,nt1)
                    enddo
                 enddo
                 do nt1 = 1,ntau
                    Xker_table(nt1,nw) = Vhelp(nt1)/sigma(nt1) !! This has changed !!
                 enddo
              enddo
              Deallocate (U,sigma)
           else
              Allocate( Sigma(ntau), xqmc1(Ntau) )
              !Call Diag(cov,U,sigma)
              do nt = 1,ntau
                 sigma(nt) = 1.d0/sqrt(cov(nt,nt))
              enddo
              xqmc1 = 0.d0
              do nt1 = 1,ntau
                 xqmc1(nt1) = xqmc(nt1)*sigma(nt1)
              enddo
              ! Transform the Kernel
              allocate ( Vhelp(Ntau) )
              do nw = 1,Ndis
                 do nt1 = 1,ntau
                    Xker_table(nt1,nw) = Xker_table(nt1,nw)*sigma(nt1) !! This has changed !!
                 enddo
              enddo
              deallocate(Sigma)
           endif
           Allocate(G_Mean(Ntau))
           G_mean = 0.d0
           ! write(6,*) ' There are ', Ngamma,' delta-functions for a spectrum'
           ! Write(6,*) ' Initializing'
           Do Ns = 1,NSims
              do ng = 1,NGamma
                 Xn_tot(ng,1,ns) = ranf_wrap()
                 Xn_tot(ng,2,ns) = 1.d0/dble(Ngamma)
              enddo
           enddo
           Xn_m_tot = 0.d0
           En_m_tot = 0.d0
           Xn_e_tot = 0.d0
           En_e_tot = 0.d0
           ! D(om) = 1/(Om_en_1 - Om_st_1)
           D = 1.d0 / (Om_en_1 - Om_st_1)
           Open (Unit=44,File='Max_stoch_log', Status="unknown", position="append")
           Write(44,*) "N E W  R U N "
           nc = 0
           do Nb = 1,Nbins
              do ns = 1,NSims
                 do ng = 1,Ngamma
                    Xn(ng,1) = Xn_tot(ng,1,ns)
                    Xn(ng,2) = Xn_tot(ng,2,ns)
                 enddo
                 Alpha = Alpha_tot(ns)
                 Call MC(Xtau, Xker_table, Xn, Alpha, NSweeps, Xn_m, En, En_m, &
                      & Acc_1, Acc_2 ) ! Just one bin
                 do ng = 1,Ngamma
                    Xn_tot(ng,1,ns) = Xn(ng,1)
                    Xn_tot(ng,2,ns) = Xn(ng,2)
                 enddo
                 En_tot(ns) = En ! this is the energy of the configuration Xn_tot for simulation ns
                 Write(44,2003) 1.d0/Alpha, En_m, Acc_1, Acc_2
                 if (nb.gt.nwarm) then
                    if (ns.eq.1) nc = nc + 1
                    do nd = 1,NDis
                       Xn_m(nd) = Xn_m(nd) * D * dble(Ndis)
                       Xn_m_tot(nd,ns) = Xn_m_tot(nd,ns) + Xn_m(nd)
                       Xn_e_tot(nd,ns) = Xn_e_tot(nd,ns) + Xn_m(nd)*Xn_m(nd)
                    enddo
                    En_m_tot(ns) = En_m_tot(ns) + En_m
                    En_e_tot(ns) = En_e_tot(ns) + En_m*En_m
                 endif
              enddo
              ! Exchange
              Acc_1 = 0.d0
              Do Nex = 1, 2*NSims
                 nalp1= nint( ranf_wrap()*dble(NSims-1) + 0.5 ) ! 1..(NSims-1)
                 nalp2 = nalp1 + 1
                 DeltaE = (Alpha_tot(nalp1)*En_tot(nalp2) + Alpha_tot(nalp2)*En_tot(nalp1))&
                      & -(Alpha_tot(nalp1)*En_tot(nalp1) + Alpha_tot(nalp2)*En_tot(nalp2))
                 Ratio = exp(-DeltaE)
                 if (Ratio.gt.ranf_wrap()) Then
                    Acc_1 = Acc_1 + 1.0
                    !Switch confs an Energies.
                    do ng = 1,Ngamma
                       Xn(ng,1) = Xn_tot(ng,1,nalp1)
                       Xn(ng,2) = Xn_tot(ng,2,nalp1)
                    enddo
                    do ng = 1,Ngamma
                       Xn_tot(ng,1,nalp1) = Xn_tot(ng,1,nalp2)
                       Xn_tot(ng,2,nalp1) = Xn_tot(ng,2,nalp2)
                       Xn_tot(ng,1,nalp2) = Xn(ng,1)
                       Xn_tot(ng,2,nalp2) = Xn(ng,2)
                    enddo
                    En_m = En_tot(nalp1)
                    En_tot(nalp1) = En_tot(nalp2)
                    En_tot(nalp2) = En_m
                 endif
              enddo
              Acc_1 = Acc_1/dble(Nex)
              Write(44,*) 'Acc Exchange: ', Acc_1
           enddo
           !Open(Unit=66,File="energies",status="unknown")
           do ns = Nsims,Nsims
              En_m_tot(ns) = En_m_tot(ns) / dble(nc)
              En_e_tot(ns) = En_e_tot(ns) / dble(nc)
              En_e_tot(ns) = ( En_e_tot(ns) - En_m_tot(ns)**2)/dble(nc)
              if ( En_e_tot(ns) .gt. 0.d0) then
                 En_e_tot(ns) = sqrt(En_e_tot(ns))
              else
                 En_e_tot(ns) = 0.d0
              endif
              !write(66,*) Alpha_tot(ns), En_m_tot(ns), En_e_tot(ns)
           enddo
           !close(60)
           Chisq = En_e_tot(Nsims)
           Close(44)
           do ns = Nsims,Nsims
              do nd = 1,Ndis
                 Xn_m_tot(nd,ns) = Xn_m_tot(nd,ns) / dble(nc) ! * delta /(dble(nc)*pi)
                 Xn_e_tot(nd,ns) = Xn_e_tot(nd,ns) / dble(nc) ! * delta /(dble(nc)*pi)
                 Xn_e_tot(nd,ns) = (Xn_e_tot(nd,ns) - Xn_m_tot(nd,ns)* Xn_m_tot(nd,ns))/dble(nc)
                 if (Xn_e_tot(nd,ns).gt.0.d0) then
                    Xn_e_tot(nd,ns) = sqrt(Xn_e_tot(nd,ns))
                 else
                    Xn_e_tot(nd,ns) = 0.d0
                 endif
                 om = PhiM1(dble(nd)/dble(NDis))
                 Aom = Xn_m_tot(nd,ns) * Xmom1
                 Err = Xn_e_tot(nd,ns) * Xmom1
                 !write(66,2001) om, Aom, Err
                 if (ns.eq.Nsims) then
                    Aom_res(nd) = Aom
                    xom_res(nd) = om
                 endif
              enddo
              !Close(66)
           enddo
           ! Reset the input data
           xqmc = XMOM1* xqmc
           cov = ((XMOM1)**2)* cov
           DeAllocate (Xn_tot)
           DeAllocate (En_m_tot, En_e_tot, En_tot )
           DeAllocate (Xn_m_tot, Xn_e_tot )
           DeAllocate (Xn)
           DeAllocate (Xn_m, Xn_e)
           DeAllocate( G_Mean )
           DeAllocate( xqmc1 )
           Deallocate( Xker_table )
           DeAllocate(Iseed_vec)

2001     format(F14.7,2x,F14.7,2x,F14.7)
2004     format(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)
2005     format(F14.7,2x,F14.7,2x,F14.7,2x,F14.7,2x,F14.7)
2003     format('Alpha, En_m, Acc ', F14.7,2x,F14.7,2x,F14.7,2x,F14.7,2x,F14.7)
         end Subroutine MaxEnt_stoch_fit
 !***********
         Real (Kind=8) Function Phim1(x)
           Implicit None
           ! Flat Default with sum 1. This is the correct sum rule for the data!
           ! D(om) = 1/(Om_en_1 - Om_st_1)
           Real (Kind=8) :: x
           PhiM1 = x*(Om_en_1 - Om_st_1) + Om_st_1
         end Function Phim1
         
         Integer Function NPhim1(x)
           Implicit None
           ! Flat Default with sum 1. This is the correct sum rule for the data!
           ! D(om) = 1/(Om_en_1 - Om_st_1)
           Real (Kind=8) :: x, om
           om = x*(Om_en_1 - Om_st_1) + Om_st_1
           NPhiM1 = Nint ( (om - Om_st_1)/Dom + 0.75 )
         end Function NPhim1
!***********
         Subroutine Sum_Xn(Xn_m,Xn)
           Implicit none
           Real (Kind=8), Dimension(:,:) :: Xn
           Real (Kind=8), Dimension(:) :: Xn_m
           Real (Kind=8) :: X
           do nd = 1,NDis
              X = dble( nd )/dble( NDis )
              do ng = 1,Ngamma
                 Xn_m(nd) = Xn_m(nd) + Xn(ng,2)/( (X-Xn(ng,1))**2 + Delta2)
                 !aimag( cmplx(Xn(ng,2),0.d0)/cmplx( X-Xn(ng,1), -Delta) )
              enddo
           enddo
         end Subroutine Sum_Xn
!***********
         Subroutine Sum_Xn_Boxes(Xn_m,Xn)
           Implicit none
           Real (Kind=8), Dimension(:,:) :: Xn
           Real (Kind=8), Dimension(:) :: Xn_m
           Real (Kind=8) :: X
           do ng = 1,Ngamma
              X = Xn(ng,1)
              nd = Nint(dble(NDis)*X + 0.5 )
              Xn_m(nd) = Xn_m(nd) + Xn(ng,2)
           Enddo
         end Subroutine Sum_Xn_Boxes
!***********
         Subroutine MC(Xtau, Xker_table, Xn, Alpha, NSweeps, Xn_m, En, En_m, Acc_1,Acc_2)
           !Implicit Real (KIND=8) (A-G,O-Z)
           !Implicit Integer (H-N)
           Implicit None
           Real (Kind=8), Dimension(:,:) :: Xn, Xker_table
           Real (Kind=8), Dimension(:) :: Xtau, Xn_m
           Real (Kind=8) :: Alpha, En_m, s, ratio, A_gamma, Z_gamma, Acc_1, Acc_2
           Integer :: NSweeps, nl, Lambda_max, ng1, ng2
           !Local
           Real (Kind=8), Allocatable :: h(:), Deltah(:), A_gamma_p(:), Z_gamma_p(:), &
                & A_gamma_o(:), Z_gamma_o(:)
           Real (Kind=8), Allocatable :: XKER_stor(:,:), XKER_new(:)
           Real (Kind=8) :: X, En, DeltaE
           Integer, Allocatable :: Lambda(:)
           Integer :: nsw, Nacc_1, Nacc_2, nw
           Allocate (h(ntau), Deltah(ntau) )
           Allocate (Lambda(2), Z_gamma_p(2), A_gamma_p(2), &
                & Z_gamma_o(2), A_gamma_o(2) ) ! Max of moves of two walkers.
           Allocate ( XKer_stor(Ntau,Ngamma), XKer_New(Ntau) )
           Xn_m = 0.d0
           En_m = 0.d0
           ! Setup h(tau)
           do nt = 1,Ntau
              X = 0.d0
              do ng = 1,Ngamma
                 A_gamma = xn(ng,1)
                 Z_gamma = xn(ng,2)
                 XKer_stor( nt, ng ) = XKER_table(nt, NPhiM1(A_gamma) )
                 ! XKER(xtau(nt),PhiM1(A_gamma),beta)
                 X = X + Xker_stor(nt,ng)*Z_gamma
              enddo
              h(nt) = X - xqmc1(nt)
           enddo
              NAcc_1 = 0; NAcc_2 = 0;
              do nsw = 1,Nsweeps
                 ! Weight sharing moves.
                 do ng = 1,Ngamma
                    x = ranf_wrap()
                    if (x.gt.0.5) then
                       ! Weight sharing moves.
                       Lambda_max = 2
                       Lambda(1) = nint(ranf_wrap()*dble(Ngamma) + 0.5)
                       do
                          Lambda(2) = nint(ranf_wrap()*dble(Ngamma) + 0.5)
                          if ( Lambda(2) .ne. Lambda(1) ) exit
                       enddo
                       ng1 = Lambda(1)
                       ng2 = Lambda(2)
                       A_gamma_o(1) = Xn(ng1,1)
                       A_gamma_o(2) = Xn(ng2,1)
                       Z_gamma_o(1) = Xn(ng1,2)
                       Z_gamma_o(2) = Xn(ng2,2)
                       A_gamma_p(1) = Xn(ng1,1)
                       A_gamma_p(2) = Xn(ng2,1)
                       s = (Z_gamma_o(1) + Z_gamma_o(2))*ranf_wrap() - Z_gamma_o(1)
                       Z_gamma_p(1) = Z_gamma_o(1) + s
                       Z_gamma_p(2) = Z_gamma_o(2) - s
                       ! Kernel stays unchanged.
                       ! Compute Delta H
                       do nt = 1,ntau
                          X = Xker_stor(nt,ng1)*( Z_gamma_p(1) - Z_gamma_o(1) ) + &
                               & Xker_stor(nt,ng2)*( Z_gamma_p(2) - Z_gamma_o(2) )
                          Deltah(nt) = X
                       enddo
                    else
                       Lambda_max = 1
                       Lambda(1) = nint(ranf_wrap()*dble(Ngamma) + 0.5)
                       ng1 = Lambda(1)
                       Z_gamma_o(1) = Xn(ng1,2)
                       Z_gamma_p(1) = Xn(ng1,2)
                       A_gamma_o(1) = Xn(ng1,1)
                       A_gamma_p(1) = xpbc( Xn(ng1,1) + (ranf_wrap() - 0.5)*DeltaXMAX, 1.d0 )
                       !om = PhiM1(A_gamma_p(1))
                       nw = NPhiM1(A_gamma_p(1))
                       Xker_new = Xker_table(:, nw)
                       Deltah = ( Xker_new - Xker_stor(:, ng1) ) * Z_gamma_o(1)
                    endif
                    DeltaE = 0.d0
                    do nt = 1,ntau
                       DeltaE = DeltaE + (Deltah(nt) + 2.0 * h(nt) ) *Deltah(nt)
                    enddo
                    Ratio = exp( -alpha * DeltaE )
                    ! write(6,*) ' Ratio : ',Ratio, DeltaE
                    if (Ratio .gt. ranf_wrap()) Then
                       ! write(6,*) 'Accepted'
                       if (Lambda_max.eq.1) then
                          Nacc_1 = Nacc_1 + 1
                          ng1 = Lambda(1)
                          Xker_stor(:, ng1) = Xker_new
                       endif
                       if (Lambda_max.eq.2) Nacc_2 = Nacc_2 + 1
                       do nl = 1,Lambda_max
                          Xn(Lambda(nl),1) = A_gamma_p(nl)
                          Xn(Lambda(nl),2) = Z_gamma_p(nl)
                      enddo
                      h = h + Deltah
                    endif
                 enddo
                 En = dot_product(h, h)
                 En_m = En_m + En
                 Call Sum_Xn_Boxes( Xn_m, Xn )
              enddo
              Acc_1 = dble(Nacc_1)/dble(Ngamma*NSweeps)
              Acc_2 = dble(Nacc_2)/dble(Ngamma*NSweeps)
              En_m = En_m/dble( nsweeps )
              Xn_m = Xn_m/dble( nsweeps )
           Deallocate ( h, Deltah )
           Deallocate ( Lambda, Z_gamma_p, A_gamma_p, Z_gamma_o, A_gamma_o )
           Deallocate ( XKER_stor, XKER_new )
2005       format(I4,2x,I4,2x,F14.7,2x,F14.7,' --> ',F14.7,2x,F14.7)
2006       format(I4,2x,F14.7, ' --> ',F14.7)
         end Subroutine MC
!**********
         real (Kind=8) function xpbc(X,XL)
           real (kind=8) :: X, XL
           XPBC = X
           if (X.GT. XL ) XPBC = X - XL
           if (X.LT. 0.0) XPBC = X + XL
         end function xpbc
end Module MaxEnt_stoch_mod
