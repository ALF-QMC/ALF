!  Copyright (C) 2018-2024 The ALF project
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


Module MaxEnt_stoch_mod

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This module provides  an implementation of  the  stochastic  stochastic  analytical 
!> continuation.  It  allows to  specify  a  default  model,  and  follows  rather precicely the 
!> article  of  Kevin Beach https://arxiv.org/abs/cond-mat/0403055 
!> 
!
!--------------------------------------------------------------------
       Use runtime_error_mod
       Use MyMats
       Use Random_Wrap
       Use Files_mod
       use iso_fortran_env, only: output_unit, error_unit

       Integer, private :: NTAU, nt, Ngamma, ng, Ndis, nd,  L_seed
       Integer, private, allocatable:: Iseed_vec(:)
       Real (Kind=Kind(0.d0)), private :: Delta, Delta2, OM_st_1, Om_en_1, DeltaXMAX, Beta, Pi, Dom_table, Dom_spectral, &
                                      &   Dx_spectral, Dx_table
       Real (Kind=Kind(0.d0)), allocatable, private :: XQMC1(:)
       Integer, allocatable,  private ::  Phim1_func(:), Phi_func(:)
       
       contains

!--------------------------------------------------------------------
         Subroutine MaxEnt_stoch(XQMC, Xtau, COV,Xmom1, XKER, Back_Trans_Aom, Beta_1, Alpha_tot,&
              & Ngamma_1, OM_ST, OM_EN, Ndis_1, Nsweeps, NBins, NWarm, F, Default_provided)

           Implicit None

           Real (Kind=Kind(0.d0)), Dimension(:) :: XQMC, Xtau, Alpha_tot
           Real (Kind=Kind(0.d0)), Dimension(:,:) :: COV
           Real (Kind=Kind(0.d0)), Dimension(:),  Intent(In), allocatable,   optional :: Default_provided
           Real (Kind=Kind(0.d0)), External :: XKER, Back_trans_Aom, F
           Real (Kind=Kind(0.d0)) :: OM_ST, OM_EN, Beta_1, Xmom1, Err
           Integer :: Nsweeps, NBins, Ngamma_1, Ndis_1, nw, nt1

           ! Local
           Integer NSims, ns, nb, nc, Nwarm, nalp1, nalp2, Nex, p_star, Ndis_table, &
                & io_error, io_error1, i, n, nc1, nx
           Real (Kind=Kind(0.d0)), Allocatable :: Xn_M_tot(:,:), En_M_tot(:), Xn_E_tot(:,:), En_E_tot(:), &
                & Xn_tot(:,:,:), En_tot(:), Default(:)
           Real (Kind=Kind(0.d0)), Allocatable :: G_Mean(:), Xn_m(:), Xn_e(:), Xn(:,:), Vhelp(:), Default_table(:), A(:)
           Real (Kind=Kind(0.d0)) :: En_M, X, Alpha, Acc_1, Acc_2, En, DeltaE, Ratio
           Real (Kind=Kind(0.d0)) :: Aom, om, XMAX, tau
           Real (Kind=Kind(0.d0)) :: CPUT, F_A
           Integer :: ICPU_1, ICPU_2, N_P_SEC
           Character (64) :: File_root, File1, File_conf, File_Aom
           Real (Kind=Kind(0.d0)), allocatable :: Xker_table(:,:), U(:,:), sigma(:)
           ! Space for moments.
           Real (Kind=Kind(0.d0)), allocatable:: Mom_M_tot(:,:), Mom_E_tot(:,:)
           Real (Kind=Kind(0.d0)), allocatable ::  F_A_m(:), F_A_e(:)
           

           Pi        = acos(-1.d0)
           NDis      = Ndis_1
           DeltaXMAX = 0.01
           delta     = 0.001
           delta2    = delta*delta
           Ngamma    = Ngamma_1
           Beta      = Beta_1 ! Physical temperature for calculation of the kernel.
           Ntau      = Size(xqmc,1)
           NSims = Size(Alpha_tot,1)
           Allocate (Xn_tot(Ngamma,2,NSims))
           Allocate (En_m_tot(NSims), En_e_tot(NSims), En_tot(NSims) )
           Allocate (F_A_m(NSims), F_A_e(NSims))
           Allocate (Xn_m_tot(NDis,NSims), Xn_e_tot(NDis,NSims) )
           Allocate (Mom_M_tot(4,Nsims), Mom_E_tot(4,Nsims) )
           Allocate (Xn(Ngamma,2))
           Allocate (Xn_m(NDis), Xn_e(NDis), A(Ndis) )
           Om_st_1 = OM_st; Om_en_1 = OM_en
           ! Setup table for the Kernel
           Ndis_table   = 50000
           Dx_table     =  1.d0 /dble(Ndis_table)
           Dom_table    = (OM_EN_1 - OM_ST_1)/dble(Ndis_table)
           Dom_spectral = (OM_EN_1 - OM_ST_1)/dble(Ndis)
           Dx_spectral  = 1.d0/dble(Ndis)
           Allocate ( Xker_table(Ntau, Ndis_table) )
           Allocate ( Default_table(Ndis_table))
           do nt = 1,Ntau
              do nw = 1,Ndis_table
                 tau = xtau(nt)
                 Om = OM_st + dble(nw-1)*Dom_table
                 Xker_table(nt,nw) = Xker(tau,om,beta)
              enddo
           enddo
           ! Normalize  the  data
           xqmc = xqmc / XMOM1
           cov  = cov / ((XMOM1)**2)
           ! Diagonalize the covariance
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
           
           ! Setup  the  Default.
           !D = Xmom1 / (Om_en_1 - Om_st_1)  !  Flat default with correct  sum-rule
           allocate(Default(Ndis))
           IF (.not.Present(Default_provided))  then 
              Default  =  Xmom1 / (Om_en_1 - Om_st_1)  !  Flat default with correct  sum-rule
           else
              If (size(Default_provided,1) .ne. size(Default,1) )  then 
                 write(error_unit,*) 'Default_provided in MaxEnt_stoch has wrong dimensions'
                 CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
              endif
              Default =  Default_provided
           endif
           Call Set_default_table(Default, Default_table, Xmom1)
           Call Get_seed_Len(L_seed)
           Allocate(Iseed_vec(L_seed))
           Iseed_vec = 0
           Call Ranset(Iseed_vec)
           File_conf = "dump_conf"
           File_Aom  = "dump_Aom"
           Open(unit=41,file=File_conf,status='old',action='read', iostat=io_error)
           Open(unit=42,file=File_Aom, status='old',action='read', iostat=io_error1)
           If (io_error == 0 .and. io_error1 == 0 ) then
              Nwarm = 0
              read(41,*) (Iseed_vec(I), I = 1,size(Iseed_vec,1))
              do ns = 1,Nsims
                 do ng = 1,Ngamma
                    read(41,*) Xn_tot(ng,1,ns), Xn_tot(ng,2,ns)
                 enddo
                 read(41,*) En_m_tot(ns), En_e_tot(ns), F_A_m(ns), F_A_e(ns)
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
              F_A_m    = 0.d0
              F_A_e    = 0.d0
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
                    
                    !  Given   Xn_m  =  n(x) dx   with  x =  phi(om)  the  following  loop  sets:
                    !  A(om)   =  n(phi(om)) * D (om )
                    do nd = 1,NDis
                       om = om_st_1 + dble(nd -1)*Dom_spectral
                       nw = Int((om - OM_st_1)/dom_table) + 1    ! Index  for  default_table
                       x  = Phi(om)
                       nx = Int(x/Dx_spectral) +  1              ! Index  for box  n(x)
                       A(nd)    = Xn_m(nx)*Default_table(nw)/Dx_spectral
                    enddo
                    Xn_m = A
                    do  nd  =  1,Ndis
                       Xn_m_tot(nd,ns) = Xn_m_tot(nd,ns) + Xn_m(nd)
                       Xn_e_tot(nd,ns) = Xn_e_tot(nd,ns) + Xn_m(nd)*Xn_m(nd)
                    enddo
                    En_m_tot(ns) = En_m_tot(ns) + En_m
                    En_e_tot(ns) = En_e_tot(ns) + En_m*En_m
                    ! Compute  (F,A)
                    F_A =  0.d0 
                    do ng = 1,Ngamma
                       F_A = F_A +   F(Phim1(Xn_tot(ng,1,ns)),beta) * Xn_tot(ng,2,Ns)
                    enddo  
                    F_A =  F_A *Xmom1
                    F_A_m(ns) = F_A_m(ns) + F_A
                    F_A_e(ns) = F_A_e(ns) + F_A*F_A
                    ! End compute  (F,A)
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
              write(41,*) En_m_tot(ns), En_e_tot(ns), F_A_m(ns), F_A_e(ns)
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
           Open(Unit=67,File="(F,A).dat",status="unknown")
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

              F_A_m(ns) = F_A_m(ns) / dble(nc)
              F_A_e(ns) = F_A_e(ns) / dble(nc)
              F_A_e(ns) = ( F_A_e(ns) -F_A_m(ns)**2)/dble(nc)
              if ( F_A_e(ns) .gt. 0.d0) then
                 F_A_e(ns) = sqrt(F_A_e(ns))
              else
                 F_A_e(ns) = 0.d0
              endif
              write(67,*) Alpha_tot(ns), F_A_m(ns), F_A_e(ns)
           enddo
           close(66)
           close(67)
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
              write(66,"(F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,F12.6,2x,F12.6)") &
                   & Alpha_tot(ns), Mom_m_tot(1,ns), Mom_e_tot(1,ns), &
                   & Mom_m_tot(2,ns), Mom_e_tot(2,ns),Mom_m_tot(3,ns), Mom_e_tot(3,ns),  &
                   & Mom_m_tot(4,ns), Mom_e_tot(4,ns) 
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
                 om =  Om_st_1 + dble(nd-1)*Dom_spectral ! PhiM1(dble(nd)/dble(NDis)) HERE
                 Aom = Xn_m_tot(nd,ns) ! * Xmom1
                 Err = Xn_e_tot(nd,ns) ! * Xmom1
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
                 om =  Om_st_1 +  dble(nd-1)*Dom_spectral  !  PhiM1(dble(nd)/dble(NDis)) HERE
                 Aom = Xn_m(nd) ! * Xmom1
                 Err = Xn_e(nd) ! * Xmom1
                 Xn_m(nd) = Back_Trans_Aom(Aom,Beta,om)
                 Xn_e(nd) = Back_Trans_Aom(Err,Beta,om)
                 IF (Xn_m(nd) .gt. XMAX ) XMAX = Xn_m(nd)
              enddo
              do nd = 1,Ndis
                 om =  Om_st_1 +  dble(nd-1)*Dom_spectral  !  PhiM1(dble(nd)/dble(NDis)) HERE
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
           DeAllocate (F_A_e, F_A_m)
           DeAllocate (Xn)
           DeAllocate (Xn_m, Xn_e)
           DeAllocate( G_Mean )
           DeAllocate( xqmc1 )
           Deallocate( Xker_table )
2001       format(F14.7,2x,F14.7,2x,F14.7)
!2004       format(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)
2005       format(F14.7,2x,F14.7,2x,F14.7,2x,F14.7,2x,F14.7)
2003       format('Alpha, En_m, Acc ', F14.7,2x,F24.12,2x,F14.7,2x,F14.7,2x,F14.7)
           
         end Subroutine MaxEnt_stoch
!------------------------------------------------------------------------------------------------
!        Sets various  tables  for  the  Default  model 
!------------------------------------------------------------------------------------------------
         Subroutine Set_default_table(Default, Default_table, Xmom1)

            Implicit none 
            Real (Kind=Kind(0.d0)), allocatable, Dimension(:), Intent(In) ::  Default          
            Real (Kind=Kind(0.d0)), allocatable, Dimension(:), Intent(InOut) ::  Default_table
            Real (Kind=Kind(0.d0)), Intent(In)::  Xmom1

            Integer :: nw, nw1, nw_d, nx
            Real (Kind= Kind(0.d0)) :: om,  dom, a, b, x, x1, f1,f2
            Logical ::  Test=.false.
         

            dom = (Om_en_1 -  Om_st_1)/dble(Ndis)
            nw_d = 1
            f1 = Default(nw_d); f2 = Default(nw_d + 1) 
            x1 = Om_st_1 +  dble(nw_d) * dom
            b = (f2 - f1 ) / dom ;  a = f1 - b *x1
            Do nw = 1,Size(Default_table,1)
               om = Om_st_1 + dble(nw)*Dom_table
               Default_table(nw) = a + b*om
               if (om > x1 + dom .and. nw_d < Ndis-1 )   then
                  nw_d = nw_d + 1 
                  f1 = Default(nw_d); f2 = Default(nw_d + 1) 
                  x1 = Om_st_1 +  dble(nw_d) * dom
                  b = (f2 - f1 ) / dom ;  a = f1 - b *x1
               endif
            enddo

            !Reimpose  the  sum-rule  since  the linear  interpolation will   violate it.
            x1 = 0.d0
            Do nw = 1,Size(Default_table,1)
               x1 =  x1 + Default_table(nw)
            enddo
            x1 = x1 * Dom_table 
            Default_table  =   Default_table * Xmom1/X1

            ! Setup  Phimi_func
            Allocate(Phim1_func(Size(Default_table,1)), Phi_func(Size(Default_table,1)))
            x = 0.D0 
            Do nx = 1,Size(Default_table,1)
               x = dble(nx)*Dx_table
               x1 = 0.d0
               do nw1 = 1,Size(Default_table,1) 
                  x1 = x1 + dom_table*Default_table(nw1)/Xmom1
                  if (x1 > x ) exit
               enddo
               Phim1_func(nx) = nw1
            enddo
            Do nw = 1,Size(Default_table,1)
               x1 = 0.d0
               do nw1 = 1,nw
                  x1 = x1 + dom_table*Default_table(nw1)/Xmom1
               enddo
               nx = Int(x1/Dx_table) + 1
               Phi_func(nw) = nx
            enddo
            If (Test)   then 
               ! Check the Phim1   function
               do nx = 1,Size(Default_table,1)
                  x = dble(nx)*dx_table
                  om = Om_st_1 + dble(nx-1)*dom_table
                  write(13,*) x, Phim1(x), om, Phi(om), Phi(Phim1(x))
               enddo
            endif
         end Subroutine Set_default_table

!-------------------
!         Subroutine MaxEnt_stoch_fit(XQMC, Xtau, COV, Lcov, XKER, Xmom1, Beta_1, Alpha_tot,&
!                 & Ngamma_1, OM_ST, OM_EN, Nsweeps, NBins, NWarm, Aom_res,&
!                 & xom_res, Chisq )
!              
!           Implicit None
!           
!           Real (Kind=Kind(0.d0)), Dimension(:) :: XQMC, Xtau, Alpha_tot, Aom_res, Xom_res
!           Real (Kind=Kind(0.d0)), Dimension(:,:) :: COV
!           Real (Kind=Kind(0.d0)), external :: XKER
!           Real (Kind=Kind(0.d0)) :: CHISQ, OM_ST, OM_EN, Beta_1, Xmom1, Err
!           Integer :: Nsweeps, NBins, Ngamma_1, nw, nt1, Lcov

!           ! Reset the input data
!           xqmc = XMOM1* xqmc
!           cov = ((XMOM1)**2)* cov
!         end Subroutine MaxEnt_stoch_fit
!-------------------
 
!--------------------------------------------------------------------------------------
!        Uses   the  Phi_func  table  to   generate   
!        Phi(om) =  (1/Xmom1) \int_{-inf}^{om} D(om).  Here  D(om) is  normalized to Xmom1
!--------------------------------------------------------------------------------------
         Real (Kind=Kind(0.d0)) Function Phi(om)

            Implicit None
            ! Flat Default with sum xmom1. 
            ! D(om) = xmom1/(Om_en_1 - Om_st_1)
            Real (Kind=Kind(0.d0)) :: om
            Integer ::  nw

            nw = int((om-OM_st_1)/Dom_table) + 1
            If (nw > size(Phi_func,1) ) nw = size(Phi_func,1) 
            Phi =  dble(Phi_func(nw)-1)*Dx_table
         end Function Phi
!--------------------------------------------------------------------------------------
!        Uses   the  Phmi_func  table  to   generate   
!        Phim1(x) = om  for    Phi(x) =  om
!--------------------------------------------------------------------------------------

         Real (Kind=Kind(0.d0)) Function Phim1(x)

            Implicit None
            ! Flat Default with sum xmom1. 
            ! D(om) = xmom1/(Om_en_1 - Om_st_1)
            Real (Kind=Kind(0.d0)) :: x,  test
            Integer ::  nx

            nx = int(x/Dx_table) + 1
            If (nx > size(Phim1_func,1) ) nx = size(Phim1_func,1) 
            PhiM1 = Om_st_1 + dble(Phim1_func(nx)-1)*Dom_table
            !PhiM1 = x*(Om_en_1 - Om_st_1) + Om_st_1  !  This is  the special  case  for a  flat  default
         end Function Phim1
         
         Integer Function NPhim1(x)
           Implicit None
           ! Flat Default with sum  xmom1
           ! D(om) = xmom1/(Om_en_1 - Om_st_1)
           Real (Kind=Kind(0.d0)) :: x, om
           om = Phim1(x)
           NPhiM1 = int ( (om - Om_st_1)/Dom_table ) + 1 
         end Function NPhim1

!--------------------------------------------------------------------------------------
!        This   function is  outdated
!--------------------------------------------------------------------------------------
         Subroutine Sum_Xn(Xn_m,Xn)
           Implicit none
           Real (Kind=Kind(0.d0)), Dimension(:,:) :: Xn
           Real (Kind=Kind(0.d0)), Dimension(:) :: Xn_m
           Real (Kind=Kind(0.d0)) :: X
           do nd = 1,NDis
              X = dble( nd )/dble( NDis )
              do ng = 1,Ngamma
                 Xn_m(nd) = Xn_m(nd) + Xn(ng,2)/( (X-Xn(ng,1))**2 + Delta2)
                 !aimag( cmplx(Xn(ng,2),0.d0)/cmplx( X-Xn(ng,1), -Delta) )
              enddo
           enddo
         end Subroutine Sum_Xn

!--------------------------------------------------------------------------------------
!        Given        n(x)   =  sum_{ng=1}^{\gamma} xn(ng,2) \delta[x - xn(ng,2)] 
!        Cummulates   Xn_m(n) =  n(x) dx      with  box  distributions
!--------------------------------------------------------------------------------------
         Subroutine Sum_Xn_Boxes(Xn_m,Xn)
           Implicit none
           Real (Kind=Kind(0.d0)), Dimension(:,:) :: Xn
           Real (Kind=Kind(0.d0)), Dimension(:) :: Xn_m
           Real (Kind=Kind(0.d0)) :: om  
           Integer :: nd, ng

           do ng = 1,Ngamma
             nd  = int (xn(ng,1)/Dx_spectral) + 1
             Xn_m(nd) = Xn_m(nd) + Xn(ng,2)
           enddo


         end Subroutine Sum_Xn_Boxes

!--------------------------------------------------------------------------------------
!        Carries out  local  updates.
!--------------------------------------------------------------------------------------

         Subroutine MC(Xtau, Xker_table, Xn, Alpha, NSweeps, Xn_m, En, En_m, Acc_1,Acc_2)
           Implicit None
           Real (Kind=Kind(0.d0)), Dimension(:,:) :: Xn, Xker_table
           Real (Kind=Kind(0.d0)), Dimension(:) :: Xtau, Xn_m
           Real (Kind=Kind(0.d0)) :: Alpha, En_m, s, ratio, A_gamma, Z_gamma, Acc_1, Acc_2
           Integer :: NSweeps, nl, Lambda_max, ng1, ng2
           !Local
           Real (Kind=Kind(0.d0)), Allocatable :: h(:), Deltah(:), A_gamma_p(:), Z_gamma_p(:), &
                & A_gamma_o(:), Z_gamma_o(:)
           Real (Kind=Kind(0.d0)), Allocatable :: XKER_stor(:,:), XKER_new(:)
           Real (Kind=Kind(0.d0)) :: X, En, DeltaE
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
                    ! Shifting  delta_functions
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
         end Subroutine MC
!**********
         real (Kind=Kind(0.d0)) function xpbc(X,XL)
           real (Kind=Kind(0.d0)) :: X, XL
           XPBC = X
           if (X.GT. XL ) XPBC = X - XL
           if (X.LT. 0.0) XPBC = X + XL
         end function xpbc
end Module MaxEnt_stoch_mod
