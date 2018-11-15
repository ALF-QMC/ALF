     Program Test

       Use MaxEnt_stoch_mod
       use Files_mod
       
       Implicit none
       
       Real (Kind=Kind(0.d0)), Dimension(:), allocatable :: XQMC, XTAU, Alpha_tot, xom, A, &
            &    ERROR, ARES
       Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: Xcov
       
       Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: Xcov_st
       Real (Kind=Kind(0.d0)), Dimension(:), allocatable :: Xqmc_st, Xtau_st
       Real (Kind=Kind(0.d0)), External :: XKER, Back_trans_Aom, F_Fit
       Character (len=64) :: File1, File2
       Character (len=1)  :: Fermion_type
       Complex (Kind=Kind(0.d0)) :: Z
       Character (Len=64) :: command
       Logical            :: skip

       Integer ::  N_alpha, nc, nw, Ndis, nwp, nbins, NSweeps, NWarm, Ntau, Ngamma, nt, N_alpha_1, io_error
       Integer :: L_cov, Lt, nt_st, nt_en, nt1, nbins_qmc
       Real (Kind=Kind(0.d0)) :: alpha_st, R, Xmom1, X,  X1, X2, om, omp, pi, Om_st, Om_en, err
       Real (Kind=Kind(0.d0)) :: Dom, delta, beta, dtau

       ! Manually set number of OpenMM threads to 1
       ! !$ essentially asks the compiler, if OpenMP is available/activated
       !    -> if yes, then call function
       !    -> if no,  then no need to call function (linker issue if still called)
       !$ Call omp_set_num_threads(1)

       
       open(unit=30,file='paramSAC',status='old',action='read', iostat=io_error) 
       Open(unit=50,File='Info',Status="unknown")
       if (io_error.eq.0) then
          write(50,*) 'Reading in from paramSAC'
          read(30,*)  Ngamma, OM_st, OM_en, Ndis,  NBins, NSweeps, Nwarm
          read(30,*)  N_alpha, alpha_st, R
          read(30,*)  L_cov
       else
          write(6,*) 'No file paramSAC! ' 
          stop
       endif
       close(30)
       

       ! Read the data
       open (unit=10,File="g_dat", status="unknown")
       Read(10,*)  Lt, nbins_qmc
       Allocate ( Xtau_st(LT), Xqmc_st(Lt), Xcov_st(Lt,Lt) )
       Xcov_st = 0.d0
       Do nt = 1,Lt
          Read(10,*) Xtau_st(nt), Xqmc_st(nt), err
          xcov_st(nt,nt) = err*err
       Enddo
       If (L_cov == 1 ) then
          Do nt = 1,Lt
             Do nt1 = 1,Lt
                Read(10,*) xcov_st(nt1,nt) 
             Enddo
          Enddo
       endif
       close(10)

       !Set the range 
       nt_st = 1
       if (  sqrt(xcov_st(1,1)) < 10.D-8) nt_st=2
       do nt = 1,lt
          X1 = Xqmc_st(nt)
          X2 = sqrt(xcov_st(nt,nt))
          nt_en = nt -1 
          if ( abs(X1) - abs(X2)  < 0.d0 .or. X1 < 0.d0 ) exit
       enddo
       

       Xmom1 =  Xqmc_st(1)

       Ntau = nt_en - (nt_st-1)
       Allocate ( Xcov(ntau,ntau), Xqmc(ntau), Xtau(ntau) )
       do nt = 1, Ntau
          xqmc(nt)  =  Xqmc_st(nt + (nt_st-1) )
          xtau(nt)  =  xtau_st(nt + (nt_st-1) )
       enddo
       do nt = 1,Ntau
          do nt1 = 1,Ntau
             xcov(nt1,nt)  =  xcov_st(nt1 + (nt_st-1), nt + (nt_st-1) )
          enddo
       enddo
       deallocate ( Xtau_st, Xqmc_st, Xcov_st )
       Allocate (Xqmc_st(ntau))
       Xqmc_st = Xqmc
       dtau = Xtau(2) - Xtau(1)
          
       Allocate (Alpha_tot(N_alpha) )
       do nt = 1,N_alpha
          alpha_tot(nt) = alpha_st*(R**(nt-1))
       enddo
       write(50,*) 'First Moment  ', Xmom1
       close(50)

       if (nt_en >= nbins_qmc/3) nt_en =  nbins_qmc/3
       Write(6,*) nt_st, nt_en
       
          
       If (L_cov == 1 ) then
          Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER, Back_Trans_Aom, Beta, &
               &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, L_cov) 
       else
          Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER, Back_Trans_Aom, Beta, &
               &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm)
       endif
       
       N_alpha_1 = N_alpha - 10  
       File1 ="Aom_ps"
       file2 =File_i(File1,N_alpha_1)
       Write(50,*) file2
       Open(Unit=66,File=file2,status="unknown")
       Allocate (xom(Ndis), A(Ndis))
       do nw = 1,Ndis
          read(66,*) xom(nw), A(nw), x, x1, x2
       enddo
       close(66)
       Dom = xom(2) - xom(1)
       
       
       ! Compute the real frequency Green function.
       delta = 4.0*Dom
       Open (Unit=43,File="Green", Status="unknown", action="write")
       pi = acos(-1.d0)
       do nw = 1,Ndis
          Z = cmplx(0.d0,0.d0,Kind(0.d0))
          om = xom(nw)
          do nwp = 1,Ndis
             omp = xom(nwp)
             Z = Z + A(nwp)/cmplx( om -  omp, delta, kind(0.d0))
          enddo
          Z = Z * dom
          write(43,"('X'2x,F14.7,2x,F16.8,2x,F16.8)" ) xom(nw), dble(Z), -Aimag(Z)/pi
       enddo
       close(43)
       
       
       Open (Unit=70,File="data_out", status="unknown")
       do nt = 1,Ntau
          X = 0.d0
          do  nw = 1, Ndis
             X = X + Xker(Xtau(nt),xom(nw), beta)*A(nw)
          enddo
          X = X *dom
          write(70,2005) xtau(nt), xqmc_st(nt), sqrt(xcov(nt,nt)), X
       enddo
2005   format(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)
       Close(70)
       
       Command = "rm dump*"
       Call System (Command)
       Command = "ls"
       Call System (Command)
       
          
     end Program Test
     
     
     
     Real (Kind=Kind(0.d0)) function XKER(tau,om, beta)
       
       Implicit None
       real (Kind=Kind(0.d0)) :: tau, om, pi, beta

       !pi = 3.1415927
       XKER = exp(-tau*om) !/ ( 1.d0 + exp(-Beta*om) ) ! /pi

     end function XKER

     Real (Kind=Kind(0.d0)) function Back_trans_Aom(Aom, om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) ::  Aom, om, beta

       Back_trans_Aom = Aom 

     end function BACK_TRANS_AOM

     
     
