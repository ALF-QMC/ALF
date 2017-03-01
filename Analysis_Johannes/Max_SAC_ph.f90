     Program Test

       Use MaxEnt_stoch_mod
  
       Implicit Real (KIND=8) (A-G,O-Z)
       Implicit Integer (H-N)
       
       Real (Kind=8), Dimension(:), allocatable :: XQMC, XTAU, Alpha_tot, om_bf, alp_bf
       Real (Kind=8), Dimension(:,:), allocatable :: XCOV, Xn
       Real (Kind=8) :: X_moments(2), Xerr_moments(2)
       Real (Kind=8), External :: XKER, Back_trans_Aom
       Character (Len=64) :: command


       Open(unit=50,File='Info',Status="unknown")
       
       open(unit=30,file='paramSAC_ph',status='old',action='read', iostat=io_error) 
       if (io_error.eq.0) then
          write(50,*) 'Reading in from paramSAC_ph'
          read(30,*)  Ngamma, OM_st, OM_en, Ndis,  NBins, NSweeps, Nwarm
          read(30,*)  N_alpha, alpha_st, R
          read(30,*) L_cov
       else
          write(6,*) 'No file paramSAC_ph! ' 
          stop
       endif


       open (unit=10,File="g_dat", status="unknown") 
       read(10,*)  ntau
       Allocate ( XCOV(NTAU,NTAU), XQMC(NTAU),XTAU(NTAU) )
       XCOV  = 0.d0
       Do nt = 1,NTAU
          read(10,*)  xtau(nt), xqmc(nt), err 
          xcov(nt,nt) = err*err
       Enddo
       if (L_cov.eq.1) then
          do nt = 1,ntau
             do nt1 = 1,ntau
                read(10,*) xcov(nt,nt1)
             enddo
          enddo
       endif
       close(10)

       dtau = Xtau(2) - Xtau(1)
       Beta = Xtau(Ntau)

       pi = acos(-1.d0)
       xmom1 = pi * xqmc(ntau) 
       Allocate (Alpha_tot(N_alpha) )
       do nt = 1,N_alpha
          alpha_tot(nt) = alpha_st*(R**(nt-1))
       enddo
       write(50,*) 'First Moment, Beta  ', Xmom1, Beta
       close(50)

       
       Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER, Back_Trans_Aom, Beta, &
            &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm)

       Command = "rm dump*"
       Call System (Command)
       Command = "ls"
       Call System (Command)

       open (unit=10,File="g_dat", status="unknown") 
       read(10,*)  ntau
       XCOV  = 0.d0
       Do nt = 1,NTAU
          read(10,*)  xtau(nt), xqmc(nt), err 
          xcov(nt,nt) = err*err
       Enddo
       close(10)
       pi = acos(-1.d0)
       xmom1 = pi * xqmc(ntau) 

       Open (Unit = 10,File="Best_fit", Status ="unknown") 
       Allocate (om_bf(Ngamma), alp_bf(Ngamma) )
       DO i = 1, Ngamma
          read(10,*)  om_bf(i), alp_bf(i)
       Enddo
       close(10) 

       Open (Unit = 11,File="Data_out", Status ="unknown") 
       do nt = 1,Ntau
          X = 0.d0
          tau = xtau(nt)
          do i = 1,Ngamma
             X = X + alp_bf(i)*XKER(tau,om_bf(i), beta)
          enddo
          Write(11,"(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)")  xtau(nt), xqmc(nt),  sqrt(xcov(nt,nt)), xmom1*X
       enddo
       close(11)

     end Program Test
     


     Real (Kind=8) function XKER(tau,om, beta)

       Implicit None
       real (Kind=8) :: tau, om, pi, beta

       pi = 3.1415927

       XKER = (exp(-tau*om) + exp(-( beta - tau )*om ) )/( pi*(1.d0 + exp( - beta * om ) ) )

     end function XKER

     Real (Kind=8) function Back_trans_Aom(Aom, om, beta)

       Implicit None
       real (Kind=8) ::  Aom, om, beta

       !Back_trans_Aom = Aom/(1.d0 + exp(-beta*om) )   ! This gives S(q,om) = chi(q,om)/(1 - e^(-beta om))
       Back_trans_Aom = Aom

     end function BACK_TRANS_AOM

     
