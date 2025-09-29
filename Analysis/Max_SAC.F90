!  Copyright (C) 2016-2024 The ALF project
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


Program MaxEnt_Wrapper
!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> General wrapper for maxent.  Handles particle, partile-hole, particle-particle,
!> channels as well as zero temperature.  See documentation for details.
!>
!
!--------------------------------------------------------------------
       Use runtime_error_mod
       Use MaxEnt_stoch_mod
       Use MaxEnt_mod
       use MaxEnt_Wrapper_mod
       use iso_fortran_env, only: output_unit, error_unit
#ifdef _OPENMP
       use check_omp_num_threads_mod
#endif

       Implicit  None

       Real (Kind=Kind(0.d0)), Dimension(:)  , allocatable :: XQMC, XQMC_st, XTAU, Xtau_st, &
            &                                                 Alpha_tot, om_bf, alp_bf, xom, A
       Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: XCOV, XCOV_st
       Real (Kind=Kind(0.d0))                              :: X_moments(2), Xerr_moments(2), ChiSq
       Character (Len=64)                                  :: command, File1, File2
       Complex (Kind=Kind(0.d0))                           :: Z
       Logical                                             :: Test =.false.
       

       Integer                :: Ngamma, Ndis,  NBins, NSweeps, Nwarm, N_alpha, N_cov
       Integer                :: N_skip, N_rebin, N_Back, N_auto, Norb
       Real (Kind=Kind(0.d0)) :: OM_st, OM_en,  alpha_st, R, Tolerance
       Logical                :: Checkpoint,  Stochastic, Default_model_exists, Particle_channel_PH
       Character (Len=:), allocatable :: Channel
       Character (Len=1)      :: Char, Char1
       Character (len=64)     :: str_temp
       ! Space  for classic MaxEnt
       Real (Kind=Kind(0.d0)), allocatable ::  Xker_classic(:,:),  A_classic(:),  Default(:)

       Integer                :: nt, nt1, io_error, n,nw, nwp, ntau, N_alpha_1, i,  nbin_qmc
       Integer                :: ntau_st, ntau_en, ntau_new, Ntau_old
       Real (Kind=Kind(0.d0)) :: dtau, pi, xmom1, x,x1,x2, tau, omp, om, Beta,err, delta, Dom
       Real (Kind=Kind(0.d0)) :: Zero, Alpha_classic_st=100000.d0
       Integer ::  N_BZ_Zones     =  1 
       Logical ::  Extended_Zone = .false.

       NAMELIST /VAR_Max_Stoch/ Ngamma, Ndis,  NBins, NSweeps, Nwarm, N_alpha, &
            &                   OM_st, OM_en,  alpha_st, R,  Checkpoint, Tolerance, &
            &                   Stochastic

       NAMELIST /VAR_errors/    N_skip, N_rebin, N_cov,  N_Back, N_auto,  N_BZ_Zones,  Extended_Zone

#ifdef _OPENMP
       call check_omp_num_threads()
#endif

       Stochastic = .true. !  This is  the  default
       open(unit=30,file='parameters',status='old',action='read', iostat=io_error)
       if (io_error.eq.0) then
          READ(30,NML=VAR_errors)
          READ(30,NML=VAR_Max_Stoch)
       else
          write(error_unit,*) 'No file parameters '
          CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
       endif
       close(30)
       
       INQUIRE(FILE="Default", EXIST=Default_model_exists)

       open (unit=10,File="g_dat", status="unknown")
       read(10,*)  ntau, nbin_qmc, Beta, Norb,  str_temp
       Channel  = trim(str_temp)
       Allocate ( XCOV(NTAU,NTAU), XQMC(NTAU),XTAU(NTAU) )
       XCOV  = 0.d0
       Do nt = 1,NTAU
          read(10,*)  xtau(nt), xqmc(nt), err
          xcov(nt,nt) = err*err
       Enddo
       if (N_cov.eq.1) then
          do nt = 1,ntau
             do nt1 = 1,ntau
                read(10,*) xcov(nt,nt1)
             enddo
          enddo
       endif
       close(10)

       dtau = Xtau(2) - Xtau(1)
       11 format(A20, ': ', A)
       12 format(A20, ': ', I10)
       13 format(A20, ': ', *(F14.7))
       14 format(A20, ': ', (L1))
       15 format((E26.17E3))

       If (Stochastic) then 
          Open(unit=50,File='Info_MaxEnt',Status="unknown")
       else
          Open(unit=50,File='Info_MaxEnt_cl',Status="unknown")
       endif
       write(50,11) 'Channel', Channel
       If (str_to_upper(Channel) == "PH" .or. str_to_upper(Channel) == "P_PH" )  then
          Write(50,"(A72)")  'Om_start is set to zero. PH  and  P_PH channels corresponds to symmetric data'
          Om_st = 0.d0
       endif
       Write(50, 12) "Covariance", N_cov
       Write(50, 13) "Om_st", Om_st
       Write(50, 13) "Om_en", Om_en
       Write(50, 13) "Delta Om",  (Om_en - Om_st)/real(Ndis,kind(0.d0))
       Write(50, 14) "Default model exists",Default_model_exists 
       If (Stochastic) then
         Write(50, 14)  'Checkpoint' , Checkpoint
         Write(50, 12) "Bins",    NBins
         Write(50, 12) "Sweeps",  NSweeps
         Write(50, 12) "Warm",    Nwarm
         If (N_alpha <= 10 ) then
            Write(error_unit,*) 'Not enough temperatures: N_alpha has to be bigger than 10'
            CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
         Endif
         Write(50,'(A54)')  "Tempertaure  set:  1/T(n=1..N_alpha) = Alpha_st*R**(n)"
         Write(50,12) "N_Alpha", N_alpha
         Write(50,13) "Alpha_st",alpha_st
         Write(50,13) "R", R
       endif
       Zero= 1.D-10
       pi = acos(-1.d0)
       Ntau_st = 1
       Ntau_en = Ntau
       Select Case (str_to_upper(Channel))
       Case ("PH")
          xmom1 = pi * xqmc(1)
       Case ("PP")
          xmom1 =  pi * (xqmc(1) + xqmc(ntau) )
       Case ("P")
          xmom1 =  pi * ( xqmc(1) + xqmc(ntau) )
          !  Remove the tau = beta point from the data since it is  correlated
          !  due to the sum rule with  the tau=0 data point. Also if the tau = 0
          !  data point has no fluctations (due to particle-hole symmetry for instance)
          !  it will be removed.
          Ntau_en = Ntau - 1
          Ntau_st = 1
          if ( xcov(1,1) < zero )  ntau_st = 2
       Case ("P_PH")
          xmom1 =  pi*xqmc(1)
          if ( xcov(1,1) < zero )  ntau_st = 2
          Particle_channel_PH = .true.
       Case ("T0")
          xmom1 =  pi*xqmc(1)
          Ntau_st = 1
          if ( xcov(1,1) < zero )  ntau_st = 2
       Case default
          Write(error_unit,*) "Channel '" // Channel // "' not yet implemented"
          CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
       end Select
       Ntau_old = Ntau
       Call Rescale ( XCOV, XQMC,XTAU, Ntau_st, Ntau_en, Tolerance, NTAU)
       Write(50,"(A32, I4,A4, I4)") trim('Data has been rescaled from Ntau'), NTAU_old,trim("to"), Ntau
       If ( Ntau <= 4 ) then
          write(error_unit,*) 'Not enough data!'
          CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
       Endif
       If (  nbin_qmc > 2*Ntau .and. N_cov == 0  )   Write(50,"(A72)") 'Consider using the covariance. You seem to have enough bins'
       If (  nbin_qmc < 2*Ntau .and. N_cov == 1  )   Write(50,"(A72)") 'Not enough bins for a reliable estimate of the covariance '

       !Store
       Allocate ( XCOV_st(NTAU,NTAU), XQMC_st(NTAU),XTAU_st(NTAU) )
       XCOV_st = XCOV
       XQMC_st = XQMC
       XTAU_st = XTAU
       Allocate (A_classic(Ndis), Default(Ndis), XKer_classic(size(Xqmc,1),Ndis))
       If (Default_model_exists) then
          Open(Unit=10,file="Default",status="unknown") 
          read (10,*) Char 
          rewind(10) 
          If (Char == "X" )   then
             do nw = 1,Ndis
               read(10,*) Char1,X, X1, Default(nw)
             enddo
          else 
             do nw = 1,Ndis
               read(10,*) X,Default(nw)
             enddo
          endif
          close(10)
       endif
       Call Set_default(Default,beta,Channel, OM_st, Om_en, xmom1,Default_model_exists,Stochastic)
      
       Allocate (Alpha_tot(N_alpha) )
       do nt = 1,N_alpha
          alpha_tot(nt) = alpha_st*(R**(nt-1))
       enddo
       write(50,13) "First Moment",  Xmom1
       write(50,13) "Beta", Beta

       Select Case (str_to_upper(Channel))
       Case ("PH")
          If  (Stochastic)  then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_ph, Back_Trans_ph, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, F_QFI_ph, Default)
             ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
          else
             Call Set_Ker_classic(Xker_ph,Xker_classic,Om_st,Om_en,beta,xtau_st)
             Call  MaxEnt( XQMC, XCOV, A_classic, XKER_classic, Alpha_classic_st, CHISQ ,DEFAULT)
          endif 
       Case ("PP")
          If  (Stochastic) then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_pp, Back_Trans_pp, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, F , Default)
             ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
          else
             Call Set_Ker_classic(Xker_pp,Xker_classic,Om_st,Om_en,beta,xtau_st)
             Call  MaxEnt( XQMC, XCOV, A_classic, XKER_classic, Alpha_classic_st, CHISQ ,DEFAULT)
          endif
       Case ("P")
          If  (Stochastic)  then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_p, Back_Trans_p, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm , F, Default)
             ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
          else  ! Classic
             Call Set_Ker_classic(Xker_p,Xker_classic,Om_st,Om_en,beta,xtau_st)
             Call  MaxEnt( XQMC, XCOV, A_classic, XKER_classic, Alpha_classic_st, CHISQ ,DEFAULT)
          endif  
       Case ("P_PH")
          If  (Stochastic)  then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_p_ph, Back_Trans_p, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm ,F, Default)
          else  ! Classic
             Call Set_Ker_classic(Xker_p_ph,Xker_classic,Om_st,Om_en,beta,xtau_st)
             Call  MaxEnt( XQMC, XCOV, A_classic, XKER_classic, Alpha_classic_st, CHISQ ,DEFAULT)
          endif  
       Case ("T0")
          If (Stochastic)  then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_T0, Back_Trans_T0, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm,F,Default)
             ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
          else
             Call Set_Ker_classic(Xker_T0,Xker_classic,Om_st,Om_en,beta,xtau_st)
             Call  MaxEnt( XQMC, XCOV, A_classic, XKER_classic, Alpha_classic_st, CHISQ ,DEFAULT)
          endif
       Case default
          Write(error_unit,*) "Channel '" // Channel // "' not yet implemented"
          CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
       end Select

       Allocate (xom(Ndis), A(Ndis))
       If  ( Stochastic )   then
          If ( .not.  Checkpoint  ) then
            Command = "rm dump*"
            Call EXECUTE_COMMAND_LINE(Command)
            Command = "ls"
            Call EXECUTE_COMMAND_LINE(Command)
          endif
          Open (Unit=10,File="energies",status="unknown")

          Do n = 1,N_alpha
             Read(10,*) X,X1,X2
          enddo
          Write(50,13) "Best Chisq", X1
          close(50)

          Open (Unit = 10,File="Best_fit", Status ="unknown")
          Allocate (om_bf(Ngamma), alp_bf(Ngamma) )
          DO i = 1, Ngamma
            read(10,*)  om_bf(i), alp_bf(i)
          Enddo
          close(10)

          Open (Unit = 11,File="Data_out", Status ="unknown")
          do nt = 1,Ntau
             X = 0.d0
             tau = xtau_st(nt)
             Select Case (str_to_upper(Channel))
                Case ("PH")
                   do i = 1,Ngamma
                      X = X + alp_bf(i)*Xker_ph(tau,om_bf(i), beta)
                   enddo
                Case ("PP")
                   do i = 1,Ngamma
                      X = X + alp_bf(i)*Xker_pp(tau,om_bf(i), beta)
                   enddo
                Case ("P_PH")
                   do i = 1,Ngamma
                      X = X + alp_bf(i)*Xker_p_ph(tau,om_bf(i), beta)
                   enddo
                Case ("P")
                   do i = 1,Ngamma
                      X = X + alp_bf(i)*Xker_p(tau,om_bf(i), beta)
                   enddo
                Case ("T0")
                   do i = 1,Ngamma
                      X = X + alp_bf(i)*Xker_T0(tau,om_bf(i), beta)
                   enddo
                Case default
                   Write(error_unit,*) "Channel '" // Channel // "' not yet implemented"
                   CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
             end Select
             Write(11,"(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)")  xtau_st(nt), xqmc_st(nt),  sqrt(xcov_st(nt,nt)), xmom1*X
          enddo
          close(11)
          N_alpha_1 = N_alpha - 10
          File1 ="Aom_ps"
          file2 =File_i(File1,N_alpha_1)
          Open(Unit=66,File=file2,status="unknown")
          do nw = 1,Ndis
             read(66,*) xom(nw), A(nw), x, x1, x2
          enddo
          close(66)
          Dom = xom(2) - xom(1)
          Open (Unit=43,File="Green", Status="unknown", action="write")
       else 
          DOM  =  (OM_En -  OM_St)/dble(Ndis)
          A =  A_classic/Dom
          do  nw  = 1,Ndis
             xom(nw) =  OM_St +  dble(nw-1)*dom
          enddo
          Open (Unit = 11,File="Data_out_cl", Status ="unknown")
          Do Nt = 1,Ntau
             X = 0.d0
             Do nw = 1,Ndis
                X = X  +  A_classic(nw)*Xker_classic(nt,nw)
             enddo
             Write(11,"(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)")  xtau_st(nt), xqmc_st(nt),  sqrt(xcov_st(nt,nt)), X
          enddo
          close(11)
          Write(50,13) "Final value of  alpha ", 1.d0/Alpha_classic_st
          Write(50,13) "CHISQ" , CHISQ
          close(50)
          Select Case (str_to_upper(Channel))
             Case ("PH")
                do  nw  = 1,Ndis
                   A(nw) =  Back_trans_ph(A(nw), xom(nw), beta)
                enddo
             Case ("PP")
                do  nw  = 1,Ndis
                   A(nw) =  Back_trans_pp(A(nw), xom(nw), beta)
                enddo
             Case ("P")
                do  nw  = 1,Ndis
                   A(nw) =  Back_trans_p(A(nw), xom(nw), beta)
                enddo
             Case ("P_PH")
                do  nw  = 1,Ndis
                   A(nw) =  Back_trans_p(A(nw), xom(nw), beta)
                enddo
             Case ("T0")
                do  nw  = 1,Ndis
                   A(nw) =  Back_trans_T0(A(nw), xom(nw), beta)
                enddo
             Case default
                Write(error_unit,*) "Channel '" // Channel // "' not yet implemented"
                CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
             end Select
          Open (Unit=43,File="Green_cl", Status="unknown", action="write")
       endif 

       ! Compute the real frequency Green function.
       delta = Dom
       pi = acos(-1.d0)
       x  = 0.d0
       x1 = 0.d0
       x2 = 0.d0
       do nw = 1,Ndis
          Z = cmplx(0.d0,0.d0,Kind(0.d0))
          om = xom(nw)
          do nwp = 1,Ndis
             omp = xom(nwp)
             If  (str_to_upper(Channel) == "P_PH" .and.  omp > 0.00001d0  )  then 
               Z = Z + A(nwp)/cmplx(  om -  omp, delta, kind(0.d0)) &
                   & + A(nwp)/cmplx(  om +  omp, delta, kind(0.d0)) 
            else
               Z = Z + A(nwp)/cmplx( om -  omp, delta, kind(0.d0))
            endif
          enddo
          Z = Z * dom
          x  = x  - Aimag(Z)/pi
          x1 = x1 - om*Aimag(Z)/pi
          x2 = x2 - om*om*Aimag(Z)/pi
          If (Test)   then 
            write(43,"('X',2x,F14.7,2x,F16.8,2x,F16.8,2x,F14.7,2x,F14.7,2x,F14.7)")  & 
                & xom(nw), dble(Z), -Aimag(Z)/pi,  X*dom, x1*dom, x2*dom
          else
            write(43,"('X',2x,F14.7,2x,F16.8,2x,F16.8)")  & 
                & xom(nw), dble(Z), -Aimag(Z)/pi
          endif   
       enddo
       close(43)

     end Program MaxEnt_Wrapper
