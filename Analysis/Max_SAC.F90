!  Copyright (C) 2016-2020 The ALF project
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
       Use MaxEnt_stoch_mod

       use iso_fortran_env, only: output_unit, error_unit
       Implicit  None

       Interface
          Subroutine  Rescale ( XCOV, XQMC,XTAU, Ntau_st, Ntau_en, Tolerance,  Ntau)
            Implicit none
            Real (Kind=Kind(0.d0)), INTENT(INOUT),allocatable ::  XCOV(:,:), XQMC(:), XTAU(:)
            Real (Kind=Kind(0.d0)), INTENT (IN) :: Tolerance
            Integer,  INTENT(IN)    ::  Ntau_st, Ntau_en
            Integer,  INTENT(INOUT) ::  Ntau
          end Subroutine Rescale
       end Interface

       Real (Kind=Kind(0.d0)), Dimension(:)  , allocatable :: XQMC, XQMC_st, XTAU, Xtau_st, &
            &                                                 Alpha_tot, om_bf, alp_bf, xom, A
       Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: XCOV, XCOV_st
       Real (Kind=Kind(0.d0))                              :: X_moments(2), Xerr_moments(2)
       Real (Kind=Kind(0.d0)), External                    :: XKER_ph, Back_trans_ph, XKER_pp, Back_trans_pp, &
            &                                                 XKER_p, Back_trans_p, XKER_T0, Back_trans_T0
       Character (Len=64)                                  :: command, File1, File2
       Complex (Kind=Kind(0.d0))                           :: Z


       Integer                :: Ngamma, Ndis,  NBins, NSweeps, Nwarm, N_alpha, N_cov
       Integer                :: N_skip, N_rebin, N_Back, N_auto, Norb
       Real (Kind=Kind(0.d0)) :: OM_st, OM_en,  alpha_st, R, Tolerance
       Logical                :: Checkpoint
       Character (Len=2)      :: Channel

       Integer                :: nt, nt1, io_error, n,nw, nwp, ntau, N_alpha_1, i,  nbin_qmc
       Integer                :: ntau_st, ntau_en, ntau_new, Ntau_old
       Real (Kind=Kind(0.d0)) :: dtau, pi, xmom1, x,x1,x2, tau, omp, om, Beta,err, delta, Dom
       Real (Kind=Kind(0.d0)) :: Zero

       NAMELIST /VAR_Max_Stoch/ Ngamma, Ndis,  NBins, NSweeps, Nwarm, N_alpha, &
            &                   OM_st, OM_en,  alpha_st, R,  Checkpoint, Tolerance

       NAMELIST /VAR_errors/    N_skip, N_rebin, N_cov,  N_Back, N_auto

       open(unit=30,file='parameters',status='old',action='read', iostat=io_error)
       if (io_error.eq.0) then
          READ(30,NML=VAR_errors)
          READ(30,NML=VAR_Max_Stoch)
       else
          write(error_unit,*) 'No file parameters '
          error stop 1
       endif
       close(30)


       open (unit=10,File="g_dat", status="unknown")
       read(10,*)  ntau, nbin_qmc, Beta, Norb, Channel
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

       Open(unit=50,File='Info_MaxEnt',Status="unknown")
       write(50,*) 'Channel      :: ', Channel
       If (str_to_upper(Channel) == "PH" )  then
          Write(50,*)  'Om_start is set to zero. PH channel corresponds to symmetric data '
          Om_st = 0.d0
       endif
       Write(50, "('Covariance         :: ',I2)")  N_cov
       Write(50, "('Checkpoint         :: ',L1)")  Checkpoint
       Write(50, "('Om_st, Om_en       :: ',2x,F12.6,2x,F12.6)") Om_st, Om_en
       Write(50, "('Delta Om           :: ',2x,F12.6)")  (Om_en - Om_st)/real(Ndis,kind(0.d0))
       Write(50, "('Bins, Sweeps, Warm :: ',2x,I4,2x,I4,2x,I4)") NBins, NSweeps, Nwarm
       If (N_alpha <= 10 ) then
          Write(error_unit,*) 'Not enough temperatures: N_alpha has to be bigger than 10'
          error stop 1
       Endif
       Write(50, "('N_Alpha, Alpha_st,R:: ',2x,I4,F12.6,2x,F12.6)") N_alpha, alpha_st, R

       Zero = 1.D-10
       pi = acos(-1.d0)
       Ntau_st = 1
       Ntau_en = Ntau
       Select Case (str_to_upper(Channel))
       Case ("PH")
          xmom1 = pi * xqmc(1)
       Case ("PP")
          xmom1 = 2.d0* pi * xqmc(1)
       Case ("P")
          xmom1 =  pi * ( xqmc(1) + xqmc(ntau) )
          !  Remove the tau = beta point from the data since it is  correlated
          !  due to the sum rule with  the tau=0 data point. Also if the tau = 0
          !  data point has no fluctations (due to particle-hole symmetry for instance)
          !  it will be removed.
          Ntau_en = Ntau - 1
          Ntau_st = 1
          if ( xcov(1,1) < zero )  ntau_st = 2
       Case ("T0")
          xmom1 =  pi*xqmc(1)
          Ntau_st = 1
          if ( xcov(1,1) < zero )  ntau_st = 2
       Case default
          Write(error_unit,*) "Channel not yet implemented"
          error stop 1
       end Select
       Ntau_old = Ntau
       Call Rescale ( XCOV, XQMC,XTAU, Ntau_st, Ntau_en, Tolerance, NTAU)
       Write(50,"('Data has been rescaled from Ntau  ',  I4,' to ', I4)")  NTAU_old, Ntau
       If ( Ntau <= 4 ) then
          write(error_unit,*) 'Not enough data!'
          error stop 1
       Endif
       If (  nbin_qmc > 2*Ntau .and. N_cov == 0  )   Write(50,*) 'Consider using the covariance. You seem to have enough bins'
       If (  nbin_qmc < 2*Ntau .and. N_cov == 1  )   Write(50,*) 'You do not seem to have enough bins for a reliable estimate of the covariance '


       ! Store
       Allocate ( XCOV_st(NTAU,NTAU), XQMC_st(NTAU),XTAU_st(NTAU) )
       XCOV_st = XCOV
       XQMC_st = XQMC
       XTAU_st = XTAU


       Allocate (Alpha_tot(N_alpha) )
       do nt = 1,N_alpha
          alpha_tot(nt) = alpha_st*(R**(nt-1))
       enddo
       write(50,"('First Moment, Beta  ',2x,F12.6,2x,F12.6)")  Xmom1, Beta

       Select Case (str_to_upper(Channel))
       Case ("PH")
          If (N_cov == 1 ) then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_ph, Back_Trans_ph, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, N_Cov)
          else
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_ph, Back_Trans_ph, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm)
          endif
          ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
       Case ("PP")
          If (N_cov == 1 ) then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_pp, Back_Trans_pp, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, N_Cov)
          else
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_pp, Back_Trans_pp, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm)
          endif
          ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
       Case ("P")
          If (N_cov == 1 ) then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_p, Back_Trans_p, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, N_Cov)
          else
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_p, Back_Trans_p, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm)
          endif
          ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
       Case ("T0")
          If (N_cov == 1 ) then
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_T0, Back_Trans_T0, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm, N_Cov)
          else
             Call MaxEnt_stoch(XQMC, Xtau, Xcov, Xmom1, XKER_T0, Back_Trans_T0, Beta, &
                  &            Alpha_tot, Ngamma, OM_ST, OM_EN, Ndis, Nsweeps, NBins, NWarm)
          endif
          ! Beware: Xqmc and cov are modified in the MaxEnt_stoch call.
       Case default
          Write(error_unit,*) "Channel not yet implemented"
          error stop 1
       end Select

       If ( .not.  Checkpoint ) then
          Command = "rm dump*"
          Call EXECUTE_COMMAND_LINE(Command)
          Command = "ls"
          Call EXECUTE_COMMAND_LINE(Command)
       endif
       Open (Unit=10,File="energies",status="unknown")

       Do n = 1,N_alpha
          Read(10,*) X,X1,X2
       enddo
       Write(50,"('Chisq at lowest temperature: 1/T, Chi^2, Error',2x,F12.6,2x,F12.6,2x,F12.6)") X,X1,X2
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
          Case ("P")
             do i = 1,Ngamma
                X = X + alp_bf(i)*Xker_p(tau,om_bf(i), beta)
             enddo
          Case ("T0")
             do i = 1,Ngamma
                X = X + alp_bf(i)*Xker_T0(tau,om_bf(i), beta)
             enddo
          Case default
             Write(error_unit,*) "Channel not yet implemented"
             error stop 1
          end Select
          Write(11,"(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)")  xtau_st(nt), xqmc_st(nt),  sqrt(xcov_st(nt,nt)), xmom1*X
       enddo
       close(11)

       N_alpha_1 = N_alpha - 10
       File1 ="Aom_ps"
       file2 =File_i(File1,N_alpha_1)
       Open(Unit=66,File=file2,status="unknown")
       Allocate (xom(Ndis), A(Ndis))
       do nw = 1,Ndis
          read(66,*) xom(nw), A(nw), x, x1, x2
       enddo
       close(66)
       Dom = xom(2) - xom(1)

       ! Compute the real frequency Green function.
       delta = Dom
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
          write(43,"('X',2x,F14.7,2x,F16.8,2x,F16.8)") xom(nw), dble(Z), -Aimag(Z)/pi
       enddo
       close(43)


     end Program MaxEnt_Wrapper



     Real (Kind=Kind(0.d0)) function XKER_ph(tau,om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) :: tau, om, pi, beta

       pi = 3.1415927

       XKER_ph = (exp(-tau*om) + exp(-( beta - tau )*om ) )/( pi*(1.d0 + exp( - beta * om ) ) )

     end function XKER_ph

     Real (Kind=Kind(0.d0)) function XKER_pp(tau,om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) :: tau, om, pi, beta

       pi = 3.1415927

       XKER_pp = exp(-tau*om) / ( pi*(1.d0 + exp( - beta * om ) ) )

     end function XKER_pp

     Real (Kind=Kind(0.d0)) function XKER_p(tau,om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) :: tau, om, pi, beta

       pi = 3.1415927

       XKER_p  = exp(-tau*om) / ( pi*(1.d0 + exp( - beta * om ) ) )

     end function XKER_p

     Real (Kind=Kind(0.d0)) function XKER_T0(tau,om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) :: tau, om, pi, beta

       pi = 3.1415927

       XKER_T0  = exp(-tau*om) / pi

     end function XKER_T0


     Real (Kind=Kind(0.d0)) function Back_trans_ph(Aom, om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) ::  Aom, om, beta

       Back_trans_ph = Aom/(1.d0 + exp(-beta*om) )
       ! This gives S(q,om) = chi(q,om)/(1 - e^(-beta om))

     end function BACK_TRANS_PH

     Real (Kind=Kind(0.d0)) function Back_trans_pp(Aom, om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) ::  Aom, om, beta
       real (Kind=Kind(0.d0)) :: Zero

       Zero = 1.D-8

       if ( abs(om) < zero ) then
          Back_trans_pp = beta * Aom/2.d0
       else
          Back_trans_pp = Aom * (1.d0 - exp(-beta*om) ) / (om *( 1.d0 + exp(-beta*om) ) )
       endif
       ! This gives  = chi(q,om)/omega

     end function BACK_TRANS_PP

     Real (Kind=Kind(0.d0)) function Back_trans_p(Aom, om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) ::  Aom, om, beta

       Back_trans_p =  Aom

     end function BACK_TRANS_P

     Real (Kind=Kind(0.d0)) function Back_trans_T0(Aom, om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) ::  Aom, om, beta

       Back_trans_T0 =  Aom

     end function BACK_TRANS_T0


     Subroutine  Rescale ( XCOV, XQMC,XTAU, Ntau_st, Ntau_en, Tolerance,  Ntau)

       Implicit none

       Real (Kind=Kind(0.d0)), INTENT(INOUT), allocatable ::  XCOV(:,:), XQMC(:), XTAU(:)
       Real (Kind=Kind(0.d0)), INTENT (IN) :: Tolerance

       Integer,  INTENT(IN)    ::  Ntau_st, Ntau_en
       Integer,  INTENT(INOUT) ::  Ntau

       !Local
       Integer :: Ntau_new, nt, nt1
       Integer, allocatable :: List(:)
       Real (Kind=Kind(0.d0)), dimension(:,:), allocatable  ::  XCOV_st
       Real (Kind=Kind(0.d0)), dimension(:)  , allocatable  ::  XQMC_st, XTAU_st



       ! Count the number of elements
       ntau_new = 0
       Do nt = ntau_st,ntau_en
          if ( sqrt(xcov(nt,nt))/xqmc(nt) < Tolerance .and. xqmc(nt) > 0.d0  ) then
             ntau_new  = ntau_new + 1
          endif
       Enddo
       Allocate ( XCOV_st(NTAU_new,NTAU_new), XQMC_st(NTAU_new),XTAU_st(NTAU_new), List(NTAU_new) )
       ntau_new = 0
       Do nt = ntau_st,ntau_en
          if ( sqrt(xcov(nt,nt))/xqmc(nt) < Tolerance .and. xqmc(nt) > 0.d0 ) then
             ntau_new  = ntau_new + 1
             List(ntau_new) = nt
          endif
       Enddo
       do nt = 1,ntau_new
          XQMC_st(nt) = XQMC( List(nt) )
          XTAU_st(nt) = XTAU( List(nt) )
       enddo
       do nt = 1,ntau_new
          do nt1 = 1,ntau_new
             Xcov_st(nt,nt1) = Xcov(List(nt), List(nt1) )
          enddo
       enddo
       NTAU = NTAU_New
       Deallocate (XCOV, XQMC,XTAU )
       allocate   (XCOV(NTAU,NTAU), XQMC(NTAU),XTAU(NTAU) )
       XCOV = XCOV_st
       XQMC = XQMC_st
       XTAU = XTAU_st
       Deallocate (XCOV_st, XQMC_st,XTAU_st, List )




       end Subroutine Rescale
