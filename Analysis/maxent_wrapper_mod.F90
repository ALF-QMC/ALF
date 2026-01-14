!  Copyright (C) 2016-2025 The ALF project
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


module MaxEnt_Wrapper_mod
   implicit none
   
contains
     Real (Kind=Kind(0.d0)) function XKER_ph(tau,om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) :: tau, om, pi, beta

       pi = 3.1415927

       XKER_ph = (exp(-tau*om) + exp(-( beta - tau )*om ) )/( pi*(1.d0 + exp( - beta * om ) ) )

     end function XKER_ph
     
     Real (Kind=Kind(0.d0)) function XKER_ph_c(tau,om, beta)
        ! Kernal for A_c(om), same as XKER_ph
       Implicit None
       real (Kind=Kind(0.d0)) :: tau, om, pi, beta

       pi = 3.1415927

       XKER_ph_c = (exp(-tau*om) + exp(-( beta - tau )*om ) )/( pi*(1.d0 + exp( - beta * om ) ) )

     end function XKER_ph_c

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

     Real (Kind=Kind(0.d0)) function F(om, beta)
      Implicit None
      real (Kind=Kind(0.d0)) ::  om, beta
      F = 1.d0 
     end function F

     Real (Kind=Kind(0.d0)) function F_QFI_ph(om, beta)
      Implicit None
      real (Kind=Kind(0.d0)) ::  om, beta
      real (Kind=Kind(0.d0)) :: pi
      pi = 3.1415927
      F_QFI_ph = (4.d0/pi) * ( (exp(beta*om) - 1.d0)/( exp(beta*om) + 1.d0 ) )**2

     end function F_QFI_ph
     
     Real (Kind=Kind(0.d0)) function F_QFI_ph_c(om, beta)
      ! will improve
      Implicit None
      real (Kind=Kind(0.d0)) ::  om, beta
      real (Kind=Kind(0.d0)) :: pi
      pi = 3.1415927
      F_QFI_ph_c = (4.d0/pi) * ( (exp(beta*om) - 1.d0)/( exp(beta*om) + 1.d0 ) )**2

     end function F_QFI_ph_c

     Real (Kind=Kind(0.d0)) function Back_trans_ph(Aom, om, beta )

       Implicit None
       real (Kind=Kind(0.d0)), intent(in) ::  Aom, beta, om

       Back_trans_ph = Aom/(1.d0 + exp(-beta*om) )
       ! This gives S(q,om) = chi(q,om)/(1 - e^(-beta om))

     end function BACK_TRANS_PH
     
     Real (Kind=Kind(0.d0)) function Back_trans_ph_c(Aom, om, beta)
       Implicit None
       real (Kind=Kind(0.d0)), intent(in) ::  Aom, beta, om
       real (Kind=Kind(0.d0)) :: Zero
       ! same as Back_trans_pp, since Back_trans_pp gives = chi(q,om)/omega
       !                                                  = A(q,om)*tanh(beta om/2)/om
       Zero = 1.D-8
       if ( abs(om) < zero ) then
          Back_trans_ph_c = beta * Aom/2.d0
       else
          Back_trans_ph_c = Aom * (1.d0 - exp(-beta*om) ) / (om *( 1.d0 + exp(-beta*om) ) )
          ! This gives sigma'(q,om) = A_c(q,om)*(1 - e^(-beta om))/(1 + e^(-beta om))/om
          !                         = A_c(q,om)*tanh(beta om/2)/om
       endif
       If (Back_trans_ph_c   < 0.d0)  then 
         Write(6,*)  Aom,om,beta
       endif

     end function BACK_TRANS_PH_C

     Real (Kind=Kind(0.d0)) function Back_trans_pp(Aom, om, beta)

       Implicit None
       real (Kind=Kind(0.d0)), intent(in) ::  Aom, beta, om
       real (Kind=Kind(0.d0)) :: Zero

       Zero = 1.D-8
       if ( abs(om) < Zero ) then
          Back_trans_pp = beta * Aom/2.d0
       else
          Back_trans_pp = Aom * (1.d0 - exp(-beta*om) ) / (om *( 1.d0 + exp(-beta*om) ) )
       endif
       If (Back_trans_pp   < 0.d0)  then 
         Write(6,*)  Aom,om,beta
       endif
       !Back_trans_pp = Aom  
       ! This gives  = chi(q,om)/omega

     end function BACK_TRANS_PP

     Real (Kind=Kind(0.d0)) function XKER_p_ph(tau,om, beta)

       Implicit None
       real (Kind=Kind(0.d0)) :: tau, om, pi, beta

       pi = 3.1415927
       XKER_p_ph  =  (exp(-tau*om)  + exp(-(beta-tau)*om)) / (pi*(1.d0 + exp( -beta * om )) )

     end function XKER_p_ph


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

    Subroutine  Set_Ker_classic(Xker, Xker_classic, Om_ST,  Om_en, beta,xtau)
       Implicit  none  

       Real (Kind=Kind(0.d0)), External :: Xker
       Real (Kind=Kind(0.d0)), INTENT(INOUT), allocatable ::  Xker_classic(:,:)
       Real (Kind=Kind(0.d0)), INTENT(IN ), allocatable ::  Xtau(:)
       Real (Kind=Kind(0.d0)), INTENT(IN) :: Om_ST, Om_en, beta

       Integer ::  nw, Ndis, ntau, nt 
       Real (Kind=Kind(0.d0))  :: Om, Dom

       Ntau =  size(Xker_classic,1)
       Ndis =  size(Xker_classic,2)
       DOM  =  (OM_En -  OM_St)/dble(Ndis)
       do  nw  = 1,Ndis
          om =  OM_St +  dble(nw)*dom
          Do nt  = 1,Ntau
               Xker_classic(nt,nw)  =  XKER(xtau(nt),om, beta)
          Enddo
       enddo
    end Subroutine Set_Ker_classic

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Sets the default model.
!--------------------------------------------------------------------
   Subroutine Set_default(Default,beta,Channel, OM_st, Om_en, xmom1,Default_model_exists,Stochastic)

       use runtime_error_mod
       use iso_fortran_env, only: output_unit, error_unit
       use Files_mod
       
       Implicit none

       Real (Kind=Kind(0.d0)), INTENT(INOUT), allocatable ::  Default(:)
       Real (Kind=Kind(0.d0)), INTENT(IN) ::  beta, xmom1,  Om_st,  Om_en
       Character (Len=*), INTENT(IN)      :: Channel
       Logical,  INTENT(IN)               :: Default_model_exists, Stochastic
       Integer :: Ndis, Nw
       Real (Kind = Kind(0.d0)) ::   Dom, X, Om,  Zero = 1.D-8

       Ndis = size(Default,1)
       Dom = (OM_en - Om_st)/dble(Ndis)
       Select Case (str_to_upper(Channel))
       case("P", "P_PH")
         If (.not. Default_model_exists ) Default = Xmom1/(Om_en - Om_st)
         Default = Default*Dom
       case("PH")
         If (.not. Default_model_exists ) Default = 1.d0/(Om_en - Om_st) ! Flat  default   
         !Compute   sum rule  for  A(om)
         X  = 0.d0
         Do  nw = 1, Ndis
             Om = Om_st + dble(nw)*dom
             Default(nw)  = (1.d0 + exp(-beta*om)) * Default(nw)
             X = X + Default(nw) 
         enddo
         X = X*dom
         Default =  Default*Xmom1/X
         Default =  Default*dom
       case("PH_C")
       ! the Back transformation is same as that in PP channel, 
         If (.not. Default_model_exists )  then 
            Default = 1.d0/(Om_en - Om_st) ! Flat  default   
            Default =  Default*Xmom1
            Default =  Default*dom
         else
            !Compute   sum rule  for  A_c(om)
            X  = 0.d0
            Do nw = 1, Ndis
            ! Default(om) : sigma'(om) -> A_c(om)
            ! See Back_trans_ph_c/Back_trans_pp
            ! Default(om) = (1 - exp(-beta*om))/(1 + exp(-beta*om))*A(om)/om 
               Om = Om_st + dble(nw)*dom
               if ( abs(om) < zero ) then
                  Default(nw) = Default(nw)*2.d0/ beta 
               else
                  Default(nw) = Default(nw) * (om *( 1.d0 + exp(-beta*om) ) )/ (1.d0 - exp(-beta*om) ) 
               endif
               !Default(nw)  = (1.d0 + exp(-beta*om)) * Default(nw)
               X = X + Default(nw) 
            enddo
            X = X*dom
            Default =  Default*Xmom1/X
            Default =  Default*dom
         endif
       case("T0")
         If (.not. Default_model_exists ) Default = Xmom1/(Om_en - Om_st)
         Default = Default*Dom
       case("PP")
         If (.not. Default_model_exists )  then 
            Default = 1.d0/(Om_en - Om_st) ! Flat  default   
            Default =  Default*Xmom1
            Default =  Default*dom
         else
            !Compute   sum rule  for  A(om)
            X  = 0.d0
            Do nw = 1, Ndis
               Om = Om_st + dble(nw)*dom
               if ( abs(om) < zero ) then
                  Default(nw) = Default(nw)*2.d0/ beta 
               else
                  Default(nw) = Default(nw) * (om *( 1.d0 + exp(-beta*om) ) )/ (1.d0 - exp(-beta*om) ) 
               endif
               !Default(nw)  = (1.d0 + exp(-beta*om)) * Default(nw)
               X = X + Default(nw) 
            enddo
            X = X*dom
            Default =  Default*Xmom1/X
            Default =  Default*dom
         endif
       case  default
         Write(error_unit,*) "Channel '" // Channel // "' for  default model not yet implemented"
         CALL Terminate_on_error(ERROR_MAXENT,__FILE__,__LINE__)
       end Select
       Open (Unit=10,File="Default_used", status="Unknown")
       X = 0.d0
       Do  nw = 1, Ndis
          Om = Om_st + dble(nw)*Dom
          X = X + Default(nw)
          Write(10,"(F14.7,2x,F14.7)") Om, Default(nw)/dom
       enddo
       Write(10,'("# Testing  sum rule for  default : ", F14.7,2x,F14.7)' )  X, Xmom1
       close(10)
       If  (Stochastic) Default = Default/dom

   end subroutine Set_default
   
end module MaxEnt_Wrapper_mod
