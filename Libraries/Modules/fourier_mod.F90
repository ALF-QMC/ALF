Module  Fourier
  Use MaxEnt_mod
  Use MaxEnt_stoch_mod
  Use Matrix
  use iso_fortran_env, only: output_unit, error_unit

  interface Matz_tau
     module procedure Matz_tau_T, Matz_tau_T0, Matz_tau_T0_all, Matz_tau_T_all, Matz_tau_T_all_C, &
          &           Matz_tau_T_cdmft
  end interface

  interface Matz_tau_Bose
     module procedure Matz_tau_T_Bose
  end interface

  interface Tau_Matz
     module procedure Tau_Matz_T, Tau_Matz_T0, Tau_Matz_T0_all, Tau_Matz_T_all,&
          &   tau_matz_spline, tau_matz_spline_all, Tau_Matz_T_stoch, Tau_Matz_T_all_stoch, &
          &   Tau_Matz_T_all_stoch_C, Tau_Matz_T0_stoch , Tau_Matz_T_all_stoch_cdmft
  end interface

  interface Tau_Matz_Bose
     module procedure Tau_Matz_T_Bose
  end interface

  contains

!********
    subroutine  Matz_tau_T(griom, xiom, grtau, xtau, beta)
      implicit none
      ! Given the  G(i omega) calculates G(tau).
      real    (Kind=Kind(0.d0)), Dimension(:) :: xiom,  xtau
      real    (Kind=Kind(0.d0)) :: beta
      complex (Kind=Kind(0.d0)), Dimension(:) :: griom
      real    (Kind=Kind(0.d0)), Dimension(:) :: grtau


      Integer          :: Niom, Ntau, nt, niw, Ntail
      Real    (Kind=Kind(0.d0)) :: a,b, x
      complex (Kind=Kind(0.d0)), Dimension(:), allocatable :: griom1

      Niom = size( xiom ,1 )
      Ntau = size( xtau, 1 )

      allocate ( griom1(Niom) )

      a = 0.d0
      b = 0.d0
      Ntail = 10
      do niw = Niom - Ntail, Niom
         a = a + dble( griom(niw) * cmplx(0.d0,xiom(niw), kind(0.D0) ) )
         b = b + dble( griom(niw) * ( -xiom(niw) * xiom(niw)))
      enddo
      a = a/dble(Ntail + 1)
      b = b/dble(Ntail + 1)
      write(6,*) 'Fourier: a, b ', a, b
      a = 1.d0
      do niw = 1,Niom
         griom1(niw)  = griom(niw) - a/cmplx(0.d0,xiom(niw), kind(0.D0)) &
                &                 +b/ ( xiom(niw)*xiom(niw) )
      enddo

      do  nt = 1,Ntau
         x = 0.d0
         do niw = 1,Niom
            x = x + dble(exp( cmplx(0.d0,-xiom(niw)*xtau(nt), kind(0.D0)) ) *griom1(niw))*2.d0
         enddo
         grtau(nt) = x/beta   - a/2.d0  + b*xtau(nt)/2.d0 - beta * b/4.d0
      enddo

      deallocate ( griom1 )

    end subroutine Matz_tau_T

!------------------
    subroutine Matz_Tau_T_Bose(griom, xiom, grtau, xtau, beta)
      !Working on this
      implicit none
      ! Given the  G(i omega) calculates G(tau) ! for bosons.
      real    (Kind=Kind(0.d0)), Dimension(:) :: xiom,  xtau
      real    (Kind=Kind(0.d0)) :: beta
      real    (Kind=Kind(0.d0)), Dimension(:) :: griom
      real    (Kind=Kind(0.d0)), Dimension(:) :: grtau


      Integer          :: Niom, Ntau, nt, niw
      Real    (Kind=Kind(0.d0)) :: x
      complex (Kind=Kind(0.d0)), Dimension(:), allocatable :: griom1

      Niom = size( xiom ,1 )
      Ntau = size( xtau, 1 )

      allocate ( griom1(Niom) )

      ! No tail really necessary since decays as 1/Om**2
      !a = 0.d0
      !b = 0.d0
      !Ntail = 10
      !do niw = Niom - Ntail, Niom
      !   a = a + dble( griom(niw) * cmplx(0.d0,xiom(niw) ) )
      !   b = b + dble( griom(niw) * ( cmplx(0.d0,xiom(niw)) *cmplx(0.d0,xiom(niw))  ) )
      !enddo
      !a = a/dble(Ntail + 1)
      !b = b/dble(Ntail + 1)
      !!write(6,*) 'Fourier: a, b ', a, b
      !!a = 1.d0
      !do niw = 1,Niom
      !   griom1(niw)   = griom(niw) - cmplx(a,0.d0)/   cmplx( 0.d0,xiom(niw) ) &
      !          &                   - cmplx(b,0.d0)/ ( cmplx( 0.d0,xiom(niw) ) * cmplx(0.d0,xiom(niw)) )
      !enddo

      do  nt = 1,Ntau
         x = 0.d0
         do niw = 1,Niom
            if ( xiom(niw).gt.0.d0) then
               x = x + dble(exp( cmplx(0.d0,-xiom(niw)*xtau(nt), kind(0.D0)) ) * griom(niw))*2.d0
            else
               x = x + dble(exp( cmplx(0.d0,-xiom(niw)*xtau(nt), kind(0.D0)) ) * griom(niw))
            endif
         enddo
         grtau(nt) = x/beta   ! - a/2.d0  + b*xtau(nt)/2.d0 - beta * b/4.d0
      enddo

      deallocate ( griom1 )

    end subroutine Matz_Tau_T_Bose



!********
    subroutine  Matz_tau_T0(griom, xiom, grt0, gr0t,  xtau, beta)
      implicit none
      ! Given the  G(i omega) calculates G(tau).
      real    (Kind=Kind(0.d0)), Dimension(:) :: xiom,  xtau
      real    (Kind=Kind(0.d0)) :: beta
      complex (Kind=Kind(0.d0)), Dimension(:) :: griom
      real    (Kind=Kind(0.d0)), Dimension(:) :: grt0, gr0t



      Integer          :: Niom, Ntau, nt, niw, Ntail
      Real    (Kind=Kind(0.d0)) :: a,b, xp, xm
      complex (Kind=Kind(0.d0)), Dimension(:), allocatable :: griom1

      Niom = size( xiom ,1 )
      Ntau = size( xtau, 1 )

      allocate ( griom1(Niom) )

      a = 0.d0
      b = 0.d0
      Ntail = 10
      do niw = Niom - Ntail, Niom
         a = a + dble( griom(niw) * cmplx(0.d0, xiom(niw), kind(0.D0)) )
         b = b + dble( griom(niw) * ( -xiom(niw) * xiom(niw)) )
      enddo
      a = a/dble(Ntail + 1)
      b = b/dble(Ntail + 1)
      write(6,*) 'Fourier: a, b ', a, b
      a = 1.d0
      do niw = 1,Niom
         griom1(niw)   = griom(niw) - a/  cmplx(0.d0,xiom(niw), kind(0.D0)) &
                &                   + b/( xiom(niw) * xiom(niw) )
      enddo

      do nt = 1,Ntau
         xp = 0.d0
         xm = 0.d0
         do niw = 1,Niom
            xp = xp + dble(exp( cmplx(0.d0,-xiom(niw)*xtau(nt), kind(0.D0)) ) *griom1(niw))*2.d0
            xm = xm + dble(exp( cmplx(0.d0, xiom(niw)*xtau(nt), kind(0.D0)) ) *griom1(niw))*2.d0
         enddo
         grt0(nt) = xp/beta   - a/2.d0  + b*xtau(nt)   /2.d0 - beta * b/4.d0
         gr0t(nt) = xm/beta   + a/2.d0  - b*(-xtau(nt))/2.d0 - beta * b/4.d0
      enddo

      deallocate ( griom1 )

    end subroutine Matz_tau_T0

!**********
    subroutine  Matz_tau_T0_all(g_iom, xiom, g_t0, g_0t,  xtau, beta)
      implicit none
      ! Given the  G(i omega) calculates G(tau).
      real    (Kind=Kind(0.d0)), Dimension(:) :: xiom,  xtau
      real    (Kind=Kind(0.d0)) :: beta
      Type    (Mat_C), Dimension(:,:) :: g_iom
      Type    (Mat_R), Dimension(:,:) :: g_t0, g_0t

      Complex (Kind=Kind(0.d0)), Dimension(:), allocatable :: giom
      Real    (Kind=Kind(0.d0)), Dimension(:), allocatable :: gt0, g0t

      Integer :: Niom, Ntau, nt, niw, Norb, LQ_C
      Integer :: nk, no1, no2
      Real (Kind=Kind(0.D0)) :: Z1, Z2

      Write (6,*) "Size of griom: ", size(g_iom,1), size(g_iom,2)
      Write (6,*) "Size of grt0 : ", size(g_t0,1), size(g_t0,2)
      Write (6,*) "# of orbitals: ", Size(g_t0(1,1)%el,1), Size(g_t0(1,1)%el,2)
      Ntau = size(g_t0,2)
      If ( Ntau.ne.size(g_0t,2) .OR. Ntau.ne.size(xtau,1)  ) Then
         write(6,*) 'Error in ntau! (Fourier, Matz_tau_T0_all) '
      endif
      LQ_c = size(g_t0,1)
      If ( LQ_c.ne.size(g_0t,1) .OR. LQ_C.ne.size(g_iom,1)  ) Then
         write(6,*) 'Error in LQ_C! (Fourier, Matz_tau_T0_all) '
      endif
      Niom = size(g_iom,2)
      If ( Niom.ne.size(xiom,1)  ) Then
         write(6,*) 'Error in Niom! (Fourier, Matz_tau_T0_all) '
      endif

      Norb = Size(g_t0(1,1)%el,1)
      Allocate (giom(Niom), gt0(Ntau), g0t(Ntau) )


      Do nk = 1,LQ_C
         Do no1 = 1,Norb
            Do no2 = 1,Norb
               If (no1.eq.no2)  then
                  do niw = 1,Niom
                     giom(niw) = g_iom(nk,niw)%el(no1,no1)
                  enddo
               elseif (no2.gt.no1) then
                  ! Build Gamma
                  do niw  = 1,Niom
                     giom(niw) = ( g_iom(nk,niw)%el(no1,no1) + &
                          &        g_iom(nk,niw)%el(no2,no2) + &
                          &        g_iom(nk,niw)%el(no1,no2) + &
                          &        g_iom(nk,niw)%el(no2,no1)  ) / 2.D0
                  enddo
               else
                  ! Build eta
                  do niw = 1,Niom
                     giom(niw) = ( g_iom(nk,niw)%el(no1,no1) + &
                          &        g_iom(nk,niw)%el(no2,no2) - &
                          &        g_iom(nk,niw)%el(no1,no2) - &
                          &        g_iom(nk,niw)%el(no2,no1)  ) / 2.D0
                  enddo
               endif
               Call Matz_tau_T0(giom, xiom, gt0, g0t,  xtau, beta)
               do nt = 1,ntau
                  g_0t(nk,nt)%el(no1,no2) = g0t(nt)
                  g_t0(nk,nt)%el(no1,no2) = gt0(nt)
               enddo
            enddo
         enddo
         do nt = 1,ntau
            do no1 = 1,Norb
               do no2 = no1+1, Norb
                  Z1 = g_0t(nk,nt)%el(no1,no2)
                  Z2 = g_0t(nk,nt)%el(no2,no1)
                  g_0t(nk,nt)%el(no1,no2)  = (Z1 - Z2 )/2.D0
                  g_0t(nk,nt)%el(no2,no1)  = (Z1 - Z2 )/2.D0

                  Z1 = g_t0(nk,nt)%el(no1,no2)
                  Z2 = g_t0(nk,nt)%el(no2,no1)
                  g_t0(nk,nt)%el(no1,no2)  = (Z1 - Z2 )/2.D0
                  g_t0(nk,nt)%el(no2,no1)  = (Z1 - Z2 )/2.D0
               enddo
            enddo
         enddo

      enddo

      Deallocate (giom, gt0, g0t )

    end subroutine Matz_tau_T0_all

!**********
    subroutine  Matz_tau_T_all(g_iom, xiom, g_t0,  xtau, beta)
      implicit none
      ! Given the  G(i omega) calculates G(tau).
      real    (Kind=Kind(0.d0)), Dimension(:) :: xiom,  xtau
      real    (Kind=Kind(0.d0)) :: beta
      Type    (Mat_C), Dimension(:,:) :: g_iom
      Type    (Mat_R), Dimension(:,:) :: g_t0

      Complex (Kind=Kind(0.d0)), Dimension(:), allocatable :: giom
      Real    (Kind=Kind(0.d0)), Dimension(:), allocatable :: gt0

      Integer :: Niom, Ntau, nt, niw, Norb, LQ_C
      Integer :: nk, no1, no2
      Real (Kind=Kind(0.D0)) :: Z1, Z2

      Write (6,*) "Size of griom: ", size(g_iom,1), size(g_iom,2)
      Write (6,*) "Size of grt0 : ", size(g_t0,1), size(g_t0,2)
      Write (6,*) "# of orbitals: ", Size(g_t0(1,1)%el,1), Size(g_t0(1,1)%el,2)
      Ntau = size(g_t0,2)
      If ( Ntau.ne.size(xtau,1)  ) Then
         write(6,*) 'Error in ntau! (Fourier, Matz_tau_T0_all) '
      endif
      LQ_c = size(g_t0,1)
      If ( LQ_C.ne.size(g_iom,1)  ) Then
         write(6,*) 'Error in LQ_C! (Fourier, Matz_tau_T0_all) '
      endif
      Niom = size(g_iom,2)
      If ( Niom.ne.size(xiom,1)  ) Then
         write(6,*) 'Error in Niom! (Fourier, Matz_tau_T0_all) '
      endif

      Norb = Size(g_t0(1,1)%el,1)
      Allocate (giom(Niom), gt0(Ntau) )


      Do nk = 1,LQ_C
         Do no1 = 1,Norb
            Do no2 = 1,Norb
               If (no1.eq.no2)  then
                  do niw = 1,Niom
                     giom(niw) = g_iom(nk,niw)%el(no1,no1)
                  enddo
               elseif (no2.gt.no1) then
                  ! Build Gamma
                  do niw = 1,Niom
                     giom(niw) = ( g_iom(nk,niw)%el(no1,no1) + &
                          &        g_iom(nk,niw)%el(no2,no2) + &
                          &        g_iom(nk,niw)%el(no1,no2) + &
                          &        g_iom(nk,niw)%el(no2,no1)  ) / 2.D0
                  enddo
               else
                  ! Build eta
                  do niw = 1,Niom
                     giom(niw) = ( g_iom(nk,niw)%el(no1,no1) + &
                          &        g_iom(nk,niw)%el(no2,no2) - &
                          &        g_iom(nk,niw)%el(no1,no2) - &
                          &        g_iom(nk,niw)%el(no2,no1)  ) / 2.D0
                  enddo
               endif
               Call Matz_tau_T(giom, xiom, gt0, xtau, beta)
               !write(6,*) 'Back in Matz_tau_T_all'
               do nt = 1,ntau
                  g_t0(nk,nt)%el(no1,no2) = gt0(nt)
               enddo
            enddo
         enddo
         do nt = 1,ntau
            do no1 = 1,Norb
               do no2 = no1+1, Norb
                  Z1 = g_t0(nk,nt)%el(no1,no2)
                  Z2 = g_t0(nk,nt)%el(no2,no1)
                  g_t0(nk,nt)%el(no1,no2)  = (Z1 - Z2 )/2.D0
                  g_t0(nk,nt)%el(no2,no1)  = (Z1 - Z2 )/2.D0
               enddo
            enddo
         enddo

      enddo

      Deallocate (giom, gt0 )

    end subroutine Matz_tau_T_all
!**********

!**********
    subroutine  Matz_tau_T_cdmft(g_iom, xiom, g_t0,  xtau, beta)
      implicit none
      ! Given the  G(i omega) calculates G(tau).
      real    (Kind=Kind(0.d0)), Dimension(:) :: xiom,  xtau
      real    (Kind=Kind(0.d0)) :: beta
      Type    (Mat_C), Dimension(:) :: g_iom
      Type    (Mat_R), Dimension(:) :: g_t0

      Complex (Kind=Kind(0.d0)), Dimension(:), allocatable :: giom
      Real    (Kind=Kind(0.d0)), Dimension(:), allocatable :: gt0

      Integer :: Niom, Ntau, nt, niw, Norb
      Integer :: no1, no2
      Real (Kind=Kind(0.D0)) :: Z1, Z2

      Write (6,*) "Size of griom: ", size(g_iom,1)
      Write (6,*) "Size of grt0 : ", size(g_t0,1)
      Write (6,*) "# of orbitals: ", Size(g_t0(1)%el,1)
      Ntau = size(g_t0,1)
      If ( Ntau.ne.size(xtau,1)  ) Then
         write(6,*) 'Error in ntau! (Fourier, Matz_tau_T0_all) '
      endif
      Niom = size(g_iom,1)
      If ( Niom.ne.size(xiom,1)  ) Then
         write(6,*) 'Error in Niom! (Fourier, Matz_tau_T0_all) '
      endif

      Norb = Size(g_t0(1)%el,1)
      Allocate ( giom(Niom), gt0(Ntau) )


      Do no1 = 1,Norb
         Do no2 = 1,Norb
            If (no1.eq.no2)  then
               do niw = 1,Niom
                  giom(niw) = g_iom(niw)%el(no1,no1)
               enddo
            elseif (no2.gt.no1) then
               ! Build Gamma
               do niw = 1,Niom
                  giom(niw) = ( g_iom(niw)%el(no1,no1) + &
                       &        g_iom(niw)%el(no2,no2) + &
                       &        g_iom(niw)%el(no1,no2) + &
                       &        g_iom(niw)%el(no2,no1)  ) / 2.D0
               enddo
            else
               ! Build eta
               do niw = 1,Niom
                  giom(niw) = ( g_iom(niw)%el(no1,no1) + &
                       &        g_iom(niw)%el(no2,no2) - &
                       &        g_iom(niw)%el(no1,no2) - &
                       &        g_iom(niw)%el(no2,no1)  ) / 2.D0
               enddo
            endif
            Call Matz_tau_T(giom, xiom, gt0, xtau, beta)
            !write(6,*) 'Back in Matz_tau_T_all'
            do nt = 1,ntau
               g_t0(nt)%el(no1,no2) = gt0(nt)
            enddo
         enddo
      enddo
      do nt = 1,ntau
         do no1 = 1,Norb
            do no2 = no1+1, Norb
               Z1 = g_t0(nt)%el(no1,no2)
               Z2 = g_t0(nt)%el(no2,no1)
               g_t0(nt)%el(no1,no2)  = (Z1 - Z2 )/2.D0
               g_t0(nt)%el(no2,no1)  = (Z1 - Z2 )/2.D0
            enddo
         enddo
      enddo

      Deallocate (giom, gt0 )

    end subroutine Matz_tau_T_cdmft
!**********


!----------
    subroutine  Matz_tau_T_all_C(g_iom, xiom, g_t0,  xtau, beta)
      implicit none
      ! Given the  G(i omega) calculates G(tau).
      real    (Kind=Kind(0.d0)), Dimension(:) :: xiom,  xtau
      real    (Kind=Kind(0.d0)) :: beta
      Type    (Mat_C), Dimension(:,:) :: g_iom
      Type    (Mat_C), Dimension(:,:) :: g_t0

      Complex (Kind=Kind(0.d0)), Dimension(:), allocatable :: giom
      Real    (Kind=Kind(0.d0)), Dimension(:), allocatable :: gt0

      Integer :: Niom, Ntau, nt, niw, Norb, LQ_C
      Integer :: nk, no1, no2
      Complex (Kind=Kind(0.d0)) :: Z1, Z2

      Write (6,*) "In Matz_tau_T_all_C"
      Write (6,*) "Size of griom: ", size(g_iom,1), size(g_iom,2)
      Write (6,*) "Size of grt0 : ", size(g_t0,1), size(g_t0,2)
      Write (6,*) "# of orbitals: ", Size(g_t0(1,1)%el,1), Size(g_t0(1,1)%el,2)
      Ntau = size(g_t0,2)
      If ( Ntau.ne.size(xtau,1)  ) Then
         write(6,*) 'Error in ntau! (Fourier, Matz_tau_T0_all) '
      endif
      LQ_c = size(g_t0,1)
      If ( LQ_C.ne.size(g_iom,1)  ) Then
         write(6,*) 'Error in LQ_C! (Fourier, Matz_tau_T0_all) '
      endif
      Niom = size(g_iom,2)
      If ( Niom.ne.size(xiom,1)  ) Then
         write(6,*) 'Error in Niom! (Fourier, Matz_tau_T0_all) '
      endif

      Norb = Size(g_t0(1,1)%el,1)
      Allocate (giom(Niom), gt0(Ntau) )


      Do nk = 1,LQ_C
         Do no1 = 1,Norb
            Do no2 = 1,Norb
               If (no1.eq.no2)  then
                  do niw = 1,Niom
                     giom(niw) = g_iom(nk,niw)%el(no1,no1)
                  enddo
               elseif (no2.gt.no1) then
                  ! Build Gamma
                  do niw = 1,Niom
                     giom(niw) = ( g_iom(nk,niw)%el(no1,no1) + &
                          &        g_iom(nk,niw)%el(no2,no2) + &
                          &        g_iom(nk,niw)%el(no1,no2) + &
                          &        g_iom(nk,niw)%el(no2,no1)  ) / 2.D0
                  enddo
               else
                  ! Build eta
                  do niw = 1,Niom
                     giom(niw) = ( g_iom(nk,niw)%el(no1,no1) +                                &
                          &        g_iom(nk,niw)%el(no2,no2) + cmplx(0.d0,1.d0, kind(0.D0))*(    &
                          &        g_iom(nk,niw)%el(no2,no1) - g_iom(nk,niw)%el(no1,no2)))/2.d0
                  enddo
               endif
               Call Matz_tau_T(giom, xiom, gt0, xtau, beta)
               !write(6,*) 'Back in Matz_tau_T_all'
               do nt = 1,ntau
                  g_t0(nk,nt)%el(no1,no2) = cmplx(gt0(nt), 0.d0, kind(0.D0))
               enddo
            enddo
         enddo
         do nt = 1,ntau
            do no1 = 1,Norb
               do no2 = no1+1, Norb
                  Z1 = g_t0(nk,nt)%el(no1,no2) - &
                       &      (g_t0(nk,nt)%el(no1,no1) + g_t0(nk,nt)%el(no2,no2) )/2.d0
                  Z2 = g_t0(nk,nt)%el(no2,no1) - &
                       &      (g_t0(nk,nt)%el(no1,no1) + g_t0(nk,nt)%el(no2,no2) )/2.d0
                  g_t0(nk,nt)%el(no1,no2)  =  Z1 + cmplx(0.0,1.d0, kind(0.D0)) * Z2
                  g_t0(nk,nt)%el(no2,no1)  =  Z1 - cmplx(0.0,1.d0, kind(0.D0)) * Z2
               enddo
            enddo
         enddo

      enddo

      Deallocate (giom, gt0 )

    end subroutine Matz_tau_T_all_C

!------------



    subroutine Tau_Matz_T(griom, xiom, grtau, xtau, beta, A, xom, cov)
      Implicit none

      !Arguments
      Complex (Kind=Kind(0.d0)), Dimension(:)   :: griom
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: xiom, xom, grtau, xtau, A
      Real    (Kind=Kind(0.d0)), Dimension(:,:) :: cov
      Real    (Kind=Kind(0.d0)) :: Beta

      ! Local
      Real (Kind=Kind(0.d0)), Dimension(:), allocatable :: xqmc
      Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: xker

      Integer :: Nom, Ntau, Niom, Niw, Nt, Nw
      Real (Kind=Kind(0.d0)) :: Alpha_st,  Chisq, x

      Nom   = Size(Xom ,1)
      Niom  = Size(Xiom,1)
      Ntau  = Size(Xtau,1)
      Allocate (Xqmc(Ntau), Xker(Ntau,Nom) )
      xqmc = -grtau
      ! Setup data for MaxEnt.
      do nt = 1,ntau
         do nw = 1,Nom
            XKer(nt,nw) = EXP(-xtau(nt)*xom(nw) ) / ( 1.d0 + EXP( -BETA*xom(nw) ) )
         Enddo
      Enddo


      Alpha_st = 1000000.D0
      Chisq    = 0.d0
      Call MaxEnt(XQMC, COV,  A, XKER, ALPHA_ST, CHISQ )

      do niw = 1,niom
         griom(niw) = sum(A/cmplx(-xom, xiom(niw), kind(0.D0)))
      enddo

      open (unit=60,file='data_out', status='unknown', position='append')
      do nt = 1,ntau
         x  = dot_product(xker(nt, :), a)
         write(60,2004) xtau(nt), xqmc(nt), sqrt(cov(nt,nt)), x
      enddo
      close(60)
2004  format(f16.8,2x,f16.8,2x,f16.8,2x,f16.8)

      deallocate (Xqmc, Xker)

    end subroutine Tau_Matz_T

!--------------------

     subroutine Tau_Matz_T_stoch(griom, xiom, grtau, xtau, beta, cov,  &
          &                Alpha_tot, OM_ST, OM_EN, Nsweeps, NBins, NWarm, Xker_func )
      Implicit none

      !Arguments
      Complex (Kind=Kind(0.d0)), Dimension(:)   :: griom
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: xiom, grtau, xtau
      Real    (Kind=Kind(0.d0)), Dimension(:,:) :: cov
      Real    (Kind=Kind(0.d0))                 :: Beta, OM_ST, OM_EN
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: Alpha_tot
      Real    (Kind=Kind(0.d0)), external       :: xker_func
      Integer                          :: Nsweeps, NBins, NWarm

      ! Local
      Real (Kind=Kind(0.d0)), Dimension(:  ), allocatable :: xqmc, A, xom

      Integer ::  Ntau, Niom, Niw, Nt, Nw, Ndis, Ngamma, Lcov
      Real (Kind=Kind(0.d0)) ::  Chisq, x, dom, xmom1

      Ndis   =  5000
      Allocate ( A(ndis),xom(ndis) )
      Niom   = Size(Xiom,1)
      Ntau   = Size(Xtau,1)
      Allocate ( Xqmc(Ntau) )
      Ngamma = Nint(dble(Ntau)*1.5)
      If (Ngamma.lt. 200 ) Ngamma = 200
      Lcov   = 0
      xqmc   = -grtau
      xmom1  = 1.d0
      Call MaxEnt_stoch_fit(xqmc, xtau, cov, Lcov, xker_func, Xmom1, Beta, Alpha_tot,&
           &                Ngamma, OM_ST, OM_EN, Nsweeps, NBins, NWarm, A, &
           &                xom , Chisq )


      dom = xom(2) - xom(1)
      do niw = 1,niom
         griom(niw) = dom * sum(A/cmplx(-xom, xiom(niw), kind(0.D0)))
      enddo

      open (unit=60,file='data_out', status='unknown', position='append')
      do nt = 1,ntau
         x  = 0.d0
         do nw = 1,ndis
            x = x + Xker_func(Xtau(nt),xom(nw), beta)*a(nw)
         enddo
         x = x*dom
         write(60,2004) xtau(nt), xqmc(nt), sqrt(cov(nt,nt)), x
      enddo
      close(60)
2004  format(f16.8,2x,f16.8,2x,f16.8,2x,f16.8)


      deallocate (Xqmc)
      deallocate (A,xom )

    end subroutine Tau_Matz_T_stoch

!--------------------

    subroutine Tau_Matz_T_Bose(griom, xiom, grtau, xtau, beta, A, xom, cov)
      ! Working on this.
      implicit none
      ! Arguments
      Real   ( Kind=Kind(0.d0) ) , Dimension(:)   :: griom
      Real   ( Kind=Kind(0.d0) ) , Dimension(:)   :: xiom, xom, grtau, xtau, A
      Real   ( Kind=Kind(0.d0) ) , Dimension(:,:) :: cov
      Real   ( Kind=Kind(0.d0) )  :: Beta

      ! Local
      Real (Kind=Kind(0.d0)), Dimension(:  ), allocatable :: xqmc
      Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: xker

      Integer :: Nom, Ntau, Niom, Niw, Nt, Nw
      Real (Kind=Kind(0.d0)) :: Alpha_st,  Chisq, x, Zero

      Nom   = Size(Xom ,1)
      Zero = 1.D-10
      Do Nw = 1,Nom
         if ( xom(Nw) .lt. -Zero ) then
            Write(error_unit,*) 'Tau_Matz_T_Bose: Frequencies should be larger than zero'
            error stop 1
         endif
      enddo
      Niom  = Size(Xiom,1)
      Ntau  = Size(Xtau,1)
      Allocate ( Xqmc(Ntau), Xker(Ntau,Nom) )
      ! Setup data for MaxEnt.

      xqmc =  grtau
      do nt = 1,ntau
         !write(6,*) xtau(nt), xqmc(nt), sqrt(cov(nt,nt)), Beta
         do nw = 1,Nom
            if (Xom(nw).gt.Zero) then
               XKer(nt,nw) =  xom(nw)*(EXP(-xtau(nt)*xom(nw))/(1.d0-EXP( -BETA*xom(nw) ) ) - &
                    &                  EXP( xtau(nt)*xom(nw))/(1.d0-EXP(  BETA*xom(nw) ) )    )
            else
               Xker(nt,nw) = 2.d0/Beta
            endif
         Enddo
      Enddo


      Alpha_st = 1000000.0
      Chisq    = 0.d0
      Call MaxEnt(XQMC, COV,  A, XKER, ALPHA_ST, CHISQ )

      do niw = 1,niom
         x = 0.d0
         If ( abs(xiom(niw)).gt.Zero) then
            do nw = 1,nom
               x = x +  2.d0*A(nw) * xom(nw)* xom(nw)/( xom(nw)**2 +  xiom(niw)**2)
            enddo
         else
            do nw = 1,nom
               x = x +  2.d0*A(nw)
            enddo
         endif
         griom(niw) = x
      enddo

      !  A( nw )  =  A(w)*Dom
      !  A(w)   = (1/pi)*chi''(w)/w
      open (unit=60,file='data_out', status='unknown', position='append')
      do nt = 1,ntau
         x  = 0.d0
         do nw = 1,nom
            x = x + xker(nt,nw)*a(nw)
         enddo
         write(60,2004) xtau(nt), xqmc(nt), sqrt(cov(nt,nt)), x
      enddo
      close(60)
2004  format(f16.8,2x,f16.8,2x,f16.8,2x,f16.8)


      deallocate (Xqmc, Xker)

    end subroutine Tau_Matz_T_Bose


!--------------------
!!!!!! To be tested !!!!!
     subroutine Tau_Matz_T0_stoch(griom, xiom, g_t0, cov_t0, g_0t, cov_0t, xtau, beta, Rel_Err,  &
          &                Alpha_tot, OM_ST, OM_EN, Nsweeps, NBins, NWarm, Xker_func )
      Implicit none

      !Arguments
      Complex (Kind=Kind(0.d0)), Dimension(:)   :: griom
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: xiom,g_t0, g_0t, xtau
      Real    (Kind=Kind(0.d0)), Dimension(:,:) :: cov_t0, cov_0t
      Real    (Kind=Kind(0.d0))                 :: Beta, OM_ST, OM_EN, Rel_Err
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: Alpha_tot
      Real    (Kind=Kind(0.d0)), external       :: xker_func
      Integer                          :: Nsweeps, NBins, NWarm

      ! Local
      Real (Kind=Kind(0.d0)), Dimension(:  ), allocatable :: xqmc, A_t0, A_0t, xom

      Integer ::  Ntau, Niom, Niw, Nt, Nw, Ndis, Ngamma, Lcov
      Real (Kind=Kind(0.d0)) ::  Chisq, x, dom, xmom1
      Complex (Kind=Kind(0.d0)) :: z

      Ndis   =  5000
      Allocate ( A_0t(ndis),A_t0(ndis), xom(ndis) )
      Niom   = Size(Xiom,1)
      Ntau   = Size(Xtau,1)
      Allocate ( Xqmc(Ntau) )
      Ngamma = Nint(dble(Ntau)*1.5)
      If (Ngamma.lt. 200 ) Ngamma = 200
      Lcov   = 0
      xqmc   = -g_t0
      xmom1  = xqmc(1)
      Call MaxEnt_stoch_fit(xqmc, xtau, cov_t0, Lcov, xker_func, Xmom1, Beta, Alpha_tot,&
           &                Ngamma, OM_ST, OM_EN, Nsweeps, NBins, NWarm, A_t0, &
           &                xom , Chisq )

      Lcov   = 0
      xqmc   = g_0t
      xmom1  = xqmc(1)
      Call MaxEnt_stoch_fit(xqmc, xtau, cov_0t, Lcov, xker_func, Xmom1, Beta, Alpha_tot,&
           &                Ngamma, OM_ST, OM_EN, Nsweeps, NBins, NWarm, A_0t, &
           &                xom , Chisq )

      dom = xom(2) - xom(1)
      do niw = 1,niom
         z = cmplx(0.d0, 0.d0, kind(0.D0))
         do nw = 1, ndis
            z = z + A_t0(nw)/cmplx(-xom(nw), xiom(niw), kind(0.D0)) + &
                &   A_0t(nw)/cmplx( xom(nw), xiom(niw), kind(0.D0))
         enddo
         griom(niw) = z*dom
      enddo


      open (unit=60,file='data_out', status='unknown', position='append')
      do nt = ntau,1,-1
         x  = 0.d0
         do nw = 1,ndis
            x = x + Xker_func(Xtau(nt),xom(nw), beta)*A_0t(nw)
         enddo
         x = x*dom
         write(60,2004) -xtau(nt), g_0t(nt), sqrt(cov_0t(nt,nt)), x
      enddo
      do nt = 1,ntau
         x  = 0.d0
         do nw = 1,ndis
            x = x + Xker_func(Xtau(nt),xom(nw), beta)*A_t0(nw)
         enddo
         x = x*dom
         write(60,2004) xtau(nt), -g_t0(nt), sqrt(cov_t0(nt,nt)), x
      enddo
      close(60)
2004  format(f16.8,2x,f16.8,2x,f16.8,2x,f16.8)

      deallocate (Xqmc)
      deallocate (A_0t,A_t0,xom )

    end subroutine Tau_Matz_T0_stoch

!--------------------

    subroutine Tau_Matz_T0(griom, xiom, g_t0, cov_t0, g_0t, cov_0t,  xtau, A_0t, A_t0, xom, &
         &                 Rel_Err, Beta)

      Implicit none

      !Arguments
      Complex (Kind=Kind(0.d0)), Dimension(:)   :: griom
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: xiom, g_t0, g_0t, xtau, A_0t, A_t0, xom
      Real    (Kind=Kind(0.d0)), Dimension(:,:) :: cov_t0, cov_0t
      Real    (Kind=Kind(0.d0))  :: Rel_Err
      Real    (Kind=Kind(0.d0)), optional      :: Beta

      ! Local
      Real (Kind=Kind(0.d0)), Dimension(:), allocatable :: xqmc
      Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: xker

      Integer :: Nom, Ntau, Niom, Niw, Nt, Nw
      Real (Kind=Kind(0.d0)) :: Alpha_st,  Chisq, x

      Complex (Kind=Kind(0.d0)) :: z

      Nom   = Size(Xom ,1)
      Niom  = Size(Xiom,1)
      Ntau  = Size(Xtau,1)

      Allocate (Xqmc(Ntau), Xker(Ntau,Nom))
      Write(6,*) ' Calling Max_Ent from T=0 routine. '
      ! t > 0
      ! Setup data for MaxEnt.
      If (Present(Beta)) Then
         do nt = 1,ntau
            do nw = 1,Nom
               XKer(nt,nw) = EXP(-xtau(nt)*xom(nw) )/ ( 1.d0 + EXP( -BETA*xom(nw)) )
            Enddo
         Enddo
      else
         do nt = 1,ntau
            do nw = 1,Nom
               XKer(nt,nw) = EXP(-xtau(nt)*xom(nw) )
            Enddo
         Enddo
      endif
      Alpha_st = 100000.0
      Chisq    = 0.d0
      xqmc = -g_t0
      Open (Unit=13,file='In_MaxEnt_T0', status='unknown', position='append')
      do nt = 1,ntau
         write(13,2001) xtau(nt), xqmc(nt), sqrt(cov_t0(nt,nt))
      enddo
      write(13,*)
      close(13)
      Call MaxEnt(XQMC, COV_t0,  A_t0, XKER, ALPHA_ST, CHISQ, Rel_err )

      Alpha_st = 100000.0
      Chisq    = 0.d0
      xqmc =  g_0t
      Open (Unit=13,file='In_MaxEnt_T0', status='unknown', position='append')
      do nt = 1,ntau
         write(13,2001) xtau(nt), xqmc(nt), sqrt(cov_0t(nt,nt))
      enddo
      write(13,*)
      close(13)
      Call MaxEnt(XQMC, COV_0t,  A_0t, XKER, ALPHA_ST, CHISQ, Rel_err)

!      do niw = 1,niom
!         z = cmplx(0.d0,0.d0)
!         do nw = 1,nom
!            z = z + cmplx(A(nw),0.d0)/cmplx(-xom(nw), xiom(niw))
!         enddo
!         griom(niw) = z
!      enddo


      do niw = 1,niom
         z = cmplx(0.d0, 0.d0, kind(0.D0))
         do nw = 1,nom
            z = z + A_t0(nw)/cmplx(-xom(nw), xiom(niw), kind(0.D0)) + &
                &   A_0t(nw)/cmplx( xom(nw), xiom(niw), kind(0.D0))
         enddo
         griom(niw) = z
      enddo


      open (unit=60,file='data_out', status='unknown', position='append')
      do nt = 1,ntau
         x  = 0.d0
         do nw = 1,nom
            x = x + xker(nt,nw)*a_t0(nw)
         enddo
         write(60,2004) xtau(nt), -g_t0(nt), sqrt(cov_t0(nt,nt)), x
      enddo
      write(60,*)
      do nt = 1,ntau
         x  = 0.d0
         do nw = 1,nom
            x = x + xker(nt,nw)*a_0t(nw)
         enddo
         write(60,2004) xtau(nt), g_0t(nt), sqrt(cov_0t (nt,nt)), x
      enddo
      write(60,*)

      close(60)
2004  format(f16.8,2x,f16.8,2x,f16.8,2x,f16.8)

2001  format(F16.8,2x,F16.8,2x,F16.8)

      deallocate (Xqmc, Xker)



    end subroutine Tau_Matz_T0

!************
    subroutine Tau_Matz_T0_all( g_iom_mat, xiom, g_t0_mat, error_t0_mat, g_0t_mat, error_0t_mat, &
         &                      xtau, xom, Rel_err )

      Implicit none

      !Arguments
      Type    (Mat_C),  Dimension(:,:) :: g_iom_mat
      Type    (Mat_R),  Dimension(:,:) :: g_t0_mat, g_0t_mat
      Type    (Mat_R),  Dimension(:,:) :: error_t0_mat, error_0t_mat
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: xiom,  xtau, xom
      Real    (Kind=Kind(0.d0)) :: Rel_err

      Complex (Kind=Kind(0.d0)), Dimension(:),   allocatable :: giom
      Real    (Kind=Kind(0.d0)), Dimension(:),   allocatable :: gt0, g0t, A0t, At0
      Real    (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: covt0, cov0t

      Integer ::   Nom, Ntau, Niom, Nt, Nw

      !   ******* Input  is | g_11                           , (g_11 +  g_22 + g_12 + g_21)/2.0 |
      !                     |(g_11 +  g_22 - g_12 - g_21)/2.0,  g_22                            |
      ! As a function of tau.
      !
      !   ******* Output is | g_11 , g_12 |
      !                     | g_21 , g_22 |
      ! As a funtion of omega_m

      ! Local
      Integer ::   LQ_c,   Norb
      Integer ::   nk, no1,no2
      Complex (Kind=Kind(0.d0)) :: Zp

      Ntau = size(g_t0_mat,2)
      If ( Ntau.ne.size(g_0t_mat,2) .OR. Ntau.ne.size(xtau,1)  ) Then
         write(6,*) 'Error in ntau! (Fourier, Tau_Matz_T0_all) '
      endif
      LQ_c = size(g_t0_mat,1)
      If ( LQ_c.ne.size(g_0t_mat,1) .OR. LQ_C.ne.size(g_iom_mat,1)  ) Then
         write(6,*) 'Error in LQ_C! (Fourier, Tau_Matz_T0_all) '
      endif
      Niom = size(g_iom_mat,2)
      If ( Niom.ne.size(xiom,1)  ) Then
         write(6,*) 'Error in Niom! (Fourier, Tau_Matz_T0_all) '
      endif

      Norb = Size(g_t0_mat(1,1)%el,1)
      Nom  = Size(xom,1)

      allocate(giom(Niom), gt0(Ntau), g0t(Ntau), A0t(Nom), At0(Nom), &
           &   covt0(Ntau,Ntau), cov0t(Ntau,Ntau)                   )


      do nk = 1,LQ_C
         do no1 = 1,Norb
            do no2 = 1,Norb
               do nt = 1,Ntau
                  gt0(nt) =  g_t0_mat(nk,nt)%el(no1,no2)
                  g0t(nt) =  g_0t_mat(nk,nt)%el(no1,no2)
               enddo
               covt0 = 0.D0; cov0t = 0.D0
               do nt = 1,Ntau
                  covt0(nt,nt) = (error_t0_mat(nk,nt)%el(no1,no2))**2
                  cov0t(nt,nt) = (error_0t_mat(nk,nt)%el(no1,no2))**2
               enddo
               Write(6,* ) ' Nk is : ',  nk
               call Tau_Matz_T0(giom, xiom, gt0, covt0, g0t, cov0t,  xtau,  A0t, At0, xom, Rel_err)
               do nw = 1,Niom
                  g_iom_mat(nk,nw)%el(no1,no2) = giom(nw)
               enddo
            enddo
         enddo
         do no1 = 1,Norb
            do no2 = no1 +  1, Norb
               do nw = 1,Niom
                  Zp =  g_iom_mat(nk,nw)%el(no1,no2) - g_iom_mat(nk,nw)%el(no2,no1)
                  g_iom_mat(nk,nw)%el(no1,no2) = Zp/2.D0
                  g_iom_mat(nk,nw)%el(no2,no1) = Zp/2.D0
               enddo
            enddo
         enddo
      enddo
      deallocate( giom, gt0, g0t, A0t, At0, covt0, cov0t )

    end subroutine Tau_Matz_T0_all


!-----------------
    subroutine Tau_Matz_T_all_stoch_C(g_iom_mat, xiom, g_t0_mat, error_t0_mat,  xtau, Beta,&
         &                            Alpha_tot, OM_ST, OM_EN, Nsweeps, NBins, NWarm, Xker_func )

      Implicit none

      !Arguments
      Type    (Mat_C),  Dimension(:,:) :: g_iom_mat
      Type    (Mat_C),  Dimension(:,:) :: g_t0_mat
      Type    (Mat_R),  Dimension(:,:) :: error_t0_mat
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: xiom,  xtau
      Real    (Kind=Kind(0.d0))                 :: Beta, OM_St, OM_EN
      Real    (Kind=Kind(0.d0)), external       :: Xker_func
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: Alpha_tot
      Integer                          :: Nsweeps, NBins, NWarm

      Complex (Kind=Kind(0.d0)), Dimension(:),   allocatable :: giom
      Real    (Kind=Kind(0.d0)), Dimension(:),   allocatable :: gt0
      Real    (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: covt0

      Integer ::   Ntau, Niom, Nt, Nw

      !   ******* Input  is | g_11                               , (g_11 +  g_22 + g_12 + g_21)/2.0 |
      !                     |(g_11 +  g_22 - i[g_12 - g_21]) /2.0,  g_22                            |
      ! As a function of tau.  Note that input is real since.
      !   With gamma =  (c + d)/sqrt(2)  and eta = (c + i d)/sqrt(2)
      !   ******* Input  is | cc*       , gamma gamma* |
      !                     | eta eta*  ,  dd*         |

      !
      !   ******* Output is | g_11 , g_12 |
      !                     | g_21 , g_22 |
      ! As a funtion of omega_m

      ! Local
      Integer ::   LQ_c,   Norb
      Integer ::   nk, no1,no2
      Complex (Kind=Kind(0.d0)) :: Z1, Z2

      Ntau = size(g_t0_mat,2)
      If (  Ntau.ne.size(xtau,1)  ) Then
         write(6,*) 'Error in ntau! (Fourier, Tau_Matz_T_all) '
      endif
      LQ_c = size(g_t0_mat,1)
      If (  LQ_C.ne.size(g_iom_mat,1)  ) Then
         write(6,*) 'Error in LQ_C! (Fourier, Tau_Matz_T_all) '
      endif
      Niom = size(g_iom_mat,2)
      If ( Niom.ne.size(xiom,1)  ) Then
         write(6,*) 'Error in Niom! (Fourier, Tau_Matz_T_all) '
      endif

      Norb = Size(g_t0_mat(1,1)%el,1)

      allocate(giom(Niom), gt0(Ntau),  covt0(Ntau,Ntau) )


      do nk = 1,LQ_C
         do no1 = 1,Norb
            do no2 = 1,Norb
               do nt = 1,Ntau
                  gt0(nt) =  dble(g_t0_mat(nk,nt)%el(no1,no2))
               enddo
               covt0 = 0.D0
               do nt = 1,Ntau
                  covt0(nt,nt) = (error_t0_mat(nk,nt)%el(no1,no2))**2
               enddo
               Write(6,* ) ' Nk is : ',  nk
               Call Tau_Matz_T_stoch(giom, xiom, gt0, xtau, beta, covt0, &
                    &                Alpha_tot, OM_ST, OM_EN, Nsweeps, NBins, NWarm, Xker_func )
               do nw = 1,Niom
                  g_iom_mat(nk,nw)%el(no1,no2) = giom(nw)
               enddo
            enddo
         enddo
         do no1 = 1,Norb
            do no2 = no1 +  1, Norb
               do nw = 1,Niom
                  Z1 = g_iom_mat(nk,nw)%el(no1,no2)   - &
                       & (g_iom_mat(nk,nw)%el(no1,no1)+g_iom_mat(nk,nw)%el(no2,no2))/2.d0
                  Z2 = g_iom_mat(nk,nw)%el(no2,no1)   - &
                       & (g_iom_mat(nk,nw)%el(no1,no1)+g_iom_mat(nk,nw)%el(no2,no2))/2.d0
                  g_iom_mat(nk,nw)%el(no1,no2) = Z1 + cmplx(0.0,1.d0, kind(0.D0))*Z2
                  g_iom_mat(nk,nw)%el(no2,no1) = Z1 - cmplx(0.0,1.d0, kind(0.D0))*Z2
               enddo
            enddo
         enddo
      enddo
      deallocate( giom, gt0,  covt0)
    end subroutine Tau_Matz_T_all_stoch_C


!-----------------
    subroutine Tau_Matz_T_all_stoch(g_iom_mat, xiom, g_t0_mat, error_t0_mat,  xtau, Beta,&
         &                          Alpha_tot, OM_ST, OM_EN, Nsweeps, NBins, NWarm, Xker_func )

      Implicit none

      !Arguments
      Type    (Mat_C),  Dimension(:,:) :: g_iom_mat
      Type    (Mat_R),  Dimension(:,:) :: g_t0_mat
      Type    (Mat_R),  Dimension(:,:) :: error_t0_mat
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: xiom,  xtau
      Real    (Kind=Kind(0.d0))                 :: Beta, OM_St, OM_EN
      Real    (Kind=Kind(0.d0)), external       :: Xker_func
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: Alpha_tot
      Integer                          :: Nsweeps, NBins, NWarm

      Complex (Kind=Kind(0.d0)), Dimension(:),   allocatable :: giom
      Real    (Kind=Kind(0.d0)), Dimension(:),   allocatable :: gt0
      Real    (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: covt0

      Integer ::   Ntau, Niom, Nt, Nw

      !   ******* Input  is | g_11                           , (g_11 +  g_22 + g_12 + g_21)/2.0 |
      !                     |(g_11 +  g_22 - g_12 - g_21)/2.0,  g_22                            |
      ! As a function of tau.
      !
      !   ******* Output is | g_11 , g_12 |
      !                     | g_21 , g_22 |
      ! As a funtion of omega_m

      ! Local
      Integer ::   LQ_c,   Norb
      Integer ::   nk, no1,no2
      Complex (Kind=Kind(0.d0)) :: Zp

      Ntau = size(g_t0_mat,2)
      If (  Ntau.ne.size(xtau,1)  ) Then
         write(6,*) 'Error in ntau! (Fourier, Tau_Matz_T_all) '
      endif
      LQ_c = size(g_t0_mat,1)
      If (  LQ_C.ne.size(g_iom_mat,1)  ) Then
         write(6,*) 'Error in LQ_C! (Fourier, Tau_Matz_T_all) '
      endif
      Niom = size(g_iom_mat,2)
      If ( Niom.ne.size(xiom,1)  ) Then
         write(6,*) 'Error in Niom! (Fourier, Tau_Matz_T_all) '
      endif

      Norb = Size(g_t0_mat(1,1)%el,1)

      allocate(giom(Niom), gt0(Ntau),  covt0(Ntau,Ntau) )


      do nk = 1,LQ_C
         do no1 = 1,Norb
            do no2 = 1,Norb
               do nt = 1,Ntau
                  gt0(nt) =  g_t0_mat(nk,nt)%el(no1,no2)
               enddo
               covt0 = 0.D0
               do nt = 1,Ntau
                  covt0(nt,nt) = (error_t0_mat(nk,nt)%el(no1,no2))**2
               enddo
               Write(6,* ) ' Nk is : ',  nk
               Call Tau_Matz_T_stoch(giom, xiom, gt0, xtau, beta, covt0, &
                    &                Alpha_tot, OM_ST, OM_EN, Nsweeps, NBins, NWarm, Xker_func )
               do nw = 1,Niom
                  g_iom_mat(nk,nw)%el(no1,no2) = giom(nw)
               enddo
            enddo
         enddo
         do no1 = 1,Norb
            do no2 = no1 +  1, Norb
               do nw = 1,Niom
                  Zp =  g_iom_mat(nk,nw)%el(no1,no2) - g_iom_mat(nk,nw)%el(no2,no1)
                  g_iom_mat(nk,nw)%el(no1,no2) = Zp/2.D0
                  g_iom_mat(nk,nw)%el(no2,no1) = Zp/2.D0
               enddo
            enddo
         enddo
      enddo
      deallocate( giom, gt0,  covt0)
    end subroutine Tau_Matz_T_all_stoch

!-----------------------
    subroutine Tau_Matz_T_all_stoch_cdmft(g_iom_mat, xiom, g_t0_mat, error_t0_mat,  xtau, Beta,&
         &                          Alpha_tot, OM_ST, OM_EN, Nsweeps, NBins, NWarm, Xker_func )

      Implicit none

      !Arguments
      Type    (Mat_C),  Dimension(:) :: g_iom_mat
      Type    (Mat_R),  Dimension(:) :: g_t0_mat
      Type    (Mat_R),  Dimension(:) :: error_t0_mat
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: xiom,  xtau
      Real    (Kind=Kind(0.d0))                 :: Beta, OM_St, OM_EN
      Real    (Kind=Kind(0.d0)), external       :: Xker_func
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: Alpha_tot
      Integer                          :: Nsweeps, NBins, NWarm

      Complex (Kind=Kind(0.d0)), Dimension(:),   allocatable :: giom
      Real    (Kind=Kind(0.d0)), Dimension(:),   allocatable :: gt0
      Real    (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: covt0

      Integer ::   Ntau, Niom, Nt, Nw

      !   ******* Input  is | g_11                           , (g_11 +  g_22 + g_12 + g_21)/2.0 |
      !                     |(g_11 +  g_22 - g_12 - g_21)/2.0,  g_22                            |
      ! As a function of tau.  and generalization thereof for larger matrices.
      !
      !   ******* Output is | g_11 , g_12 |
      !                     | g_21 , g_22 |
      ! As a funtion of omega_m

      ! Local
      Integer ::   Norb
      Integer ::   no1,no2
      Complex (Kind=Kind(0.d0)) :: Zp

      Ntau = size(g_t0_mat,1)
      If (  Ntau.ne.size(xtau,1)  ) Then
         write(6,*) 'Error in ntau! (Fourier, Tau_Matz_T_all) '
      endif
      Niom = size(g_iom_mat,1)
      If ( Niom.ne.size(xiom,1)  ) Then
         write(6,*) 'Error in Niom! (Fourier, Tau_Matz_T_all) '
      endif

      Norb = Size(g_t0_mat(1)%el,1)

      allocate(giom(Niom), gt0(Ntau),  covt0(Ntau,Ntau) )


      do no1 = 1,Norb
         do no2 = 1,Norb

            do nt = 1,Ntau
               gt0(nt) =  g_t0_mat(nt)%el(no1,no2)
            enddo
            covt0 = 0.D0
            do nt = 1,Ntau
               covt0(nt,nt) = (error_t0_mat(nt)%el(no1,no2))**2
            enddo
            Call Tau_Matz_T_stoch(giom, xiom, gt0, xtau, beta, covt0, &
                 &                Alpha_tot, OM_ST, OM_EN, Nsweeps, NBins, NWarm, Xker_func )
            do nw = 1,Niom
               g_iom_mat(nw)%el(no1,no2) = giom(nw)
            enddo
         enddo
      enddo
      do no1 = 1,Norb
         do no2 = no1 +  1, Norb
            do nw = 1,Niom
               Zp =  g_iom_mat(nw)%el(no1,no2) - g_iom_mat(nw)%el(no2,no1)
               g_iom_mat(nw)%el(no1,no2) = Zp/cmplx(2.0,0.0)
               g_iom_mat(nw)%el(no2,no1) = Zp/cmplx(2.0,0.0)
            enddo
         enddo
      enddo
      deallocate( giom, gt0,  covt0)
    end subroutine Tau_Matz_T_all_stoch_cdmft

!-----------------
    subroutine Tau_Matz_T_all(g_iom_mat, xiom, g_t0_mat, error_t0_mat,  xtau, xom, Beta )

      Implicit none

      !Arguments
      Type    (Mat_C),  Dimension(:,:) :: g_iom_mat
      Type    (Mat_R),  Dimension(:,:) :: g_t0_mat
      Type    (Mat_R),  Dimension(:,:) :: error_t0_mat
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: xiom,  xtau, xom
      Real    (Kind=Kind(0.d0)) :: Beta

      Complex (Kind=Kind(0.d0)), Dimension(:),   allocatable :: giom
      Real    (Kind=Kind(0.d0)), Dimension(:),   allocatable :: gt0, At0
      Real    (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: covt0

      Integer ::   Nom, Ntau, Niom, Nt, Nw

      !   ******* Input  is | g_11                           , (g_11 +  g_22 + g_12 + g_21)/2.0 |
      !                     |(g_11 +  g_22 - g_12 - g_21)/2.0,  g_22                            |
      ! As a function of tau.
      !
      !   ******* Output is | g_11 , g_12 |
      !                     | g_21 , g_22 |
      ! As a funtion of omega_m

      ! Local
      Integer ::   LQ_c,   Norb
      Integer ::   nk, no1,no2
      Complex (Kind=Kind(0.d0)) :: Zp

      Ntau = size(g_t0_mat,2)
      If (  Ntau.ne.size(xtau,1)  ) Then
         write(6,*) 'Error in ntau! (Fourier, Tau_Matz_T_all) '
      endif
      LQ_c = size(g_t0_mat,1)
      If (  LQ_C.ne.size(g_iom_mat,1)  ) Then
         write(6,*) 'Error in LQ_C! (Fourier, Tau_Matz_T_all) '
      endif
      Niom = size(g_iom_mat,2)
      If ( Niom.ne.size(xiom,1)  ) Then
         write(6,*) 'Error in Niom! (Fourier, Tau_Matz_T_all) '
      endif

      Norb = Size(g_t0_mat(1,1)%el,1)
      Nom  = Size(xom,1)

      allocate(giom(Niom), gt0(Ntau),  At0(Nom), covt0(Ntau,Ntau) )


      do nk = 1,LQ_C
         do no1 = 1,Norb
            do no2 = 1,Norb
               do nt = 1,Ntau
                  gt0(nt) =  g_t0_mat(nk,nt)%el(no1,no2)
               enddo
               covt0 = 0.D0
               do nt = 1,Ntau
                  covt0(nt,nt) = (error_t0_mat(nk,nt)%el(no1,no2))**2
               enddo
               Write(6,* ) ' Nk is : ',  nk
               Call Tau_Matz_T(giom, xiom, gt0, xtau, beta, At0, xom, covt0)
               do nw = 1,Niom
                  g_iom_mat(nk,nw)%el(no1,no2) = giom(nw)
               enddo
            enddo
         enddo
         do no1 = 1,Norb
            do no2 = no1 +  1, Norb
               do nw = 1,Niom
                  Zp =  g_iom_mat(nk,nw)%el(no1,no2) - g_iom_mat(nk,nw)%el(no2,no1)
                  g_iom_mat(nk,nw)%el(no1,no2) = Zp/cmplx(2.0,0.0)
                  g_iom_mat(nk,nw)%el(no2,no1) = Zp/cmplx(2.0,0.0)
               enddo
            enddo
         enddo
      enddo
      deallocate( giom, gt0, At0, covt0)

    end subroutine Tau_Matz_T_all


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------
  subroutine tau_matz_spline(nspl,griom, xiom, grtau, xtau)
    implicit none

    integer, intent(in) :: nspl
    real(Kind=Kind(0.d0)), dimension(:), intent(in) :: xiom,xtau,grtau
    complex(Kind=Kind(0.d0)), dimension(size(xiom)), intent(out)  :: griom

    integer :: itau,iom,ntau,niom
    real(Kind=Kind(0.d0)) :: dx
    real(Kind=Kind(0.d0)), dimension(:), allocatable :: xtau_spl,grtau_spl

    ntau = size(xtau)
    niom = size(xiom)

    allocate(xtau_spl(0:nspl),grtau_spl(0:nspl))

    dx = xtau(ntau) / dble(nspl)
    do itau = 0,nspl
       xtau_spl(itau) = dx * dble(itau)
    enddo

    call aspline(xtau,grtau,xtau_spl,grtau_spl)

!!$    open(10,file='spline.dat',position='append')
!!$    do itau = 0,nspl
!!$       write(10,*) xtau_spl(itau),grtau_spl(itau)
!!$    enddo
!!$    write(10,*)
!!$    write(10,*)
!!$    close(10)

    griom = (0.d0,0.d0)
    do iom = 1,niom
       do itau = 0,nspl
          griom(iom) = griom(iom) &
               +  exp(cmplx(0.d0,xiom(iom)*xtau_spl(itau), kind(0.D0))) * grtau_spl(itau) * dx
       enddo
    enddo

    deallocate(xtau_spl,grtau_spl)

  end subroutine tau_matz_spline

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------
 subroutine tau_matz_spline_all(nspl, g_iom_mat, xiom, g_t0_mat, xtau)

      Implicit none

      !Arguments
      integer :: nspl
      Type    (Mat_C),  Dimension(:,:) :: g_iom_mat
      Type    (Mat_R),  Dimension(:,:) :: g_t0_mat
      Real    (Kind=Kind(0.d0)), Dimension(:)   :: xiom,  xtau

      Complex (Kind=Kind(0.d0)), Dimension(:),   allocatable :: giom
      Real    (Kind=Kind(0.d0)), Dimension(:),   allocatable :: gt0

      Integer ::   Ntau, Niom, Nt, Nw

      !   ******* Input  is | g_11                           , (g_11 +  g_22 + g_12 + g_21)/2.0 |
      !                     |(g_11 +  g_22 - g_12 - g_21)/2.0,  g_22                            |
      ! As a function of tau.
      !
      !   ******* Output is | g_11 , g_12 |
      !                     | g_21 , g_22 |
      ! As a funtion of omega_m

      ! Local
      Integer ::   LQ_c,   Norb
      Integer ::   nk, no1,no2
      Complex (Kind=Kind(0.d0)) :: Zp

      Ntau = size(g_t0_mat,2)
      If (  Ntau.ne.size(xtau,1)  ) Then
         write(6,*) 'Error in ntau! (Fourier, Tau_Matz_T_all) '
      endif
      LQ_c = size(g_t0_mat,1)
      If (  LQ_C.ne.size(g_iom_mat,1)  ) Then
         write(6,*) 'Error in LQ_C! (Fourier, Tau_Matz_T_all) '
      endif
      Niom = size(g_iom_mat,2)
      If ( Niom.ne.size(xiom,1)  ) Then
         write(6,*) 'Error in Niom! (Fourier, Tau_Matz_T_all) '
      endif

      Norb = Size(g_t0_mat(1,1)%el,1)

      allocate(giom(Niom), gt0(Ntau))


      do nk = 1,LQ_C
         do no1 = 1,Norb
            do no2 = 1,Norb
               do nt = 1,Ntau
                  gt0(nt) =  g_t0_mat(nk,nt)%el(no1,no2)
               enddo
               Write(6,* ) ' Nk is : ',  nk
               Call tau_matz_spline(nspl, giom, xiom, gt0, xtau)
               do nw = 1,Niom
                  g_iom_mat(nk,nw)%el(no1,no2) = giom(nw)
               enddo
            enddo
         enddo
         do no1 = 1,Norb
            do no2 = no1 +  1, Norb
               do nw = 1,Niom
                  Zp =  g_iom_mat(nk,nw)%el(no1,no2) - g_iom_mat(nk,nw)%el(no2,no1)
                  g_iom_mat(nk,nw)%el(no1,no2) = Zp/2.D0
                  g_iom_mat(nk,nw)%el(no2,no1) = Zp/2.D0
               enddo
            enddo
         enddo
      enddo
      deallocate( giom, gt0)

    end subroutine Tau_matz_spline_all

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------
! equidistant x-axis values
  subroutine aspline(x,y,x_new,y_new)
    implicit none

    real(Kind=Kind(0.d0)), dimension(:), intent(in) :: x,y,x_new
    real(Kind=Kind(0.d0)), dimension(:), intent(out) :: y_new

    integer :: i,j,n1,n2
    real(Kind=Kind(0.d0)), dimension(:), allocatable:: x_tmp,y_tmp,t
    real(Kind=Kind(0.d0)) :: dx,a,b,m1,m2,m3,m4

    n1 = size(x)
    n2 = size(x_new)

    allocate(x_tmp(n1+4),y_tmp(n1+4)) ! add two points at both sides

    dx = x(2)-x(1)
    x_tmp         = 0.d0
    y_tmp         = 0.d0
    x_tmp(3:n1+2) = x(:)
    y_tmp(3:n1+2) = y(:)

!Corner points
    x_tmp(1)    = x(1)  - 2.d0 * dx
    x_tmp(2)    = x(1)  - dx
    x_tmp(n1+3) = x(n1) + dx
    x_tmp(n1+4) = x(n1) + 2.d0 * dx

    y_tmp(n1+3) = yup(n1+3,x_tmp,y_tmp)
    y_tmp(n1+4) = yup(n1+4,x_tmp,y_tmp)
    y_tmp(2)    = ydn(2,x_tmp,y_tmp)
    y_tmp(1)    = ydn(1,x_tmp,y_tmp)

! Slopes
    allocate(t(n1))
    do i = 1,n1
       j = i + 2
       m1 = slope(dx,y_tmp(j-2),y_tmp(j-1))
       m2 = slope(dx,y_tmp(j-1),y_tmp(j))
       m3 = slope(dx,y_tmp(j),y_tmp(j+1))
       m4 = slope(dx,y_tmp(j+1),y_tmp(j+2))
       a = dabs(m4-m3) * m2 + dabs(m2-m1) * m3
       b = dabs(m4-m3) + dabs(m2-m1)
       if (b /= 0.d0) then
          t(i) = a / b
       else
          t(i) = 0.5d0 * (m2+m3)
       end if
    enddo

! Interpolate
    do i = 1,n2
       do j = 1,n1-1
          if (x_new(i) >= x(j) .and. x_new(i) <= x(j+1) ) &
               y_new(i) = poly(x(j),x(j+1),y(j),y(j+1),t(j),t(j+1),x_new(i))
       enddo
    enddo

    deallocate(x_tmp,y_tmp,t)

  end subroutine aspline

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------
  real(Kind=Kind(0.d0)) function yup(n,x,y)
    implicit none

    integer, intent(in) :: n
    real(Kind=Kind(0.d0)), dimension(:), intent(in) :: x,y

    yup = (2.d0 &
         * (y(n-1)-y(n-2))/(x(n-1)-x(n-2)) - (y(n-2)-y(n-3))/(x(n-2)-x(n-3)))  &
         * (x(n)-x(n-1)) + y(n-1)

  end function yup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------
  real(Kind=Kind(0.d0)) function ydn(n,x,y)
    implicit none

    integer, intent(in) :: n
    real(Kind=Kind(0.d0)), dimension(:), intent(in) :: x,y

    ydn = (-2.d0 &
         * (y(n+2)-y(n+1))/(x(n+2)-x(n+1)) + (y(n+3)-y(n+2))/(x(n+3)-x(n+2)))  &
         * (x(n+1)-x(n)) + y(n+1)

  end function ydn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------
  real(Kind=Kind(0.d0)) function slope(dx,y_dn,y_up)
    implicit none

    real(Kind=Kind(0.d0)), intent(in) :: dx,y_dn,y_up

    slope = (y_up - y_dn) / dx

  end function slope

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------
  real(Kind=Kind(0.d0)) function poly(x1,x2,y1,y2,t1,t2,x)
    implicit none

    real(Kind=Kind(0.d0)), intent(in) :: x1,x2,y1,y2,t1,t2,x
    real(Kind=Kind(0.d0)) :: p0,p1,p2,p3

    p0 = y1
    p1 = t1
    p2 = (3.d0*(y2-y1)/(x2-x1)-2.d0*t1-t2)/(x2-x1)
    p3 = (t1+t2-2.d0*(y2-y1)/(x2-x1))/(x2-x1)**2

    poly = p0 + p1 * (x-x1) + p2 * (x-x1)**2 + p3 * (x-x1)**3

  end function poly

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------

end Module Fourier
