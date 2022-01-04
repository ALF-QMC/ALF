

function E_kin(GRC)
   Implicit none
   Complex (Kind=Kind(0.d0)), intent(in) :: GRC(Ndim,Ndim,N_FL)

   Complex (Kind=Kind(0.d0)) :: E_kin

   Integer :: n, nf, I, J

   E_kin = cmplx(0.d0, 0.d0, kind(0.D0))
   Do n  = 1,Size(Op_T,1)
      Do nf = 1,N_FL
         Do I = 1,Size(Op_T(n,nf)%O,1)
            Do J = 1,Size(Op_T(n,nf)%O,2)
               E_kin = E_kin + Op_T(n,nf)%O(i, j)*Grc( Op_T(n,nf)%P(I), Op_T(n,nf)%P(J), nf )
            ENddo
         Enddo
      Enddo
   Enddo
   E_kin = E_kin * dble(N_SUN)
end function E_kin


Subroutine Obser(GR,Phase,Ntau, Mc_step_weight)
  Implicit none

  Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
  Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
  Integer, INTENT(IN)          :: Ntau
  Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight
  !Local
  Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZP, ZS, Z

  Complex (Kind=Kind(0.d0)) :: Z_z_ising1, Z_x_ising1, Z_m1, Z_chi
  Complex (Kind=Kind(0.d0)) :: Z_z_ising2, Z_x_ising2, Z_m2
  Complex (Kind=Kind(0.d0)) :: Zkin, Zpot, Zrho
  Complex (Kind=Kind(0.d0)) :: d_Z, d_X
  Integer :: nf, I, J, I1, nc1, imj, J1, no_I, no_J, nt1, nt, dnt, Ntau1, N_ising, n

  ZP = PHASE/Real(Phase, kind(0.D0))
  ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

!   Write(6,*) "Ntau =", Ntau
!   Write(6,*) "i_sweep =", i_sweep, ZP, ZS

  Do nf = 1,N_FL
    Do I = 1,Ndim
      Do J = 1,Ndim
        GRC(I, J, nf) = -GR(J, I, nf)
      Enddo
      GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
    Enddo
  Enddo
  ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >

  Z_z_ising1 = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_x_ising1 = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_m1       = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_z_ising2 = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_x_ising2 = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_m2       = cmplx(0.d0, 0.d0, kind(0.D0))

  nt1 = Ntau+1; if ( Ntau == Ltrot ) nt1 = 1

  Zkin = E_kin(GRC)

  Zpot = cmplx(0.d0, 0.d0, kind(0.D0))
  do I = 1, size(Op_V, 1)
    I1 = Op_V(I,1)%P(1)
    do nc1 = 2, size( Op_V(I,1)%O, 1 )
      Zpot = Zpot + nsigma%f(I,Ntau) * GRC(I1, Op_V(I,1)%P(nc1) ,1) * Op_V(I,1)%O(1,nc1)
      Zpot = Zpot + nsigma%f(I,Ntau) * GRC(Op_V(I,1)%P(nc1), I1 ,1) * Op_V(I,1)%O(nc1,1)
    enddo
  enddo
  Zpot = Zpot * dble(N_SUN) * ham_t * ham_xi

  call Obs_scal(3)%measure( [Zkin, Zpot, Zkin+Zpot], Phase )

  Zrho = cmplx(0.d0, 0.d0, kind(0.D0))
  do I = 1, Ndim
    Zrho = Zrho + Grc(i, i, 1)
  enddo
  Zrho = Zrho * dble(N_SUN)

  call Obs_scal(1)%measure( [Zrho], Phase )

  do I = 1,Latt%N
    Z_z_ising1 = Z_z_ising1 + nsigma%f(I,Ntau)
    Z_m1 = Z_m1 + nsigma%f(1,Ntau)*nsigma%f(I,Ntau)

    if ( nsigma%f(I,nt1) == nsigma%f(I,Ntau) ) then
      Z_x_ising1 = Z_x_ising1 + eq_x_ising
    else
      Z_x_ising1 = Z_x_ising1 + neq_x_ising
    endif
  enddo
  Z_z_ising1 = Z_z_ising1/Latt%N
  Z_x_ising1 = Z_x_ising1/Latt%N
  Z_m1 = Z_m1/Latt%N

  call Obs_scal(2)%measure( [Z_z_ising1], Phase )
  call Obs_scal(4)%measure( [Z_x_ising1], Phase )
  call Obs_scal(5)%measure( [Z_m1, Z_m1**2, Z_m1**4], Phase )

  If ( Model_vers == 3 ) then
    do I = 1+Latt%N, 2*Latt%N
      Z_z_ising2 = Z_z_ising2 + nsigma%f(I,Ntau)
      Z_m2 = Z_m2 + nsigma%f(1+Latt%N,Ntau)*nsigma%f(I,Ntau)

      if ( nsigma%f(I,nt1) == nsigma%f(I,Ntau) ) then
        Z_x_ising2 = Z_x_ising2 + eq_x_ising
      else
        Z_x_ising2 = Z_x_ising2 + neq_x_ising
      endif
    enddo
    Z_z_ising2 = Z_z_ising2/Latt%N
    Z_x_ising2 = Z_x_ising2/Latt%N
    Z_m2 = Z_m2/Latt%N

    call Obs_scal(2)%measure( [Z_z_ising2], Phase )
    call Obs_scal(4)%measure( [Z_x_ising2], Phase )
    call Obs_scal(5)%measure( [Z_m2, Z_m2**2, Z_m2**4], Phase )
  endif


  ! counting up correlation functions
!   DO I = 1,Size(Obs_eq,1)
  DO I = 1,3
    Obs_eq(I)%N        = Obs_eq(I)%N + 1
    Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
  ENDDO

  ! Compute Ising X-X and Z-Z correlation functions
  nt1 = Ntau+1; if ( Ntau == Ltrot ) nt1 = 1
  Do I = 1,Latt%N
    do no_I = 1, size(Obs_eq(1)%Obs_Latt0,1)
      I1 = I + (no_I-1)*Latt%N
      if ( nsigma%f(I1,nt1) == nsigma%f(I1,Ntau) ) then
        Z_x_ising1 = eq_x_ising
      else
        Z_x_ising1 = neq_x_ising
      endif
      Obs_eq(1)%Obs_Latt0(no_I) = Obs_eq(1)%Obs_Latt0(no_I) + Z_x_ising1      * ZP*ZS
  !     Obs_eq(2)%Obs_Latt0(1) = Obs_eq(2)%Obs_Latt0(1) + nsigma%f(I,Ntau) * ZP*ZS
      Do J = 1,Latt%N
        do no_J = 1, size(Obs_eq(1)%Obs_Latt0,1)
          J1 = J + (no_J-1)*Latt%N
          imj = Latt%imj(I,J)
          Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + nsigma%f(I1,Ntau) * nsigma%f(J1,Ntau) * ZP*ZS
          if ( nsigma%f(J1,nt1) == nsigma%f(J1,Ntau) ) then
            Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + Z_x_ising1 * eq_x_ising  * ZP*ZS
          else
            Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + Z_x_ising1 * neq_x_ising * ZP*ZS
          endif
        enddo
      enddo
    enddo
  enddo

!   If ( Model_vers == 3 ) then
!     Do I = 1+Latt%N,2*Latt%N
!       if ( nsigma%f(I,nt1) == nsigma%f(I,Ntau) ) then
!         Z_x_ising1 = eq_x_ising
!       else
!         Z_x_ising1 = neq_x_ising
!       endif
!       Obs_eq(1)%Obs_Latt0(1) = Obs_eq(1)%Obs_Latt0(1) + Z_x_ising1      * ZP*ZS
! !       Obs_eq(2)%Obs_Latt0(1) = Obs_eq(2)%Obs_Latt0(1) + nsigma%f(I,Ntau) * ZP*ZS
!       Do J = 1+Latt%N,2*Latt%N
!         imj = Latt%imj(I-Latt%N,J-Latt%N)
!         Obs_eq(2)%Obs_Latt(imj,1,1,1) = Obs_eq(2)%Obs_Latt(imj,1,1,1) + nsigma%f(I,Ntau) * nsigma%f(J,Ntau) * ZP*ZS
!         if ( nsigma%f(J,nt1) == nsigma%f(J,Ntau) ) then
!           Obs_eq(1)%Obs_Latt(imj,1,1,1) = Obs_eq(1)%Obs_Latt(imj,1,1,1) + Z_x_ising1 * eq_x_ising  * ZP*ZS
!         else
!           Obs_eq(1)%Obs_Latt(imj,1,1,1) = Obs_eq(1)%Obs_Latt(imj,1,1,1) + Z_x_ising1 * neq_x_ising * ZP*ZS
!         endif
!       enddo
!     enddo
!   endif

  ! Computing time-displaced X-X and Z-Z correlation functions and chi

  nBlub2 = nBlub2 + 1
  if ( nBlub2 > 48 ) then
    nBlub2 = 0

  DO I = 4,5
    Obs_eq(I)%N        = Obs_eq(I)%N + 1
    Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
  ENDDO

  Z_chi     = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_x_ising1 = cmplx(0.d0, 0.d0, kind(0.D0))
  Ntau1 = Ntau+1; if ( Ntau == Ltrot ) Ntau1 = 1

  do dnt= 0, Ltrot
    nt = Ntau + dnt
    if ( nt > Ltrot ) nt = nt - Ltrot
    nt1 = nt+1; if ( nt == Ltrot ) nt1 = 1
    do I = 1,Latt%N
      if ( nsigma%f(I,nt1) == nsigma%f(I,nt) ) then
        Z_x_ising1 = eq_x_ising
      else
        Z_x_ising1 = neq_x_ising
      endif
      Z_x_ising1 = Z_x_ising1 + Z_x_ising1
      Obs_eq(5)%Obs_Latt0(1) = Obs_eq(5)%Obs_Latt0(1) + Z_x_ising1    * ZP*ZS
!       Obs_eq(4)%Obs_Latt0(1) = Obs_eq(4)%Obs_Latt0(1) + nsigma%f(I,nt) * ZP*ZS
      do J = 1,Latt%N
        imj = Latt%imj(I,J)
        Obs_eq(4)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(4)%Obs_Latt(imj,dnt+1,1,1) + nsigma%f(I,Ntau) * nsigma%f(J,nt) * ZP*ZS
        if ( nt == Ntau .and. I == J ) then
          Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + 1  * ZP*ZS
          Z_chi = Z_chi + 1
        elseif ( nsigma%f(J,Ntau1) == nsigma%f(J,Ntau) ) then
          Z_chi = Z_chi + Z_x_ising1 * eq_x_ising
          Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + Z_x_ising1 * eq_x_ising  * ZP*ZS
        else
          Z_chi = Z_chi + Z_x_ising1 * neq_x_ising
          Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + Z_x_ising1 * neq_x_ising * ZP*ZS
        endif
      enddo
    enddo
  enddo

  If ( Model_vers == 3 ) then
    do dnt= 0, Ltrot
      nt = Ntau + dnt
      if ( nt > Ltrot ) nt = nt - Ltrot
      nt1 = nt+1; if ( nt == Ltrot ) nt1 = 1
      do I = 1+Latt%N,2*Latt%N
        if ( nsigma%f(I,Ntau1) == nsigma%f(I,Ntau) ) then
          Z_x_ising1 = eq_x_ising
        else
          Z_x_ising1 = neq_x_ising
        endif
        Z_x_ising1 = Z_x_ising1 + Z_x_ising1
        Obs_eq(5)%Obs_Latt0(1) = Obs_eq(5)%Obs_Latt0(1) + Z_x_ising1    * ZP*ZS
!         Obs_eq(4)%Obs_Latt0(1) = Obs_eq(4)%Obs_Latt0(1) + nsigma%f(I,nt) * ZP*ZS
        do J = 1+Latt%N,2*Latt%N
          imj = Latt%imj(I-Latt%N,J-Latt%N)
          Obs_eq(4)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(4)%Obs_Latt(imj,dnt+1,1,1) + nsigma%f(I,Ntau) * nsigma%f(J,nt) * ZP*ZS
          if ( nt == Ntau .and. I == J ) then
            Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + 1  * ZP*ZS
            Z_chi = Z_chi + 1
          elseif ( nsigma%f(J,nt) == nsigma%f(J,Ntau) ) then
            Z_chi = Z_chi + Z_x_ising1 * eq_x_ising
            Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + Z_x_ising1 * eq_x_ising  * ZP*ZS
          else
            Z_chi = Z_chi + Z_x_ising1 * neq_x_ising
            Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + Z_x_ising1 * neq_x_ising * ZP*ZS
          endif
        enddo
      enddo
    enddo
  endif

  If ( Model_vers == 3 ) then
    N_ising = 2*Latt%N
  else
    N_ising = Latt%N
  endif

  call Obs_scal(6)%measure( [ Z_chi/(N_ising**2), Z_x_ising1/N_ising/Ltrot ], Phase )

  endif

  ! Compute Green-function
  Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
  Do I1 = 1,Ndim
    I    = List(I1,1)
    no_I = List(I1,2)
    Do J1 = 1,Ndim
      J    = List(J1,1)
      no_J = List(J1,2)
      imj = Latt%imj(I,J)
      Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + Z * GRC(I1,J1,1) *  ZP*ZS
    enddo
  enddo

end Subroutine Obser


Subroutine measure_hist(ham, Phase)
  Implicit none
  class(ham_Nematic_Dirac) , intent(inout) :: ham
  Complex (Kind=Kind(0.d0)), Intent(IN)    :: Phase

  Complex (Kind=Kind(0.d0)) :: ZP, ZS

  Complex (Kind=Kind(0.d0)) :: Z_z_ising1, Z_x_ising1, Z_m1, Z_m1_sq, Z_m1_quad
  Complex (Kind=Kind(0.d0)) :: Z_z_ising2, Z_x_ising2, Z_m2, Z_m2_sq, Z_m2_quad
  Complex (Kind=Kind(0.d0)) :: Z_m_temp
  real    (Kind=Kind(0.d0)) :: B
  Complex (Kind=Kind(0.d0)) :: d_Z, d_X
  Integer :: I, nt1, nt


  ZP = PHASE/Real(Phase, kind(0.D0))
  ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))


  Z_z_ising1 = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_x_ising1 = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_m1       = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_m1_sq    = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_m1_quad  = cmplx(0.d0, 0.d0, kind(0.D0))

  do nt = 1, Ltrot
    nt1 = nt+1; if ( nt == Ltrot ) nt1 = 1
    Z_m_temp = cmplx(0.d0, 0.d0, kind(0.D0))
    do I = 1,Latt%N
      Z_z_ising1 = Z_z_ising1 + nsigma%f(I,nt)
      Z_m_temp = Z_m_temp + nsigma%f(1,nt)*nsigma%f(I,nt)

      if ( nsigma%f(I,nt1) == nsigma%f(I,nt) ) then
        Z_x_ising1 = Z_x_ising1 + eq_x_ising
      else
        Z_x_ising1 = Z_x_ising1 + neq_x_ising
      endif
    enddo
    Z_m_temp = Z_m_temp/Latt%N
    Z_m1 = Z_m1 + Z_m_temp
    Z_m1_sq = Z_m1_sq + Z_m_temp**2
    Z_m1_quad = Z_m1_quad + Z_m_temp**4
  enddo

  Z_z_ising1 = Z_z_ising1/(Latt%N*Ltrot)
  Z_x_ising1 = Z_x_ising1/(Latt%N*Ltrot)
  Z_m1 = Z_m1/Ltrot
  Z_m1_sq = Z_m1_sq/Ltrot
  Z_m1_quad = Z_m1_quad/Ltrot

  call Obs_scal(7)%measure( [Z_z_ising1], Phase )
  call Obs_scal(8)%measure( [Z_x_ising1], Phase )
  call Obs_scal(9)%measure( [Z_m1, Z_m1**2, Z_m1**4], Phase )

  call Obs_hist(1)%measure( dble(Z_x_ising1), dble(ZS) )
  call Obs_hist(2)%measure( dble(Z_m1), dble(ZS) )
  B = (3.d0 - dble(Z_m1_quad/Z_m1_sq**2))/2.d0
  call Obs_hist(3)%measure( B, dble(ZS) )

  If ( Model_vers == 3 ) then
    Z_z_ising2 = cmplx(0.d0, 0.d0, kind(0.D0))
    Z_x_ising2 = cmplx(0.d0, 0.d0, kind(0.D0))
    Z_m2       = cmplx(0.d0, 0.d0, kind(0.D0))
    Z_m2_sq    = cmplx(0.d0, 0.d0, kind(0.D0))
    Z_m2_quad  = cmplx(0.d0, 0.d0, kind(0.D0))

    do nt = 1, Ltrot
      nt1 = nt+1; if ( nt == Ltrot ) nt1 = 1
      Z_m_temp = cmplx(0.d0, 0.d0, kind(0.D0))
      do I = 1+Latt%N, 2*Latt%N
        Z_z_ising2 = Z_z_ising2 + nsigma%f(I,nt)
        Z_m_temp = Z_m_temp + nsigma%f(1+Latt%N,nt)*nsigma%f(I,nt)

        if ( nsigma%f(I,nt1) == nsigma%f(I,nt) ) then
          Z_x_ising2 = Z_x_ising2 + eq_x_ising
        else
          Z_x_ising2 = Z_x_ising2 + neq_x_ising
        endif
      enddo
      Z_m_temp = Z_m_temp/Latt%N
      Z_m2 = Z_m2 + Z_m_temp
      Z_m2_sq = Z_m2_sq + Z_m_temp**2
      Z_m2_quad = Z_m2_quad + Z_m_temp**4
    enddo
    
    Z_z_ising2 = Z_z_ising2/(Latt%N*Ltrot)
    Z_x_ising2 = Z_x_ising2/(Latt%N*Ltrot)
    Z_m2 = Z_m2/Ltrot
    Z_m2_sq = Z_m2_sq/Ltrot
    Z_m2_quad = Z_m2_quad/Ltrot

    d_Z = (Z_z_ising1 - Z_z_ising2)**2
    d_X = (Z_x_ising1 - Z_x_ising2)**2

    call Obs_scal(7)%measure( [Z_z_ising2], Phase )
    call Obs_scal(8)%measure( [Z_x_ising2], Phase )
    call Obs_scal(9)%measure( [Z_m2, Z_m2**2, Z_m2**4], Phase )

    call Obs_scal(10)%measure( [d_Z], Phase )
    call Obs_scal(11)%measure( [d_X], Phase )

    call Obs_hist(1)%measure( dble(Z_x_ising2), dble(ZS) )
    call Obs_hist(2)%measure( dble(Z_m2), dble(ZS) )
    call Obs_hist(3)%measure( (3.d0 - dble(Z_m2_quad/Z_m2_sq**2))/2.d0, dble(ZS) )

    call Obs_hist(4)%measure( dble(d_Z), dble(ZS) )
    call Obs_hist(5)%measure( dble(d_X), dble(ZS) )
  endif

end Subroutine measure_hist
