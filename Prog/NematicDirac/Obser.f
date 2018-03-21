Subroutine Obser(GR,Phase,Ntau)
  Implicit none
  
  Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
  Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
  Integer, INTENT(IN)          :: Ntau
  !Local 
  Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZP, ZS, Z

  Complex (Kind=Kind(0.d0)) :: Zkin, Z_F_by_xi, Z_z_ising, Z_x_ising, Z_m, Z_chi, Z_x_ising1, Z_m1, Z_mt
  Integer :: nf, I, J, I1, nc1, imj, J1, no_I, no_J, nt1, nt, dnt, Ntau1
  
  ZP = PHASE/Real(Phase, kind(0.D0))
  ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
  
!   Write(6,*) "Ntau =", Ntau
!   Write(6,*) "i_sweep =", i_sweep, ZP, ZS
  
IF ( (Lattice_type =="BipartiteSquare" .and. Model == "NematicDirac") .or. &
   & (Lattice_type =="Square"          .and. Model == "NematicDirac2")  .or. &
   & (Lattice_type =="BipartiteSquare" .and. Model == "NematicDirac3") ) then

  Do nf = 1,N_FL
    Do I = 1,Ndim
      Do J = 1,Ndim
        GRC(I, J, nf) = -GR(J, I, nf)
      Enddo
      GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
    Enddo
  Enddo
  ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >

  ! Compute scalar observables. 
!   Do I = 1,Size(Obs_scal,1)
  Do I = 1,5
    Obs_scal(I)%N         =  Obs_scal(I)%N + 1
    Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
  Enddo
  

!           Zkin = cmplx(0.d0, 0.d0, kind(0.D0))conjg(
!           Do nf = 1,N_FL
!             Do J = 1,Ndim
!               Zkin = Zkin + sum(Op_T(1,nf)%O(:, j)*Grc(:, j, nf))
!             ENddo
!           Enddo
!           Zkin = Zkin * dble(N_SUN)
!           Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin * ZP*ZS

  Z_F_by_xi = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_z_ising = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_x_ising = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_m       = cmplx(0.d0, 0.d0, kind(0.D0))
  nt1 = Ntau+1; if ( Ntau == Ltrot ) nt1 = 1
  
  do I = 1,Latt%N
    Z_z_ising = Z_z_ising + nsigma(I,Ntau)
    Z_m = Z_m + nsigma(1,Ntau)*nsigma(I,Ntau)
    
    if ( nsigma(I,nt1) == nsigma(I,Ntau) ) then
      Z_x_ising = Z_x_ising + eq_x_ising
    else
      Z_x_ising = Z_x_ising + neq_x_ising
    endif
      
    I1 = Op_V(I,1)%P(1)
!     Do nc1 = 1,N_coord
!       Z_F_by_xi = Z_F_by_xi + nsigma(I,Ntau) * GRC(I1, Op_V(I,1)%P(nc1+1) ,1) * Op_V(I,1)%O(1,nc1+1)
!       Z_F_by_xi = Z_F_by_xi + nsigma(I,Ntau) * GRC(Op_V(I,1)%P(nc1+1), I1 ,1) * Op_V(I,1)%O(nc1+1,1)
!     enddo
  enddo
  Obs_scal(2)%Obs_vec(1) = Obs_scal(2)%Obs_vec(1) + Z_z_ising/Latt%N * ZP*ZS
  Obs_scal(3)%Obs_vec(1) = Obs_scal(3)%Obs_vec(1) + Z_F_by_xi/Latt%N * cmplx(ham_t, 0, kind(0.D0)) * ZP*ZS
  Obs_scal(4)%Obs_vec(1) = Obs_scal(4)%Obs_vec(1) + Z_x_ising/Latt%N * ZP*ZS 
  Z_m = Z_m/Latt%N
  Obs_scal(5)%Obs_vec(1) = Obs_scal(5)%Obs_vec(1) + Z_m    * ZP*ZS
  Obs_scal(5)%Obs_vec(2) = Obs_scal(5)%Obs_vec(2) + Z_m**2 * ZP*ZS
  Obs_scal(5)%Obs_vec(3) = Obs_scal(5)%Obs_vec(3) + Z_m**4 * ZP*ZS


  ! counting up correlation functions
!   DO I = 1,Size(Obs_eq,1)
  DO I = 1,3
    Obs_eq(I)%N        = Obs_eq(I)%N + 1
    Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
  ENDDO
      Obs_eq(6)%N        = Obs_eq(6)%N + 1
      Obs_eq(6)%Ave_sign = Obs_eq(6)%Ave_sign + real(ZS,kind(0.d0))
            
  ! Compute Ising X-X and Z-Z correlation functions
  nt1 = Ntau+1; if ( Ntau == Ltrot ) nt1 = 1
  Do I = 1,Latt%N
    if ( nsigma(I,nt1) == nsigma(I,Ntau) ) then
      Z_x_ising = eq_x_ising
    else
      Z_x_ising = neq_x_ising
    endif
    Obs_eq(1)%Obs_Latt0(1) = Obs_eq(1)%Obs_Latt0(1) + Z_x_ising      * ZP*ZS /Latt%N
    Obs_eq(2)%Obs_Latt0(1) = Obs_eq(2)%Obs_Latt0(1) + nsigma(I,Ntau) * ZP*ZS /Latt%N
    Do J = 1,Latt%N
      imj = latt%imj(I,J)
      Obs_eq(2)%Obs_Latt(imj,1,1,1) = Obs_eq(2)%Obs_Latt(imj,1,1,1) + nsigma(I,Ntau) * nsigma(J,Ntau) * ZP*ZS /Latt%N
      if ( nsigma(J,nt1) == nsigma(J,Ntau) ) then
        Obs_eq(1)%Obs_Latt(imj,1,1,1) = Obs_eq(1)%Obs_Latt(imj,1,1,1) + Z_x_ising * eq_x_ising  * ZP*ZS /Latt%N
      else
        Obs_eq(1)%Obs_Latt(imj,1,1,1) = Obs_eq(1)%Obs_Latt(imj,1,1,1) + Z_x_ising * neq_x_ising * ZP*ZS /Latt%N
      endif
    enddo
  enddo
  
  ! Computing time-displaced X-X and Z-Z correlation functions and chi
!   if ( Ntau == 1 ) then
!     nBlub = nBlub + 1
!     if ( nBlub == 2 ) then !trigers after every full up-down sweep, if LOBS_ST = 1
!       nBlub = 0
!       i_sweep = i_sweep+1
  
  nBlub2 = nBlub2 + 1
  if ( nBlub2 > 48 ) then
    nBlub2 = 0
      
  DO I = 4,5
    Obs_eq(I)%N        = Obs_eq(I)%N + 1
    Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
  ENDDO
  Do I = 6,6
    Obs_scal(I)%N         =  Obs_scal(I)%N + 1
    Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
  Enddo
      
  Z_chi     = cmplx(0.d0, 0.d0, kind(0.D0))
  Z_x_ising = cmplx(0.d0, 0.d0, kind(0.D0))
  Ntau1 = Ntau+1; if ( Ntau == Ltrot ) Ntau1 = 1
  
  do dnt= 0, Ltrot
    nt = Ntau + dnt
    if ( nt > Ltrot ) nt = nt - Ltrot
    nt1 = nt+1; if ( nt == Ltrot ) nt1 = 1
    do I = 1,Latt%N
      if ( nsigma(I,Ntau1) == nsigma(I,Ntau) ) then
        Z_x_ising1 = eq_x_ising
      else
        Z_x_ising1 = neq_x_ising
      endif
      Z_x_ising = Z_x_ising + Z_x_ising1
      Obs_eq(5)%Obs_Latt0(1) = Obs_eq(5)%Obs_Latt0(1) + Z_x_ising    * ZP*ZS /Latt%N
      Obs_eq(4)%Obs_Latt0(1) = Obs_eq(4)%Obs_Latt0(1) + nsigma(I,nt) * ZP*ZS /Latt%N
      do J = 1,Latt%N
        imj = latt%imj(I,J)
        Obs_eq(4)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(4)%Obs_Latt(imj,dnt+1,1,1) + nsigma(I,Ntau) * nsigma(J,nt) * ZP*ZS /Latt%N
        if ( nt == Ntau .and. I == J ) then
          Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + 1  * ZP*ZS /Latt%N
          Z_chi = Z_chi + 1
        elseif ( nsigma(J,nt) == nsigma(J,Ntau) ) then
          Z_chi = Z_chi + Z_x_ising1 * eq_x_ising
          Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + Z_x_ising1 * eq_x_ising  * ZP*ZS /Latt%N
        else
          Z_chi = Z_chi + Z_x_ising1 * neq_x_ising
          Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) = Obs_eq(5)%Obs_Latt(imj,dnt+1,1,1) + Z_x_ising1 * neq_x_ising * ZP*ZS /Latt%N
        endif
      enddo
    enddo
  enddo
  Obs_scal(6)%Obs_vec(1) = Obs_scal(6)%Obs_vec(1) + Z_chi     * ZP*ZS /(Latt%N**2)
  Obs_scal(6)%Obs_vec(2) = Obs_scal(6)%Obs_vec(2) + Z_x_ising * ZP*ZS /Latt%N/Ltrot
  
  endif
          
  ! Compute Green-function and spinZ correlations
  Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
  Do I1 = 1,Ndim
    I    = List(I1,1)
    no_I = List(I1,2)
    Do J1 = 1,Ndim
      J    = List(J1,1)
      no_J = List(J1,2)
      imj = latt%imj(I,J)
      Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + Z * GRC(I1,J1,1) *  ZP*ZS
      ! SpinZ
      Obs_eq(6)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(6)%Obs_Latt(imj,1,no_I,no_J) + &
      &               Z * GRC(I1,J1,1) * GR(I1,J1,1) * ZP*ZS
    enddo
  enddo



      
! elseIF ( Lattice_type =="Square" .and. Model == "yyhe" ) then
!   Do nf = 1,N_FL
!     Do I = 1,Ndim
!       Do J = 1,Ndim
!         GRC(I, J, nf) = -GR(J, I, nf)
!       Enddo
!       GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
!     Enddo
!   Enddo
!   ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >
! 
!   ! Compute scalar observables. 
!   Do I = 1,Size(Obs_scal,1)
!     Obs_scal(I)%N         =  Obs_scal(I)%N + 1
!     Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
!   Enddo
! 
!   Z_F_by_xi = cmplx(0.d0, 0.d0, kind(0.D0))
!   Z_z_ising = cmplx(0.d0, 0.d0, kind(0.D0))
!   Z_x_ising = cmplx(0.d0, 0.d0, kind(0.D0))
!   Z_m       = cmplx(0.d0, 0.d0, kind(0.D0))
!   nt1 = Ntau+1; if ( Ntau == Ltrot ) nt1 = 1
!   
!   do I = 1,Ndim
!   Z_z_ising = Z_z_ising + nsigma(I,Ntau)
!   Z_m = Z_m + nsigma(1,Ntau)*nsigma(I,Ntau)
!     
!     if ( nsigma(I,nt1) == nsigma(I,Ntau) ) then
!       Z_x_ising = Z_x_ising + eq_x_ising
!     else
!       Z_x_ising = Z_x_ising + neq_x_ising
!     endif
!       
! !     I1 = Op_V(I,1)%P(1)
! !     Do nc1 = 1,N_coord
! !       Z_F_by_xi = Z_F_by_xi + nsigma(I,Ntau) * GRC(I1, Op_V(I,1)%P(nc1+1) ,1) * Op_V(I,1)%O(1,nc1+1)
! !       Z_F_by_xi = Z_F_by_xi + nsigma(I,Ntau) * GRC(Op_V(I,1)%P(nc1+1), I1 ,1) * Op_V(I,1)%O(nc1+1,1)
! !     enddo
!   enddo
!   Obs_scal(2)%Obs_vec(1) = Obs_scal(2)%Obs_vec(1) + Z_z_ising/Ndim * ZP*ZS
! !   Obs_scal(3)%Obs_vec(1) = Obs_scal(3)%Obs_vec(1) + Z_F_by_xi/Latt%N * cmplx(ham_t, 0, kind(0.D0)) * ZP*ZS
!   Obs_scal(4)%Obs_vec(1) = Obs_scal(4)%Obs_vec(1) + Z_x_ising/Ndim * ZP*ZS 
!   Z_m = Z_m/Ndim
!   Obs_scal(5)%Obs_vec(1) = Obs_scal(5)%Obs_vec(1) + Z_m    * ZP*ZS
!   Obs_scal(5)%Obs_vec(2) = Obs_scal(5)%Obs_vec(2) + Z_m**2 * ZP*ZS
!   Obs_scal(5)%Obs_vec(3) = Obs_scal(5)%Obs_vec(3) + Z_m**4 * ZP*ZS
! 
!   ! counting up correlation functions
!   DO I = 1,Size(Obs_eq,1)
! !   DO I = 1,4
!     Obs_eq(I)%N        = Obs_eq(I)%N + 1
!     Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
!   ENDDO
!             
!   ! Compute Ising X-X and Z-Z correlation functions
!   nt1 = Ntau+1; if ( Ntau == Ltrot ) nt1 = 1
!   Do I = 1,Latt%N
!     if ( nsigma(I,nt1) == nsigma(I,Ntau) ) then
!       Z_x_ising = eq_x_ising
!     else
!       Z_x_ising = neq_x_ising
!     endif
!     Obs_eq(1)%Obs_Latt0(1) = Obs_eq(1)%Obs_Latt0(1) + Z_x_ising      * ZP*ZS /Latt%N
!     Obs_eq(2)%Obs_Latt0(1) = Obs_eq(2)%Obs_Latt0(1) + nsigma(I,Ntau) * ZP*ZS /Latt%N
!     Do J = 1,Latt%N
!       imj = latt%imj(I,J)
!       Obs_eq(2)%Obs_Latt(imj,1,1,1) = Obs_eq(2)%Obs_Latt(imj,1,1,1) + nsigma(I,Ntau) * nsigma(J,Ntau) * ZP*ZS /Latt%N
!       if ( nsigma(J,nt1) == nsigma(J,Ntau) ) then
!         Obs_eq(1)%Obs_Latt(imj,1,1,1) = Obs_eq(1)%Obs_Latt(imj,1,1,1) + Z_x_ising * eq_x_ising  * ZP*ZS /Latt%N
!       else
!         Obs_eq(1)%Obs_Latt(imj,1,1,1) = Obs_eq(1)%Obs_Latt(imj,1,1,1) + Z_x_ising * neq_x_ising * ZP*ZS /Latt%N
!       endif
!     enddo
!   enddo
!           
!   ! Compute Green-function
!   Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
!   Do I1 = 1,Ndim
!     I    = List(I1,1)
!     no_I = List(I1,2)
!     if ( nsigma(I1,nt1) == nsigma(I1,Ntau) ) then
!       Z_x_ising = eq_x_ising
!     else
!       Z_x_ising = neq_x_ising
!     endif
!     Obs_eq(1)%Obs_Latt0(no_I) = Obs_eq(1)%Obs_Latt0(no_I) + Z_x_ising       * ZP*ZS /Ndim
!     Obs_eq(2)%Obs_Latt0(no_I) = Obs_eq(2)%Obs_Latt0(no_I) + nsigma(I1,Ntau) * ZP*ZS /Ndim
!     Do J1 = 1,Ndim
!       J    = List(J1,1)
!       no_J = List(J1,2)
!       imj = latt%imj(I,J)
!       if ( nsigma(J,nt1) == nsigma(J,Ntau) ) then
!         Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + Z_x_ising * eq_x_ising  * ZP*ZS /Ndim
!       else
!         Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + Z_x_ising * neq_x_ising * ZP*ZS /Ndim
!       endif
!       Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) + nsigma(I1,Ntau) * nsigma(J1,Ntau) * ZP*ZS /Ndim
!       ! Green
!       Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) + ( GRC(I1,J1,1) + GRC(I1,J1,2)) *  ZP*ZS 
!     enddo
!   enddo
  
endif

end Subroutine Obser
