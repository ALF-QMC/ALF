        Subroutine ObserT(NT, GT0,G0T,G00,GTT, PHASE, Mc_step_weight)
          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS
          Integer :: IMJ, I, J, I1, J1, no_I, no_J

          real (Kind=Kind(0.d0)) :: Z_z_ising
          Integer                   :: Ntau

!           Write(6,*) "ObserT NT =" , NT

          Z_z_ising = 0.d0
          do I = 1,Latt%N
            do Ntau = 1, Ltrot
              Z_z_ising = Z_z_ising + nsigma%f(I,Ntau)
            enddo
          enddo

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          If (NT == 0 ) then
             Obs_tau(1)%N = Obs_tau(1)%N + 1
             Obs_tau(1)%Ave_sign = Obs_tau(1)%Ave_sign + Real(ZS,kind(0.d0))
!              if ( Z_z_ising > 0 ) then
!                 Obs_tau(2)%N = Obs_tau(2)%N + 1
!                 Obs_tau(2)%Ave_sign = Obs_tau(2)%Ave_sign + Real(ZS,kind(0.d0))
!              else
!                 Obs_tau(3)%N = Obs_tau(3)%N + 1
!                 Obs_tau(3)%Ave_sign = Obs_tau(3)%Ave_sign + Real(ZS,kind(0.d0))
!              endif
          endif
             Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & +  Z * GT0(I1,J1,1) * ZP* ZS
!                    if ( Z_z_ising > 0 ) then
!                       Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
!                            & +  Z * GT0(I1,J1,1) * ZP* ZS
!                    else
!                       Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
!                            & +  Z * GT0(I1,J1,1) * ZP* ZS
!                    endif
                Enddo
             Enddo

        end Subroutine OBSERT
