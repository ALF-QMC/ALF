        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE)
          Implicit none
          
          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: Z, ZP, ZS
          Integer :: IMJ, I, J, I1, J1, no_I, no_J
          
!           Write(6,*) "ObserT NT =" , NT

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          If (NT == 0 ) then 
             DO I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             ENDDO
          endif
          If ( Model == "Hubbard_SU2" .or. Model == "Hubbard_SU2_Ising" .or. Model == "NematicDirac" .or. &
             & Model == "NematicDirac2" .or. Model == "NematicDirac3" ) then 
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
                   
!                    ! SpinZ
!                    Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
!                         &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
!                    
!                    ! SpinXY
!                    Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
!                         &      - Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
!                    
!                    ! Den
!                    Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)  &
!                         & + ( Z*Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1))*       &
!                         &         (cmplx(1.d0,0.d0,kind(0.d0)) - G00(J1,J1,1))  -     &
!                         &     Z * GT0(I1,J1,1)*G0T(J1,I1,1)                                ) * ZP * ZS
                Enddo
!                 Obs_tau(4)%Obs_Latt0(no_I) = Obs_tau(4)%Obs_Latt0(no_I) + &
!                      &         Z*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1)) * ZP * ZS
             Enddo
          Elseif ( Model == "Hubbard_Mz"  ) then 
             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   !Green
                   Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &   +   ( GT0(I1,J1,1) + GT0(I1,J1,2) ) * ZP* ZS

                   !SpinZ
                   Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                       & +  ( &
                       &    (GTT(I1,I1,1) -  GTT(I1,I1,2) ) * ( G00(J1,J1,1)  -  G00(J1,J1,2) )   &
                       &  - (G0T(J1,I1,1) * GT0(I1,J1,1)  +  G0T(J1,I1,2) * GT0(I1,J1,2) )    )*ZP*ZS

                   !SpinXY
                   Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        &  - &
                        &   (G0T(J1,I1,1) * GT0(I1,J1,2)  +  G0T(J1,I1,2) * GT0(I1,J1,1))*ZP*ZS
                   !Den
                   Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & +  (                                        &  
                        &    (cmplx(2.D0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2) ) * &
                        &    (cmplx(2.D0,0.d0,kind(0.d0)) - G00(J1,J1,1) - G00(J1,J1,2) )   &
                        & -  ( G0T(J1,I1,1) * GT0(I1,J1,1) + G0T(J1,I1,2) * GT0(I1,J1,2) )  )*ZP*ZS     

                enddo
             
                Obs_tau(4)%Obs_Latt0(no_I) =  Obs_tau(4)%Obs_Latt0(no_I) + &
                     &       (cmplx(2.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1) - GTT(I1,I1,2)) * ZP * ZS
             Enddo
          Endif
          
        end Subroutine OBSERT
