! compile with
!gfortran -std=f2003 -I ../../Prog/ -I ../../Libraries/Modules/ -L ../../Libraries/Modules/ 9-Op-Phase.F90  ../../Prog/Operator_mod.o ../../Prog/Fields_mod.o ../../Libraries/Modules/modules_90.a -llapack -lblas


Program Testop_PHASE

  Use Operator_mod
  Use Fields_mod

  implicit none

  Complex(Kind = Kind(0.D0)) :: Phasenew, Phaseold, diff
  Integer :: i,n, N_SUN,  nf,nt
  Type (Fields)  :: nsigma
  Type(Operator) :: Op(3,3)


  Call nsigma%make(3,3)
  do nf = 1, 3
     do nt = 1, 3
        N_SUN = 3
        Call Op_make(Op(nf,nt), 3)
        do i = 1, Op(nf, nt)%N
           Op(nf, nt)%O(i,i) = real(2*i-3,kind(0.d0))
           Op(nf, nt)%P(i)   = i
        enddo
        Op(nf, nt)%g = 2.d0
        Op(nf ,nt)%type  = nf !  The type of an operator depends only on nf
        Op(nf ,nt)%alpha = cmplx(dble(nf), dble(nt), kind(0.D0))
        Call Op_Set(Op(nf,nt))
        nsigma%f(nf, nt) = real((2*Mod(nf, 2) - 1)*(1 + Mod(nf+nt, 2)),kind(0.d0))
        nsigma%t(nf    ) = nf 
     enddo
  enddo
  Phasenew = 1.D0
  Phaseold = 1.D0
  
  Do nf = 1, 3
     Call Op_Phase(Phasenew, Op, NSigma, nf)
  enddo
  Phasenew=Phasenew**N_SUN
  
  ! check against old version
  do nf = 1, Size(Op,2)
     do n = 1, size(Op,1)
        do nt = 1, size(nsigma%f,2)
           Phaseold=Phaseold*exp(cmplx(0.d0,Aimag(Op(n,nf)%g*Op(n,nf)%alpha) ,kind(0.D0)) * nsigma%Phi(n,nt)  )
        enddo
     enddo
  enddo
  Phaseold = Phaseold**N_SUN
  
  diff = Phaseold - Phasenew
  
  if (abs(DBLE(diff)) > MAX(ABS(DBLE(Phaseold)), ABS(DBLE(Phasenew)))*1D-14) then
     write (*,*) "ERROR in real", Phaseold, Phasenew, DBLE(diff), MAX(ABS(DBLE(Phaseold)), ABS(DBLE(Phasenew)))*1D-14
     STOP 2
  endif
  if (abs(AIMAG(diff)) > MAX(ABS(AIMAG(Phaseold)), ABS(AIMAG(Phasenew)))*2D-14) then
     write (*,*) "ERROR in imag", Phaseold, Phasenew, diff, MAX(ABS(AIMAG(Phaseold)), ABS(AIMAG(Phasenew)))*2D-14
     STOP 3
  endif

  Call nsigma%clear() 
  
  write (*,*) "success"
  
end Program TESTOP_PHASE
