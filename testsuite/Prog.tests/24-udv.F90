Program  test

  Use Random_Wrap
  Use UDV_State_mod
  Use MyMats 

  Implicit none

  Interface
     Subroutine Check_udv(Tmp, UDV, Diff)
       
       Use UDV_State_mod
       Implicit none
       
       Type (UDV_state), Intent(IN) :: udv
       Complex (Kind=Kind(0.d0)), allocatable  :: Tmp(:,:)
       Real (Kind=Kind(0.d0)), INTENT(out) :: Diff
     End Subroutine Check_udv
  end Interface
  Type (UDV_state) :: udvr, udvl
  Integer :: Ndim = 16, Seed_Vec(1), I, I1, nth
  Complex(Kind =Kind(0.d0)), allocatable :: TMP(:,:),  TMP1(:,:),  TmpR(:,:), TmpL(:,:)
  Real(Kind =Kind(0.d0)) :: Diff

  allocate (TMP(Ndim,Ndim), TMP1(Ndim,Ndim), TmpL(Ndim,Ndim), TmpR(Ndim,Ndim) )
  call  udvr%alloc(Ndim)
  call  udvl%alloc(Ndim)
  call  udvr%reset("r")
  call  udvl%reset("l")

  Seed_vec(1) = 4782347
  Call Ranset(Seed_vec)

  TmpR = cmplx(0.d0,0.d0,Kind=Kind(0.d0))
  TmpL = cmplx(0.d0,0.d0,Kind=Kind(0.d0))
  Do I = 1,Ndim
     TmpR(I,I) = cmplx(1.d0,0.d0,Kind=Kind(0.d0))
     TmpL(I,I) = cmplx(1.d0,0.d0,Kind=Kind(0.d0))
  Enddo
  Do nth = 1,10
     Do I = 1,Ndim
        Do I1 = 1,Ndim
           Tmp(I1,I) = cmplx(4.d0*(ranf_wrap() -0.5D0) , 2.d0*(ranf_wrap() -0.5D0), Kind=Kind(0.d0))
        Enddo
     Enddo
     Call MMULT(Tmp1,Tmp,udvr%U)
     udvr%U = Tmp1
     call udvr%decompose()
     Call Mmult(Tmp1, Tmp,TmpR)
     TmpR = Tmp1
     Call Check_udv(TmpR, udvr, diff)
     Write(6,*) 'Right: ', Diff
     
     Call MMULT(Tmp1,Tmp,udvl%U)
     udvl%U = Tmp1
     call udvl%decompose()
     Call Mmult(Tmp1, Tmp,TmpL)
     TmpL = Tmp1
     Call Check_udv(TmpL, udvl, diff)
     Write(6,*) 'Left: ', Diff
  Enddo
  !call udvr%print()

end Program test


Subroutine Check_udv(Tmp, UDV, Diff)
  
  Use UDV_State_mod
  Implicit none

  Type (UDV_state), Intent(IN) :: udv
  Complex (Kind=Kind(0.d0)), allocatable  :: Tmp(:,:)
  Real (Kind=Kind(0.d0)), INTENT(out) :: Diff

  Integer :: I, I1, M
  Real (Kind = Kind(0.d0)) :: X
  Complex (Kind = Kind(0.d0)) :: Z 
  Complex (Kind = Kind(0.d0)) :: DM   ! scale M; getscale also works in the STABLOG build
  
  If (udv%side == "r" .or. udv%side == "R" ) Then
     Diff = 0.d0
     Do I = 1, size(udv%U,1)
        Do I1 = 1,  size(udv%V,2)
           Z = cmplx(0.d0,0.d0,Kind = Kind(0.d0))
           Do M = 1, size(udv%U,2)
              call udv%getscale(DM, M)
              Z =  Z  +    udv%U(I,M)*DM * udv%V(M,I1)
           Enddo
           Z = Z - TMP(I,I1)
           X = abs(Z)
           If ( X > Diff ) Diff = X
        Enddo
     Enddo
  endif
  If (udv%side == "l" .or. udv%side == "L" ) Then
     Diff = 0.d0
     Do I = 1, size(udv%U,1)
        Do I1 = 1,  size(udv%V,2)
           Z = cmplx(0.d0,0.d0,Kind = Kind(0.d0))
           Do M = 1, size(udv%U,2)
              call udv%getscale(DM, M)
              Z =  Z  +    udv%U(I,M)*DM * conjg(udv%V(I1,M))
           Enddo
           Z = Z - TMP(I,I1)
           X = abs(Z)
           If ( X > Diff ) Diff = X
        Enddo
     Enddo
  endif


End Subroutine Check_Udv
