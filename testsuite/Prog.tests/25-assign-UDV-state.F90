Program  test

  Use Random_Wrap
  Use UDV_State_mod
  Use MyMats 

  Implicit none

  Interface
     Subroutine Check_assign_UDV_state(udv1, udv2, zero_in)
       
       Use UDV_State_mod
       Implicit none
       
       CLASS(UDV_state), Intent(IN) :: udv1, udv2
       Complex (Kind=Kind(0.d0)), allocatable  :: Tmp(:,:)
       Real (Kind=Kind(0.d0)), INTENT(IN), optional :: zero_in
     End Subroutine Check_assign_UDV_state
  end Interface
  CLASS(UDV_state), allocatable :: udv1(:), udv2(:)
  Integer :: Ndim = 16, Seed_Vec(1), I, I1
  Complex(Kind =Kind(0.d0)), allocatable :: TMP(:,:,:)
  
  Integer :: N_udv = 10, n

  allocate ( TMP(Ndim,Ndim,N_udv) )
  allocate ( udv1(N_udv), udv2(N_udv) )

  Seed_vec(1) = 4782347
  Call Ranset(Seed_vec)
  
  do n=1, N_udv
     call  udv1(n)%alloc(Ndim)
     call  udv2(n)%alloc(Ndim)
     Do I = 1,Ndim
        Do I1 = 1,Ndim
           Tmp(I1,I,n) = cmplx(4.d0*(ranf_wrap() -0.5D0) , 2.d0*(ranf_wrap() -0.5D0), Kind=Kind(0.d0))
        Enddo
     Enddo
     
     if ( n <= N_udv/2 ) then
        call udv1(n)%reset("r")
     else
        call udv1(n)%reset("l")
     endif
     udv1(n)%U(:,:) = Tmp(:,:,n)
     call udv1(n)%decompose()
  enddo
  
  do n=1, N_udv
     call udv2(n)%reset("r")
     udv2(n) = udv1(n)
  enddo
  
  do n=1, N_udv
     call Check_assign_UDV_state(udv1(n), udv2(n))
  enddo
  
  print*, "success"

end Program test


Subroutine Check_assign_UDV_state(udv1, udv2, zero_in)
  
  Use UDV_State_mod
  Implicit none

  CLASS(UDV_state), Intent(IN) :: udv1, udv2
  Real (Kind=Kind(0.d0)), INTENT(IN), optional :: zero_in
  
  Integer :: I, I1
  Real (Kind=Kind(0.d0)) :: zero
  Complex (Kind=Kind(0.d0)) :: D1, D2   ! scales; getscale also works in the STABLOG build
  
  if ( present(zero_in) ) then
     zero = zero_in
  else
     zero = 1.D-10
  endif
  
  if ( udv1%ndim   .ne. udv2%ndim   ) stop 1 !"ndim does not match"
  if ( udv1%n_part .ne. udv2%n_part ) stop 2 !"n_part does not match"
  if ( udv1%side   .ne. udv2%side   ) stop 3 !"side does not match"
  
  Do I = 1,udv1%Ndim
     call udv1%getscale(D1, I)
     call udv2%getscale(D2, I)
     if ( abs( D1 - D2 ) > zero ) stop 4 !"difference too big"
     Do I1 = 1,udv1%Ndim
        if ( abs( udv1%U(I1,I) - udv2%U(I1,I) ) > zero ) stop 5 !"difference too big"
        if ( abs( udv1%V(I1,I) - udv2%V(I1,I) ) > zero ) stop 6 !"difference too big"
     Enddo
  Enddo

End Subroutine Check_assign_UDV_state
