! compile with
! gfortran -I ../../Libraries/libmscbdecomp/ -L ../../Libraries/libmscbdecomp/ 1-ZeroDiag-lmult.F90 -lmscbdecomp

subroutine exectest(op_T, mys)
  Use mscbOpT_mod
  Use Exponentials_mod
  Use colorvertex_mod
  Use Operator_mod
  Use graphdata_mod
  implicit none
        Type(Operator), intent(in) :: Op_T
        class(CmplxEulermscbOpT), allocatable :: ee
        integer :: ndim
        real(kind=kind(0.D0)) :: weight, sumdiag, sumoff
        real(kind=kind(0.D0)), allocatable :: mys(:)
        COMPLEX (KIND=kind(0.D0)), DIMENSION(:, :), allocatable :: mat
        integer :: i

        ndim = size(op_t%O, 1)
        allocate (ee, mat(ndim, ndim))
        
        weight = 1.0

        ! initialize as identity matrix
        mat = 0
        do i = 1, ndim
            mat(i,i) = 1
        enddo
        call ee%init(Op_T)
        
        call ee%lmult(mat)
        call ee%lmult(mat)
        call ee%lmult(mat)
        call ee%lmultinv(mat)
        call ee%lmultinv(mat)
        call ee%lmultinv(mat)

        ! test for Trace(mat) = ndim
        sumdiag = 0
        sumoff = 0
        do i = 1, ndim
            sumdiag = sumdiag + DBLE(mat(i,i))
        enddo
        do i = 1, ndim-3
            sumoff = sumoff + DBLE(mat(i,i+2))
        enddo
        write (*,*) "lmult: ", sumoff, sumdiag
        if (abs(sumdiag - ndim) > ndim*1E-15) then
        ERROR STOP 2
        endif
        if (abs(sumoff) > 1E-14) then !FIXME: this limit is a bit scale less...
        ERROR STOP 4
        endif
        

        ! initialize as identity matrix
        mat = 0
        do i = 1, ndim
            mat(i,i) = 1
        enddo
        

        call ee%rmult(mat)
        call ee%rmult(mat)
        call ee%rmult(mat)
        call ee%rmultinv(mat)
        call ee%rmultinv(mat)
        call ee%rmultinv(mat)
        ! test for Trace(mat) = ndim
        sumdiag = 0
        sumoff = 0
        do i = 1, ndim
            sumdiag = sumdiag + DBLE(mat(i,i))
        enddo
        do i = 1, ndim-3
            sumoff = sumoff + DBLE(mat(i,i+2))
        enddo
        write(*,*) "rmult  ", sumoff, sumdiag
        if (abs(sumdiag - ndim) > ndim*1E-15) then
        ERROR STOP 3
        endif
        if (abs(sumoff) > 1E-14) then !FIXME: this limit is a bit scale less...
        ERROR STOP 6
        endif

        call ee%adjointaction(mat)
        ! test for Trace(mat) = ndim
        sumdiag = 0
        sumoff = 0
        do i = 1, ndim
            sumdiag = sumdiag + DBLE(mat(i,i))
        enddo
        do i = 1, ndim-3
            sumoff = sumoff + DBLE(mat(i,i+2))
        enddo
        write (*,*) "adjoint: ", sumoff, sumdiag
        if (abs(sumdiag - ndim) > ndim*1E-15) then
        ERROR STOP 7
        endif
        if (abs(sumoff) > 1E-14) then !FIXME: this limit is a bit scale less...
        ERROR STOP 14
        endif
        
        call ee%dealloc()
        
end subroutine 

Program EulerExpTest

  Use Exponentials_mod
  Use MvG_mod
  Use Operator_mod
  Use Lattices_v3
  implicit none 
  
  interface
  subroutine exectest(op_T, mys)
  Use Exponentials_mod
  Use colorvertex_mod
  Use Operator_mod
  implicit none
        Type(Operator), intent(in) :: Op_T
        real(kind=kind(0.D0)), allocatable :: mys(:)
    end subroutine
  end interface
        type(EulerExp) :: ee
        Type(Operator) :: Op_T
        integer :: ndim, L1, L2, Ix, Iy
        Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
        
        Type (Lattice) :: Lat
        Type (Unit_cell) :: Lat_unit
        real(kind=kind(0.D0)) :: weight, sumdiag, sumoff, ham_t
        real(kind=kind(0.D0)), allocatable :: mys(:)
        COMPLEX (KIND=kind(0.D0)), DIMENSION(:,:), allocatable :: input
        integer :: i, j

        L1 = 8
        L2 = 8
        ham_t = 1.0
        

          Lat_unit%Norb    = 1
          Lat_unit%N_coord = 2
          allocate(Lat_unit%Orb_pos_p(Lat_unit%Norb,2))
          Lat_unit%Orb_pos_p(1, :) = [0.d0, 0.d0]

          a1_p(1) =  1.0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0  ; a2_p(2) =  1.d0
          L1_p    =  dble(L1)*a1_p
          L2_p    =  dble(L2)*a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Lat )
         
          Ndim = Lat%N*Lat_unit%Norb
          write (*,*) Ndim
          Call Op_make(Op_T, Ndim)
          Do I = 1,Lat%N
            Ix = Lat%nnlist(I,1,0)
            Op_T%O(I,  Ix) = cmplx(-Ham_T,    0.d0, kind(0.D0))
            Op_T%O(Ix, I ) = cmplx(-Ham_T,    0.d0, kind(0.D0))
            If ( L2 > 1 ) then
                Iy = Lat%nnlist(I,0,1)
                Op_T%O(I,  Iy) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                Op_T%O(Iy, I ) = cmplx(-Ham_T,    0.d0, kind(0.D0))
            endif
            Op_T%O(I,  I ) = 0.D0
            Op_T%P(i) = i
         Enddo
         Op_T%g      = -0.1
         Op_T%alpha  =  cmplx(0.d0,0.d0, kind(0.D0))
         Call Op_set(Op_T)


        allocate(mys(ndim) )

        ! Start with zero diagonals
        mys = 0.0 ! initialize chemical potential to zero
        call exectest(op_t, mys)

        ! Now test homogeneous exponentials
        mys = 0.5
        call exectest(op_t, mys)
        
        ! Now test general exponentials
       do i = 1, ndim
         mys(i) = 0.1*i
       enddo
       call exectest(op_t, mys)
       call Op_clear(Op_T, ndim)
       deallocate(mys, Lat_unit%Orb_pos_p, Lat%L2_p, Lat%L1_p, Lat%a1_p, Lat%a2_p, Lat%b1_p, Lat%b2_p, Lat%BZ1_p, Lat%BZ2_p, Lat%b1_perp_p, Lat%b2_perp_p, Lat%List, Lat%Invlist, Lat%Listk, Lat%Invlistk, Lat%nnlist, Lat%imj)
       write (*,*) "success"
end Program EulerExpTest
