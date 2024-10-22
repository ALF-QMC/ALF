!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module provides a set of predefined trial wave functions.
!>
!
!--------------------------------------------------------------------

module Predefined_Trial

   use runtime_error_mod
   use Lattices_v3
   use Operator_mod
   use WaveFunction_mod
   use MyMats
   use Predefined_Hoppings
   use Random_wrap
   use iso_fortran_env, only: output_unit, error_unit

   implicit none

contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!>    Sets the trial wave function corresponding to the solution of the non-interacting
!>    tight binding Hamiltonian on the given lattice. Twisted boundary conditions (Phi_X=0.01)
!>    are implemented so as to generate a non-degenerate trial wave functions.
!> @param [in]  Lattice_type
!>    Character(64)
!> \verbatim
!>    Square,  Honeycomb, Pi_Flux
!> \endverbatim
!> @param [in]  Latt_unit
!>    Type(Unit_cell)
!> \verbatim
!>     Contains number of orbitals per unit cell and positions, as well as coordination number
!> \endverbatim
!> @param [in]  Ndim
!>    Integer
!> \verbatim
!>     Number of orbitals
!> \endverbatim
!> @param [in]  List, Invlist
!>    Integer(:,:)
!> \verbatim
!>    List(I=1.. Ndim,1)    =   Unit cell of site I
!>    List(I=1.. Ndim,2)    =   Orbital index  of site I
!>    Invlist(Unit_cell,Orbital) = site I
!> \endverbatim
!> @param [in]    Latt
!>    Type(Lattice)
!> \verbatim
!>    The Lattice
!> \endverbatim
!> @param [in]  N_part
!>    Integer
!> \verbatim
!>    Particle number for each flavor
!> \endverbatim
!> @param [in]  N_FL
!>    Integer
!> \verbatim
!>    Flavor
!> \endverbatim
!> @param [out]  WF_L, WF_R
!>    Type(Wavefunction)(N_FL)
!> \verbatim
!>    Wavefunction
!>    Also sets the degeneracy:  E(N_part + 1) - E(N_part). Energy eigenvalues are ordered in ascending order.
!> \endverbatim
!>
!------------------------------------------------------------------
   subroutine Predefined_TrialWaveFunction(Lattice_type, Ndim, List, Invlist, Latt, Latt_unit, &
        &                                  N_part, alpha, N_FL, N_slat, WF_L, WF_R)

      implicit none
      character(len=64), intent(IN)                :: Lattice_type
      integer, intent(IN)                           :: Ndim, N_FL, N_slat, N_part
      integer, intent(IN), dimension(:, :)           :: List, Invlist
      real(Kind=kind(0.d0)), intent(in)             :: alpha
      type(Lattice), intent(in)                   :: Latt
      type(Unit_cell), intent(in)                   :: Latt_Unit
      type(WaveFunction), intent(out), dimension(:, :), allocatable :: WF_L, WF_R

      type(operator), dimension(:, :), allocatable  :: OP_tmp
      type(Hopping_Matrix_type), allocatable       :: Hopping_Matrix_tmp(:)
      real(Kind=kind(0.d0))                        :: Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, Dimer, mass
      real(Kind=kind(0.d0))                        :: ham_U, sgn_i, sgn_updn, rmu
      logical                                       :: Checkerboard, Symm, Kekule_Trial, hatree_fock, l_cmplx_trial

      type(Lattice)                                :: Latt_Kekule
      real(Kind=kind(0.d0))  :: A1_p(2), A2_p(2), L1_p(2), L2_p(2), x_p(2), x1_p(2), hop(3), del_p(2)
      real(Kind=kind(0.d0))  :: delta = 0.01, Ham_T1, Ham_T2, Ham_Tperp

      integer :: N, nf, I, I1, I2, nc, nc1, IK_u, I_u, J1, lp, J, N_Phi, ns, no, k1, N_part_tot, Np_arr(2), i0, j0
      logical :: Test = .false., Bulk = .true.
      complex(Kind=kind(0.d0)) :: Z_norm

      real(Kind=kind(0.d0)), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
             & ham_tx_vec(:), ham_ty_vec(:), Ham_T2_vec(:), Ham_lambda_vec(:)
      integer, allocatable ::   N_Phi_vec(:)
      real(Kind=kind(0.d0)), allocatable :: eig_sort_arr(:,:)

      allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
           &    ham_tx_vec(N_FL), ham_ty_vec(N_FL), Ham_lambda_vec(N_FL), N_Phi_vec(N_FL))

      Checkerboard = .false.
      Kekule_Trial = .false.
      hatree_fock = .false.
      l_cmplx_trial = .false.
      Symm = .false.

      N_Phi = 0
      Phi_X = 0.d0
      Phi_Y = 0.d0
      Bulk = .false.
      Ham_T = 1.d0
      Ham_T2 = 0.d0
      Ham_Tperp = 0.d0
      Ham_Chem = 0.d0
      Dtau = 1.d0

      N_Phi_vec = N_Phi
      Phi_X_vec = Phi_X
      Phi_Y_vec = Phi_Y
      Ham_T_vec = Ham_T
      Ham_Tperp_vec = Ham_Tperp
      Ham_T2_vec = Ham_T2
      Ham_Chem_vec = Ham_Chem
      Ham_lambda_vec = 0.d0

      select case (Lattice_type)

      case ('N_leg_ladder')
         if ( l_cmplx_trial ) then
            ham_tx_vec(1) = ham_t; 
            ham_tx_vec(2) = ham_t*alpha; 
            ham_ty_vec(1) = ham_t*alpha; 
            ham_ty_vec(2) = ham_t; 
            Phi_X_vec = 0.002
            call set_hopping_parameters_n_ladder_anisotropic(Hopping_Matrix_tmp, ham_tx_vec, ham_ty_vec, Ham_Chem_vec, &
                   & Phi_X_vec, Phi_Y_vec, Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit)
         else

            allocate (op_tmp(1, n_fl))
            do n = 1, n_fl
               call op_make(op_tmp(1, n), ndim)
               do i = 1, ndim
                  op_tmp(1, n)%P(i) = i
               end do
               op_tmp(1, n)%g = cmplx(1.d0, 0.d0, kind(0.d0))
               op_tmp(1, n)%alpha = cmplx(0.d0, 0.d0, kind(0.d0))
            end do

            do i = 1, latt%n
               do no = 1, Latt_unit%Norb - 1
                  delta = 0.005d0*ranf_wrap()
                  ham_t = 1.d0 + delta

                  i1 = invlist(i, no)
                  j1 = invlist(latt%nnlist(i, 1, 0), no)
                  k1 = invlist(i, no + 1)

                  ham_tx_vec(1) = ham_t; 
                  ham_tx_vec(2) = ham_t*alpha; 
                  ham_ty_vec(1) = ham_t*alpha; 
                  ham_ty_vec(2) = ham_t; 
                  do n = 1, n_fl
                     op_tmp(1, n)%o(i1, j1) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
                     op_tmp(1, n)%o(j1, i1) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
                     op_tmp(1, n)%o(i1, k1) = cmplx(-ham_ty_vec(n), 0.d0, kind(0.d0))
                     op_tmp(1, n)%o(k1, i1) = cmplx(-ham_ty_vec(n), 0.d0, kind(0.d0))
                  end do
               end do
               delta = 0.005d0*ranf_wrap()
               ham_t = 1.d0 + delta
               no = Latt_unit%Norb
               i1 = invlist(i, no)
               j1 = invlist(latt%nnlist(i, 1, 0), no)
               ham_tx_vec(1) = ham_t; 
               ham_tx_vec(2) = ham_t*alpha; 
               do n = 1, n_fl
                  op_tmp(1, n)%o(i1, j1) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
                  op_tmp(1, n)%o(j1, i1) = cmplx(-ham_tx_vec(n), 0.d0, kind(0.d0))
               end do
            end do

            do n = 1, n_fl
               call op_set(op_tmp(1, n))
            end do

        endif

      case default
         write (error_unit, *) 'No predefined trial wave function for this lattice.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end select

      if (l_cmplx_trial )   &
      &     call Predefined_Hoppings_set_OPT(Hopping_Matrix_tmp, List, Invlist, Latt, Latt_unit, Dtau, Checkerboard, Symm, OP_tmp)

      do nf = 1, N_FL
         call Diag(op_tmp(1, nf)%o, op_tmp(1, nf)%u, op_tmp(1, nf)%e)
      end do

      !! sort eigen state
      allocate(eig_sort_arr(3,n_fl*Ndim))
      eig_sort_arr(:,:) = 0.d0
      nc = 0
      do nf = 1, N_fl
          do i1 = 1, ndim
              nc = nc + 1
              eig_sort_arr(1,nc) = op_tmp(1, nf)%e(i1)
              eig_sort_arr(1+nf,nc) = i1
          enddo
      enddo
      call s_heapsort(n_fl*Ndim,3,eig_sort_arr)
      
      N_part_tot = 2*N_part
      i0 = 0; j0 = 0
      do nc = 1, N_part_tot
          i1 = int(eig_sort_arr(nc,2))
          if (i1 .ne. 0) then
              i0 = i0 + 1
          else
              j0 = j0 + 1
          endif
      enddo
      np_arr(1) = i0
      np_arr(2) = j0

      !! allocate wf
      allocate (wf_l(n_fl, n_slat), wf_r(n_fl, n_slat))
      do ns = 1, n_slat
         do n=1,N_FL
            call wf_alloc(wf_l(n, ns), ndim, np_arr(n))
            call wf_alloc(wf_r(n, ns), ndim, np_arr(n))
         enddo
      end do

      do nf = 1, N_FL
         do i2 = 1, np_arr(nf)
            do i1 = 1, ndim
               wf_l(nf, 1)%p(i1,i2) = op_tmp(1, nf)%u(i1,i2)
               wf_r(nf, 1)%p(i1,i2) = op_tmp(1, nf)%u(i1,i2)
            end do
         end do
         wf_l(nf, 1)%degen = op_tmp(1, nf)%e(np_arr(nf) + 1) - op_tmp(1, nf)%e(np_arr(nf))
         wf_r(nf, 1)%degen = op_tmp(1, nf)%e(np_arr(nf) + 1) - op_tmp(1, nf)%e(np_arr(nf))
      end do

      do nf = 1, N_FL
         call WF_overlap(WF_L(nf, 1), WF_R(nf, 1), Z_norm)
         WF_L(nf, 1)%P(:, :) = WF_L(nf, 1)%P(:, :)/sqrt(dble(N_slat))
         WF_R(nf, 1)%P(:, :) = WF_R(nf, 1)%P(:, :)/sqrt(dble(N_slat))
      end do

      do nf = 1, N_FL
      do ns = 2, n_slat
         WF_L(nf, ns)%Degen = WF_L(nf, 1)%Degen
         WF_R(nf, ns)%Degen = WF_R(nf, 1)%Degen
         WF_L(nf, ns)%P(:, :) = WF_L(nf, 1)%P(:, :)
         WF_R(nf, ns)%P(:, :) = WF_R(nf, 1)%P(:, :)
      end do
      end do

      do nf = 1, N_FL
         call Op_clear(OP_tmp(1, nf), Ndim)
      end do
      deallocate (OP_tmp)
      call Predefined_hoppings_clear(Hopping_Matrix_tmp)

      deallocate (Ham_T_vec, Ham_Tperp_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, N_Phi_vec)
      deallocate (ham_tx_vec, ham_ty_vec)

   end subroutine Predefined_TrialWaveFunction
  
   subroutine s_heapsort(n,width,brxyz)
! modified from BSTATE code by Xiaoyan Xu 2013.9
!########################################################################
! Former discription
!C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
!C-----CALL HPSORT(KG,IG1,IG2,IG3,GX,GY,GZ,GR)
!      subroutine hpsort(kn,n,bxyz,br)
!c                           @(#)hpsort.f 9.1 97/05/08 14:48:33 
!C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!C     THIS SUBROUTINE SORTS THE ELEMENTS OF ARRAY BR IN ASCENDING
!C     ORDER. THE CODE IS BASED UPON HEAPSORT ALGORITHM WHICH HAS
!C     N*LOG(N) SCALING BEHAVIOUR, PARTICULARLY FAVOURABLE FOR LARGE
!C     VECTORS BR
!C                              STEFAN BL"UGEL, ISSP, NOV 1989
!C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!########################################################################


! description
! This subroutine sorts the elements of array brxyz(width,n) in ascending order
! based on brxyz(1,:).

      IMPLICIT LOGICAL(A-Z)
      integer n, width
      real(kind=kind(0.d0)), dimension(width, n) :: brxyz
      integer :: l, i, ii, iheap
      real(kind=kind(0.d0)) :: brr
      real(kind=kind(0.d0)), dimension(width-1) :: bxyz
!     integer :: kn, n,l,i,ii,iheap
!     real*8    bxyz(kn,3),br(*),brr,bxx,byy,bzz
!C
!C
!C=====> BUILD-UP OF THE HEAP
!C
!C-----> LOOP OVER ALL HIERACHICAL LEVELS OF THE HEAP TREE
!C
      DO 10 L = N/2 , 1 , -1
!         brr = br(l)
         brr = brxyz(1,l)
         bxyz = brxyz(2:width,l)
!         bxx = bxyz(l,1)
!         byy = bxyz(l,2)
!         bzz = bxyz(l,3)
         I   = L
         II  = L + L
!C
!C-----> GO DOWN ALL THE HIERACHICAL LEVELS OF THE HEAP TREE
!C
  20    IF ( II .LE. N ) THEN
!C
!C-----> COMPARE NEIGHBOURING ELEMENTS
           IF ( II .LT. N ) THEN
              IF ( BRXYZ(1,II) .LT. BRXYZ(1,II+1) ) II = II + 1
           ENDIF
!C
!C-----> COMPARE THE ELEMENTS OF TWO HIRACHICAL LEVELS
!C       PROMOTE LARGER ONE, DEMOTE SMALLER ONE
           IF ( BRR .LT. BRXYZ(1,II) ) THEN
              BRXYZ(1,I) = BRXYZ(1,II)
              brxyz(2:width,i) = brxyz(2:width,ii)
!              bxyz(i,1) = bxyz(ii,1)
!              bxyz(i,2) = bxyz(ii,2)
!              bxyz(i,3) = bxyz(ii,3)
              I     = II
              II    = II + II
           ELSE
!C
!C-----> THIS PART OF THE TREE IS ORDERED , STOP BY PUTTING II=N+1
              II    = N + 1
           END IF
                                                GO TO 20
        END IF
!C-----> PUT ELEMENTS IN THE PROPER SLOT
        brxyz(1,i) = brr
        brxyz(2:width,i) = bxyz(:)
!        bxyz(i,1) = bxx
!        bxyz(i,2) = byy
!        bxyz(i,3) = bzz
  10  continue
!C
!C=====> NOW COLLECT ALL ELEMENTS FROM THE HEAP
!C
      DO 30 IHEAP = N , 1 , -1
!C
         brr = brxyz(1,iheap)
         bxyz = brxyz(2:width,iheap)
!         bxx = bxyz(iheap,1)
!         byy = bxyz(iheap,2)
!         bzz = bxyz(iheap,3)
!C
!C-----> THE FIRST ELEMENT IS ALWAYS THE LARGEST
              brxyz(1,iheap) = brxyz(1,1)
              brxyz(2:width,iheap) = brxyz(2:width,1)
!              bxyz(iheap,1) = bxyz(1,1)
!              bxyz(iheap,2) = bxyz(1,2)
!              bxyz(iheap,3) = bxyz(1,3)
!C-----> NOW LARGEST ELEMENT OF ALL BR(I) WITH 1<=I<=IHEAP IS STORED
!C
         I   = 1
         II  = 2
!C
!C-----> NOW GENERATE LARGEST ELEMENT OF BR(I) WITH 1<=I<=IHEAP-1
!C
  40    IF ( II .LE. IHEAP - 1 ) THEN
!C
!C-----> COMPARE NEIGHBOURING ELEMENTS
           IF ( II .LT. IHEAP - 1 ) THEN
              IF ( BRXYZ(1,II) .LT. BRXYZ(1,II+1) ) II = II + 1
           ENDIF
!C
!C-----> PROMOTE EACH ELEMENT OF THE TWIG OF BR UNTIL BRR > BR(I)
           IF ( BRR .LT. BRXYZ(1,II) ) THEN
              brxyz(1,i) = brxyz(1,ii)
              brxyz(2:width,i) = brxyz(2:width,ii)
!              bxyz(i,1) = bxyz(ii,1)
!              bxyz(i,2) = bxyz(ii,2)
!              bxyz(i,3) = bxyz(ii,3)
              I     = II
              II    = II + II
           ELSE
!C
!C-----> THIS PART OF THE TREE IS PROMOTED , STOP BY PUTTING II=IHEAP+1
              II    = IHEAP + 1
           END IF
                                                GO TO 40
        END IF
!C-----> PUT ELEMENTS IN THE PROPER SLOT
        brxyz(1,i) = brr
        brxyz(2:width,i) = bxyz(:)
!        bxyz(i,1) = bxx
!        bxyz(i,2) = byy
!        bxyz(i,3) = bzz
  30  continue
  end subroutine s_heapsort

end module Predefined_Trial
