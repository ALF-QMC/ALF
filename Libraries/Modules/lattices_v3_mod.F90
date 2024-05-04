!  Copyright (C) 2018-2023 The ALF project
!
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.


     Module Lattices_v3

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> This module generates one and two dimensional Bravais lattices and the unit cell
!
!--------------------------------------------------------------------
         Use Matrix
         Use runtime_error_mod
         use iso_fortran_env, only: output_unit, error_unit
         implicit none
         private
         public :: Unit_cell, Lattice, Make_lattice, Fourier_K_to_R, &
                   Fourier_R_to_K, Inv_K, Inv_R, Iscalar, npbc, Xnorm, Print_latt, &
                   clear_lattice

         Type Unit_cell
            !> Number of orbitals
            Integer :: Norb
            !> Coordination number
            Integer :: N_coord
            !> Positions of orbitals: Orb_pos(1..Norb,2)
            Real (Kind=Kind(0.d0)), pointer :: Orb_pos_p(:,:)
         end type Unit_cell

         Type Lattice
            Integer          :: N, Ns
            Integer, pointer :: list(:,:), invlist(:,:), nnlist(:,:,:), nnlistk(:,:,:), listk(:,:), &
                 &              invlistk(:,:),  imj(:,:)
            Real (Kind=Kind(0.d0)), pointer :: a1_p(:), a2_p(:), a3_p(:), b1_p(:), b2_p(:), b3_p(:), &
                 &                    BZ1_p(:), BZ2_p(:), BZ3_p(:), &
                 &                    L1_p(:), L2_p(:), L3_p(:), b1_perp_p(:), b2_perp_p(:), b3_perp_p(:)
         end Type Lattice

         Interface Iscalar
            module procedure Iscalar_II, Iscalar_IR, Iscalar_RR
         end Interface
         Interface npbc
            module procedure npbc_I2D, npbc_R2D, npbc_I3D, npbc_R3D, npbc_R_B2D, npbc_R_B3D
         end Interface
         Interface Inv_K
            module procedure Inv_K2D, Inv_K3D
         end Interface
         Interface Inv_R
            module procedure Inv_R2D, Inv_R3D
         end Interface
         Interface Xnorm
            module procedure Xnorm_I, Xnorm_R
         end Interface
         Interface Fourier_K_to_R
            module procedure FT_K_to_R, FT_K_to_R_Mat, FT_K_to_R_C,  FT_K_to_R_Mat_C
         end Interface
         Interface Fourier_R_to_K
            module procedure FT_R_to_K, FT_R_to_K_mat, FT_R_to_K_C
         end Interface
         Interface Make_lattice
            module procedure Make_lattice3D, Make_lattice2D
         end Interface

       Contains

       subroutine Make_lattice3D(L1_p, L2_p, L3_p, a1_p, a2_p, a3_p, Latt)

         ! This is for a general tilted square lattice defined by the vector a1, a2, a3
         ! L1_p,  L2_p, L3_p define cluster topology. ( Tilted etc.)
         ! L1_p = n*a1_p + m *a2_p + o * a3_p


         Implicit none

         Real (Kind=Kind(0.d0)),  dimension(:) :: L1_p, L2_p, L3_p, a1_p, a2_p, a3_p
         Type (Lattice) :: Latt

         Real (Kind=Kind(0.d0)), dimension(:), allocatable :: xk_p, b1_p, b2_p, b3_p, BZ1_p, BZ2_p, BZ3_p, b_p
         Real (Kind=Kind(0.d0)), dimension(:), allocatable :: x_p, x1_p, a_p,d_p
         Real (Kind=Kind(0.d0)), allocatable :: Mat(:,:), Mat_inv(:,:)

         Integer :: ndim, L, L1, nc, i, i1,i2, i3, L_f, LQ, o, n,m, nd1,nd2,nd3,nr, nnr1, nnr2, nnr3, nnr, nr1, imj_1, imj_2, imj_3
         Integer :: imj
         Real    (Kind=Kind(0.d0)) :: Zero,pi, X

         ndim = size(L1_p)
         allocate (Latt%L3_p(ndim), Latt%L2_p(ndim), Latt%L1_p(ndim), Latt%a1_p(ndim) , Latt%a2_p(ndim), &
              &    Latt%a3_p(ndim), Latt%b1_p(ndim), Latt%b2_p(ndim), Latt%b3_p(ndim) Latt%BZ1_p(ndim), &
              &    Latt%BZ2_p(ndim), Latt%BZ3_p(ndim) )
         allocate (Latt%b1_perp_p(ndim), Latt%b2_perp_p(ndim), Latt%b3_perp_p(ndim) )
         Zero = 1.D-5
         Latt%L1_p = L1_p
         Latt%L2_p = L2_p
         Latt%L3_p = L3_p
         Latt%a1_p = a1_p
         Latt%a2_p = a2_p
         Latt%a3_p = a3_p


         !Compute the Reciprocal lattice vectors.
         Allocate ( b1_p(ndim), b2_p(ndim), b3_p(ndim), xk_p(ndim), b_p(ndim) )
         Allocate ( BZ1_p(ndim), BZ2_p(ndim),  BZ3_p(ndim))
         Allocate ( x_p(ndim),  x1_p(ndim), d_p(ndim),  a_p(ndim) )


         pi   = acos(-1.d0)

         ! Setup the 3X3 matrix to determine  BZ1_p, BZ2_p, BZ3_p
         Allocate ( Mat(3 , 3), Mat_inv( 3 , 3 ) )
         Mat(1,1) = dble(a1_p(1))
         Mat(1,2) = dble(a1_p(2))
         Mat(1,3) = dble(a1_p(3))
         Mat(2,1) = dble(a2_p(1))
         Mat(2,2) = dble(a2_p(2))
         Mat(2,3) = dble(a2_p(3))
         Mat(3,1) = dble(a3_p(1))
         Mat(3,2) = dble(a3_p(2))
         Mat(3,3) = dble(a3_p(3))
         X = Mat(1,1)*(Mat(2,2)*Mat(3,3) - Mat(2,3)*Mat(3,2)) -&
            &(Mat(1,2)*(Mat(2,1)*Mat(3,3) - Mat(2,3)*Mat(3,1)))+&
            &(Mat(1,3)*(Mat(2,1)*Mat(3,2) - Mat(2,2)*Mat(3,1)))
         Mat_inv(1,1) =  (Mat(2,2)*Mat(3,3) - Mat(2,3)*Mat(3,2))/X
         Mat_inv(2,2) =  (Mat(1,1)*Mat(3,3) - Mat(1,3)*Mat(3,1))/X
         Mat_inv(3,3) =  (Mat(1,1)*Mat(2,2) - Mat(1,2)*Mat(2,1))/X
         Mat_inv(1,2) = -(Mat(1,2)*Mat(3,3) - Mat(1,3)*Mat(3,2))/X
         Mat_inv(2,1) = -(Mat(2,1)*Mat(3,3) - Mat(2,3)*Mat(3,1))/X
         Mat_inv(2,3) = -(Mat(1,1)*Mat(2,3) - Mat(1,3)*Mat(2,1))/X
         Mat_inv(3,2) = -(Mat(1,1)*Mat(3,2) - Mat(1,2)*Mat(3,1))/X
         Mat_inv(3,1) =  (Mat(2,1)*Mat(3,2) - Mat(2,2)*Mat(3,1))/X
         Mat_inv(1,3) =  (Mat(1,2)*Mat(2,3) - Mat(1,3)*Mat(2,2))/X
         BZ1_p(1)      = 2.d0*pi*Mat_inv(1,1)
         BZ1_p(2)      = 2.d0*pi*Mat_inv(2,1)
         BZ1_p(3)      = 2.d0*pi*Mat_inv(3,1)
         BZ2_p(1)      = 2.d0*pi*Mat_inv(1,2)
         BZ2_p(2)      = 2.d0*pi*Mat_inv(2,2)
         BZ2_p(3)      = 2.d0*pi*Mat_inv(3,2)
         BZ3_p(1)      = 2.d0*pi*Mat_inv(1,3)
         BZ3_p(2)      = 2.d0*pi*Mat_inv(2,3)
         BZ3_p(3)      = 2.d0*pi*Mat_inv(3,3)
         Latt%BZ1_p = BZ1_p
         Latt%BZ2_p = BZ2_p
         Latt%BZ3_p = BZ3_p




         ! K-space Quantization  from periodicity in L1_p and L2_p
         X =  2.d0*pi / ( Iscalar(BZ1_p,L1_p) * ( Iscalar(BZ2_p,L2_p) * Iscalar(BZ3_p,L3_p) -   &
               &          Iscalar(BZ2_p,L1_p) * Iscalar(BZ1_p,L2_p) ) - &
               &          Iscalar(BZ1_p,L2_p) * ( Iscalar(BZ2_p,L1_p) * Iscalar(BZ3_p,L3_p) -   &
               &          Iscalar(BZ2_p,L3_p) * Iscalar(BZ3_p,L1_p) ) + &
               &          Iscalar(BZ1_p,L3_p) * ( Iscalar(BZ2_p,L1_p) * Iscalar(BZ3_p,L2_p) -   &
               &          Iscalar(BZ2_p,L2_p) * Iscalar(BZ3_p,L1_p) ) )
         X = abs(X)
         !(L2g3 L3g2 n - L2g2 L3g3 n ) g1 + (-L2g3 L3g1 n + 
         !L2g1 L3g3 n ) g2 + (L2g2 L3g1 n - L2g1 L3g2 n ) g3
         b3_p = X*( (Iscalar(BZ2_p,L2_p) * Iscalar(BZ3_p,L3_p)-Iscalar(BZ3_p,L2_p) * Iscalar(BZ2_p,L3_p)) * BZ1_p &
               & +  (Iscalar(BZ3_p,L2_p) * Iscalar(BZ1_p,L3_p)-Iscalar(BZ1_p,L2_p) * Iscalar(BZ3_p,L3_p)) * BZ2_p &
               & +  (Iscalar(BZ1_p,L2_p) * Iscalar(BZ2_p,L3_p)-Iscalar(BZ2_p,L2_p) * Iscalar(BZ1_p,L3_p)) * BZ3_p )
         !(-L1g3 L3g2 m + L1g2 L3g3 m ) g1 + (L1g3 L3g1 m  - 
         ! L1g1 L3g3 m ) g2 + (-L1g2 L3g1 m  + 
         ! L1g1 L3g2 m ) g3
         b2_p = X*( (Iscalar(BZ1_p,L1_p) * Iscalar(BZ3_p,L3_p)-Iscalar(BZ3_p,L1_p) * Iscalar(BZ1_p,L3_p)) * BZ2_p &
               & +  (Iscalar(BZ2_p,L1_p) * Iscalar(BZ1_p,L3_p)-Iscalar(BZ1_p,L1_p) * Iscalar(BZ2_p,L3_p)) * BZ3_p &
               & +  (Iscalar(BZ3_p,L1_p) * Iscalar(BZ2_p,L3_p)-Iscalar(BZ2_p,L1_p) * Iscalar(BZ3_p,L3_p)) * BZ1_p )
         ! (l L1g3 L2g2 - l L1g2 L2g3) g1 + (-l L1g3 L2g1 + 
         !l L1g1 L2g3) g2 + (l L1g2 L2g1 - l L1g1 L2g2 ) g3
         b3_p = X*( (Iscalar(BZ1_p,L1_p) * Iscalar(BZ2_p,L2_p)-Iscalar(BZ2_p,L1_p) * Iscalar(BZ1_p,L2_p)) * BZ3_p &
               & +  (Iscalar(BZ2_p,L1_p) * Iscalar(BZ3_p,L2_p)-Iscalar(BZ3_p,L1_p) * Iscalar(BZ2_p,L2_p)) * BZ1_p &
               & +  (Iscalar(BZ3_p,L1_p) * Iscalar(BZ1_p,L2_p)-Iscalar(BZ1_p,L1_p) * Iscalar(BZ3_p,L2_p)) * BZ2_p ) 
         Latt%b1_p  = b1_p
         Latt%b2_p  = b2_p
         Latt%b3_p  = b3_p


         ! Setup the 3X3 matrix to determine  b1_perp_p, b2_perp_p, b3_perp_p
         Mat(1,1) = dble(b1_p(1))
         Mat(1,2) = dble(b1_p(2))
         Mat(1,3) = dble(b1_p(3))
         Mat(2,1) = dble(b2_p(1))
         Mat(2,2) = dble(b2_p(2))
         Mat(2,3) = dble(b2_p(3))
         Mat(3,1) = dble(b3_p(1))
         Mat(3,2) = dble(b3_p(2))
         Mat(3,3) = dble(b3_p(3))
         X = Mat(1,1)*(Mat(2,2)*Mat(3,3) - Mat(2,3)*Mat(3,2)) -&
               &(Mat(1,2)*(Mat(2,1)*Mat(3,3) - Mat(2,3)*Mat(3,1)))+&
               &(Mat(1,3)*(Mat(2,1)*Mat(3,2) - Mat(2,2)*Mat(3,1)))
         Mat_inv(1,1) =  (Mat(2,2)*Mat(3,3) - Mat(2,3)*Mat(3,2))/X
         Mat_inv(2,2) =  (Mat(1,1)*Mat(3,3) - Mat(1,3)*Mat(3,1))/X
         Mat_inv(3,3) =  (Mat(1,1)*Mat(2,2) - Mat(1,2)*Mat(2,1))/X
         Mat_inv(1,2) = -(Mat(1,2)*Mat(3,3) - Mat(1,3)*Mat(3,2))/X
         Mat_inv(2,1) = -(Mat(2,1)*Mat(3,3) - Mat(2,3)*Mat(3,1))/X
         Mat_inv(2,3) = -(Mat(1,1)*Mat(2,3) - Mat(1,3)*Mat(2,1))/X
         Mat_inv(3,2) = -(Mat(1,1)*Mat(3,2) - Mat(1,2)*Mat(3,1))/X
         Mat_inv(3,1) =  (Mat(2,1)*Mat(3,2) - Mat(2,2)*Mat(3,1))/X
         Mat_inv(1,3) =  (Mat(1,2)*Mat(2,3) - Mat(1,3)*Mat(2,2))/X
         Latt%b1_perp_p(1)      = Mat_inv(1,1)
         Latt%b1_perp_p(2)      = Mat_inv(2,1)
         Latt%b1_perp_p(3)      = Mat_inv(3,1)
         Latt%b2_perp_p(1)      = Mat_inv(1,2)
         Latt%b2_perp_p(2)      = Mat_inv(2,2)
         Latt%b2_perp_p(3)      = Mat_inv(3,2)
         Latt%b3_perp_p(1)      = Mat_inv(1,3)
         Latt%b3_perp_p(2)      = Mat_inv(2,3)
         Latt%b3_perp_p(3)      = Mat_inv(3,3)

         Deallocate ( Mat,  Mat_inv )



         ! Count the number of lattice points
         L      =   abs(nint ( Iscalar(Latt%BZ1_p,L1_p) / (2.d0*pi) ))
         L1     =   abs(nint ( Iscalar(Latt%BZ2_p,L1_p) / (2.d0*pi) ))
         if (L1 .gt. L) L = L1
         L1     =   abs(nint ( Iscalar(Latt%BZ3_p,L1_p) / (2.d0*pi) ))
         if (L1 .gt. L) L = L1
         L1     =   abs(nint ( Iscalar(Latt%BZ1_p,L2_p) / (2.d0*pi) ))
         if (L1 .gt. L) L = L1
         L1     =   abs(nint ( Iscalar(Latt%BZ2_p,L2_p) / (2.d0*pi) ))
         if (L1 .gt. L) L = L1
         L1     =   abs(nint ( Iscalar(Latt%BZ3_p,L2_p) / (2.d0*pi) ))
         if (L1 .gt. L) L = L1
         L1     =   abs(nint ( Iscalar(Latt%BZ1_p,L3_p) / (2.d0*pi) ))
         if (L1 .gt. L) L = L1
         L1     =   abs(nint ( Iscalar(Latt%BZ2_p,L3_p) / (2.d0*pi) ))
         if (L1 .gt. L) L = L1
         L1     =   abs(nint ( Iscalar(Latt%BZ3_p,L3_p) / (2.d0*pi) ))
         if (L1 .gt. L) L = L1
         nc = 0
         do i1 = -L,L
            do i2 = -L,L
               do i3 = -L,L
                  x_p  = dble(i1)*a1_p + dble(i2)*a2_p + dble(i3)*a3_p
                  L_f = 1
                  ! CHECK: not sure how many compinations necessary
                  do i = 1,13
                     if (i.eq.1) a_p =  L3_p
                     if (i.eq.2) a_p =  L2_p
                     if (i.eq.3) a_p =  L1_p
                     if (i.eq.4) a_p =  L2_p - L1_p
                     if (i.eq.5) a_p =  L2_p + L1_p
                     if (i.eq.6) a_p =  L3_p - L1_p
                     if (i.eq.7) a_p =  L3_p + L1_p
                     if (i.eq.8) a_p =  L3_p - L2_p
                     if (i.eq.9) a_p =  L3_p + L2_p
                     if (i.eq.10) a_p =  L3_p + L2_p + L1_p
                     if (i.eq.11) a_p =  L3_p + L2_p - L1_p
                     if (i.eq.12) a_p =  L3_p - L2_p + L1_p
                     if (i.eq.13) a_p =  -L3_p + L2_p + L1_p
                     if  (  Iscalar(x_p, a_p)  .le.  xnorm(a_p)**2/2.d0 + Zero   .and.   &
                          & Iscalar(x_p, a_p)  .ge. -xnorm(a_p)**2/2.d0 + Zero    ) then
                        L_f = L_f * 1
                     else
                        L_f = 0
                     endif
                  enddo
                  if (L_f .eq. 1) then
                     nc = nc + 1
                  endif
               enddo
            enddo
         enddo
         LQ = nc
         Latt%Ns = LQ
         Latt%N  = LQ
         !Write(6,*) L, LQ


         Allocate ( Latt%List(LQ,ndim), Latt%Invlist(-L:L, -L:L, -L:L ) )
         Latt%List = 0
         Latt%Invlist = 0
         !Setting up real space lattice
         nc = 0
         do i1 = -L,L
            do i2 = -L,L
               do i3 = -L,L
                  x_p  = dble(i1)*a1_p + dble(i2)*a2_p + dble(i3)*a3_p
                  L_f = 1
                  do i = 1,4
                     if (i.eq.1) a_p =  L3_p
                     if (i.eq.2) a_p =  L2_p
                     if (i.eq.3) a_p =  L1_p
                     if (i.eq.4) a_p =  L2_p - L1_p
                     if (i.eq.5) a_p =  L2_p + L1_p
                     if (i.eq.6) a_p =  L3_p - L1_p
                     if (i.eq.7) a_p =  L3_p + L1_p
                     if (i.eq.8) a_p =  L3_p - L2_p
                     if (i.eq.9) a_p =  L3_p + L2_p
                     if (i.eq.10) a_p =  L3_p + L2_p + L1_p
                     if (i.eq.11) a_p =  L3_p + L2_p - L1_p
                     if (i.eq.12) a_p =  L3_p - L2_p + L1_p
                     if (i.eq.13) a_p =  -L3_p + L2_p + L1_p
                     if  (  Iscalar( x_p, a_p )    .le.  xnorm(a_p)**2/2.d0 + Zero   .and.   &
                          & Iscalar( x_p, a_p )    .ge. -xnorm(a_p)**2/2.d0 + Zero    ) then
                        L_f = L_f * 1
                     else
                        L_f = 0
                     endif
                  enddo
                  if (L_f .eq. 1) then
                     nc = nc + 1
                     Latt%list(nc,1) = i1
                     Latt%list(nc,2) = i2
                     Latt%list(nc,3) = i3
                     Latt%invlist(i1, i2, i3 ) = nc
                  endif
               enddo
            enddo
         enddo


         Allocate ( Latt%Listk(LQ,ndim), Latt%Invlistk(-L:L, -L:L, -L:L) )
         Latt%Listk = 0
         Latt%Invlistk = 0
         nc = 0
         do m = -L,L
            do n = -L,L
               do o = -L,L
                  xk_p = dble(m) * b1_p + dble(n) *  b2_p + dble(o) * b3_p
               L_f = 1
               do i = 1,13
                  if (i.eq.1) b_p = BZ3_p
                  if (i.eq.2) b_p = BZ2_p
                  if (i.eq.3) b_p = BZ1_p
                  if (i.eq.4) b_p = BZ2_p - BZ1_p
                  if (i.eq.5) b_p = BZ2_p + BZ1_p
                  if (i.eq.6) b_p = BZ2_p - BZ3_p
                  if (i.eq.7) b_p = BZ2_p + BZ3_p
                  if (i.eq.8) b_p = BZ3_p - BZ1_p
                  if (i.eq.9) b_p = BZ3_p + BZ1_p
                  if (i.eq.10) b_p = BZ1_p + BZ2_p + BZ3_p
                  if (i.eq.11) b_p = BZ1_p + BZ2_p - BZ3_p
                  if (i.eq.12) b_p = -BZ1_p + BZ2_p + BZ3_p
                  if (i.eq.13) b_p = BZ1_p - BZ2_p + BZ3_p
                  if  (  Iscalar( xk_p, b_p )    .le.  xnorm(b_p)**2/2.d0 + Zero   .and.   &
                       & Iscalar( xk_p, b_p )    .ge. -xnorm(b_p)**2/2.d0 + Zero    ) then
                     L_f = L_f * 1
                  else
                     L_f = 0
                  endif
               enddo
               if (L_f .eq. 1) then
                  !write(11,"(F14.7,2x,F14.7)")  xk_p(1), xk_p(2)
                  nc = nc + 1
                  Latt%listk(nc,1) = m
                  Latt%listk(nc,2) = n
                  Latt%listk(nc,3) = o
                  Latt%invlistk(m,n,o) = nc
               endif
               enddo
            enddo
         enddo
         If (nc.ne.Latt%N) Then
            write(error_unit,*) 'Make_lattice: Error ', nc, Latt%N
            Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
         endif

         !Setup nnlist
         Allocate ( Latt%nnlist(LQ,-1:1,-1:1,-1:1) )

         do nr = 1, Latt%N
            do nd1 = -1,1
               do nd2 = -1,1
                  do nd3 = -1,1
                     d_p = dble(nd1)*a1_p + dble(nd2)*a2_p + dble(nd3)*a3_p
                     x_p  = dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p &
                     &       + dble(Latt%list(nr,3))*Latt%a3_p  + d_p
                     ! Why is this repeated four times?
                     call npbc(x1_p, x_p , Latt%L1_p, Latt%L2_p, Latt%L3_p)
                     call npbc(x_p , x1_p, Latt%L1_p, Latt%L2_p, Latt%L3_p)
                     call npbc(x1_p, x_p , Latt%L1_p, Latt%L2_p, Latt%L3_p)
                     call npbc(x_p , x1_p, Latt%L1_p, Latt%L2_p, Latt%L3_p)
                     nnr1 =  nint ( Iscalar(Latt%BZ1_p,x_p) / (2.d0*pi) )
                     nnr2 =  nint ( Iscalar(Latt%BZ2_p,x_p) / (2.d0*pi) )
                     nnr3 =  nint ( Iscalar(Latt%BZ3_p,x_p) / (2.d0*pi) )
                     nnr  = Latt%invlist(nnr1,nnr2,nnr3)
                     Latt%nnlist(nr,nd1,nd2,nd3) = nnr
                     if ( nnr < 1  .or.  nnr > Latt%N ) then
                         write(6,*) "Error in nnlist ", nnr
                        x1_p =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p +&
                                    & dble(Latt%list(nr,3))*Latt%a3_p
                        !Write(91,"(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)") x1_p(1), x1_p(2), d_p(1), d_p(2)
                        Write(91,"(F14.7,2x,F14.7)") x1_p(1) , x1_p(2)
                        Write(91,*)
                     endif
                  enddo
               enddo
            enddo
         enddo

         !Setup nnlistk
         Allocate ( Latt%nnlistk(LQ,-1:1,-1:1,-1:1) )
         do nr = 1, Latt%N
            do nd1 = -1,1
               do nd2 = -1,1
                  do nd3 = -1,1
                     d_p = dble(nd1)*b1_p + dble(nd2)*b2_p + dble(nd3)*n3_p
                     x_p  = dble(Latt%listk(nr,1))*Latt%b1_p + dble(Latt%listk(nr,2))*Latt%b2_p &
                           &+ dble(Latt%listk(nr,3))*Latt%b3_p + d_p
                     Latt%nnlistk(nr,nd1,nd2,nd3)  = Inv_K(x_p,Latt)
                  enddo
               enddo
            enddo
         enddo

         !Setup imj
         Allocate ( Latt%imj(LQ,LQ) )
         do nr = 1, Latt%N
            x_p = dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*a2_p + dble(Latt%list(nr,3))*a3_p
            do nr1 = 1,Latt%N
               x1_p = dble(Latt%list(nr1,1))*Latt%a1_p + dble(Latt%list(nr1,2))*a2_p + dble(Latt%list(nr1,3))*a3_p
               d_p = x_p - x1_p
               call npbc(x1_p  , d_p , Latt%L1_p, Latt%L2_p, Latt%L3_p)
               call npbc(d_p , x1_p, Latt%L1_p, Latt%L2_p, Latt%L3_p)
               imj_1 =  nint ( Iscalar(Latt%BZ1_p,d_p) / (2.d0*pi) )
               imj_2 =  nint ( Iscalar(Latt%BZ2_p,d_p) / (2.d0*pi) )
               imj_3 =  nint ( Iscalar(Latt%BZ3_p,d_p) / (2.d0*pi) )
               imj   = Latt%invlist(imj_1,imj_2,imj_3)
               Latt%imj(nr,nr1) = imj
            enddo
         enddo

         deallocate ( b1_p, b2_p, xk_p, b_p )
         deallocate ( BZ1_p, BZ2_p )
         deallocate ( x_p,  x1_p, d_p,  a_p )

       end subroutine MAKE_LATTICE3D
        

         subroutine Make_lattice2D(L1_p, L2_p, a1_p, a2_p, Latt)

           ! This is for a general tilted square lattice defined by the vector a1, a2
           ! L1_p,  L2_p define cluster topology. ( Tilted etc.)
           ! L1_p = n*a1_p + m *a2_p


           Implicit none

           Real (Kind=Kind(0.d0)),  dimension(:) :: L1_p, L2_p, a1_p, a2_p
           Type (Lattice) :: Latt

           Real (Kind=Kind(0.d0)), dimension(:), allocatable :: xk_p, b1_p, b2_p, BZ1_p, BZ2_p, b_p
           Real (Kind=Kind(0.d0)), dimension(:), allocatable :: x_p, x1_p, a_p,d_p
           Real (Kind=Kind(0.d0)), allocatable :: Mat(:,:), Mat_inv(:,:)

           Integer :: ndim, L, L1, nc, i, i1,i2, L_f, LQ, n,m, nd1,nd2,nr, nnr1, nnr2, nnr, nr1, imj_1, imj_2
           Integer :: imj
           Real    (Kind=Kind(0.d0)) :: Zero,pi, X

           ndim = size(L1_p)
           allocate (Latt%L2_p(ndim), Latt%L1_p(ndim), Latt%a1_p(ndim) , Latt%a2_p(ndim), &
                &    Latt%b1_p(ndim), Latt%b2_p(ndim), Latt%BZ1_p(ndim), Latt%BZ2_p(ndim) )
           allocate (Latt%b1_perp_p(ndim), Latt%b2_perp_p(ndim) )
           Zero = 1.D-5
           Latt%L1_p = L1_p
           Latt%L2_p = L2_p
           Latt%a1_p = a1_p
           Latt%a2_p = a2_p


           !Compute the Reciprocal lattice vectors.
           Allocate ( b1_p(ndim), b2_p(ndim), xk_p(ndim), b_p(ndim) )
           Allocate ( BZ1_p(ndim), BZ2_p(ndim) )
           Allocate ( x_p(ndim),  x1_p(ndim), d_p(ndim),  a_p(ndim) )


           pi   = acos(-1.d0)

           ! Setup the 2X2 matrix to determine  BZ1_p, BZ2_p
           Allocate ( Mat(2 , 2), Mat_inv( 2 , 2 ) )
           Mat(1,1) = dble(a1_p(1))
           Mat(1,2) = dble(a1_p(2))
           Mat(2,1) = dble(a2_p(1))
           Mat(2,2) = dble(a2_p(2))
           X = Mat(1,1)*Mat(2,2) - Mat(2,1)*Mat(1,2)
           Mat_inv(1,1) =  Mat(2,2)/X
           Mat_inv(2,2) =  Mat(1,1)/X
           Mat_inv(1,2) = -Mat(1,2)/X
           Mat_inv(2,1) = -Mat(2,1)/X
           BZ1_p(1)      = 2.d0*pi*Mat_inv(1,1)
           BZ1_p(2)      = 2.d0*pi*Mat_inv(2,1)
           BZ2_p(1)      = 2.d0*pi*Mat_inv(1,2)
           BZ2_p(2)      = 2.d0*pi*Mat_inv(2,2)
           Latt%BZ1_p = BZ1_p
           Latt%BZ2_p = BZ2_p




           ! K-space Quantization  from periodicity in L1_p and L2_p
           X =  2.d0*pi / ( Iscalar(BZ1_p,L1_p) * Iscalar(BZ2_p,L2_p) -   &
                &           Iscalar(BZ2_p,L1_p) * Iscalar(BZ1_p,L2_p)   )
           X = abs(X)
           b1_p = X*( Iscalar(BZ2_p,L2_p) * BZ1_p - Iscalar(BZ1_p,L2_p) * BZ2_p )
           b2_p = X*( Iscalar(BZ1_p,L1_p) * BZ2_p - Iscalar(BZ2_p,L1_p) * BZ1_p )
           Latt%b1_p  = b1_p
           Latt%b2_p  = b2_p


           ! Setup the 2X2 matrix to determine  b1_perp_p, b2_perp_p
           Mat(1,1) = dble(b1_p(1))
           Mat(1,2) = dble(b1_p(2))
           Mat(2,1) = dble(b2_p(1))
           Mat(2,2) = dble(b2_p(2))
           X = Mat(1,1)*Mat(2,2) - Mat(2,1)*Mat(1,2)
           Mat_inv(1,1) =  Mat(2,2)/X
           Mat_inv(2,2) =  Mat(1,1)/X
           Mat_inv(1,2) = -Mat(1,2)/X
           Mat_inv(2,1) = -Mat(2,1)/X
           Latt%b1_perp_p(1)      = Mat_inv(1,1)
           Latt%b1_perp_p(2)      = Mat_inv(2,1)
           Latt%b2_perp_p(1)      = Mat_inv(1,2)
           Latt%b2_perp_p(2)      = Mat_inv(2,2)

           Deallocate ( Mat,  Mat_inv )



           ! Count the number of lattice points.
           L      =   abs(nint ( Iscalar(Latt%BZ1_p,L1_p) / (2.d0*pi) ))
           L1     =   abs(nint ( Iscalar(Latt%BZ2_p,L1_p) / (2.d0*pi) ))
           if (L1 .gt. L) L = L1
           L1     =   abs(nint ( Iscalar(Latt%BZ1_p,L2_p) / (2.d0*pi) ))
           if (L1 .gt. L) L = L1
           L1     =   abs(nint ( Iscalar(Latt%BZ2_p,L2_p) / (2.d0*pi) ))
           if (L1 .gt. L) L = L1
           nc = 0
           do i1 = -L,L
              do i2 = -L,L
                 x_p  = dble(i1)*a1_p + dble(i2)*a2_p
                 L_f = 1
                 do i = 1,4
                    if (i.eq.1) a_p =  L2_p
                    if (i.eq.2) a_p =  L1_p
                    if (i.eq.3) a_p =  L2_p - L1_p
                    if (i.eq.4) a_p =  L2_p + L1_p
                    if  (  Iscalar(x_p, a_p)  .le.  xnorm(a_p)**2/2.d0 + Zero   .and.   &
                         & Iscalar(x_p, a_p)  .ge. -xnorm(a_p)**2/2.d0 + Zero    ) then
                       L_f = L_f * 1
                    else
                       L_f = 0
                    endif
                 enddo
                 if (L_f .eq. 1) then
                    nc = nc + 1
                 endif
              enddo
           enddo
           LQ = nc
           Latt%Ns = LQ
           Latt%N  = LQ
           !Write(6,*) L, LQ


           Allocate ( Latt%List(LQ,ndim), Latt%Invlist(-L:L, -L:L ) )
           Latt%List = 0
           Latt%Invlist = 0
           !Setting up real space lattice
           nc = 0
           do i1 = -L,L
              do i2 = -L,L
                 x_p  = dble(i1)*a1_p + dble(i2)*a2_p
                 L_f = 1
                 do i = 1,4
                    if (i.eq.1) a_p =  L2_p
                    if (i.eq.2) a_p =  L1_p
                    if (i.eq.3) a_p =  L2_p - L1_p
                    if (i.eq.4) a_p =  L2_p + L1_p
                    if  (  Iscalar( x_p, a_p )    .le.  xnorm(a_p)**2/2.d0 + Zero   .and.   &
                         & Iscalar( x_p, a_p )    .ge. -xnorm(a_p)**2/2.d0 + Zero    ) then
                       L_f = L_f * 1
                    else
                       L_f = 0
                    endif
                 enddo
                 if (L_f .eq. 1) then
                    nc = nc + 1
                    Latt%list(nc,1) = i1
                    Latt%list(nc,2) = i2
                    Latt%invlist(i1, i2 ) = nc
                 endif
              enddo
           enddo


           Allocate ( Latt%Listk(LQ,ndim), Latt%Invlistk(-L:L, -L:L) )
           Latt%Listk = 0
           Latt%Invlistk = 0
           nc = 0
           do m = -L,L
              do n = -L,L
                 xk_p = dble(m) * b1_p + dble(n) *  b2_p
                 L_f = 1
                 do i = 1,4
                    if (i.eq.1) b_p = BZ2_p
                    if (i.eq.2) b_p = BZ1_p
                    if (i.eq.3) b_p = BZ2_p - BZ1_p
                    if (i.eq.4) b_p = BZ2_p + BZ1_p
                    if  (  Iscalar( xk_p, b_p )    .le.  xnorm(b_p)**2/2.d0 + Zero   .and.   &
                         & Iscalar( xk_p, b_p )    .ge. -xnorm(b_p)**2/2.d0 + Zero    ) then
                       L_f = L_f * 1
                    else
                       L_f = 0
                    endif
                 enddo
                 if (L_f .eq. 1) then
                    !write(11,"(F14.7,2x,F14.7)")  xk_p(1), xk_p(2)
                    nc = nc + 1
                    Latt%listk(nc,1) = m
                    Latt%listk(nc,2) = n
                    Latt%invlistk(m,n) = nc
                 endif
              enddo
           enddo
           If (nc.ne.Latt%N) Then
              write(error_unit,*) 'Make_lattice: Error ', nc, Latt%N
              Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif

           !Setup nnlist
           Allocate ( Latt%nnlist(LQ,-1:1,-1:1) )

           do nr = 1, Latt%N
              do nd1 = -1,1
                 do nd2 = -1,1
                    d_p = dble(nd1)*a1_p + dble(nd2)*a2_p
                    x_p  = dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p  + d_p
                    call npbc(x1_p, x_p , Latt%L1_p, Latt%L2_p)
                    call npbc(x_p , x1_p, Latt%L1_p, Latt%L2_p)
                    call npbc(x1_p, x_p , Latt%L1_p, Latt%L2_p)
                    call npbc(x_p , x1_p, Latt%L1_p, Latt%L2_p)
                    nnr1 =  nint ( Iscalar(Latt%BZ1_p,x_p) / (2.d0*pi) )
                    nnr2 =  nint ( Iscalar(Latt%BZ2_p,x_p) / (2.d0*pi) )
                    nnr  = Latt%invlist(nnr1,nnr2)
                    Latt%nnlist(nr,nd1,nd2) = nnr
                    if ( nnr < 1  .or.  nnr > Latt%N ) then
                        write(6,*) "Error in nnlist ", nnr
                        x1_p =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p
                        !Write(91,"(F14.7,2x,F14.7,2x,F14.7,2x,F14.7)") x1_p(1), x1_p(2), d_p(1), d_p(2)
                        Write(91,"(F14.7,2x,F14.7)") x1_p(1) , x1_p(2)
                        Write(91,*)
                     endif
                 enddo
              enddo
           enddo

           !Setup nnlistk
           Allocate ( Latt%nnlistk(LQ,-1:1,-1:1) )
           do nr = 1, Latt%N
              do nd1 = -1,1
                 do nd2 = -1,1
                    d_p = dble(nd1)*b1_p + dble(nd2)*b2_p
                    x_p  = dble(Latt%listk(nr,1))*Latt%b1_p + dble(Latt%listk(nr,2))*Latt%b2_p  + d_p
                    Latt%nnlistk(nr,nd1,nd2)  = Inv_K(x_p,Latt)
                 enddo
              enddo
           enddo

           !Setup imj
           Allocate ( Latt%imj(LQ,LQ) )
           do nr = 1, Latt%N
              x_p = dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*a2_p
              do nr1 = 1,Latt%N
                 x1_p = dble(Latt%list(nr1,1))*Latt%a1_p + dble(Latt%list(nr1,2))*a2_p
                 d_p = x_p - x1_p
                 call npbc(x1_p  , d_p , Latt%L1_p, Latt%L2_p)
                 call npbc(d_p , x1_p, Latt%L1_p, Latt%L2_p)
                 imj_1 =  nint ( Iscalar(Latt%BZ1_p,d_p) / (2.d0*pi) )
                 imj_2 =  nint ( Iscalar(Latt%BZ2_p,d_p) / (2.d0*pi) )
                 imj   = Latt%invlist(imj_1,imj_2)
                 Latt%imj(nr,nr1) = imj
              enddo
           enddo

           deallocate ( b1_p, b2_p, xk_p, b_p )
           deallocate ( BZ1_p, BZ2_p )
           deallocate ( x_p,  x1_p, d_p,  a_p )

         end subroutine MAKE_LATTICE2D

         subroutine Clear_Lattice(Latt)
           Implicit none
           Type (Lattice) :: Latt

           if(size(Latt%a1_p).eq.3) then
            deallocate (Latt%L3_p, Latt%a3_p, Latt%b3_p, Latt%BZ3_p)
            deallocate (Latt%b3_perp_p)
           endif

           deallocate (Latt%L2_p, Latt%L1_p, Latt%a1_p , Latt%a2_p, &
                &    Latt%b1_p, Latt%b2_p, Latt%BZ1_p, Latt%BZ2_p )
           deallocate (Latt%b1_perp_p, Latt%b2_perp_p )
           deallocate ( Latt%Listk, Latt%Invlistk, Latt%imj, Latt%nnlistk, &
                &    Latt%nnlist,Latt%List, Latt%Invlist)
            

         end subroutine Clear_Lattice


 !********
         subroutine npbc_I2D(nr_p, n_p, L1_p, L2_p)

           Implicit none

           integer, dimension(:) ::  n_p, L1_p, L2_p
           integer, dimension(:), intent(out) :: nr_p

           integer, dimension(:), allocatable :: x_p
           Real (Kind=Kind(0.d0)) :: Zero, X
           Integer :: Ndim, i

           Zero = 1.D-8
           nr_p = n_p
           ndim = size(n_p)

           allocate (x_p(ndim))

           do  i = 1,4
              if (i.eq.1) x_p = L2_p
              if (i.eq.2) x_p = L1_p
              if (i.eq.3) x_p = L2_p - L1_p
              if (i.eq.4) x_p = L2_p + L1_p

              X = real(Iscalar(nr_p,x_p),kind(0.d0))/(Xnorm(x_p)**2)
              if (X .ge.   0.5D0 + Zero  ) nr_p = nr_p - x_p
              if (X .le.  -0.5D0 + Zero  ) nr_p = nr_p + x_p
           enddo

           deallocate(x_p)

         end subroutine npbc_I2D

         subroutine npbc_I3D(nr_p, n_p, L1_p, L2_p, L3_p)

            Implicit none
 
            integer, dimension(:) ::  n_p, L1_p, L2_p, L3_p
            integer, dimension(:), intent(out) :: nr_p
 
            integer, dimension(:), allocatable :: x_p
            Real (Kind=Kind(0.d0)) :: Zero, X
            Integer :: Ndim, i
 
            Zero = 1.D-8
            nr_p = n_p
            ndim = size(n_p)
 
            allocate (x_p(ndim))
 
            do  i = 1,13
               if (i.eq.1) x_p = L3_p
               if (i.eq.2) x_p = L2_p
               if (i.eq.3) x_p = L1_p
               if (i.eq.4) x_p = L2_p - L1_p
               if (i.eq.5) x_p = L2_p + L1_p
               if (i.eq.6) x_p = L2_p - L3_p
               if (i.eq.7) x_p = L2_p + L3_p
               if (i.eq.8) x_p = L3_p - L1_p
               if (i.eq.9) x_p = L3_p + L1_p
               if (i.eq.10) x_p = L3_p + L2_p + L1_p
               if (i.eq.11) x_p = L3_p + L2_p - L1_p
               if (i.eq.12) x_p = L3_p - L2_p + L1_p
               if (i.eq.13) x_p = - L3_p + L2_p + L1_p
 
               X = real(Iscalar(nr_p,x_p),kind(0.d0))/(Xnorm(x_p)**2)
               if (X .ge.   0.5D0 + Zero  ) nr_p = nr_p - x_p
               if (X .le.  -0.5D0 + Zero  ) nr_p = nr_p + x_p
            enddo
 
            deallocate(x_p)
 
         end subroutine npbc_I3D

         subroutine npbc_R2D(nr_p, n_p, L1_p, L2_p)

            Implicit none
            Real (Kind=Kind(0.d0)), dimension(:) :: n_p, L1_p, L2_p
            Real (Kind=Kind(0.d0)), dimension(:), intent(out) :: nr_p
 
            Real (Kind=Kind(0.d0)), dimension(:), allocatable :: x_p
 
            Real (Kind=Kind(0.d0)) :: Zero, X
            Integer :: ndim, i
 
            Zero = 1.D-8
            nr_p = n_p
            ndim = size(n_p)
            allocate(x_p(ndim))
            do i = 1,4
               if (i.eq.1) x_p = L2_p
               if (i.eq.2) x_p = L1_p
               if (i.eq.3) x_p = L2_p - L1_p
               if (i.eq.4) x_p = L2_p + L1_p
               X =  Iscalar(nr_p,x_p)/(Xnorm(x_p)**2)
               if (X .ge.   0.5D0 + Zero  ) nr_p = nr_p - x_p
               if (X .le.  -0.5D0 + Zero  ) nr_p = nr_p + x_p
            enddo
 
            deallocate(x_p)
 
          end subroutine npbc_R2D

         subroutine npbc_R3D(nr_p, n_p, L1_p, L2_p, L3_p)

           Implicit none
           Real (Kind=Kind(0.d0)), dimension(:) :: n_p, L1_p, L2_p, L3_p
           Real (Kind=Kind(0.d0)), dimension(:), intent(out) :: nr_p

           Real (Kind=Kind(0.d0)), dimension(:), allocatable :: x_p

           Real (Kind=Kind(0.d0)) :: Zero, X
           Integer :: ndim, i

           Zero = 1.D-8
           nr_p = n_p
           ndim = size(n_p)
           allocate(x_p(ndim))
           do i = 1,13
              if (i.eq.1) x_p = L3_p
              if (i.eq.2) x_p = L2_p
              if (i.eq.3) x_p = L1_p
              if (i.eq.4) x_p = L2_p - L1_p
              if (i.eq.5) x_p = L2_p + L1_p
              if (i.eq.6) x_p = L2_p - L3_p
              if (i.eq.7) x_p = L2_p + L3_p
              if (i.eq.8) x_p = L3_p - L1_p
              if (i.eq.9) x_p = L3_p + L1_p
              if (i.eq.10) x_p = L3_p + L2_p + L1_p
              if (i.eq.11) x_p = L3_p + L2_p - L1_p
              if (i.eq.12) x_p = L3_p - L2_p + L1_p
              if (i.eq.13) x_p = - L3_p + L2_p + L1_p
              X =  Iscalar(nr_p,x_p)/(Xnorm(x_p)**2)
              if (X .ge.   0.5D0 + Zero  ) nr_p = nr_p - x_p
              if (X .le.  -0.5D0 + Zero  ) nr_p = nr_p + x_p
           enddo

           deallocate(x_p)

         end subroutine npbc_R3D

!********
         subroutine npbc_R_B2D(nr_p, n_p, L1_p, L2_p, N1, N2 )

           !n_p = nr_p + N1* L1_p  + N2 * L2_p

           Implicit none
           Real (Kind=Kind(0.d0)), dimension(:), intent(in) :: n_p, L1_p, L2_p
           Real (Kind=Kind(0.d0)), dimension(:), intent(out) :: nr_p
           Integer, intent(inout) :: N1, N2

           Real (Kind=Kind(0.d0)), dimension(:), allocatable :: x_p

           Real (Kind=Kind(0.d0)) :: Zero, X
           Integer :: ndim, i,  Del_N1, Del_N2

           Zero = 1.D-8
           nr_p = n_p
           ndim = size(n_p)
           allocate(x_p(ndim))
           do i = 1,4
              Del_N2 = 0;  Del_N1 = 0
              if (i.eq.1) Del_N2 =1
              if (i.eq.2) Del_N1 =1
              if (i.eq.3) then
                 Del_N1 = -1
                 Del_N2 =  1
              endif
              if (i.eq.4) then
                 Del_N1 =  1
                 Del_N2 =  1
              endif
              x_p = real(Del_N2,kind(0.d0)) * L2_p   + real(Del_N1,kind(0.d0))*L1_p
              X =  Iscalar(nr_p,x_p)/(Xnorm(x_p)**2)
              if (X .ge.   0.5D0 + Zero  ) then
                 nr_p = nr_p - x_p
                 N2 = N2 + Del_N2
                 N1 = N1 + Del_N1
              endif
              if (X .le.  -0.5D0 + Zero  ) then
                 nr_p = nr_p + x_p
                 N2 = N2 - Del_N2
                 N1 = N1 - Del_N1
              endif
           enddo

           deallocate(x_p)

         end subroutine npbc_R_B2D

         subroutine npbc_R_B3D(nr_p, n_p, L1_p, L2_p, L3_p, N1, N2, N3 )

            !n_p = nr_p + N1* L1_p  + N2 * L2_p
 
            Implicit none
            Real (Kind=Kind(0.d0)), dimension(:), intent(in) :: n_p, L1_p, L2_p, L3_p
            Real (Kind=Kind(0.d0)), dimension(:), intent(out) :: nr_p
            Integer, intent(inout) :: N1, N2, N3
 
            Real (Kind=Kind(0.d0)), dimension(:), allocatable :: x_p
 
            Real (Kind=Kind(0.d0)) :: Zero, X
            Integer :: ndim, i,  Del_N1, Del_N2, Del_N3
 
            Zero = 1.D-8
            nr_p = n_p
            ndim = size(n_p)
            allocate(x_p(ndim))
            do i = 1,4
               Del_N2 = 0;  Del_N1 = 0; Del_N3 = 0
               if (i.eq.1) Del_N3 =1
               if (i.eq.2) Del_N2 =1
               if (i.eq.3) Del_N1 =1
               if (i.eq.4) then
                  Del_N1 = -1
                  Del_N2 =  1
               endif
               if (i.eq.5) then
                  Del_N1 =  1
                  Del_N2 =  1
               endif
               if (i.eq.6) then
                  Del_N3 = -1
                  Del_N2 =  1
               endif
               if (i.eq.7) then
                  Del_N3 =  1
                  Del_N2 =  1
               endif
               if (i.eq.8) then
                  Del_N1 = -1
                  Del_N3 =  1
               endif
               if (i.eq.9) then
                  Del_N1 =  1
                  Del_N3 =  1
               endif
               if (i.eq.10) then
                  Del_N1 =  1
                  Del_N2 =  1
                  Del_N3 =  1
               endif
               if (i.eq.11) then
                  Del_N1 =  1
                  Del_N2 =  1
                  Del_N3 =  -1
               endif
               if (i.eq.12) then
                  Del_N1 =  1
                  Del_N2 =  -1
                  Del_N3 =  1
               endif
               if (i.eq.13) then
                  Del_N1 =  -1
                  Del_N2 =  1
                  Del_N3 =  1
               endif
               x_p = real(Del_N2,kind(0.d0)) * L2_p   + real(Del_N1,kind(0.d0))*L1_p + real(Del_N3,kind(0.d0))*L3_p
               X =  Iscalar(nr_p,x_p)/(Xnorm(x_p)**2)
               if (X .ge.   0.5D0 + Zero  ) then
                  nr_p = nr_p - x_p
                  N3 = N3 + Del_N3
                  N2 = N2 + Del_N2
                  N1 = N1 + Del_N1
               endif
               if (X .le.  -0.5D0 + Zero  ) then
                  nr_p = nr_p + x_p
                  N3 = N3 - Del_N3
                  N2 = N2 - Del_N2
                  N1 = N1 - Del_N1
               endif
            enddo
 
            deallocate(x_p)
 
          end subroutine npbc_R_B3D

 !********
         integer Function Inv_K2D(XK_P,Latt)

           Implicit None
           Real (Kind=Kind(0.d0))  :: XK_P(2)
           Type (Lattice) :: Latt

           Integer :: nkx, nky, nk
           Real (Kind=Kind(0.d0)) :: XK1_P(2), XK2_P(2), Zero

           call npbc(xk1_p, xk_p , Latt%BZ1_p, Latt%BZ2_p)
           call npbc(xk2_p, xk1_p, Latt%BZ1_p, Latt%BZ2_p)

           nkx = nint (Iscalar(XK2_P,Latt%b1_perp_p) )
           nky = nint (Iscalar(XK2_P,Latt%b2_perp_p) )
           nk = Latt%Invlistk(nkx,nky)

           !Test
           Zero  = 1.D-10
           XK1_P = Latt%listk(nk,1)*latt%b1_p + Latt%listk(nk,2)*latt%b2_p
           XK1_P = XK1_P - XK2_P
           if (Xnorm(XK1_P)  < Zero ) then
              Inv_K = nk
           else
              write(error_unit,*) 'Lattice: Error in Inv_K'
              Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif

 !!$          nk = 1
 !!$          do
 !!$             XK1_P = Latt%listk(nk,1)*latt%b1_p + Latt%listk(nk,2)*latt%b2_p
 !!$             if (Xnorm(XK1_P - XK_P)  < Zero ) then
 !!$                Inv_K = nk
 !!$                exit
 !!$             elseif (nk < Latt%N) then
 !!$                nk = nk + 1
 !!$             else
 !!$                write(6,*) 'Error in Inv_K Lattice_new'
 !!$                Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
 !!$             endif
 !!$          enddo

         end Function Inv_K2D


         integer Function Inv_K3D(XK_P,Latt)

           Implicit None
           Real (Kind=Kind(0.d0))  :: XK_P(3)
           Type (Lattice) :: Latt

           Integer :: nkx, nky, nkz, nk
           Real (Kind=Kind(0.d0)) :: XK1_P(3), XK2_P(3), Zero

           call npbc(xk1_p, xk_p , Latt%BZ1_p, Latt%BZ2_p, Latt%BZ3_p)
           call npbc(xk2_p, xk1_p, Latt%BZ1_p, Latt%BZ2_p, Latt%BZ3_p)

           nkx = nint (Iscalar(XK2_P,Latt%b1_perp_p) )
           nky = nint (Iscalar(XK2_P,Latt%b2_perp_p) )
           nkz = nint (Iscalar(XK2_P,Latt%b3_perp_p) )
           nk = Latt%Invlistk(nkx,nky,nkz)

           !Test
           Zero  = 1.D-10
           XK1_P = Latt%listk(nk,1)*latt%b1_p + Latt%listk(nk,2)*latt%b2_p + Latt%listk(nk,3)*latt%b3_p
           XK1_P = XK1_P - XK2_P
           if (Xnorm(XK1_P)  < Zero ) then
              Inv_K = nk
           else
              write(error_unit,*) 'Lattice: Error in Inv_K'
              Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
           endif

 !!$          nk = 1
 !!$          do
 !!$             XK1_P = Latt%listk(nk,1)*latt%b1_p + Latt%listk(nk,2)*latt%b2_p
 !!$             if (Xnorm(XK1_P - XK_P)  < Zero ) then
 !!$                Inv_K = nk
 !!$                exit
 !!$             elseif (nk < Latt%N) then
 !!$                nk = nk + 1
 !!$             else
 !!$                write(6,*) 'Error in Inv_K Lattice_new'
 !!$                Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
 !!$             endif
 !!$          enddo

         end Function Inv_K3D



 !********
         integer Function Inv_R2D(XR_P,Latt)

           Implicit None
           Real (Kind=Kind(0.d0))  :: XR_P(2)
           Type (Lattice) :: Latt

           Real (Kind=Kind(0.d0)) :: XR1_P(2), XR2_P(2)

           Integer :: n_1, n_2
           Real (Kind=Kind(0.d0)) :: pi

           pi = acos(-1.d0)
           call npbc(xr1_p, xr_p , Latt%L1_p, Latt%L2_p)
           call npbc(xr2_p, xr1_p, Latt%L1_p, Latt%L2_p)

           n_1 =  nint ( Iscalar(Latt%BZ1_p,XR2_p) / (2.d0*pi) )
           n_2 =  nint ( Iscalar(Latt%BZ2_p,XR2_p) / (2.d0*pi) )
           Inv_R  = Latt%invlist(n_1,n_2)

         end Function Inv_R2D

         integer Function Inv_R3D(XR_P,Latt)

           Implicit None
           Real (Kind=Kind(0.d0))  :: XR_P(3)
           Type (Lattice) :: Latt

           Real (Kind=Kind(0.d0)) :: XR1_P(3), XR2_P(3), XR3_P(3)

           Integer :: n_1, n_2, n_3
           Real (Kind=Kind(0.d0)) :: pi

           pi = acos(-1.d0)
           call npbc(xr1_p, xr_p , Latt%L1_p, Latt%L2_p, Latt%L3_p)
           call npbc(xr2_p, xr1_p, Latt%L1_p, Latt%L2_p, Latt%L3_p)

           n_1 =  nint ( Iscalar(Latt%BZ1_p,XR2_p) / (2.d0*pi) )
           n_2 =  nint ( Iscalar(Latt%BZ2_p,XR2_p) / (2.d0*pi) )
           n_3 =  nint ( Iscalar(Latt%BZ3_p,XR3_p) / (2.d0*pi) )
           Inv_R  = Latt%invlist(n_1,n_2,n_3)

         end Function Inv_R3D
 !********

         integer function Iscalar_II(i_p, j_p)
           Implicit none
           integer, dimension(:), intent(in) :: i_p, j_p
           integer i

           Iscalar_II = 0
           !write(6,*) size(i_p)
           do i = 1,  size(i_p)
             ! write(6,*) i
              Iscalar_II = Iscalar_II + i_p(i)*j_p(i)
           enddo
         end function Iscalar_II

 !********
         Real (Kind=Kind(0.d0))  function Iscalar_IR(x_p, j_p)
           Implicit none
           Real (Kind=Kind(0.d0)), dimension(:), intent(in) ::  x_p
           integer, dimension(:), intent(in) ::  j_p
           integer i

           Iscalar_IR = 0.d0
           !write(6,*) size(i_p)
           do i = 1,  size(x_p)
             ! write(6,*) i
              Iscalar_IR = Iscalar_IR + x_p(i)*dble(j_p(i))
           enddo
         end function Iscalar_IR
 !********

         pure Real (Kind=Kind(0.d0))  function Iscalar_RR(x_p, y_p)
           Implicit none
           Real (Kind = Kind(0.D0)), dimension(:), intent(in) ::  x_p, y_p
           Iscalar_RR = dot_product(x_p, y_p)
         end function Iscalar_RR

 !********
         Real (Kind=Kind(0.d0)) function Xnorm_I(i_p)
           Implicit none
           integer, dimension(:), intent(in) :: i_p
           integer :: i

           Xnorm_I = 0.d0
           do i = 1,  size(i_p)
              Xnorm_I = Xnorm_I + dble(i_p(i)*i_p(i))
           enddo
           Xnorm_I = sqrt(Xnorm_I)
         end function Xnorm_I

 !********
         Real (Kind=Kind(0.d0)) function Xnorm_R(x_p)
           Implicit none
           Real (Kind=Kind(0.d0)), dimension(:), intent(in) :: x_p
           integer :: i

           Xnorm_R = 0.d0
           do i = 1,  size(x_p)
              Xnorm_R = Xnorm_R + x_p(i)*x_p(i)
           enddo
           Xnorm_R = sqrt(Xnorm_R)
         end function Xnorm_R

 !********
         subroutine Print_latt(Latt)

           Implicit Real (Kind=Kind(0.d0)) (A-G,O-Z)
           Implicit Integer (H-N)

           Type (Lattice) :: Latt
           Real (Kind=Kind(0.d0))  :: i_p(2),nd_p(2), i_p3(3), nd_p3(3)
           Real    (Kind=Kind(0.d0)) :: x_p(2), x_p3(3)

           Open (Unit=55,file="Latt_info", status = "unknown")
           if (size(Latt%BZ1_p).eq.2) then
               write(55,*) ' Reciprocal vector 1: ', Latt%BZ1_p(1), Latt%BZ1_p(2)
               write(55,*) ' Reciprocal vector 2: ', Latt%BZ2_p(1), Latt%BZ2_p(2)
           else
               write(55,*) ' Reciprocal vector 1: ', Latt%BZ1_p(1), Latt%BZ1_p(2), Latt%BZ1_p(3)
               write(55,*) ' Reciprocal vector 2: ', Latt%BZ2_p(1), Latt%BZ2_p(2), Latt%BZ2_p(3)
               write(55,*) ' Reciprocal vector 3: ', Latt%BZ3_p(1), Latt%BZ3_p(2), Latt%BZ3_p(3)
           endif
           if (size(Latt%BZ1_p).eq.2) then
               write(55,*) ' Latt       vector 1: ', Latt%a1_p(1), Latt%a1_p(2)
               write(55,*) ' Latt       vector 2: ', Latt%a2_p(1), Latt%a2_p(2)
           else
               write(55,*) ' Latt       vector 1: ', Latt%a1_p(1), Latt%a1_p(2), Latt%a1_p(3)
               write(55,*) ' Latt       vector 2: ', Latt%a2_p(1), Latt%a2_p(2), Latt%a2_p(3)
               write(55,*) ' Latt       vector 3: ', Latt%a3_p(1), Latt%a3_p(2), Latt%a3_p(3)
           endif
           close(55)
           Open (Unit=56,file="Real_space_latt", status = "unknown")
           Open (Unit=57,file="K_space_latt", status = "unknown")
           Open (Unit=58,file="nn_latt", status = "unknown")
           if (size(Latt%BZ1+p).eq.2) then
               do n = 1, Latt%n
                  i_p = dble(Latt%list(n,1))*Latt%a1_p + dble(Latt%list(n,2))*Latt%a2_p
                  write(56,"(F14.7,2x,F14.7)") i_p(1), i_p(2)
                  x_p = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p
                  write(57,"(F14.7,2x,F14.7)") x_p(1), x_p(2)
                  write(58,*)
                  write(58,"('I :',F14.7,2x,F14.7)") i_p(1), i_p(2)
                  do nd1 = -1,1
                     do nd2 = -1,1
                        nd_p =   dble(nd1)*Latt%a1_p + dble(nd2)*Latt%a2_p
                        nnr = Latt%nnlist(n,nd1,nd2)
                        !Write(6,*) 'nnr : ', nnr
                        i_p = dble(Latt%list(nnr,1))*Latt%a1_p + dble(Latt%list(nnr,2))*Latt%a2_p
                        write(58,"('I+(',F12.6,',',F12.6,')=',2x,F14.7,2x,F14.7)") nd_p(1),nd_p(2),i_p(1), i_p(2)
                     enddo
                  enddo
               enddo
            else
               do n = 1, Latt%n
                  i_p3 = dble(Latt%list(n,1))*Latt%a1_p + dble(Latt%list(n,2))*Latt%a2_p + dble(Latt%list(n,3))*Latt%a3_p
                  write(56,"(F14.7,2x,F14.7)") i_p3(1), i_p3(2), i_p3(3)
                  x_p3 = dble(Latt%listk(n,1))*Latt%b1_p + dble(Latt%listk(n,2))*Latt%b2_p + dble(Latt%listk(n,3))*Latt%b3_p
                  write(57,"(F14.7,2x,F14.7)") x_p3(1), x_p3(2), x_p3(3)
                  write(58,*)
                  write(58,"('I :',F14.7,2x,F14.7)") i_p3(1), i_p3(2), i_p3(3)
                  do nd1 = -1,1
                     do nd2 = -1,1
                        do nd3 = -1,1
                           nd_p3 =   dble(nd1)*Latt%a1_p + dble(nd2)*Latt%a2_p + dble(nd3)*Latt%a3_p
                           nnr = Latt%nnlist(n,nd1,nd2,nd3)
                           !Write(6,*) 'nnr : ', nnr
                           i_p3 = dble(Latt%list(nnr,1))*Latt%a1_p + dble(Latt%list(nnr,2))*Latt%a2_p + dble(Latt%list(nnr,3))*Latt%a3_p
                           write(58,"('I+(',F12.6,',',F12.6,')=',2x,F14.7,2x,F14.7)") nd_p3(1),nd_p3(2),nd_p3(3),i_p3(1), i_p3(2),i_p3(3)
                        enddo
                     enddo
                  enddo
               enddo
            endif
           close(56)
           close(57)
           close(58)
         end subroutine Print_latt

 !*******
         subroutine FT_K_to_R_Mat( Xin_K, Xout_R, Latt)

           Implicit none

           Type (Lattice), intent(in)                 :: Latt
           Type (Mat_R ), Dimension(:,:)              :: Xin_K, Xout_R
           Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: X_MAT
           Real (Kind=Kind(0.d0))                              :: XK_p(2), IR_p(2), Xk_p3(3), IR_p3(3)

           Integer :: nb, norb, LQ, nt, nr, nk
           nb      = size(Xin_K,2  )
           norb   = size(Xin_K(1,1)%el,1)
           LQ     = Latt%N

           !Write(6,*) 'Ltrot, norb ', Ltrot, norb
           !Write(6,*) Xin_K(1,1)%el(1,1)
           !Write(6,*) Xin_K(Latt%N,Ltrot)%el(1,1)

           allocate ( X_MAT(norb,norb) )

           if(size(Latt%BZ1_p).eq.2) then
               do nt = 1,nb
                  do nr = 1,LQ
                     IR_p =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p
                     X_MAT = 0.d0
                     do nk = 1,LQ
                        XK_p =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p
                        X_MAT = X_MAT + cos(Iscalar(XK_p,IR_p))*Xin_K(nk,nt)%el
                     enddo
                     Xout_R(nr,nt)%el = X_MAT/dble(LQ)
                  enddo
               enddo
           else
               do nt = 1,nb
                  do nr = 1,LQ
                     IR_p3 =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p + dble(Latt%list(nr,3))*Latt%a3_p
                     X_MAT = 0.d0
                     do nk = 1,LQ
                        XK_p3 =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p + dble(Latt%listk(nk,3))*Latt%b3_p
                        X_MAT = X_MAT + cos(Iscalar(XK_p,IR_p))*Xin_K(nk,nt)%el
                     enddo
                     Xout_R(nr,nt)%el = X_MAT/dble(LQ)
                  enddo
               enddo
            endif

           deallocate(X_Mat)
         end subroutine FT_K_to_R_Mat

 !********
         subroutine FT_K_to_R_Mat_C( Xin_K, Xout_R, Latt)

           Implicit none

           Type (Lattice), intent(in)                    :: Latt
           Type (Mat_C )   , Dimension(:,:)              :: Xin_K, Xout_R
           Complex (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: X_MAT
           Real    (Kind=Kind(0.d0))                              :: XK_p(2), IR_p(2), ang, XK_p3(3), IR_p3(3)

           Integer :: nb, norb, LQ, nt, nr, nk

           nb     = size(Xin_K,2  )
           norb   = size(Xin_K(1,1)%el,1)
           LQ     = Latt%N

           !Write(6,*) 'Ltrot, norb ', Ltrot, norb
           !Write(6,*) Xin_K(1,1)%el(1,1)
           !Write(6,*) Xin_K(Latt%N,nb)%el(1,1)

           allocate ( X_MAT(norb,norb) )

           if(size(Latt%BZ1_p).eq.2) then
               do nt = 1,nb
                  do nr = 1,LQ
                     IR_p =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p
                     X_MAT = cmplx(0.d0, 0.d0, kind(0.D0))
                     do nk = 1,LQ
                        XK_p =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p
                        ang = Iscalar(XK_p, IR_p)
                        X_MAT = X_MAT + cmplx(cos(ang), sin(ang), kind(0.D0)) * Xin_K(nk,nt)%el
                     enddo
                     Xout_R(nr,nt)%el = X_MAT/dble(LQ)
                  enddo
               enddo
            else
               do nt = 1,nb
                  do nr = 1,LQ
                     IR_p3 =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p + dble(Latt%list(nr,3))*Latt%a3_p
                     X_MAT = cmplx(0.d0, 0.d0, kind(0.D0))
                     do nk = 1,LQ
                        XK_p3 =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p + dble(Latt%listk(nk,3))*Latt%b3_p
                        ang = Iscalar(XK_p, IR_p)
                        X_MAT = X_MAT + cmplx(cos(ang), sin(ang), kind(0.D0)) * Xin_K(nk,nt)%el
                     enddo
                     Xout_R(nr,nt)%el = X_MAT/dble(LQ)
                  enddo
               enddo
            endif

           deallocate(X_Mat)

         end subroutine FT_K_to_R_Mat_C

 !********

         subroutine FT_K_to_R( Xin_K, Xout_R, Latt)

           Implicit none

           Type (Lattice), intent(in)                 :: Latt
           Real (Kind=Kind(0.d0)), Dimension(:,:)              :: Xin_K, Xout_R
           Real (Kind=Kind(0.d0))                              :: XK_p(2), IR_p(2), X_Mat, XK_p3(3), IR_p3(3)
           Integer :: LQ, nb, nt, nr, nk

           nb     = size(Xin_K,2  )
           LQ     = Latt%N

           !Write(6,*) 'Ltrot, norb ', Ltrot, norb
           !Write(6,*) Xin_K(1,1)%el(1,1)
           !Write(6,*) Xin_K(Latt%N,Ltrot)%el(1,1)

           if(size(Latt%BZ1_p).eq.2) then
               do nt = 1,nb
                  do nr = 1,LQ
                     IR_p =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p
                     X_MAT = 0.d0
                     do nk = 1,LQ
                        XK_p =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p
                        X_MAT = X_MAT + cos(Iscalar(XK_p,IR_p))*Xin_K(nk,nt)
                     enddo
                     Xout_R(nr,nt) = X_MAT/dble(LQ)
                  enddo
              enddo
            else
               do nt = 1,nb
                  do nr = 1,LQ
                     IR_p3 =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p + dble(Latt%list(nr,3))*Latt%a3_p
                     X_MAT = 0.d0
                     do nk = 1,LQ
                        XK_p3 =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p + dble(Latt%listk(nk,3))*Latt%b3_p
                        X_MAT = X_MAT + cos(Iscalar(XK_p,IR_p))*Xin_K(nk,nt)
                     enddo
                     Xout_R(nr,nt) = X_MAT/dble(LQ)
                  enddo
              enddo
            endif

         end subroutine FT_K_to_R


         subroutine FT_K_to_R_C( Xin_K, Xout_R, Latt)

           Implicit none

           Type (Lattice), intent(in)                 :: Latt
           Complex (Kind=Kind(0.d0)), Dimension(:,:)           :: Xin_K, Xout_R
           Complex (Kind=Kind(0.d0))                           :: Z
           Real    (Kind=Kind(0.d0))                           :: XK_p(2), IR_p(2), XK_p(3), IR_p(3)

           Integer :: nb, LQ, nt, nr, nk

           nb    = size(Xin_K,2  )
           LQ     = Latt%N

           !Write(6,*) 'Ltrot, norb ', Ltrot, norb
           !Write(6,*) Xin_K(1,1)%el(1,1)
           !Write(6,*) Xin_K(Latt%N,Ltrot)%el(1,1)

           if(size(Latt%BZ1_p).eq.2) then
               do nt = 1,nb
                  do nr = 1,LQ
                     IR_p =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p
                     Z = cmplx(0.d0, 0.d0, kind(0.D0))
                     do nk = 1,LQ
                        XK_p =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p
                        Z = Z + cos(Iscalar(XK_p,IR_p))*Xin_K(nk,nt)
                     enddo
                     Xout_R(nr,nt) = Z/dble(LQ)
                  enddo
               enddo
            else
               do nt = 1,nb
                  do nr = 1,LQ
                     IR_p3 =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p + dble(Latt%list(nr,3))*Latt%a3_p
                     Z = cmplx(0.d0, 0.d0, kind(0.D0))
                     do nk = 1,LQ
                        XK_p3 =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p + dble(Latt%listk(nk,3))*Latt%b3_p
                        Z = Z + cos(Iscalar(XK_p,IR_p))*Xin_K(nk,nt)
                     enddo
                     Xout_R(nr,nt) = Z/dble(LQ)
                  enddo
               enddo
            endif

         end subroutine FT_K_to_R_C


         subroutine FT_R_to_K_mat( Xin_R, Xout_K, Latt)

           Implicit none

           Type (Lattice), intent(in)                 :: Latt
           Type (Mat_R ), Dimension(:,:)              :: Xin_R, Xout_K
           Real (Kind=Kind(0.d0)), Dimension(:,:), allocatable :: X_MAT
           Real (Kind=Kind(0.d0))                              :: XK_p(2), IR_p(2), XK_p(3), IR_p(3)

           Integer :: nb, norb, nk, nt, LQ, nr

           nb     = size(Xin_R,2  )
           norb   = size(Xin_R(1,1)%el,1)
           LQ     = Latt%N

           !Write(6,*) 'Ltrot, norb ', Ltrot, norb
           !Write(6,*) Xin_R(1,1)%el(1,1)
           !Write(6,*) Xin_R(Latt%N,Ltrot)%el(1,1)

           allocate ( X_MAT(norb,norb) )

           if(size(Latt%BZ1_p).eq.2) then
               do nt = 1,nb
                  do nk = 1,LQ
                     XK_p =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p
                     X_MAT = 0.d0
                     do nr = 1,LQ
                        IR_p =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p
                        X_MAT = X_MAT + cos(Iscalar(XK_p,IR_p))*Xin_R(nr,nt)%el
                     enddo
                     Xout_K(nk,nt)%el = X_MAT/dble(LQ)
                  enddo
               enddo
            else
               do nt = 1,nb
                  do nk = 1,LQ
                     XK_p3 =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p + dble(Latt%listk(nk,3))*Latt%b3_p
                     X_MAT = 0.d0
                     do nr = 1,LQ
                        IR_p3 =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p + dble(Latt%list(nr,3))*Latt%a3_p
                        X_MAT = X_MAT + cos(Iscalar(XK_p,IR_p))*Xin_R(nr,nt)%el
                     enddo
                     Xout_K(nk,nt)%el = X_MAT/dble(LQ)
                  enddo
               enddo
            endif

           deallocate(X_Mat)
         end subroutine FT_R_to_K_mat

         subroutine FT_R_to_K( Xin_R, Xout_K, Latt)

           Implicit none

           Type (Lattice), intent(in)                 :: Latt
           Real (Kind=Kind(0.d0)),   Dimension(:)              :: Xin_R, Xout_K

           Real (Kind=Kind(0.d0))                              :: XK_p(2), IR_p(2), X_mat, XK_p3(3), IR_p3(3)

           Integer :: nk, LQ, nr

           LQ     = Latt%N

           !Write(6,*) 'Ltrot, norb ', Ltrot, norb
           !Write(6,*) Xin_R(1,1)%el(1,1)
           !Write(6,*) Xin_R(Latt%N,Ltrot)%el(1,1)
           if(size(Latt%BZ1_p).eq.2) then
               do nk = 1,LQ
                  XK_p =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p
                  X_MAT = 0.d0
                  do nr = 1,LQ
                     IR_p =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p
                     X_MAT = X_MAT + cos(Iscalar(XK_p,IR_p))*Xin_R(nr)
                  enddo
                  Xout_K(nk) = X_MAT/dble(LQ)
               enddo
            else
               do nk = 1,LQ
                  XK_p3 =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p + dble(Latt%listk(nk,3))*Latt%b3_p
                  X_MAT = 0.d0
                  do nr = 1,LQ
                     IR_p3 =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p + dble(Latt%list(nr,3))*Latt%a3_p
                     X_MAT = X_MAT + cos(Iscalar(XK_p,IR_p))*Xin_R(nr)
                  enddo
                  Xout_K(nk) = X_MAT/dble(LQ)
               enddo
            endif

         end subroutine FT_R_to_K

 !********
         subroutine FT_R_to_K_C( Xin_R, Xout_K, Latt)

           Implicit none

           Type (Lattice), intent(in)                  :: Latt
           Complex (Kind=Kind(0.d0)), Dimension(:)              :: Xin_R, Xout_K
           Complex (Kind=Kind(0.d0))                            :: X_MAT
           Real    (Kind=Kind(0.d0))                            :: XK_p(2), IR_p(2), ang, XK_p3(3), IR_p3(3)

           Integer :: LQ, nr, nk

           LQ     = Latt%N

           !Write(6,*) 'Ltrot, norb ', Ltrot, norb
           !Write(6,*) Xin_K(1,1)%el(1,1)
           !Write(6,*) Xin_K(Latt%N,nb)%el(1,1)
           if(size(Latt%BZ1_p).eq.2) then
               do nk = 1,LQ
                  XK_p =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p
                  X_MAT = cmplx(0.d0,0.d0, kind(0.D0))
                  do nr = 1,LQ
                     IR_p =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p
                     ang = -Iscalar(XK_p,IR_p)
 !                   X_MAT = X_MAT + exp( cmplx(0.d0,-(Iscalar(XK_p,IR_p)), kind(0.D0)) ) *Xin_R(nr)
                     X_MAT = X_MAT + cmplx(cos(ang), sin(ang), kind(0.D0)) * Xin_R(nr)
                  enddo
                  Xout_K(nk) = X_MAT/dble(LQ)
               enddo
            else
               do nk = 1,LQ
                  XK_p3 =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p + dble(Latt%listk(nk,3))*Latt%b3_p
                  X_MAT = cmplx(0.d0,0.d0, kind(0.D0))
                  do nr = 1,LQ
                     IR_p3 =  dble(Latt%list(nr,1))*Latt%a1_p + dble(Latt%list(nr,2))*Latt%a2_p + dble(Latt%list(nr,3))*Latt%a3_p
                     ang = -Iscalar(XK_p,IR_p)
 !                   X_MAT = X_MAT + exp( cmplx(0.d0,-(Iscalar(XK_p,IR_p)), kind(0.D0)) ) *Xin_R(nr)
                     X_MAT = X_MAT + cmplx(cos(ang), sin(ang), kind(0.D0)) * Xin_R(nr)
                  enddo
                  Xout_K(nk) = X_MAT/dble(LQ)
               enddo
            endif

        end subroutine FT_R_to_K_C

      end Module Lattices_v3
