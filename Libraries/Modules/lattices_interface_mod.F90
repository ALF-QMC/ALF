!  Copyright (C) 2023 The ALF project
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

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> C Interface for creating lattice. For use in pyALF
!>
!
!--------------------------------------------------------------------
    module lattices_interface
        use iso_c_binding, only: c_double, c_int

        use lattices_v3, only: lattice, make_lattice, clear_lattice
        
        implicit none
        Type (Lattice) :: Latt
        
        contains

        subroutine make_lattice_c(ndim, L1_p, L2_p, a1_p, a2_p) bind(c)
            integer(c_int), intent(in), value :: ndim
            real(c_double), intent(in) :: L1_p(ndim), L2_p(ndim), a1_p(ndim), a2_p(ndim)

            call make_lattice(L1_p, L2_p, a1_p, a2_p, Latt)
        end subroutine make_lattice_c

        subroutine get_l(L) bind(c)
            integer(c_int), intent(out) :: L
            L = (size(Latt%Invlist,1) - 1)/2
        end subroutine get_l

        subroutine get_n(N) bind(c)
            integer(c_int), intent(out) :: N
            N = Latt%N
        end subroutine get_n

        subroutine get_arrays( &
                ndim, L, N, &
                BZ1, BZ2, b1, b2, b1_perp, b2_perp, &
                listr, invlistr, nnlistr, &
                listk, invlistk, nnlistk, imj) bind(c)
            integer(c_int), intent(in), value :: ndim, L, N

            real(c_double), intent(out) :: BZ1(ndim), BZ2(ndim), b1(ndim), b2(ndim)
            real(c_double), intent(out) :: b1_perp(ndim), b2_perp(ndim)
            integer(c_int), intent(out) :: listr(N, ndim), invlistr(2*L+1, 2*L+1), nnlistr(N, 3, 3)
            integer(c_int), intent(out) :: listk(N, ndim), invlistk(2*L+1, 2*L+1), nnlistk(N, 3, 3)
            integer(c_int), intent(out) :: imj(N, N)

            BZ1      = Latt%BZ1_p
            BZ2      = Latt%BZ2_p
            b1       = Latt%b1_p
            b2       = Latt%b2_p
            b1_perp  = Latt%b1_perp_p
            b2_perp  = Latt%b2_perp_p
            listr    = Latt%list
            listk    = Latt%listk

            invlistr(1:L+1, 1:L+1) = Latt%invlist(0:L, 0:L) - 1
            invlistr(1:L+1, L+2:2*L+1) = Latt%invlist(0:L, -L:-1) - 1
            invlistr(L+2:2*L+1, 1:L+1) = Latt%invlist(-L:-1, 0:L) - 1
            invlistr(L+2:2*L+1, L+2:2*L+1) = Latt%invlist(-L:-1, -L:-1) - 1

            invlistk(1:L+1, 1:L+1) = Latt%invlistk(0:L, 0:L) - 1
            invlistk(1:L+1, L+2:2*L+1) = Latt%invlistk(0:L, -L:-1) - 1
            invlistk(L+2:2*L+1, 1:L+1) = Latt%invlistk(-L:-1, 0:L) - 1
            invlistk(L+2:2*L+1, L+2:2*L+1) = Latt%invlistk(-L:-1, -L:-1) - 1

            nnlistr(:, 1:2, 1:2)  = Latt%nnlist(:, 0:1, 0:1) - 1
            nnlistr(:, 1:2, 3)  = Latt%nnlist(:, 0:1, -1) - 1
            nnlistr(:, 3, 1:2)  = Latt%nnlist(:, -1, 0:1) - 1
            nnlistr(:, 3, 3)  = Latt%nnlist(:, -1, -1) - 1

            nnlistk(:, 1:2, 1:2)  = Latt%nnlistk(:, 0:1, 0:1) - 1
            nnlistk(:, 1:2, 3)  = Latt%nnlistk(:, 0:1, -1) - 1
            nnlistk(:, 3, 1:2)  = Latt%nnlistk(:, -1, 0:1) - 1
            nnlistk(:, 3, 3)  = Latt%nnlistk(:, -1, -1) - 1

            imj      = Latt%imj - 1
        end subroutine get_arrays
        
        subroutine clear_lattice_c() bind(c)
            Implicit none
            
            call clear_Lattice(Latt)
        end subroutine clear_lattice_c
        
    end module lattices_interface
