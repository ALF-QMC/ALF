!  Copyright (C) 2016 - 2023 The ALF project
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
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
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
!>
!> @brief
!> This module provides a set of predefined observables
!>
!
!--------------------------------------------------------------------

module Predefined_Obs

   use runtime_error_mod
   use Operator_mod
   use Observables
   use Lattices_v3
   use entanglement_mod
   use iso_fortran_env, only: output_unit, error_unit

   implicit none

   interface Predefined_Obs_scal_Renyi_Ent
      module procedure Predefined_Obs_scal_Renyi_Ent_gen_all, Predefined_Obs_scal_Renyi_Ent_indep, &
      & Predefined_Obs_scal_Renyi_Ent_gen_fl
   end interface
   interface Predefined_Obs_scal_Mutual_Inf
      module procedure Predefined_Obs_scal_Mutual_Inf_indep, Predefined_Obs_scal_Mutual_Inf_gen_fl, &
      & Predefined_Obs_scal_Mutual_Inf_gen_all
   end interface

contains
!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure SU(N) spin-spin coorelations
!>       If N_FL = 1 then  this routine returns
!>       \frac{2N}{N^2-1}\sum_{a=1}{N^2 - 1 }  <c^{dag}_i T^a c_i  c^{dag}_j T^a  c_j>  where  T^a are th generators of SU(N)
!>       satisfying the normalization condition:  Tr [ T^a  T^b ]= \delta_{a,b}/2
!>       Note that for SU(N) symmetry:
!>       \sum_{a=1}{N^2 - 1 }  <c^{dag}_i T^a c_i  c^{dag}_j T^a  c_j>   = \sum_{a} Tr{T^a T^a} < c^{dag}_i c_j> < c_i c^{dag}_j> =
!>       (N^2 -1)/2 < c^{dag}_i c_j> < c_i c^{dag}_j>
!--------------------------------------------------------------------
   subroutine Predefined_Obs_eq_SpinSUN_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :)
      complex(Kind=kind(0.d0)), intent(In) :: GR(:, :, :), GRC(:, :, :), ZS, ZP
      type(Obser_Latt), intent(inout)   :: Obs

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj
      complex(Kind=kind(0.d0)) :: ZZ

      if (size(List, 1) .ne. size(GR, 1) .or. size(List, 2) .ne. 2) then
         write (error_unit, *) 'List in  Predefined_Obs_eq_SpinSUN_measure has  wrong  dim.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if
      if (Obs%File_Latt .ne. "SpinZ") then
         write (error_unit, *) 'Predefined_Obs_eq_SpinSUN_measure: Wrong filename'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if
      ! Count and average sign
      Obs%N = Obs%N + 1
      Obs%Ave_sign = Obs%Ave_sign + real(ZS, kind(0.d0))

      ! Measure
      N_FL = size(GR, 3)
      if (N_FL == 1) then
         do I1 = 1, size(List, 1)
            I = List(I1, 1)
            if (I > 0) then
               no_I = List(I1, 2)
               do J1 = 1, size(List, 1)
                  J = List(J1, 1)
                  if (J > 0) then
                     no_J = List(J1, 2)
                     imj = latt%imj(I, J)
                     ZZ = GRC(I1, J1, 1)*GR(I1, J1, 1)*cmplx(dble(N_SUN), 0.d0, kind(0.d0))
                     Obs%Obs_Latt(imj, 1, no_I, no_J) = Obs%Obs_Latt(imj, 1, no_I, no_J) + ZZ*ZP*ZS
                  end if
               end do
            end if
         end do
      end if

   end subroutine Predefined_Obs_eq_SpinSUN_measure

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!>  Measure spin-spin coorelations.  SpinXY correlations.
!>   a)  If N_FL = 2 and N_SUN = 1  the routine returns:
!>
!>     SpinZ =  4 * <c^{dag}_i S^z c_i  c^{dag}_j S^z  c_j>
!>     SpinXY=  2 ( <c^{dag}_i S^x c_i  c^{dag}_j S^x  c_j> +   <c^{dag}_i S^y c_i  c^{dag}_j S^y  c_j>)
!>     SpinT =  (2 * SpinXY  +  SpinZ)/3
!--------------------------------------------------------------------
   subroutine Predefined_Obs_eq_SpinMz_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, ObsZ, ObsXY, ObsXYZ)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :)
      complex(Kind=kind(0.d0)), intent(In) :: GR(:, :, :), GRC(:, :, :), ZS, ZP
      type(Obser_Latt), intent(inout)   :: ObsZ, ObsXY, ObsXYZ

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj
      complex(Kind=kind(0.d0)) :: ZXY, ZZ

      if (size(List, 1) .ne. size(GR, 1) .or. size(List, 2) .ne. 2) then
         write (error_unit, *) 'List in Predefined_Obs_eq_SpinMz_measure  has  wrong  dim.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      if (ObsZ%File_Latt .ne. "SpinZ" .and. ObsXY%File_Latt .ne. "SpinXY" .and.  &
         & ObsXYZ%File_Latt .ne. "SpinT") then
         write (error_unit, *) 'Predefined_Obs_eq_SpinMz_measure: Wrong filename'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if
      ! Count and average sign
      ObsZ%N = ObsZ%N + 1
      ObsZ%Ave_sign = ObsZ%Ave_sign + real(ZS, kind(0.d0))
      ObsXY%N = ObsXY%N + 1
      ObsXY%Ave_sign = ObsXY%Ave_sign + real(ZS, kind(0.d0))
      ObsXYZ%N = ObsXYZ%N + 1
      ObsXYZ%Ave_sign = ObsXYZ%Ave_sign + real(ZS, kind(0.d0))

      ! Measure
      N_FL = size(GR, 3)
      if (N_FL == 2 .and. N_SUN == 1) then
         do I1 = 1, size(List, 1)
            I = List(I1, 1)
            if (I > 0) then
               no_I = List(I1, 2)
               do J1 = 1, size(List, 1)
                  J = List(J1, 1)
                  if (J > 0) then
                     no_J = List(J1, 2)
                     imj = latt%imj(I, J)
                     ZXY = GRC(I1, J1, 1)*GR(I1, J1, 2) + GRC(I1, J1, 2)*GR(I1, J1, 1)
                     ZZ = GRC(I1, J1, 1)*GR(I1, J1, 1) + GRC(I1, J1, 2)*GR(I1, J1, 2) + &
                          (GRC(I1, I1, 2) - GRC(I1, I1, 1))*(GRC(J1, J1, 2) - GRC(J1, J1, 1))

                     ObsZ%Obs_Latt(imj, 1, no_I, no_J) = ObsZ%Obs_Latt(imj, 1, no_I, no_J) + ZZ*ZP*ZS
                     ObsXY%Obs_Latt(imj, 1, no_I, no_J) = ObsXY%Obs_Latt(imj, 1, no_I, no_J) + ZXY*ZP*ZS
                     ObsXYZ%Obs_Latt(imj, 1, no_I, no_J) = ObsXYZ%Obs_Latt(imj, 1, no_I, no_J) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                  end if
               end do
               ObsZ%Obs_Latt0(no_I) = ObsZ%Obs_Latt0(no_I) + (GRC(I1, I1, 2) - GRC(I1, I1, 1))*ZP*ZS
               ObsXYZ%Obs_Latt0(no_I) = ObsXYZ%Obs_Latt0(no_I) + (GRC(I1, I1, 2) - GRC(I1, I1, 1))*ZP*ZS/sqrt(3.d0)
            end if
         end do
      end if
   end subroutine Predefined_Obs_eq_SpinMz_measure

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure equal-time Green function
!>  \sum_{s=1}^{N_SUN} \sum_{nf=1}^{N_FL} < c^{dag}_{i,s,nf} c_{j,s,nf} >
!>
!--------------------------------------------------------------------
   subroutine Predefined_Obs_eq_Green_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :)
      complex(Kind=kind(0.d0)), intent(In) :: GR(:, :, :), GRC(:, :, :), ZS, ZP
      type(Obser_Latt), intent(inout)   :: Obs

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj, nf
      complex(Kind=kind(0.d0)) :: Z

      if (size(List, 1) .ne. size(GR, 1) .or. size(List, 2) .ne. 2) then
         write (error_unit, *) 'List in  Predefined_Obs_eq_Green_measure  has  wrong  dim.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      if (Obs%File_Latt .ne. "Green") then
         write (error_unit, *) 'Predefined_Obs_eq_Green_measure: Wrong filename'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if
      ! Count and average sign
      Obs%N = Obs%N + 1
      Obs%Ave_sign = Obs%Ave_sign + real(ZS, kind(0.d0))

      ! Measure
      N_FL = size(GR, 3)
      do I1 = 1, size(List, 1)
         I = List(I1, 1)
         if (I > 0) then
            no_I = List(I1, 2)
            do J1 = 1, size(List, 1)
               J = List(J1, 1)
               if (J > 0) then
                  no_J = List(J1, 2)
                  imj = latt%imj(I, J)
                  Z = cmplx(0.d0, 0.d0, kind(0.d0))
                  do nf = 1, N_FL
                     Z = Z + GRC(I1, J1, nf)
                  end do
                  Z = Z*cmplx(dble(N_SUN), 0.d0, kind(0.d0))
                  Obs%Obs_Latt(imj, 1, no_I, no_J) = Obs%Obs_Latt(imj, 1, no_I, no_J) + Z*ZP*ZS
               end if
            end do
         end if
      end do

   end subroutine Predefined_Obs_eq_Green_measure

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure density-density.
!>  Let N_i = \sum_{s=1}^{N_SUN} \sum_{nf=1}^{N_FL}  c^{dag}_{i,s,nf} c_{i,s,nf}
!>  Routine returns:
!>  <N_i  N_j >  -  <N_i> < N_j>
!--------------------------------------------------------------------
   subroutine Predefined_Obs_eq_Den_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, ZS, ZP, Obs)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :)
      complex(Kind=kind(0.d0)), intent(In) :: GR(:, :, :), GRC(:, :, :), ZS, ZP
      type(Obser_Latt), intent(inout)   :: Obs

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj, nf
      complex(Kind=kind(0.d0)) :: ZI, ZJ, Z

      if (size(List, 1) .ne. size(GR, 1) .or. size(List, 2) .ne. 2) then
         write (error_unit, *) 'List in  Predefined_Obs_eq_Den_measure  has  wrong  dim.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      if (Obs%File_Latt .ne. "Den") then
         write (error_unit, *) 'Predefined_Obs_eq_Den_measure: Wrong filename'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if
      ! Count and average sign
      Obs%N = Obs%N + 1
      Obs%Ave_sign = Obs%Ave_sign + real(ZS, kind(0.d0))

      ! Measure

      N_FL = size(GR, 3)
      do I1 = 1, size(List, 1)
         I = List(I1, 1)
         if (I > 0) then
            no_I = List(I1, 2)
            ZI = cmplx(0.d0, 0.d0, kind(0.d0))
            do nf = 1, N_FL
               ZI = ZI + GRC(I1, I1, nf)
            end do
            ZI = ZI*cmplx(dble(N_SUN), 0.d0, kind(0.d0))
            do J1 = 1, size(List, 1)
               J = List(J1, 1)
               if (J > 0) then
                  no_J = List(J1, 2)
                  imj = latt%imj(I, J)
                  ZJ = cmplx(0.d0, 0.d0, kind(0.d0))
                  do nf = 1, N_FL
                     ZJ = ZJ + GRC(J1, J1, nf)
                  end do
                  ZJ = ZJ*cmplx(dble(N_SUN), 0.d0, kind(0.d0))
                  Z = cmplx(0.d0, 0.d0, kind(0.d0))
                  do nf = 1, N_FL
                     Z = Z + GRC(I1, J1, nf)*GR(I1, J1, nf)
                  end do
                  Z = Z*cmplx(dble(N_SUN), 0.d0, kind(0.d0))
                  Obs%Obs_Latt(imj, 1, no_I, no_J) = Obs%Obs_Latt(imj, 1, no_I, no_J) + (ZI*ZJ + Z)*ZP*ZS
               end if
            end do
            Obs%Obs_Latt0(no_I) = Obs%Obs_Latt0(no_I) + ZI*ZP*ZS
         end if
      end do
   end subroutine Predefined_Obs_eq_Den_measure

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure time displaced Green function
!>  \sum_{s=1}^{N_SUN} \sum_{nf=1}^{N_FL} < c^{dag}_{i,s,nf}(tau) c_{j,s,nf} >
!>
!--------------------------------------------------------------------
   subroutine Predefined_Obs_tau_Green_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :), NT
      complex(Kind=kind(0.d0)), intent(In) :: GT0(:, :, :), G0T(:, :, :), G00(:, :, :), GTT(:, :, :), ZS, ZP
      type(Obser_Latt), intent(inout)   :: Obs

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj, nf
      complex(Kind=kind(0.d0)) :: Z

      if (size(List, 1) .ne. size(GT0, 1) .or. size(List, 2) .ne. 2) then
         write (error_unit, *) 'List in  Predefined_Obs_tau_Green_measure  has  wrong  dim.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      if (Obs%File_Latt .ne. "Green") then
         write (error_unit, *) 'Predefined_Obs_tau_Green_measure: Wrong filename'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      ! Count and average sign
      if (NT == 0) then
         Obs%N = Obs%N + 1
         Obs%Ave_sign = Obs%Ave_sign + real(ZS, kind(0.d0))
      end if

      ! Measure
      N_FL = size(GT0, 3)
      do I1 = 1, size(List, 1)
         I = List(I1, 1)
         if (I > 0) then
            no_I = List(I1, 2)
            do J1 = 1, size(List, 1)
               J = List(J1, 1)
               if (J > 0) then
                  no_J = List(J1, 2)
                  imj = latt%imj(I, J)
                  Z = cmplx(0.d0, 0.d0, kind(0.d0))
                  do nf = 1, N_FL
                     Z = Z + GT0(I1, J1, nf)
                  end do
                  Z = Z*cmplx(1.d0/dble(N_FL), 0.d0, kind(0.d0))
                  Obs%Obs_Latt(imj, NT + 1, no_I, no_J) = Obs%Obs_Latt(imj, NT + 1, no_I, no_J) + Z*ZP*ZS
               end if
            end do
         end if
      end do

   end subroutine Predefined_Obs_tau_Green_measure

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure time displaced Spin-Spin correlations  for SU(N) models (N_FL = 1)
!>
!>  Returns:  \frac{2N}{N^2-1}\sum_{a=1}{N^2 - 1 }  <c^{dag}_i(tau) T^a c_i (tau)  c^{dag}_j T^a  c_j>  where  T^a are th generators of SU(N)
!>
!--------------------------------------------------------------------
   subroutine Predefined_Obs_tau_SpinSUN_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :), NT
      complex(Kind=kind(0.d0)), intent(In) :: GT0(:, :, :), G0T(:, :, :), G00(:, :, :), GTT(:, :, :), ZS, ZP
      type(Obser_Latt), intent(inout)   :: Obs

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj, nf
      complex(Kind=kind(0.d0)) :: Z

      if (Obs%File_Latt .ne. "SpinZ") then
         write (error_unit, *) 'Predefined_Obs_tau_SpinSUN_measure: Wrong filename'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if
      if (size(List, 1) .ne. size(GT0, 1) .or. size(List, 2) .ne. 2) then
         write (error_unit, *) 'List in  Predefined_Obs_tau_SpinSUN_measure  has  wrong  dim.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      ! Count and average sign
      if (NT == 0) then
         Obs%N = Obs%N + 1
         Obs%Ave_sign = Obs%Ave_sign + real(ZS, kind(0.d0))
      end if

      ! Measure
      N_FL = size(GT0, 3)
      Z = cmplx(dble(N_SUN), 0.d0, kind(0.d0))
      do I1 = 1, size(List, 1)
         I = List(I1, 1)
         if (I > 0) then
            no_I = List(I1, 2)
            do J1 = 1, size(List, 1)
               J = List(J1, 1)
               if (J > 0) then
                  no_J = List(J1, 2)
                  imj = latt%imj(I, J)
                  Obs%Obs_Latt(imj, NT + 1, no_I, no_J) = Obs%Obs_Latt(imj, NT + 1, no_I, no_J) &
                       &  - Z*G0T(J1, I1, 1)*GT0(I1, J1, 1)*ZP*ZS
               end if
            end do
         end if
      end do

   end subroutine Predefined_Obs_tau_SpinSUN_measure

!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure time displaced Spin-Spin correlations  for Mz models (N_FL = 2, N_SUN = 1)
!>  Returns:
!>     SpinZ =  4 * <c^{dag}_i(tau) S^z c_i(tau)  c^{dag}_j S^z  c_j >
!>     SpinXY=  2 ( <c^{dag}_i(tau) S^x c_i(tau)  c^{dag}_j S^x  c_j > +   <c^{dag}_i(tau) S^y c_i(tau)  c^{dag}_j S^y  c_j>)
!>     SpinT =  (2 * SpinXY  +  SpinZ)/3
!>
!--------------------------------------------------------------------
   subroutine Predefined_Obs_tau_SpinMz_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, ObsZ, ObsXY, ObsXYZ)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :), NT
      complex(Kind=kind(0.d0)), intent(In) :: GT0(:, :, :), G0T(:, :, :), G00(:, :, :), GTT(:, :, :), ZS, ZP
      type(Obser_Latt), intent(inout)   :: ObsZ, ObsXY, ObsXYZ

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj, nf
      complex(Kind=kind(0.d0)) :: ZZ, ZXY

      if (size(List, 1) .ne. size(GT0, 1) .or. size(List, 2) .ne. 2) then
         write (error_unit, *) 'List in Predefined_Obs_tau_SpinMz_measure has  wrong  dim.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      ! Count and average sign
      if (ObsZ%File_Latt .ne. "SpinZ" .and. ObsXY%File_Latt .ne. "SpinXY" .and.  &
         & ObsXYZ%File_Latt .ne. "SpinT") then
         write (error_unit, *) 'Predefined_Obs_tau_SpinMz_measure: Wrong filename'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if
      if (NT == 0) then
         ! Count and average sign
         ObsZ%N = ObsZ%N + 1
         ObsZ%Ave_sign = ObsZ%Ave_sign + real(ZS, kind(0.d0))
         ObsXY%N = ObsXY%N + 1
         ObsXY%Ave_sign = ObsXY%Ave_sign + real(ZS, kind(0.d0))
         ObsXYZ%N = ObsXYZ%N + 1
         ObsXYZ%Ave_sign = ObsXYZ%Ave_sign + real(ZS, kind(0.d0))
      end if

      ! Measure
      N_FL = size(GT0, 3)
      do I1 = 1, size(List, 1)
         I = List(I1, 1)
         if (I > 0) then
            no_I = List(I1, 2)
            do J1 = 1, size(List, 1)
               J = List(J1, 1)
               if (J > 0) then
                  no_J = List(J1, 2)
                  imj = latt%imj(I, J)
                  ZZ = (GTT(I1, I1, 1) - GTT(I1, I1, 2))*(G00(J1, J1, 1) - G00(J1, J1, 2))   &
                       &    - G0T(J1, I1, 1)*GT0(I1, J1, 1) - G0T(J1, I1, 2)*GT0(I1, J1, 2)
                  ZXY = -G0T(J1, I1, 1)*GT0(I1, J1, 2) - G0T(J1, I1, 2)*GT0(I1, J1, 1)
                  ObsZ%Obs_Latt(imj, NT + 1, no_I, no_J) = ObsZ%Obs_Latt(imj, NT + 1, no_I, no_J) + ZZ*ZP*ZS
                  ObsXY%Obs_Latt(imj, NT + 1, no_I, no_J) = ObsXY%Obs_Latt(imj, NT + 1, no_I, no_J) + ZXY*ZP*ZS
                  ObsXYZ%Obs_Latt(imj, NT + 1, no_I, no_J) = ObsXYZ%Obs_Latt(imj, NT + 1, no_I, no_J) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
               end if
            end do
            ObsZ%Obs_Latt0(no_I) = ObsZ%Obs_Latt0(no_I) + (GTT(I1, I1, 2) - GTT(I1, I1, 1))*ZP*ZS
            ObsXYZ%Obs_Latt0(no_I) = ObsXYZ%Obs_Latt0(no_I) + (GTT(I1, I1, 2) - GTT(I1, I1, 1))*ZP*ZS/sqrt(3.d0)
         end if
      end do

   end subroutine Predefined_Obs_tau_SpinMz_measure

!-------------------------------------------------------------------
!> @Author
!> ALF-project
!
!>  @brief
!>  Measure time displaced Den-Den correlations  for general  SU(N) models
!>  Let N_i = \sum_{s=1}^{N_SUN} \sum_{nf=1}^{N_FL}  c^{dag}_{i,s,nf} c_{i,s,nf}
!>  Routine returns:
!>        <N_i (tau)  N_j >  -  <N_i> < N_j>
!>
!--------------------------------------------------------------------
   subroutine Predefined_Obs_tau_Den_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, ZS, ZP, Obs)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :), NT
      complex(Kind=kind(0.d0)), intent(In) :: GT0(:, :, :), G0T(:, :, :), G00(:, :, :), GTT(:, :, :), ZS, ZP
      type(Obser_Latt), intent(inout)   :: Obs

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj, nf
      complex(Kind=kind(0.d0)) :: Z, ZI, ZJ

      if (size(List, 1) .ne. size(GT0, 1) .or. size(List, 2) .ne. 2) then
         write (error_unit, *) 'List in Predefined_Obs_tau_Den_measure has  wrong  dim.'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      if (Obs%File_Latt .ne. "Den") then
         write (error_unit, *) 'Predefined_Obs_tau_Den_measure: Wrong filename'
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      ! Count and average sign
      if (NT == 0) then
         Obs%N = Obs%N + 1
         Obs%Ave_sign = Obs%Ave_sign + real(ZS, kind(0.d0))
      end if
      N_FL = size(GT0, 3)
      !Measure
      do I1 = 1, size(List, 1)
         I = List(I1, 1)
         if (I > 0) then
            no_I = List(I1, 2)
            ZI = cmplx(0.d0, 0.d0, kind(0.d0))
            do nf = 1, N_FL
               ZI = ZI + cmplx(1.d0, 0.d0, kind(0.d0)) - GTT(I1, I1, nf)
            end do
            ZI = ZI*cmplx(dble(N_SUN), 0.d0, kind(0.d0))
            do J1 = 1, size(List, 1)
               J = List(J1, 1)
               if (J > 0) then
                  no_J = List(J1, 2)
                  imj = latt%imj(I, J)
                  ZJ = cmplx(0.d0, 0.d0, kind(0.d0))
                  do nf = 1, N_FL
                     ZJ = ZJ + cmplx(1.d0, 0.d0, kind(0.d0)) - G00(J1, J1, nf)
                  end do
                  ZJ = ZJ*cmplx(dble(N_SUN), 0.d0, kind(0.d0))
                  Z = cmplx(0.d0, 0.d0, kind(0.d0))
                  do nf = 1, N_FL
                     Z = Z - G0T(J1, I1, nf)*GT0(I1, J1, nf)
                  end do
                  Z = Z*cmplx(dble(N_SUN), 0.d0, kind(0.d0))
                  Obs%Obs_Latt(imj, NT + 1, no_I, no_J) = Obs%Obs_Latt(imj, NT + 1, no_I, no_J) + (ZI*ZJ + Z)*ZP*ZS
               end if
            end do
            Obs%Obs_Latt0(no_I) = Obs%Obs_Latt0(no_I) + ZI*ZP*ZS
         end if
      end do

   end subroutine Predefined_Obs_tau_Den_measure

#include  "Cotunneling_dimer_obs.F90"
!-------------------------------------------------------------------
!> @Author
!> ALF-project
!
!>  @brief
!>  Let  OP_V  be  a  type  two  operator.  Then  this  function  returns:
!>
!>  Routine returns:
!>        << ( \sum_{s=1}^{N_FL} \sum_{sigma=1}^{N_SUN} [\sum_{x,y}( c^dag_{x,s,sig} V^{s}_{x,y} c^dag_{y,s,sig} ) + \alpha_s] )^2 >>
!>
!--------------------------------------------------------------------
   complex(Kind=kind(0.d0)) function Predefined_Obs_V_Int(OP_Vint, GR, GRC, N_SUN)

      implicit none
      type(operator), intent(In) :: OP_Vint(:)
      integer, intent(In) :: N_SUN
      complex(Kind=kind(0.d0)), intent(In) :: GR(:, :, :), GRC(:, :, :)

      ! Local
      complex(Kind=kind(0.d0))  ::  Z, Z_tmp, ZC
      integer  ::  N_FL, N, I, J, I1, J1, nf
      real(Kind=kind(0.d0))  ::   Zero = 1.0d-16

      if (OP_Vint(1)%type .ne. 2) then
         write (error_unit, *) 'Predefined_Obs_V_Int  routine is   defined  fro  tpye  2  vertices.  '
         call Terminate_on_error(ERROR_GENERIC, __FILE__, __LINE__)
      end if

      N_FL = size(GR, 3)
      N = size(OP_Vint(1)%O, 1)

      Z = cmplx(0.d0, 0.d0, kind(0.d0))
      do nf = 1, N_FL
         Z_tmp = cmplx(0.d0, 0.d0, kind(0.d0))
         do J = 1, N
            do I = 1, N
               Z_tmp = Z_tmp + OP_Vint(nf)%O(i, j)*GRC(OP_Vint(nf)%P(I), OP_Vint(nf)%P(J), nf)
            end do
         end do
         Z = Z + Z_tmp + OP_Vint(nf)%alpha
      end do
      Z = Z*dble(N_SUN)

      ZC = cmplx(0.d0, 0.d0, kind(0.d0))
      do nf = 1, N_FL
         do J = 1, N
            do I = 1, N
               Z_tmp = OP_Vint(nf)%O(I, J)
               if (real(Z_tmp*conjg(Z_tmp), kind(0.d0)) > Zero) then
                  do J1 = 1, N
                     do I1 = 1, N
                        ZC = ZC + Z_tmp*OP_Vint(nf)%O(I1, J1)* &
                             &        GRC(OP_Vint(nf)%P(I), OP_Vint(nf)%P(J1), nf)*&
                             &         GR(OP_Vint(nf)%P(J), OP_Vint(nf)%P(I1), nf)
                     end do
                  end do
               end if
            end do
         end do
      end do
      ZC = ZC*dble(N_SUN)
      Predefined_Obs_V_Int = Z*Z + ZC

   end function Predefined_Obs_V_Int

   subroutine Predefined_Obs_scal_Renyi_Ent_indep(GRC, List, Nsites, N_SUN, ZS, ZP, Obs)

      implicit none
      integer, dimension(:), intent(IN)     :: List
      integer, intent(IN)                   :: Nsites, N_SUN
      complex(kind=kind(0.d0)), intent(IN) :: GRC(:, :, :), ZS, ZP
      type(Obser_Vec), intent(inout)   :: Obs

      complex(kind=kind(0.d0))  :: Renyi

      Renyi = Calc_Renyi_Ent_indep(GRC, List, Nsites, N_SUN)

      Obs%N = Obs%N + 1
      Obs%Ave_sign = Obs%Ave_sign + real(ZS, kind(0.d0))
      Obs%Obs_vec(1) = Obs%Obs_vec(1) + Renyi*ZP*ZS

   end subroutine Predefined_Obs_scal_Renyi_Ent_indep

   subroutine Predefined_Obs_scal_Renyi_Ent_gen_fl(GRC, List, Nsites, N_SUN, ZS, ZP, Obs)

      implicit none
      integer, dimension(:, :), intent(IN)   :: List
      integer, dimension(:), intent(IN)     :: Nsites, N_SUN
      complex(kind=kind(0.d0)), intent(IN) :: GRC(:, :, :), ZS, ZP
      type(Obser_Vec), intent(inout)   :: Obs

      complex(kind=kind(0.d0))  :: Renyi

      Renyi = Calc_Renyi_Ent_gen_fl(GRC, List, Nsites, N_SUN)

      Obs%N = Obs%N + 1
      Obs%Ave_sign = Obs%Ave_sign + real(ZS, kind(0.d0))
      Obs%Obs_vec(1) = Obs%Obs_vec(1) + Renyi*ZP*ZS

   end subroutine Predefined_Obs_scal_Renyi_Ent_gen_fl

   subroutine Predefined_Obs_scal_Renyi_Ent_gen_all(GRC, List, Nsites, ZS, ZP, Obs)

      implicit none
      integer, dimension(:, :, :), intent(IN)   :: List
      integer, dimension(:, :), intent(IN)     :: Nsites
      complex(kind=kind(0.d0)), intent(IN) :: GRC(:, :, :), ZS, ZP
      type(Obser_Vec), intent(inout)   :: Obs

      complex(kind=kind(0.d0))  :: Renyi

      Renyi = Calc_Renyi_Ent_gen_all(GRC, List, Nsites)

      Obs%N = Obs%N + 1
      Obs%Ave_sign = Obs%Ave_sign + real(ZS, kind(0.d0))
      Obs%Obs_vec(1) = Obs%Obs_vec(1) + Renyi*ZP*ZS

   end subroutine Predefined_Obs_scal_Renyi_Ent_gen_all

   subroutine Predefined_Obs_scal_Mutual_Inf_indep(GRC, List_A, Nsites_A, List_B, Nsites_B, N_SUN, ZS, ZP, Obs)

      implicit none
      integer, dimension(:), intent(IN)     :: List_A, List_B
      integer, intent(IN)                   :: Nsites_A, Nsites_B, N_SUN
      complex(kind=kind(0.d0)), intent(IN) :: GRC(:, :, :), ZS, ZP
      type(Obser_Vec), intent(inout)   :: Obs

      complex(kind=kind(0.d0))  :: Renyi_A, Renyi_B, Renyi_AB

      call Calc_Mutual_Inf_indep(GRC, List_A, Nsites_A, List_B, Nsites_B, N_SUN, Renyi_A, Renyi_B, Renyi_AB)

      Obs%N = Obs%N + 1
      Obs%Ave_sign = Obs%Ave_sign + real(ZS, kind(0.d0))
      Obs%Obs_vec(1) = Obs%Obs_vec(1) + Renyi_A*ZP*ZS
      Obs%Obs_vec(2) = Obs%Obs_vec(2) + Renyi_B*ZP*ZS
      Obs%Obs_vec(3) = Obs%Obs_vec(3) + Renyi_AB*ZP*ZS

   end subroutine Predefined_Obs_scal_Mutual_Inf_indep

   subroutine Predefined_Obs_scal_Mutual_Inf_gen_fl(GRC, List_A, Nsites_A, List_B, Nsites_B, N_SUN, ZS, ZP, Obs)

      implicit none
      integer, dimension(:, :), intent(IN)   :: List_A, List_B
      integer, dimension(:), intent(IN)     :: Nsites_A, Nsites_B, N_SUN
      complex(kind=kind(0.d0)), intent(IN) :: GRC(:, :, :), ZS, ZP
      type(Obser_Vec), intent(inout)   :: Obs

      complex(kind=kind(0.d0))  :: Renyi_A, Renyi_B, Renyi_AB

      call Calc_Mutual_Inf_gen_fl(GRC, List_A, Nsites_A, List_B, Nsites_B, N_SUN, Renyi_A, Renyi_B, Renyi_AB)

      Obs%N = Obs%N + 1
      Obs%Ave_sign = Obs%Ave_sign + real(ZS, kind(0.d0))
      Obs%Obs_vec(1) = Obs%Obs_vec(1) + Renyi_A*ZP*ZS
      Obs%Obs_vec(2) = Obs%Obs_vec(2) + Renyi_B*ZP*ZS
      Obs%Obs_vec(3) = Obs%Obs_vec(3) + Renyi_AB*ZP*ZS

   end subroutine Predefined_Obs_scal_Mutual_Inf_gen_fl

   subroutine Predefined_Obs_scal_Mutual_Inf_gen_all(GRC, List_A, Nsites_A, List_B, Nsites_B, ZS, ZP, Obs)

      implicit none
      integer, dimension(:, :, :), intent(IN)   :: List_A, List_B
      integer, dimension(:, :), intent(IN)     :: Nsites_A, Nsites_B
      complex(kind=kind(0.d0)), intent(IN) :: GRC(:, :, :), ZS, ZP
      type(Obser_Vec), intent(inout)   :: Obs

      complex(kind=kind(0.d0))  :: Renyi_A, Renyi_B, Renyi_AB

      call Calc_Mutual_Inf_gen_all(GRC, List_A, Nsites_A, List_B, Nsites_B, Renyi_A, Renyi_B, Renyi_AB)

      Obs%N = Obs%N + 1
      Obs%Ave_sign = Obs%Ave_sign + real(ZS, kind(0.d0))
      Obs%Obs_vec(1) = Obs%Obs_vec(1) + Renyi_A*ZP*ZS
      Obs%Obs_vec(2) = Obs%Obs_vec(2) + Renyi_B*ZP*ZS
      Obs%Obs_vec(3) = Obs%Obs_vec(3) + Renyi_AB*ZP*ZS

   end subroutine Predefined_Obs_scal_Mutual_Inf_gen_all

end module Predefined_Obs
