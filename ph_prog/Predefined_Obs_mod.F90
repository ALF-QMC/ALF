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
   use Observables
   use Lattices_v3
   use entanglement_mod
   use iso_fortran_env, only: output_unit, error_unit

   implicit none

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
   subroutine Predefined_Obs_eq_SpinSUN_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, Z_fac, Obs)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :)
      complex(Kind=kind(0.d0)), intent(In) :: GR(:, :, :), GRC(:, :, :), Z_fac
      type(Obser_Latt), intent(inout)   :: Obs

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj
      complex(Kind=kind(0.d0)) :: ZZ

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
                     Obs%Obs_Latt(imj, 1, no_I, no_J) = Obs%Obs_Latt(imj, 1, no_I, no_J) + ZZ*Z_fac
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
   subroutine Predefined_Obs_eq_SpinMz_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, Z_fac, ObsZ, ObsXY, ObsXYZ)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :)
      complex(Kind=kind(0.d0)), intent(In) :: GR(:, :, :), GRC(:, :, :), Z_fac
      type(Obser_Latt), intent(inout)   :: ObsZ, ObsXY, ObsXYZ

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj
      complex(Kind=kind(0.d0)) :: ZXY, ZZ

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

                     ObsZ%Obs_Latt(imj, 1, no_I, no_J) = ObsZ%Obs_Latt(imj, 1, no_I, no_J) + ZZ*Z_fac
                     ObsXY%Obs_Latt(imj, 1, no_I, no_J) = ObsXY%Obs_Latt(imj, 1, no_I, no_J) + ZXY*Z_fac
                     ObsXYZ%Obs_Latt(imj, 1, no_I, no_J) = ObsXYZ%Obs_Latt(imj, 1, no_I, no_J) + (2.d0*ZXY + ZZ)*Z_fac/3.d0
                  end if
               end do
               ObsZ%Obs_Latt0(no_I) = ObsZ%Obs_Latt0(no_I) + (GRC(I1, I1, 2) - GRC(I1, I1, 1))*Z_fac
               ObsXYZ%Obs_Latt0(no_I) = ObsXYZ%Obs_Latt0(no_I) + (GRC(I1, I1, 2) - GRC(I1, I1, 1))*Z_fac/sqrt(3.d0)
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
   subroutine Predefined_Obs_eq_Green_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, Z_fac, Obs)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :)
      complex(Kind=kind(0.d0)), intent(In) :: GR(:, :, :), GRC(:, :, :), Z_fac
      type(Obser_Latt), intent(inout)   :: Obs

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj, nf
      complex(Kind=kind(0.d0)) :: Z

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
                  Obs%Obs_Latt(imj, 1, no_I, no_J) = Obs%Obs_Latt(imj, 1, no_I, no_J) + Z*Z_fac
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
   subroutine Predefined_Obs_eq_Den_measure(Latt, Latt_unit, List, GR, GRC, N_SUN, Z_fac, Obs)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :)
      complex(Kind=kind(0.d0)), intent(In) :: GR(:, :, :), GRC(:, :, :), Z_fac
      type(Obser_Latt), intent(inout)   :: Obs

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj, nf
      complex(Kind=kind(0.d0)) :: ZI, ZJ, Z

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
                  Obs%Obs_Latt(imj, 1, no_I, no_J) = Obs%Obs_Latt(imj, 1, no_I, no_J) + (ZI*ZJ + Z)*Z_fac
               end if
            end do
            Obs%Obs_Latt0(no_I) = Obs%Obs_Latt0(no_I) + ZI*Z_fac
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
   subroutine Predefined_Obs_tau_Green_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, Z_fac, Obs)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :), NT
      complex(Kind=kind(0.d0)), intent(In) :: GT0(:, :, :), G0T(:, :, :), G00(:, :, :), GTT(:, :, :), Z_fac
      type(Obser_Latt), intent(inout)   :: Obs

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj, nf
      complex(Kind=kind(0.d0)) :: Z

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
                  Obs%Obs_Latt(imj, NT + 1, no_I, no_J) = Obs%Obs_Latt(imj, NT + 1, no_I, no_J) + Z*Z_fac
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
!>  Measure time displaced Spin-Spin correlations  for Mz models (N_FL = 2, N_SUN = 1)
!>  Returns:
!>     SpinZ =  4 * <c^{dag}_i(tau) S^z c_i(tau)  c^{dag}_j S^z  c_j >
!>     SpinXY=  2 ( <c^{dag}_i(tau) S^x c_i(tau)  c^{dag}_j S^x  c_j > +   <c^{dag}_i(tau) S^y c_i(tau)  c^{dag}_j S^y  c_j>)
!>     SpinT =  (2 * SpinXY  +  SpinZ)/3
!>
!--------------------------------------------------------------------
   subroutine Predefined_Obs_tau_SpinMz_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, Z_fac, ObsZ, ObsXY, ObsXYZ)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :), NT
      complex(Kind=kind(0.d0)), intent(In) :: GT0(:, :, :), G0T(:, :, :), G00(:, :, :), GTT(:, :, :), Z_fac
      type(Obser_Latt), intent(inout)   :: ObsZ, ObsXY, ObsXYZ

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj, nf
      complex(Kind=kind(0.d0)) :: ZZ, ZXY

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
                  ObsZ%Obs_Latt(imj, NT + 1, no_I, no_J) = ObsZ%Obs_Latt(imj, NT + 1, no_I, no_J) + ZZ*Z_fac
                  ObsXY%Obs_Latt(imj, NT + 1, no_I, no_J) = ObsXY%Obs_Latt(imj, NT + 1, no_I, no_J) + ZXY*Z_fac
                  ObsXYZ%Obs_Latt(imj, NT + 1, no_I, no_J) = ObsXYZ%Obs_Latt(imj, NT + 1, no_I, no_J) + (2.d0*ZXY + ZZ)*Z_fac/3.d0
               end if
            end do
            ObsZ%Obs_Latt0(no_I) = ObsZ%Obs_Latt0(no_I) + (GTT(I1, I1, 2) - GTT(I1, I1, 1))*Z_fac
            ObsXYZ%Obs_Latt0(no_I) = ObsXYZ%Obs_Latt0(no_I) + (GTT(I1, I1, 2) - GTT(I1, I1, 1))*Z_fac/sqrt(3.d0)
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
   subroutine Predefined_Obs_tau_Den_measure(Latt, Latt_unit, List, NT, GT0, G0T, G00, GTT, N_SUN, Z_fac, Obs)

      type(Lattice), intent(in)      :: Latt
      type(Unit_cell), intent(in)      :: Latt_unit
      integer, intent(In)      :: N_SUN, LIST(:, :), NT
      complex(Kind=kind(0.d0)), intent(In) :: GT0(:, :, :), G0T(:, :, :), G00(:, :, :), GTT(:, :, :), Z_fac
      type(Obser_Latt), intent(inout)   :: Obs

      ! Local
      integer :: N_FL, I, I1, J, J1, no_I, no_J, imj, nf
      complex(Kind=kind(0.d0)) :: Z, ZI, ZJ

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
                  Obs%Obs_Latt(imj, NT + 1, no_I, no_J) = Obs%Obs_Latt(imj, NT + 1, no_I, no_J) + (ZI*ZJ + Z)*Z_fac
               end if
            end do
            Obs%Obs_Latt0(no_I) = Obs%Obs_Latt0(no_I) + ZI*Z_fac
         end if
      end do

   end subroutine Predefined_Obs_tau_Den_measure

end module Predefined_Obs
