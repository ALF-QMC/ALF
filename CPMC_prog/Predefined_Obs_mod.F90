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


    Module Predefined_Obs

      use runtime_error_mod
      Use Observables
      Use Lattices_v3
      Use entanglement_mod
      use iso_fortran_env, only: output_unit, error_unit

      Implicit none

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
      Subroutine Predefined_Obs_eq_SpinSUN_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZP, Re_ZW, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:)
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:), ZP, Re_ZW
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj
        Complex (Kind=Kind(0.d0)) :: ZZ

        If ( Obs%File_Latt .ne. "SpinZ" )   then
           Write(error_unit,*) 'Predefined_Obs_eq_SpinSUN_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        ! Count and average sign
        Obs%N          = Obs%N + 1
        Obs%sum_weight = Obs%sum_weight + Re_ZW

        ! Measure
        N_FL = Size(GR,3)
        If (N_FL == 1)   then
           Do I1 = 1,Size(List,1)
              I    = List(I1,1)
              If  ( I > 0 ) then 
                 no_I = List(I1,2)
                 Do J1 = 1,Size(List,1)
                    J    = List(J1,1)
                    if  (  J  > 0 )  then 
                       no_J = List(J1,2)
                       imj  = latt%imj(I,J)
                       ZZ   = GRC(I1,J1,1) * GR(I1,J1,1) * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                       Obs%Obs_Latt(imj,1,no_I,no_J) =  Obs%Obs_Latt(imj,1,no_I,no_J) + ZZ*Re_ZW
                    endif
                 enddo
              endif
           enddo
        endif


      end Subroutine Predefined_Obs_eq_SpinSUN_measure

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
      Subroutine Predefined_Obs_eq_SpinMz_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZP, Re_ZW, ObsZ, ObsXY, ObsXYZ )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:)
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:), ZP, Re_ZW
        Type (Obser_Latt),    Intent(inout)   :: ObsZ, ObsXY, ObsXYZ

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj
        Complex (Kind=Kind(0.d0)) :: ZXY, ZZ

        If ( ObsZ%File_Latt .ne. "SpinZ" .and. ObsXY%File_Latt .ne. "SpinXY" .and.  &
           & ObsXYZ%File_Latt .ne. "SpinT"  )   then
           Write(error_unit,*) 'Predefined_Obs_eq_SpinMz_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        ! Count and average sign
        ObsZ%N            = ObsZ%N + 1
        ObsZ%sum_weight   = ObsZ%sum_weight   + Re_ZW
        ObsXY%N           = ObsXY%N + 1
        ObsXY%sum_weight  = ObsXY%sum_weight  + Re_ZW
        ObsXYZ%N          = ObsXYZ%N + 1
        ObsXYZ%sum_weight = ObsXYZ%sum_weight + Re_ZW

        ! Measure
        N_FL = Size(GR,3)
        If (N_FL == 2 .and. N_SUN == 1 ) Then
           Do I1 = 1,Size(List,1) 
              I    = List(I1,1)
              If (I  > 0 )  then 
                 no_I = List(I1,2)
                 Do J1 = 1,Size(List,1) 
                    J    = List(J1,1)
                    If (J > 0 )  then 
                       no_J = List(J1,2)
                       imj  = latt%imj(I,J)
                       ZXY  = GRC(I1,J1,1) * GR(I1,J1,2) +  GRC(I1,J1,2) * GR(I1,J1,1)
                       ZZ   = GRC(I1,J1,1) * GR(I1,J1,1) +  GRC(I1,J1,2) * GR(I1,J1,2)    + &
                            (GRC(I1,I1,2) - GRC(I1,I1,1))*(GRC(J1,J1,2) - GRC(J1,J1,1))
                       
                       ObsZ  %Obs_Latt(imj,1,no_I,no_J) =  ObsZ  %Obs_Latt(imj,1,no_I,no_J) +  ZZ  *Re_ZW
                       ObsXY %Obs_Latt(imj,1,no_I,no_J) =  ObsXY %Obs_Latt(imj,1,no_I,no_J) +  ZXY *Re_ZW
                       ObsXYZ%Obs_Latt(imj,1,no_I,no_J) =  ObsXYZ%Obs_Latt(imj,1,no_I,no_J) + (2.d0*ZXY + ZZ)*Re_ZW/3.d0
                    endif
                 enddo
                 ObsZ  %Obs_Latt0(no_I) =  ObsZ  %Obs_Latt0(no_I) +   (GRC(I1,I1,2) - GRC(I1,I1,1)) * Re_ZW
                 ObsXYZ%Obs_Latt0(no_I) =  ObsXYZ%Obs_Latt0(no_I) +   (GRC(I1,I1,2) - GRC(I1,I1,1)) * Re_ZW/sqrt(3.d0)
              endif
           enddo
        endif
      End Subroutine Predefined_Obs_eq_SpinMz_measure


!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure equal-time Green function
!>  \sum_{s=1}^{N_SUN} \sum_{nf=1}^{N_FL} < c^{dag}_{i,s,nf} c_{j,s,nf} >
!>
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_eq_Green_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZP, Re_ZW, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:)
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:), ZP, Re_ZW
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: Z

        If ( Obs%File_Latt .ne. "Green" )   then
           Write(error_unit,*) 'Predefined_Obs_eq_Green_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        ! Count and average sign
        Obs%N          = Obs%N + 1
        Obs%sum_weight = Obs%sum_weight + Re_ZW

        ! Measure
        N_FL = Size(GR,3)
        Do I1 = 1,Size(List,1)
           I    = List(I1,1)
           If ( I  > 0 ) then 
              no_I = List(I1,2)
              Do J1 = 1,Size(List,1)
                 J    = List(J1,1)
                 If  ( J > 0 ) then 
                    no_J = List(J1,2)
                    imj  = latt%imj(I,J)
                    Z = cmplx(0.d0,0.d0,Kind(0.d0))
                    Do nf = 1, N_FL
                       Z = Z +  GRC(I1,J1,nf)
                    enddo
                    Z   = Z * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                    Obs%Obs_Latt(imj,1,no_I,no_J) =  Obs%Obs_Latt(imj,1,no_I,no_J) + Z*Re_ZW
                 endif
              enddo
           endif
        enddo

      end Subroutine Predefined_Obs_eq_Green_measure

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
      Subroutine Predefined_Obs_eq_Den_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZP, Re_ZW, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:)
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:), ZP, Re_ZW
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: ZI, ZJ, Z

        If ( Obs%File_Latt .ne. "Den" )   then
           Write(error_unit,*) 'Predefined_Obs_eq_Den_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        ! Count and average sign
        Obs%N          = Obs%N + 1
        Obs%sum_weight = Obs%sum_weight + Re_ZW

        ! Measure

        N_FL = Size(GR,3)
        Do I1 = 1,Size(List,1)
           I    = List(I1,1)
           If (I > 0 )  then 
              no_I = List(I1,2)
              ZI = cmplx(0.d0,0.d0,Kind(0.d0))
              Do nf = 1,N_FL
                 ZI  = ZI + GRC(I1,I1,nf)
              enddo
              ZI  = ZI * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
              Do J1 = 1,Size(List,1)
                 J    = List(J1,1)
                 If  (J > 0 ) then 
                    no_J = List(J1,2)
                    imj  = latt%imj(I,J)
                    ZJ = cmplx(0.d0,0.d0,Kind(0.d0))
                    Do nf = 1,N_FL
                       ZJ  = ZJ + GRC(J1,J1,nf)
                    enddo
                    ZJ  = ZJ * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                    Z   = cmplx(0.d0,0.d0,Kind(0.d0))
                    Do nf = 1,N_FL
                       Z = Z + GRC(I1,J1,nf) * GR(I1,J1,nf)
                    Enddo
                    Z   =  Z * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                    Obs%Obs_Latt(imj,1,no_I,no_J) =  Obs%Obs_Latt(imj,1,no_I,no_J) + (ZI*ZJ + Z)*Re_ZW
                 endif
              enddo
              Obs%Obs_Latt0(no_I) =  Obs%Obs_Latt0(no_I) +  ZI * Re_ZW
           endif
        enddo
      end Subroutine Predefined_Obs_eq_Den_measure
      
#include  "Cotunneling_dimer_obs.F90"
      
    End Module Predefined_Obs
