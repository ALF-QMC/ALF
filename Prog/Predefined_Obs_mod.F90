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
      use Operator_mod
      Use Observables
      Use Lattices_v3
      Use entanglement_mod
      use iso_fortran_env, only: output_unit, error_unit

      Implicit none
      
      INTERFACE Predefined_Obs_scal_Renyi_Ent
        MODULE PROCEDURE Predefined_Obs_scal_Renyi_Ent_gen_all, Predefined_Obs_scal_Renyi_Ent_indep, &
        & Predefined_Obs_scal_Renyi_Ent_gen_fl
      END INTERFACE
      INTERFACE Predefined_Obs_scal_Mutual_Inf
        MODULE PROCEDURE Predefined_Obs_scal_Mutual_Inf_indep, Predefined_Obs_scal_Mutual_Inf_gen_fl, &
        & Predefined_Obs_scal_Mutual_Inf_gen_all
      END INTERFACE

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
      Subroutine Predefined_Obs_eq_SpinSUN_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:)
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj
        Complex (Kind=Kind(0.d0)) :: ZZ


        If ( Size(List,1) .ne. Size(GR,1) .or. Size(List,2) .ne. 2  )   then
           Write(error_unit,*) 'List in  Predefined_Obs_eq_SpinSUN_measure has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        If ( Obs%File_Latt .ne. "SpinZ" )   then
           Write(error_unit,*) 'Predefined_Obs_eq_SpinSUN_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        ! Count and average sign
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))

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
                       Obs%Obs_Latt(imj,1,no_I,no_J) =  Obs%Obs_Latt(imj,1,no_I,no_J) + ZZ*ZP*ZS
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
      Subroutine Predefined_Obs_eq_SpinMz_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, ObsZ, ObsXY, ObsXYZ )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:)
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: ObsZ, ObsXY, ObsXYZ

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj
        Complex (Kind=Kind(0.d0)) :: ZXY, ZZ

        If ( Size(List,1) .ne. Size(GR,1) .or. Size(List,2) .ne. 2 )   then
           Write(error_unit,*) 'List in Predefined_Obs_eq_SpinMz_measure  has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        
        If ( ObsZ%File_Latt .ne. "SpinZ" .and. ObsXY%File_Latt .ne. "SpinXY" .and.  &
           & ObsXYZ%File_Latt .ne. "SpinT"  )   then
           Write(error_unit,*) 'Predefined_Obs_eq_SpinMz_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        ! Count and average sign
        ObsZ%N          = ObsZ%N + 1
        ObsZ%Ave_sign   = ObsZ%Ave_sign + real(ZS,kind(0.d0))
        ObsXY%N         = ObsXY%N + 1
        ObsXY%Ave_sign  = ObsXY%Ave_sign + real(ZS,kind(0.d0))
        ObsXYZ%N        = ObsXYZ%N + 1
        ObsXYZ%Ave_sign = ObsXYZ%Ave_sign + real(ZS,kind(0.d0))


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
                       
                       ObsZ  %Obs_Latt(imj,1,no_I,no_J) =  ObsZ  %Obs_Latt(imj,1,no_I,no_J) +  ZZ  *ZP*ZS
                       ObsXY %Obs_Latt(imj,1,no_I,no_J) =  ObsXY %Obs_Latt(imj,1,no_I,no_J) +  ZXY *ZP*ZS
                       ObsXYZ%Obs_Latt(imj,1,no_I,no_J) =  ObsXYZ%Obs_Latt(imj,1,no_I,no_J) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                    endif
                 enddo
                 ObsZ  %Obs_Latt0(no_I) =  ObsZ  %Obs_Latt0(no_I) +   (GRC(I1,I1,2) - GRC(I1,I1,1)) * ZP*ZS
                 ObsXYZ%Obs_Latt0(no_I) =  ObsXYZ%Obs_Latt0(no_I) +   (GRC(I1,I1,2) - GRC(I1,I1,1)) * ZP*ZS/sqrt(3.d0)
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
      Subroutine Predefined_Obs_eq_Green_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:)
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: Z

        If ( Size(List,1) .ne. Size(GR,1) .or. Size(List,2) .ne. 2 )   then
           Write(error_unit,*) 'List in  Predefined_Obs_eq_Green_measure  has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        If ( Obs%File_Latt .ne. "Green" )   then
           Write(error_unit,*) 'Predefined_Obs_eq_Green_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        ! Count and average sign
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))

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
                    Obs%Obs_Latt(imj,1,no_I,no_J) =  Obs%Obs_Latt(imj,1,no_I,no_J) + Z*ZP*ZS
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
      Subroutine Predefined_Obs_eq_Den_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:)
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: ZI, ZJ, Z


        If ( Size(List,1) .ne. Size(GR,1) .or. Size(List,2) .ne. 2 )   then
           Write(error_unit,*) 'List in  Predefined_Obs_eq_Den_measure  has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        If ( Obs%File_Latt .ne. "Den" )   then
           Write(error_unit,*) 'Predefined_Obs_eq_Den_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        ! Count and average sign
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))

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
                    Obs%Obs_Latt(imj,1,no_I,no_J) =  Obs%Obs_Latt(imj,1,no_I,no_J) + (ZI*ZJ + Z)*ZP*ZS
                 endif
              enddo
              Obs%Obs_Latt0(no_I) =  Obs%Obs_Latt0(no_I) +  ZI * ZP*ZS
           endif
        enddo
      end Subroutine Predefined_Obs_eq_Den_measure


!-------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure time displaced Green function
!>  \sum_{s=1}^{N_SUN} \sum_{nf=1}^{N_FL} < c^{dag}_{i,s,nf}(tau) c_{j,s,nf} >
!>
!--------------------------------------------------------------------
      Subroutine Predefined_Obs_tau_Green_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:), NT
        Complex (Kind=Kind(0.d0)), Intent(In) :: GT0(:,:,:), G0T(:,:,:), G00(:,:,:), GTT(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: Z


        If ( Size(List,1) .ne. Size(GT0,1) .or. Size(List,2) .ne. 2 )   then
           Write(error_unit,*) 'List in  Predefined_Obs_tau_Green_measure  has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        If ( Obs%File_Latt .ne. "Green" )   then
           Write(error_unit,*) 'Predefined_Obs_tau_Green_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        ! Count and average sign
        If (NT == 0 ) then
           Obs%N        = Obs%N + 1
           Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        endif

        ! Measure
        N_FL = Size(GT0,3)
        Do I1 = 1,Size(List,1) 
           I    = List(I1,1)
           If  ( I > 0 ) then 
              no_I = List(I1,2)
              Do J1 = 1,Size(List,1) 
                 J    = List(J1,1)
                 If  (J > 0 )  then 
                    no_J = List(J1,2)
                    imj  = latt%imj(I,J)
                    Z = cmplx(0.d0,0.d0,Kind(0.d0))
                    Do nf = 1, N_FL
                       Z = Z +  GT0(I1,J1,nf)
                    enddo
                    Z   = Z * cmplx(1.d0/dble(N_FL), 0.d0, kind(0.D0))
                    Obs%Obs_Latt(imj,NT+1,no_I,no_J) =  Obs%Obs_Latt(imj,NT+1,no_I,no_J) + Z*ZP*ZS
                 endif
              enddo
           endif
        enddo

      end Subroutine Predefined_Obs_tau_Green_measure

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
      Subroutine Predefined_Obs_tau_SpinSUN_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:), NT
        Complex (Kind=Kind(0.d0)), Intent(In) :: GT0(:,:,:), G0T(:,:,:), G00(:,:,:), GTT(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: Z

        If ( Obs%File_Latt .ne. "SpinZ" )   then
           Write(error_unit,*) 'Predefined_Obs_tau_SpinSUN_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        If ( Size(List,1) .ne. Size(GT0,1) .or. Size(List,2) .ne. 2  )   then
           Write(error_unit,*) 'List in  Predefined_Obs_tau_SpinSUN_measure  has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        ! Count and average sign
        If (NT == 0 ) then
           Obs%N        = Obs%N + 1
           Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        endif

        ! Measure
        N_FL = Size(GT0,3)
        Z =  cmplx(dble(N_SUN),0.d0, kind(0.D0))
        Do I1 = 1,Size(List,1)
           I   = List(I1,1)
           If (I > 0 ) then
              no_I = List(I1,2)
              Do J1 = 1,Size(List,1)
                 J    = List(J1,1)
                 If  (J > 0)  then 
                    no_J = List(J1,2)
                    imj  = latt%imj(I,J)
                    Obs%Obs_Latt(imj,NT+1,no_I,no_J) =  Obs%Obs_Latt(imj,NT+1,no_I,no_J) &
                         &  -  Z*G0T(J1,I1,1) * GT0(I1,J1,1) *ZP*ZS
                 endif
              enddo
           endif
        enddo

      end Subroutine Predefined_Obs_tau_SpinSUN_measure


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
      Subroutine Predefined_Obs_tau_SpinMz_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, ObsZ, ObsXY, ObsXYZ )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:), NT
        Complex (Kind=Kind(0.d0)), Intent(In) :: GT0(:,:,:), G0T(:,:,:), G00(:,:,:), GTT(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: ObsZ, ObsXY, ObsXYZ

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: ZZ, ZXY

        If ( Size(List,1) .ne. Size(GT0,1) .or. Size(List,2) .ne. 2 )   then
           Write(error_unit,*) 'List in Predefined_Obs_tau_SpinMz_measure has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif


        ! Count and average sign
        If ( ObsZ%File_Latt .ne. "SpinZ" .and. ObsXY%File_Latt .ne. "SpinXY" .and.  &
           & ObsXYZ%File_Latt .ne. "SpinT"  )   then
           Write(error_unit,*) 'Predefined_Obs_tau_SpinMz_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        If (NT == 0 ) then
           ! Count and average sign
           ObsZ%N          = ObsZ%N + 1
           ObsZ%Ave_sign   = ObsZ%Ave_sign + real(ZS,kind(0.d0))
           ObsXY%N         = ObsXY%N + 1
           ObsXY%Ave_sign  = ObsXY%Ave_sign + real(ZS,kind(0.d0))
           ObsXYZ%N        = ObsXYZ%N + 1
           ObsXYZ%Ave_sign = ObsXYZ%Ave_sign + real(ZS,kind(0.d0))
        endif

        ! Measure
        N_FL = Size(GT0,3)
        Do I1 = 1,Size(List,1)
           I    = List(I1,1)
           If  (I > 0 ) then
              no_I = List(I1,2)
              Do J1 = 1,Size(List,1)
                 J    = List(J1,1)
                 if (J > 0 ) then 
                    no_J = List(J1,2)
                    imj  = latt%imj(I,J)
                    ZZ   =      (GTT(I1,I1,1) -  GTT(I1,I1,2) ) * ( G00(J1,J1,1)  -  G00(J1,J1,2) )   &
                         &    -  G0T(J1,I1,1) * GT0(I1,J1,1)  -  G0T(J1,I1,2) * GT0(I1,J1,2) 
                    ZXY  =    -  G0T(J1,I1,1) * GT0(I1,J1,2)  -  G0T(J1,I1,2) * GT0(I1,J1,1) 
                    ObsZ  %Obs_Latt(imj,NT+1,no_I,no_J) =  ObsZ %Obs_Latt(imj,NT+1,no_I,no_J)  +  ZZ*ZP*ZS
                    ObsXY %Obs_Latt(imj,NT+1,no_I,no_J) =  ObsXY%Obs_Latt(imj,NT+1,no_I,no_J)  +  ZXY*ZP*ZS
                    ObsXYZ%Obs_Latt(imj,NT+1,no_I,no_J) =  ObsXYZ%Obs_Latt(imj,NT+1,no_I,no_J) + (2.d0*ZXY + ZZ)*ZP*ZS/3.d0
                 endif
              enddo
              ObsZ  %Obs_Latt0(no_I) =  ObsZ  %Obs_Latt0(no_I) +   (GTT(I1,I1,2) - GTT(I1,I1,1)) * ZP*ZS
              ObsXYZ%Obs_Latt0(no_I) =  ObsXYZ%Obs_Latt0(no_I) +   (GTT(I1,I1,2) - GTT(I1,I1,1)) * ZP*ZS/sqrt(3.d0)
           endif
        enddo

      end Subroutine Predefined_Obs_tau_SpinMz_measure


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
      Subroutine Predefined_Obs_tau_Den_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs )

        Type (Lattice),       Intent(in)      :: Latt
        Type (Unit_cell),     Intent(in)      :: Latt_unit
        Integer,              Intent(In)      :: N_SUN, LIST(:,:), NT
        Complex (Kind=Kind(0.d0)), Intent(In) :: GT0(:,:,:), G0T(:,:,:), G00(:,:,:), GTT(:,:,:), ZS, ZP
        Type (Obser_Latt),    Intent(inout)   :: Obs

        ! Local
        Integer :: N_FL, I, I1, J, J1, no_I, no_J, imj,nf
        Complex (Kind=Kind(0.d0)) :: Z, ZI, ZJ
        
        If ( Size(List,1) .ne. Size(GT0,1) .or. Size(List,2) .ne. 2 )   then
           Write(error_unit,*) 'List in Predefined_Obs_tau_Den_measure has  wrong  dim.'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        
        If ( Obs%File_Latt .ne. "Den" )   then
           Write(error_unit,*) 'Predefined_Obs_tau_Den_measure: Wrong filename'
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif

        ! Count and average sign
        If (NT == 0 ) then
           Obs%N        = Obs%N + 1
           Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        endif
        N_FL = Size(GT0,3)
        !Measure
        Do I1 = 1,Size(List,1)
           I    = List(I1,1)
           If ( I > 0 )  then 
              no_I = List(I1,2)
              ZI = cmplx(0.d0,0.d0,Kind(0.d0))
              Do nf = 1,N_FL
                 ZI  = ZI  + cmplx(1.d0,0.d0,kind(0.d0)) -  GTT(I1,I1,nf)
              enddo
              ZI  = ZI * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
              Do J1 = 1,Size(List,1)
                 J    = List(J1,1)
                 If  ( J > 0 ) then
                    no_J = List(J1,2)
                    imj  = latt%imj(I,J)
                    ZJ = cmplx(0.d0,0.d0,Kind(0.d0))
                    Do nf = 1,N_FL
                       ZJ  = ZJ + cmplx(1.d0,0.d0,kind(0.d0)) -  G00(J1,J1,nf)
                    enddo
                    ZJ  = ZJ * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                    Z   = cmplx(0.d0,0.d0,Kind(0.d0))
                    Do nf = 1,N_FL
                       Z = Z - G0T(J1,I1,nf) * GT0(I1,J1,nf)
                    Enddo
                    Z   =  Z * cmplx(dble(N_SUN), 0.d0, kind(0.D0))
                    Obs%Obs_Latt(imj,NT+1,no_I,no_J) =  Obs%Obs_Latt(imj,NT+1,no_I,no_J) + (ZI*ZJ + Z)*ZP*ZS
                 endif
              enddo
              Obs%Obs_Latt0(no_I) =  Obs%Obs_Latt0(no_I) +  ZI * ZP*ZS
           endif
        enddo

      end Subroutine Predefined_Obs_tau_Den_measure
      
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
      Complex  (Kind=Kind(0.d0))  function  Predefined_Obs_V_Int(OP_Vint, GR, GRC, N_SUN )

        Implicit none
        type (Operator)          , Intent(In) :: OP_Vint(:) 
        Integer                  , Intent(In) :: N_SUN
        Complex (Kind=Kind(0.d0)), Intent(In) :: GR(:,:,:), GRC(:,:,:)

        ! Local
        Complex  (Kind=Kind(0.d0))  ::  Z,  Z_tmp,  ZC  
        Integer  ::  N_FL, N,  I, J,  I1,  J1, nf
        Real     (Kind=Kind(0.d0))  ::   Zero=1.0D-16

        If ( OP_Vint(1)%type  .ne. 2 )   then
           Write(error_unit,*) 'Predefined_Obs_V_Int  routine is   defined  fro  tpye  2  vertices.  '
           CALL Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
        endif
        
        N_FL  =   Size(GR,3)
        N     =   Size(OP_Vint(1)%O,1)

        Z    =  cmplx(0.d0,0.d0,kind(0.d0))
        do  nf  =  1,N_FL 
           Z_tmp     =   cmplx(0.d0,0.d0,kind(0.d0))
           do  J =  1,N 
              do  I  =  1,N
                 Z_tmp  =  Z_tmp  +   OP_Vint(nf)%O(i,j)  * GRC(OP_Vint(nf)%P(I), OP_Vint(nf)%P(J), nf)
              enddo
           enddo
           Z =   Z   +  Z_tmp   +   OP_Vint(nf)%alpha
        enddo
        Z  =  Z*dble(N_SUN)

        ZC  =  cmplx(0.d0,0.d0,kind(0.d0))
        Do nf = 1, N_FL
           Do  J =  1,N
              Do I =  1,N
                 Z_tmp  =   OP_Vint(nf)%O(I,J)
                 If  (  real(Z_tmp*conjg(Z_tmp),  kind(0.d0)) >  Zero )  then
                    Do J1 = 1,N
                       Do I1  = 1,N
                          ZC  =  ZC  +  Z_tmp*OP_Vint(nf)%O(I1,J1)* &
                               &        GRC(OP_Vint(nf)%P(I),OP_Vint(nf)%P(J1),nf)*&
                               &         GR(OP_Vint(nf)%P(J),OP_Vint(nf)%P(I1),nf)
                       Enddo
                    Enddo
                 Endif
              Enddo
           Enddo
        Enddo
        ZC =  ZC * dble(N_SUN)
        Predefined_Obs_V_Int =   Z*Z   + ZC 
        
      end function Predefined_Obs_V_Int

      
      
      Subroutine Predefined_Obs_scal_Renyi_Ent_indep(GRC, List, Nsites, N_SUN, ZS, ZP, Obs )

        Implicit none
        Integer, Dimension(:), INTENT(IN)     :: List
        Integer, INTENT(IN)                   :: Nsites, N_SUN
        Complex (kind=kind(0.d0)), INTENT(IN) :: GRC(:,:,:), ZS, ZP
        Type (Obser_Vec),    Intent(inout)   :: Obs
        
        Complex (kind=kind(0.d0))  :: Renyi

        Renyi = Calc_Renyi_Ent_indep(GRC, List, Nsites, N_SUN)
        
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        Obs%Obs_vec(1)   = Obs%Obs_vec(1) + Renyi* ZP*ZS
        
      end Subroutine Predefined_Obs_scal_Renyi_Ent_indep
      
      Subroutine Predefined_Obs_scal_Renyi_Ent_gen_fl(GRC, List, Nsites, N_SUN, ZS, ZP, Obs )

        Implicit none
        Integer, Dimension(:,:), INTENT(IN)   :: List
        Integer, Dimension(:), INTENT(IN)     :: Nsites, N_SUN
        Complex (kind=kind(0.d0)), INTENT(IN) :: GRC(:,:,:), ZS, ZP
        Type (Obser_Vec),    Intent(inout)   :: Obs
        
        Complex (kind=kind(0.d0))  :: Renyi

        Renyi = Calc_Renyi_Ent_gen_fl(GRC, List, Nsites, N_SUN)
        
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        Obs%Obs_vec(1)   = Obs%Obs_vec(1) + Renyi* ZP*ZS
        
      end Subroutine Predefined_Obs_scal_Renyi_Ent_gen_fl
      
      Subroutine Predefined_Obs_scal_Renyi_Ent_gen_all(GRC, List, Nsites, ZS, ZP, Obs )

        Implicit none
        Integer, Dimension(:,:,:), INTENT(IN)   :: List
        Integer, Dimension(:,:), INTENT(IN)     :: Nsites
        Complex (kind=kind(0.d0)), INTENT(IN) :: GRC(:,:,:), ZS, ZP
        Type (Obser_Vec),    Intent(inout)   :: Obs
        
        Complex (kind=kind(0.d0))  :: Renyi

        Renyi = Calc_Renyi_Ent_gen_all(GRC, List, Nsites)
        
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        Obs%Obs_vec(1)   = Obs%Obs_vec(1) + Renyi* ZP*ZS
        
      end Subroutine Predefined_Obs_scal_Renyi_Ent_gen_all
      
      Subroutine Predefined_Obs_scal_Mutual_Inf_indep(GRC, List_A, Nsites_A, List_B, Nsites_B, N_SUN, ZS, ZP, Obs )

        Implicit none
        Integer, Dimension(:), INTENT(IN)     :: List_A, List_B
        Integer, INTENT(IN)                   :: Nsites_A ,Nsites_B, N_SUN
        Complex (kind=kind(0.d0)), INTENT(IN) :: GRC(:,:,:), ZS, ZP
        Type (Obser_Vec),    Intent(inout)   :: Obs
        
        Complex (kind=kind(0.d0))  :: Renyi_A, Renyi_B, Renyi_AB

        call Calc_Mutual_Inf_indep(GRC, List_A, Nsites_A, List_B, Nsites_B, N_SUN, Renyi_A, Renyi_B, Renyi_AB)
        
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        Obs%Obs_vec(1)   = Obs%Obs_vec(1) + Renyi_A* ZP*ZS
        Obs%Obs_vec(2)   = Obs%Obs_vec(2) + Renyi_B* ZP*ZS
        Obs%Obs_vec(3)   = Obs%Obs_vec(3) + Renyi_AB* ZP*ZS
        
      end Subroutine Predefined_Obs_scal_Mutual_Inf_indep
      
      Subroutine Predefined_Obs_scal_Mutual_Inf_gen_fl(GRC, List_A, Nsites_A, List_B, Nsites_B, N_SUN, ZS, ZP, Obs )

        Implicit none
        Integer, Dimension(:,:), INTENT(IN)   :: List_A, List_B
        Integer, Dimension(:), INTENT(IN)     :: Nsites_A ,Nsites_B, N_SUN
        Complex (kind=kind(0.d0)), INTENT(IN) :: GRC(:,:,:), ZS, ZP
        Type (Obser_Vec),    Intent(inout)   :: Obs
        
        Complex (kind=kind(0.d0))  :: Renyi_A, Renyi_B, Renyi_AB

        call Calc_Mutual_Inf_gen_fl(GRC, List_A, Nsites_A, List_B, Nsites_B, N_SUN, Renyi_A, Renyi_B, Renyi_AB)
        
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        Obs%Obs_vec(1)   = Obs%Obs_vec(1) + Renyi_A* ZP*ZS
        Obs%Obs_vec(2)   = Obs%Obs_vec(2) + Renyi_B* ZP*ZS
        Obs%Obs_vec(3)   = Obs%Obs_vec(3) + Renyi_AB* ZP*ZS
        
      end Subroutine Predefined_Obs_scal_Mutual_Inf_gen_fl
      
      Subroutine Predefined_Obs_scal_Mutual_Inf_gen_all(GRC, List_A, Nsites_A, List_B, Nsites_B, ZS, ZP, Obs )

        Implicit none
        Integer, Dimension(:,:,:), INTENT(IN)   :: List_A, List_B
        Integer, Dimension(:,:), INTENT(IN)     :: Nsites_A ,Nsites_B
        Complex (kind=kind(0.d0)), INTENT(IN) :: GRC(:,:,:), ZS, ZP
        Type (Obser_Vec),    Intent(inout)   :: Obs
        
        Complex (kind=kind(0.d0))  :: Renyi_A, Renyi_B, Renyi_AB

        call Calc_Mutual_Inf_gen_all(GRC, List_A, Nsites_A, List_B, Nsites_B, Renyi_A, Renyi_B, Renyi_AB)
        
        Obs%N        = Obs%N + 1
        Obs%Ave_sign = Obs%Ave_sign + real(ZS,kind(0.d0))
        Obs%Obs_vec(1)   = Obs%Obs_vec(1) + Renyi_A* ZP*ZS
        Obs%Obs_vec(2)   = Obs%Obs_vec(2) + Renyi_B* ZP*ZS
        Obs%Obs_vec(3)   = Obs%Obs_vec(3) + Renyi_AB* ZP*ZS
        
      end Subroutine Predefined_Obs_scal_Mutual_Inf_gen_all


      
    End Module Predefined_Obs
