!  Copyright (C) 2018-2026 The ALF project
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
!> @author ALF-project  
!> @brief Simple 1D uniform histogram for binning sampled data.
!> @details Provides histogram type with uniform bin spacing.
!> Operations: Make_Hist (allocate), Clear_Hist (deallocate),
!> Add_Hist (bin increment), Read_Hist, Write_Hist.
!--------------------------------------------------------------------
       Module Histograms

         Type Histogram 
            Real (Kind=Kind(0.d0)), pointer :: el(:)
            Real (Kind=Kind(0.d0))  :: range_st, range_en, dis
            Real (Kind=Kind(0.d0))  :: count
            Character (16) :: File
            
         end Type Histogram
       
         Interface   Make_Hist
            module procedure Construct_Hist
         end Interface Make_Hist
         Interface Clear_Hist
            module procedure Destroy_Hist
         end Interface Clear_Hist
         
         contains

           subroutine Construct_Hist(Hist, file, range_st, range_en, dis)
             Implicit none
             type (Histogram) :: Hist
             Real  (Kind=Kind(0.d0))   :: range_st, range_en, dis
             Character (16)   :: File
             
             Integer :: n
             n = nint( ( range_en -  range_st)/dis )
             allocate ( Hist%el(n) )
             Hist%el = 0.d0
             Hist%range_st = range_st
             Hist%range_en = range_en
             Hist%dis      = dis
             Hist%file     = file
             Hist%count    = 0.d0
             
           end subroutine Construct_Hist

           subroutine Destroy_Hist(Hist)
             Implicit none
             type (Histogram) :: Hist
             
             deallocate ( Hist%el )
             Hist%el = 0.d0
             Hist%range_st = 0.d0
             Hist%range_en = 0.d0
             Hist%dis      = 0.d0
             Hist%file     = ""
             Hist%count    = 0.d0
             
           end subroutine Destroy_Hist


           subroutine Read_Hist(Hist)
             Implicit none
             type (Histogram) :: Hist
             
             integer :: io_error, nv
             Real (Kind=Kind(0.d0)) :: X,Y
             

             Open ( unit=20,file=Hist%file,status='old',action='read', iostat=io_error)
             If (io_error.eq.0) then 
                read(20,*) Hist%count
                do nv = 1,size(Hist%el,1)
                   read(20,*) X, Y
                   Hist%el(nv) = Y * Hist%count * Hist%dis
                enddo
             else
                Hist%count = 0.d0
                Hist%el =  0.d0
             endif
             close(20)
           end subroutine Read_Hist


           subroutine Write_Hist(Hist)
             Implicit none
             type (Histogram) :: Hist
             Integer :: nv
             
             Open ( unit=20,file=Hist%file,status='unknown')
             write(20,*) Hist%count
             do nv = 1,size(Hist%el,1)
                write(20,*) dble(nv)*Hist%dis + Hist%range_st, Hist%el(nv)/(Hist%count * Hist%dis)
             enddo
             close(20)

           end subroutine Write_Hist


           subroutine Add_Hist(Hist,value)
             Implicit none
             type (Histogram) :: Hist
             Real (Kind=Kind(0.d0))    :: value
             Integer :: nv

             if ( value .gt.  Hist%range_en .or. value .lt.  Hist%range_st ) then
                write(6,*) 'Error in Add_Hist: ', Hist%file, value
             else
                nv = int((value   - Hist%range_st )/Hist%dis)
                if (nv < 1) nv =1
                if (nv > size(Hist%el,1) ) nv = size(Hist%el,1)
                Hist%el(nv) = Hist%el(nv) + 1.0
                Hist%count = Hist%count + 1.0
             endif
           end subroutine Add_Hist
           

         end Module Histograms
