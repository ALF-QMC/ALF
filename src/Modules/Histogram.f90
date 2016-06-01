       Module Histograms

         Type Histogram 
            Real (Kind=8), pointer :: el(:)
            Real (Kind=8)  :: range_st, range_en, dis
            Real (Kind=8)  :: count
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
             Real  (Kind=8)   :: range_st, range_en, dis
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
             Real (Kind=8) :: X,Y
             

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
             Real (Kind=8)    :: value
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
