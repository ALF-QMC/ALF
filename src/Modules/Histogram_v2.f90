       Module Histograms_v2

         Use Log_Mesh

         Type Histogram 
            Type (logmesh) :: mesh
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

           subroutine Construct_Hist(Hist, file, range, center, dis, Type, Lambda)
             Implicit none
             type (Histogram)  :: Hist
             Real  (Kind=8)    :: Range, Center, dis, Lambda
             Character (16)    :: File
             Character(len=10) :: Type
             Integer :: n, Nw_1
             
             !Local 
             

             Nw_1   = range*2.d0/dis
             
             call Make_log_mesh(Hist%Mesh,  Lambda, Center, Range, Type, Nw_1)
             write(6,*) 'In Construct_hist: ',  Size(Hist%Mesh%Xom,1)
             n =  Size(Hist%Mesh%Xom,1) 
             allocate ( Hist%el(n) )
             Hist%el = 0.d0
             Hist%file     = file
             Hist%count    = 0.d0
             
           end subroutine Construct_Hist

           subroutine Destroy_Hist(Hist)
             Implicit none
             type (Histogram) :: Hist
             
             deallocate ( Hist%el )
             Hist%el = 0.d0
             Hist%file     = ""
             Hist%count    = 0.d0
             Call Clear_log_mesh ( Hist%Mesh )
             
           end subroutine Destroy_Hist


!!$           subroutine Read_Hist(Hist)
!!$             Implicit none
!!$             type (Histogram) :: Hist
!!$             
!!$             integer :: io_error, nv
!!$             Real (Kind=8) :: X,Y
!!$             
!!$
!!$             Open ( unit=20,file=Hist%file,status='old',action='read', iostat=io_error)
!!$             If (io_error.eq.0) then 
!!$                read(20,*) Hist%count
!!$                do nv = 1,size(Hist%el,1)
!!$                   read(20,*) X, Y
!!$                   Hist%el(nv) = Y * Hist%count * Hist%dis
!!$                enddo
!!$             else
!!$                Hist%count = 0.d0
!!$                Hist%el =  0.d0
!!$             endif
!!$             close(20)
!!$           end subroutine Read_Hist


           subroutine Write_Hist(Hist)
             Implicit none
             type (Histogram) :: Hist
             Integer :: nv
             
             Open ( unit=20,file=Hist%file,status='unknown')
             write(20,*) Hist%count
             do nv = 1,size(Hist%el,1) -1
                write(20,*) Hist%Mesh%Xom(nv), Hist%el(nv)/(Hist%count * Hist%Mesh%DXom(nv))
             enddo
             close(20)

           end subroutine Write_Hist

           Real (Kind=8) function Inter_Hist(Hist)
             Implicit none
             type (Histogram) :: Hist
             Integer :: nv
             Real (Kind=8) :: X
             
             X = 0.d0
             do nv = 1,size(Hist%el,1) -1
                X = X + Hist%el(nv)  !*  Hist%Mesh%DXom(nv)
             enddo
             Inter_Hist = X
           end function Inter_Hist

           

           subroutine Add_Hist(Hist,value)
             Implicit none
             type (Histogram) :: Hist
             Real (Kind=8)    :: value
             Integer :: nv
             
             nv = m_find(Value,Hist%Mesh) 
             Hist%el(nv) = Hist%el(nv) + 1.0
             Hist%count = Hist%count + 1.0
           end subroutine Add_Hist
           

         end Module Histograms_v2
