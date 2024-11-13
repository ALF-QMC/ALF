       Module Histograms_v2

         Use Log_Mesh
         use iso_fortran_env, only: output_unit, error_unit

         Type Histogram 
            Type (logmesh), allocatable :: mesh(:)
            Real (Kind=Kind(0.d0)), pointer :: el(:,:)
            Real (Kind=Kind(0.d0))  :: range_st, range_en, dis
            Real (Kind=Kind(0.d0))  :: count
            Character (16) :: File
            integer :: dim
            
         end Type Histogram
       
         Interface   Make_Hist
            module procedure Construct_Hist_1D, Construct_Hist_2D
         end Interface Make_Hist
         Interface Clear_Hist
            module procedure Destroy_Hist
         end Interface Clear_Hist
         
         contains

!--------------------------------------------------------------------
           subroutine Construct_Hist_1D(Hist, file, range, center, dis, Type, Lambda)
             Implicit none
             type (Histogram)  :: Hist
             Real  (Kind=Kind(0.d0)), intent(in)   :: Range, Center, dis, Lambda
             Character (16)    :: File
             Character(len=10) :: Type
             Integer :: n, Nw_1
             
             Hist%dim      = 1
             allocate (hist%mesh(Hist%dim))

             Nw_1   = NINT(range*2.d0/dis)
             
             call Make_log_mesh(Hist%Mesh(1),  Lambda, Center, Range, Type, Nw_1)
             write(6,*) 'In Construct_hist: ',  Size(Hist%Mesh(1)%Xom,1)
             n =  Size(Hist%Mesh(1)%Xom,1)
             allocate ( Hist%el(n,1) )
             Hist%el = 0.d0
             Hist%file     = file
             Hist%count    = 0.d0
             
           end subroutine Construct_Hist_1D


!--------------------------------------------------------------------
           subroutine Construct_Hist_2D(Hist, file, range, center, dis, Type, Lambda)
             Implicit none
             type (Histogram) :: Hist
             Real  (Kind=Kind(0.d0)), dimension(:), Intent(in)  :: Range, Center, dis, Lambda
             Character (16)     :: File
             Character(len=10)  :: Type

             Integer :: n(2), Nw_1, nd
             
             Hist%dim      = 2
             allocate (hist%mesh(Hist%dim))

             do nd = 1, hist%dim
               Nw_1   = NINT(range(nd)*2.d0/dis(nd))
             
               call Make_log_mesh(Hist%Mesh(nd),  Lambda(nd), Center(nd), Range(nd), Type, Nw_1)
               write(6,*) 'In Construct_hist: ',  Size(Hist%Mesh(nd)%Xom,1)
               n(nd) =  Size(Hist%Mesh(nd)%Xom,1)
             enddo

             allocate ( Hist%el(n(1),n(2)) )
             Hist%el = 0.d0
             Hist%file     = file
             Hist%count    = 0.d0

             
           end subroutine Construct_Hist_2D

!--------------------------------------------------------------------
           subroutine Destroy_Hist(Hist)
             Implicit none
             type (Histogram) :: Hist

             integer :: n
             
             deallocate ( Hist%el )
!             Hist%el = 0.d0
             Hist%file     = ""
             Hist%count    = 0.d0
             do n = 1, size(hist%mesh,1)
                Call Clear_log_mesh ( Hist%Mesh(n) )
             enddo
             
           end subroutine Destroy_Hist

!--------------------------------------------------------------------
!!$           subroutine Read_Hist(Hist)
!!$             Implicit none
!!$             type (Histogram) :: Hist
!!$
!!$             integer :: io_error, nv
!!$             Real (Kind=Kind(0.d0)) :: X,Y
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


!--------------------------------------------------------------------
           subroutine Write_Hist(Hist,Group_Comm)

             Implicit none

             type (Histogram), INTENT(IN) :: Hist
             Integer,          INTENT(IN) :: Group_Comm

             ! LOCAL
             CHARACTER (LEN=64) :: file_out
             Integer :: nv1, nv2, i, no
             real  (Kind=Kind(0.d0)), allocatable :: Tmp(:,:)
             Real     (Kind=Kind(0.d0)) :: X

#if defined(MPI)
             INTEGER        :: irank_g, isize_g, igroup, ISIZE, IRANK, IERR
             CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
             CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
             call MPI_Comm_rank(Group_Comm, irank_g, ierr)
             call MPI_Comm_size(Group_Comm, isize_g, ierr)
             igroup           = irank/isize_g
       
             No = size(Hist%el,1)*size(Hist%el,2)
             Allocate (Tmp(size(Hist%el,1),size(Hist%el,2)) )
             Tmp = 0.d0
             CALL MPI_REDUCE(Hist%el,Tmp,No,MPI_Real8,MPI_SUM, 0,Group_Comm,IERR)

             I = 1
             X = 0.d0
             CALL MPI_REDUCE(Hist%count,X,I,MPI_REAL8,MPI_SUM, 0,Group_comm,IERR)

             if (Irank_g == 0 ) then
#else
             Allocate (Tmp(size(Hist%el,1),size(Hist%el,2)) )
             Tmp = Hist%el
             X = Hist%count
#endif

#if defined(TEMPERING)
             write(file_out,'(A,I0,A,A,A)') "Temp_",igroup,"/",trim(Hist%file), "_hist"
#else
             write(file_out,'(A,A)') trim(Hist%file), "_hist"
#endif

             Open ( unit=20,file=file_out,status='unknown')
             write(20,*) x
             if (hist%dim == 1) then
               do nv1 = 1,size(Hist%el,1) -1
                   write(20,*) Hist%Mesh(1)%Xom(nv1), tmp(nv1,1)/(X * Hist%Mesh(1)%DXom(nv1))
               enddo
             else
               do nv1 = 1, size(Hist%el,1) - 1
                 do nv2 = 1, size(Hist%el,2) - 1
                   write(20,*) Hist%Mesh(1)%Xom(nv1), Hist%Mesh(2)%Xom(nv2), tmp(nv1,nv2)/(X * Hist%Mesh(1)%DXom(nv1)*Hist%Mesh(2)%DXom(nv2))
                 enddo
                 write(20,*) ""
               enddo
             endif
             close(20)

#if defined(MPI)
             endif
#endif

           end subroutine Write_Hist

!--------------------------------------------------------------------
           Real (Kind=Kind(0.d0)) function Inter_Hist(Hist)
             Implicit none
             type (Histogram), intent(in) :: Hist

             Integer :: nv1, nv2
             Real (Kind=Kind(0.d0)) :: X
             
             X = 0.d0
             
             if (hist%dim == 1) then
               do nv1 = 1,size(Hist%el,1) -1
                    X = X + Hist%el(nv1,1)
               enddo
             else
               do nv1 = 1,size(Hist%el,1) -1
                    do nv2 = 1, size(Hist%el,2) -1
                        X = X + Hist%el(nv1,nv2)
                    enddo
               enddo
             endif
             Inter_Hist = X
           end function Inter_Hist
!--------------------------------------------------------------------

           subroutine Add_Hist(Hist,value1,value2)
             Implicit none
             type (Histogram), intent(inout) :: Hist
             Real (Kind=Kind(0.d0)), intent(in)    :: value1
             real (kind=kind(0.d0)), intent(in), optional :: value2

             Integer :: nv1, nv2
             
             if (hist%dim == 1 .and. present(value2) ) then
                WRITE(error_unit,*) 'Error in adding value to histogram: one value expected, but two present'
                Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
             elseif (hist%dim == 2 .and. (.not. present(value2)) ) then
                WRITE(error_unit,*) 'Error in adding value to histogram: two values expected, but one present'
                Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
             endif
             nv1 = m_find(Value1,Hist%Mesh(1))
             nv2 = 1
             if (Hist%dim == 2) nv2 = m_find(Value2,Hist%Mesh(2))
             Hist%el(nv1,nv2) = Hist%el(nv1,nv2) + 1.0
             Hist%count = Hist%count + 1.0
           end subroutine Add_Hist
           

         end Module Histograms_v2
