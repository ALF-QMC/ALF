
      Module Log_Mesh
        use iso_fortran_env, only: output_unit, error_unit
        use runtime_error_mod
        use files_mod

        Type logmesh
           Real (Kind=Kind(0.d0))  :: Lambda, Center, Log_Lambda
           Real (Kind=Kind(0.d0))  :: Range
           Real (Kind=Kind(0.d0))  :: Om_st, Om_en, dom
           Real (Kind=Kind(0.d0))  :: Precision
           Integer        :: Nom,Nw
           Real (Kind=Kind(0.d0)), pointer :: Xom(:),DXom(:)
           Character(len=10) :: Type
        end Type logmesh

        Interface Lookup_log_mesh
           module procedure Lookup_log_mesh_R, Lookup_log_mesh_C
        end Interface
        Interface Inter_log_mesh
           module procedure Inter_log_mesh_R, Inter_log_mesh_C
        end Interface

      Contains

        !< Rng The Range

        subroutine Make_log_mesh ( Mesh,  Lambda, Center, Rng, Type, Nw_1 )

          Implicit None

          Type (logmesh)      :: Mesh
          Real (Kind=Kind(0.d0))       :: Lambda, Center, Rng
          Integer, Optional   :: Nw_1
          Integer             :: N, Nw
          Character(len=10)   :: Type

          Real (Kind=Kind(0.d0))  :: Dom, Om_st, Om_en

          Mesh%Center     = Center
          Mesh%Range      = Rng
          If (str_to_upper(Type) == "LOG" ) Then
             OM_st = Center - Rng
             OM_en = Center + Rng
             Mesh%Om_st = Om_st
             Mesh%Om_en = Om_en
             Mesh%Lambda     = Lambda
             Mesh%Type       = "Log"
             if (Present(Nw_1) ) then
                Nw = Nw_1
             else
                Nw = NINT(10.D0*log(10.D0)/log(Lambda))
             endif
             Mesh%Nw         = Nw
             Mesh%Nom        = 2*Nw + 3
             Mesh%Log_Lambda = Log(Lambda)
             Allocate   ( Mesh%Xom(2*Nw + 3), Mesh%DXom(2*Nw+3) )
             Do n = 0,Nw
                Mesh%xom (n+1          ) =  Center  -   Rng * (Lambda**(-n))
             enddo
             Mesh%xom   (Nw+2         ) =  Center
             do n = Nw,0,-1
                Mesh%xom(Nw+3 +(Nw-n) ) =  Center  +   Rng * (Lambda**(-n))
             enddo
             Mesh%Precision = Mesh%Lambda**(-Mesh%Nw)
          elseif (str_to_upper(Type) == "LIN" ) then
             Mesh%Type       = "Lin"
             If ( Present(Nw_1)  ) then
                Nw   = Nw_1
                Mesh%Nw   = Nw
                Mesh%Nom  = 2*Nw + 1
                Mesh%Type = "Lin"
                Allocate   ( Mesh%Xom(2*Nw + 1), Mesh%DXom(2*Nw+1) )
                OM_st = Center - Rng
                OM_en = Center + Rng
                Dom = Rng/dble(Nw_1)
                Mesh%Dom   = Dom
                Mesh%Om_st = Om_st
                Mesh%Om_en = Om_en
                do n  = 1,Mesh%Nom
                   Mesh%xom(n)  = Om_st + dble(n-1)*dom
                enddo
             else
                Write(error_unit,*) 'Make_log_mesh: You need to include Nw for the Lin Mesh'
                Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
             endif
          else
             Write(error_unit,*) 'Make_log_mesh: Mesh has no type!!'
             Call Terminate_on_error(ERROR_GENERIC,__FILE__,__LINE__)
          endif
          do n = 1,Mesh%Nom-1
             Mesh%DXom(n) = Mesh%xom (n+1) - Mesh%xom (n )
          enddo

        end subroutine Make_log_mesh

        subroutine Clear_log_mesh ( Mesh )
          Implicit None

          Type (logmesh)      :: Mesh

          deallocate   ( Mesh%Xom,  Mesh%DXom  )

        end subroutine Clear_log_mesh

        Integer  Function m_find(X,Mesh)

          Implicit None

          Type (logmesh) :: Mesh
          Real (Kind=Kind(0.d0))  :: X
          Integer        :: m

          if ( str_to_upper(Mesh%Type)  == "LOG" ) then
             if ( X >  (Mesh%OM_en) .or.  X < (Mesh%Om_st) ) then
                m = 0
             else
                if      ( X < Mesh%Xom(Mesh%Nw+1) ) then
                   m = 2 - Int( log ( (Mesh%Center - X)/Mesh%Range )  / Mesh%Log_Lambda )
                   !Write(6,*) 'Hi 1', X
                elseif  ( X >  Mesh%Xom(Mesh%Nw+3) ) then
                   m = 2*Mesh%Nw + 3 + Int( log ( (X- Mesh%Center) /Mesh%Range )  / Mesh%Log_Lambda )
                   !Write(6,*) 'Hi 2', X, Mesh%Center +  Mesh%Range
                elseif  ( X >  Mesh%Center )      then
                   m = Mesh%Nw+3
                else
                   m = Mesh%Nw+2
                endif
             endif
             m_find = m
          else
             m_find = int((x - Mesh%Om_st)/Mesh%dom) + 2
             if (m_find > Mesh%Nom) m_find=Mesh%Nom
             if (m_find < 2 ) m_find=2
          endif


          !Write(6,*)
          !Write(6,*) 'Point: ', X
          !if (  m > 0 ) then
          !   Write(6,*) 'Your point lies inbetween ', Mesh%Xom(m-1), ' and ', Mesh%Xom(m)
          !else
          !   Write(6,*) 'Out of range '
          !endif

        end Function m_find
!*******
        Real(Kind=Kind(0.d0)) Function  Lookup_log_mesh_R(f, x,Mesh,m_1)

          Implicit None

          Type (logmesh) :: Mesh
          Real (Kind=Kind(0.d0)), dimension(:) :: f
          Real (Kind=Kind(0.d0))  :: X
          Integer      , Optional     :: m_1

          Integer ::  m
          Real (Kind=Kind(0.d0)) :: X1,X2,Y1,Y2,a,b

          m = m_find(X,Mesh)
          if (m == 0 ) then
             Lookup_log_mesh_R = 0.d0
          else
             x1 = Mesh%xom(m-1)
             x2 = Mesh%xom(m  )
             y1 = f(m-1)
             y2 = f(m)
             a = (y1-y2)/(x1-x2)
             b = (x1*y2 - x2*y1)/(x1-x2)
             Lookup_log_mesh_R = a*x + b
          endif

          If  ( Present(m_1) )  m_1 = m

        end Function Lookup_log_mesh_R



!*******
!!$        Complex (Kind=Kind(0.d0)) Function  Lookup_log_mesh_C(f, x,Mesh,m_1)
!!$
!!$          Implicit None
!!$
!!$          Type (logmesh) :: Mesh
!!$          Complex (Kind=Kind(0.d0)), dimension(:) :: f
!!$          Real    (Kind=Kind(0.d0))               :: X
!!$          Integer      , Optional     :: m_1
!!$
!!$
!!$          Integer  ::  n, m
!!$          Complex  (Kind=Kind(0.d0)) :: X1,X2,Y1,Y2,a,b
!!$
!!$          m = m_find(X,Mesh)
!!$          if (m == 0 ) then
!!$             Lookup_log_mesh_C = cmplx(0.d0,0.d0)
!!$          else
!!$             x1 = cmplx( Mesh%xom(m-1),0.d0 )
!!$             x2 = cmplx( Mesh%xom(m  ),0.d0 )
!!$             y1 = f(m-1)
!!$             y2 = f(m  )
!!$             a = (y1-y2)/(x1-x2)
!!$             b = (x1*y2 - x2*y1)/(x1-x2)
!!$             Lookup_log_mesh_C = a*cmplx( x , 0.d0 ) + b
!!$          endif
!!$
!!$          If  ( Present(m_1) )  m_1 = m
!!$
!!$        end Function Lookup_log_mesh_C

        Complex (Kind=Kind(0.d0)) Function  Lookup_log_mesh_C(f, x,Mesh,m_1)

          Implicit None

          Type (logmesh) :: Mesh
          Complex (Kind=Kind(0.d0)), dimension(:) :: f
          Real    (Kind=Kind(0.d0))               :: X
          Integer      , Optional     :: m_1


          Integer  ::  m
          Complex  (Kind=Kind(0.d0)) :: Z1,Z2, Z
          Real     (Kind=Kind(0.d0)) :: x1,x2,t

          m = m_find(X,Mesh)
          if (m == 0 ) then
             Lookup_log_mesh_C = cmplx(0.d0, 0.d0, kind(0.D0))
          else
             x1 =  Mesh%xom(m-1)
             x2 =  Mesh%xom(m  )
             t  = (x1 - X)/(x2-x1)
             Z1 = f(m-1)
             Z2 = f(m  )
             Z  = Z1 + (Z1-Z2)*t
             Lookup_log_mesh_C = Z
          endif

          If  ( Present(m_1) )  m_1 = m

        end Function Lookup_log_mesh_C


!******
        Real (Kind=Kind(0.d0)) Function  Inter_log_mesh_R(f,Mesh)

          Implicit None

          Type     (logmesh) :: Mesh
          Real     (Kind=Kind(0.d0)), dimension(:) :: f
          Real     (Kind=Kind(0.d0))  :: X
          Integer            :: n

          X = 0.d0
          do n = 1,Mesh%Nom-1
             X = X + Mesh%DXom(n) * (f(n+1) + f(n) )
          enddo
          Inter_log_mesh_R = X / 2.d0

        end Function Inter_log_mesh_R

!******
        Complex (Kind=Kind(0.d0)) Function  Inter_log_mesh_C(f,Mesh)

          Implicit None

          Type     (logmesh) :: Mesh
          Complex     (Kind=Kind(0.d0)), dimension(:) :: f
          Complex     (Kind=Kind(0.d0))  :: Z
          Integer               :: n

          Z = cmplx(0.d0, 0.d0, kind(0.D0))
          do n = 1,Mesh%Nom-1
             Z = Z + Mesh%DXom(n) * ( f(n+1) + f(n) )
          enddo
          Inter_log_mesh_C = Z /2.d0

        end Function Inter_log_mesh_C



        subroutine Print_log_mesh(Mesh)

          Implicit None

          Type (logmesh) :: Mesh

          Integer :: n

          If (str_to_upper(Mesh%Type) == "LOG" ) Then
             Open (Unit=10,File="Log_Mesh", status="unknown" )
             Write(10,*) '#    Log Mesh     : '
             Write(10,*) '#    Lambda       : ', Mesh%Lambda
             Write(10,*) '#    Range        : ', Mesh%Range
             Write(10,*) '#    Center       : ', Mesh%Center
             Write(10,*) '#    Nom          : ', Mesh%Nom
             Write(10,*) '#    Precision    : ', Mesh%Lambda**(-Mesh%Nw)
             do n = 1,Mesh%Nom
                write(10,"(F16.8)") Mesh%xom(n)
             enddo
             close(10)
          endif

          If (str_to_upper(Mesh%Type) == "LIN" ) Then
             Open (Unit=10,File="Lin_Mesh", status="unknown" )
             Write(10,*) '#    Lin Mesh     : '
             Write(10,*) '#    Range        : ', Mesh%Range
             Write(10,*) '#    Center       : ', Mesh%Center
             Write(10,*) '#    Nom          : ', Mesh%Nom
             Write(10,*) '#    Dom          : ', Mesh%dom
             do n = 1,Mesh%Nom
                write(10,"(F16.8)") Mesh%xom(n)
             enddo
             close(10)
          endif



        end subroutine Print_log_mesh


      end Module Log_Mesh
