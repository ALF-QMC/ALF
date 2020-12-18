! MIT License
! 
! Copyright (c) 2018-2020 Florian Goth
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights 
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in 
! all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.

module colorvertex_mod
    use Exponentials_mod
    use vertex_mod
    implicit none
    
    type, extends(Vertex) :: ColorVertex
        integer, allocatable, dimension(:) :: cols !< At initialization empty. Can be used for the colors
        integer, allocatable, dimension(:) :: nbrbycol !< An array that holds for each color the local array position in this Vertex
    contains
        procedure, pass(this) :: colorvertex_init
        generic, public :: init => colorvertex_init
        procedure :: destruct => colorvertex_destruct
        procedure :: any_color_available => colorvertex_any_color_available
        procedure :: erase_edge_color => colorvertex_erase_edge_color
        procedure :: set_edge_color => colorvertex_set_edge_color
        procedure :: get_edge_color => colorvertex_get_edge_color
        procedure :: find_maximal_fan => colorvertex_find_maximal_fan
        procedure :: findfreecolor => colorvertex_findfreecolor
        procedure :: iscolorfree => colorvertex_iscolorfree
    end type ColorVertex

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This type is used to store the information contained in a symmetric
!> matrix in terms of the associated graph. The information about the
!> connectivity is stored in the vert array. The actual elements
!> are stored in the elems array.
!--------------------------------------------------------------------
    type :: GraphData
        type(ColorVertex), allocatable, dimension(:) :: verts !< An array with all vertices. Each vertex knows the index of its neighbours in this array.
        complex(kind=kind(0.D0)), allocatable, dimension(:)  :: elems !< The actual weight of the connection that was stored in the matrix.
        integer :: ndim !< corresponds to the number of columns/rows of the associated matrix. Hence this is the number of vertices.
        integer :: nredges !< The number of edges in the graph.
        integer :: deltag !< the graph degree. It is determined by the vertex that has the most connections.
        integer :: usedcolors !< the number of colors that would be used in this graph decomposition.
    end type GraphData

contains

function colorvertex_iscolorfree(this, col) result(r)
    class(ColorVertex) :: this
    integer, intent(in) :: col
    logical :: r
    integer :: k
    
    r = .true.
    do k = 1, this%degree
    ! some debug check while we are here...
        if ((this%cols(k) == col) .and. (r .eqv. .false.)) then
        write (*,*) "duplicate color found!"
        STOP
        endif
        
        if (this%cols(k) == col) r = .false.
    enddo
end function

function colorvertex_findfreecolor(this)  result(col)
    class(ColorVertex) :: this
    integer :: col
    integer :: i
    logical :: usedcols(size(this%nbrbycol))

!     do i = maxcols, 1, -1
!     if(this%nbrbycol(i) == 0) col = i
!     enddo
    usedcols = .false.
    do i = 1, this%degree
        if(this%cols(i) /= 0) usedcols(this%cols(i)) = .true.
    enddo
    col = 0
    do i = size(this%nbrbycol), 1, -1
        if(usedcols(i) .eqv. .false.) col = i
    enddo
end function colorvertex_findfreecolor

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function constructs a Vizing Fan.
!
!> @param[in] this The vertex where we construct a Vizing fan around
!> @param[in] verts the total array of all vertices
!> @param[in] v0 The incident vertex that still is uncolored.
!> @param[out] f a linear array where we construct the fan
!> @param[out] fannbr An array with the positions of the vertices at verts(fanpos)
!> @param[out] fanedgecol The colors of the edges in the fan
!> @param[out] fanlen The length of the entire fan
!> @result rescol
!--------------------------------------------------------------------
function colorvertex_find_maximal_fan(this, verts, v0, maxcols, f, fannbr, fanedgecol, fanlen) result(rescol)
    class(ColorVertex) :: this
    type(ColorVertex), dimension(:), intent(in) :: verts
    integer, intent(out) :: fanlen
    integer :: rescol 
    integer, dimension(:), intent(out) :: f, fannbr, fanedgecol
    integer, intent(in) :: v0, maxcols
    integer :: col, i, j, ctr, nrfreecols
    logical, allocatable, dimension(:) :: usedcols !< in usedcols we track the colors we have chosen during construction of the fan.
    integer, allocatable, dimension(:) :: freecols ! determine the colors that are available at this
    
    allocate (usedcols(maxcols))
    fanlen = 1
    f = 0
    f(1) = v0
    ctr = 1
    rescol = 0
    usedcols = .false.
    nrfreecols = 0

    do i = 1, maxcols
        if (this%nbrbycol(i) == 0) nrfreecols = nrfreecols + 1
    enddo
    allocate(freecols(nrfreecols))
    j = 1
    do i = 1, maxcols
        if (this%nbrbycol(i) == 0) then
        freecols(j) = i
        j = j + 1
        endif
    enddo
    do while ((ctr < this%degree) .and. (rescol == 0))
        ! determine free color at end of fan that has not already been used in the construction of the fan.
        col = 0
        i = 1
        do while ((i <= maxcols) .and. (col == 0))
            if (usedcols(i) .eqv. .false.) then
                if(verts(f(fanlen))%nbrbycol(i) == 0) then
                col = i
                endif
            endif
            i = i + 1
        enddo
        if (col /= 0) then
            ! col is now a small color that is free at f(fanlen)
            ! determine incident edge that has exactly this color
            fanlen = fanlen + 1
            fannbr(fanlen) = this%nbrbycol(col)
            fanedgecol(fanlen) = col
            f(fanlen) = this%nbrs(this%nbrbycol(col))
            usedcols(col) = .true.
            ! let's see wether we can stop the fan construction since we find matching colors
            rescol = 0
            do i = nrfreecols, 1, -1
                if (verts(f(fanlen))%nbrbycol(freecols(i)) == 0) rescol = freecols(i)
            enddo
        else
            ! force termination of fan construction
#ifndef NDEBUG
            write (*,*) "fan construction stops. No color available that has NOT been already used in the fan."
#endif
            ctr = this%degree
        endif
        ctr = ctr + 1
    enddo
    deallocate(usedcols, freecols)
end function colorvertex_find_maximal_fan

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function takes a graphdata object as e.g. determined with the
!> help of the MvG_decomp function and creates a EulerExp(=product 
!> of checkerboard exponentials) object from it.
!
!> @param gd
!> @result fe
!--------------------------------------------------------------------

function createEulerExponentialfromGraphData(gd) result(fe)
    implicit none
    type(GraphData) :: gd
    type(EulerExp) :: fe
    complex(kind=kind(0.D0)) :: weight
    integer :: k, elempos, mynbr, nbr1, l, i
    logical, allocatable, dimension(:) :: usedcols
    type(node), allocatable, dimension(:) :: nodes
    
    if ((gd%usedcolors == 0) .or. (gd%nredges == 0)) then ! check that those are available
        gd%usedcolors = 0
        gd%nredges = 0
        do i = 1, gd%ndim
            gd%deltag = max(gd%deltag, gd%verts(i)%degree)
            do k = 1, gd%verts(i)%degree
                if (gd%verts(i)%nbrs(k) > i) gd%nredges = gd%nredges + 1
                if (gd%verts(i)%nbrs(k) > gd%ndim) then
                    write(*,*) "invalid nbr!!!"
                    STOP
                endif
                gd%usedcolors = max(gd%usedcolors, gd%verts(i)%cols(k))
            enddo
        enddo
    endif

    ! set up data in an edges based layout
    k = 0
    elempos = 0
    allocate( nodes(gd%nredges), usedcols(gd%usedcolors))
    do i = 1, gd%ndim-1
        ! check validity of the coloring locally
        usedcols = .false.
        do l = 1, gd%verts(i)%degree
            if(gd%verts(i)%cols(l) == 0) then
                write (*,*) "forgotten edge found!"
                STOP
            endif
            if (usedcols(gd%verts(i)%cols(l)) .eqv. .true. ) then
                write (*,*) "invalid coloring!!"
                STOP
            else
                usedcols(gd%verts(i)%cols(l)) = .true.
            endif
        enddo
        do l = 1, gd%usedcolors
            mynbr = gd%verts(i)%nbrbycol(l)
            nbr1 = gd%verts(i)%nbrs(mynbr)
            if (nbr1 > i) then ! nbr1 could be zero if there is no such edge
                k = k+1
                nodes(k)%x = i
                nodes(k)%y = nbr1
                nodes(k)%axy = gd%elems(elempos + mynbr)
                nodes(k)%col = l
            endif
        enddo
        elempos = elempos + gd%verts(i)%degree
    enddo
    weight = 1.0
    call fe%init(nodes, gd%usedcolors, weight)
    deallocate(nodes, usedcols)
end function

function createFullExponentialfromGraphData(gd, method) result(fe)
    implicit none
    type(GraphData) :: gd
    integer, intent(in) :: method
    type(FullExp) :: fe
    complex(kind=kind(0.D0)) :: weight
    integer :: k, elempos, mynbr, nbr1, l, i
    logical, allocatable, dimension(:) :: usedcols
    type(node), allocatable, dimension(:) :: nodes
    
    if ((gd%usedcolors == 0) .or. (gd%nredges == 0)) then ! check that those are available
        gd%usedcolors = 0
        gd%nredges = 0
        do i = 1, gd%ndim
            gd%deltag = max(gd%deltag, gd%verts(i)%degree)
            do k = 1, gd%verts(i)%degree
                if (gd%verts(i)%nbrs(k) > i) gd%nredges = gd%nredges + 1
                if (gd%verts(i)%nbrs(k) > gd%ndim) then
                    write(*,*) "invalid nbr!!!"
                    STOP
                endif
                gd%usedcolors = max(gd%usedcolors, gd%verts(i)%cols(k))
            enddo
        enddo
    endif

    ! set up data in an edges based layout
    k = 0
    elempos = 0
    allocate( nodes(gd%nredges), usedcols(gd%usedcolors))
    do i = 1, gd%ndim-1
        ! check validity of the coloring locally
        usedcols = .false.
        do l = 1, gd%verts(i)%degree
            if(gd%verts(i)%cols(l) == 0) then
                write (*,*) "forgotten edge found!"
                STOP
            endif
            if (usedcols(gd%verts(i)%cols(l)) .eqv. .true. ) then
                write (*,*) "invalid coloring!!"
                STOP
            else
                usedcols(gd%verts(i)%cols(l)) = .true.
            endif
        enddo
        do l = 1, gd%usedcolors
            mynbr = gd%verts(i)%nbrbycol(l)
            nbr1 = gd%verts(i)%nbrs(mynbr)
            if (nbr1 > i) then ! nbr1 could be zero if there is no such edge
                k = k+1
                nodes(k)%x = i
                nodes(k)%y = nbr1
                nodes(k)%axy = gd%elems(elempos + mynbr)
                nodes(k)%col = l
            endif
        enddo
        elempos = elempos + gd%verts(i)%degree
    enddo
    weight = 1.0
    call fe%init(nodes, gd%usedcolors, method, weight)
    deallocate(nodes, usedcols)
end function

subroutine colorvertex_init(this, deg, maxcols)
    class(ColorVertex), intent(inout) :: this
    integer, intent(in) :: deg, maxcols

    call vertex_init(this, deg)
    allocate(this%cols(deg), this%nbrbycol(maxcols))
    this%cols = 0
    this%nbrbycol = 0
end subroutine colorvertex_init

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Tidy up allocated space.
!
!> @param[in] this The vertex that we consider
!--------------------------------------------------------------------

subroutine colorvertex_destruct(this)
    class(ColorVertex) :: this
    
    call vertex_destruct(this)
    deallocate(this%cols, this%nbrbycol)
end subroutine colorvertex_destruct

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> get the color of an incident edge.
!
!> @param[in] this The vertex that we consider
!> @param[in] vert the index of the neighbour
!> @param[in] col the value of the color
!--------------------------------------------------------------------
function colorvertex_get_edge_color(this, vert) result(col)
    class(ColorVertex) :: this
    integer, intent(in) :: vert
    integer :: col
    integer :: i
    
    col  = 0
    do i = 1, this%degree
        if(this%nbrs(i) == vert) col = this%cols(i)
    enddo
end function colorvertex_get_edge_color

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> erases the color of an incident edge.
!
!> @param[in] this The vertex that we consider
!> @param[in] vert the index of the neighbour
!> @param[in] col the value of the color
!--------------------------------------------------------------------
subroutine colorvertex_erase_edge_color(this, vert)
    class(ColorVertex) :: this
    integer, intent(in) :: vert
    integer :: i, tmp
#ifndef NDEBUG
    write (*,*) "Called erase_color!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#endif
    do i = 1, this%degree
        if(this%nbrs(i) == vert) then
            tmp = this%cols(i)
#ifndef NDEBUG
            write (*,*) "erasing color ", tmp, "of pos ", i, "which is", this%nbrs(i)
#endif
            this%cols(i) = 0
            if (tmp .ne. 0) this%nbrbycol(tmp) = 0
        endif
    enddo
#ifndef NDEBUG
    ! let's check the consistency
    do i = 1, size(this%nbrbycol)
    if (this%nbrbycol(i) .ne. 0) then
        if(this%cols(this%nbrbycol(i)) .ne. i ) then
        write (*,*) "inconsistency!"
        write (*,*) this%cols
        write (*,*) this%nbrbycol
        write (*,*) this%nbrbycol(10000)
        endif
    endif
    enddo
#endif
end subroutine colorvertex_erase_edge_color

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Set the color of an incident edge.
!
!> @param[in] this The vertex that we consider
!> @param[in] vert the index of the neighbour
!> @param[in] col the value of the color
!--------------------------------------------------------------------
subroutine colorvertex_set_edge_color(this, vert, col)
    class(ColorVertex) :: this
    integer, intent(in) :: vert, col
    integer :: i, up
    
    i = 1
    up = size(this%nbrs)
    do while ((this%nbrs(i) .ne. vert) .and. (i <= up ))
    i = i+1
    enddo
#ifndef NDEBUG
        write (*,*) "setting color", col, "to idx", i, "which is", vert
#endif
            this%cols(i) = col
            this%nbrbycol(col) = i
! #ifndef NDEBUG
!     ! let's check the consistency
!     do i = 1, size(this%nbrbycol)
!     if (this%nbrbycol(i) .ne. 0) then
!         if(this%cols(this%nbrbycol(i)) .ne. i ) then
!         write (*,*) "inconsistency!"
!         write (*,*) this%cols
!         write (*,*) this%nbrbycol
!         write (*,*) this%nbrbycol(10000)
!         endif
!     endif
!     enddo
! #endif
end subroutine colorvertex_set_edge_color

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function checks wether a free color is available and returns
!> the index of the neighbour vertex
!
!> @param[in] this The vertex that we consider.
!> @return the index of the neighbour where the edge is still colorless
!--------------------------------------------------------------------

function colorvertex_any_color_available(this) result(n)
    class(ColorVertex), intent(in) :: this
    integer :: n, i

    n = 0
    do i = this%degree, 1, -1
        if (this%cols(i) == 0) n = this%nbrs(i)
    enddo
end function colorvertex_any_color_available

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> An implementation of quicksort to sort an array of integers
!
!> @param[in] a the array in which we search
!> @param[in] first where to start sorting
!> @param[in] last where to stop sorting
!--------------------------------------------------------------------
recursive subroutine quicksort(a, first, last)
  implicit none
  integer, dimension(:), intent(inout) :: a
  integer :: x, t
  integer :: first, last
  integer :: i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> A function that transforms a matrix into our internal data structure
!
!> @param[in] A the matrix that we bring into our internal datastructure
!> @result our internal data structure of graph vertices
!--------------------------------------------------------------------
function mat2verts(A) result(gd)
    implicit none
    complex (kind=kind(0.d0)), intent(in) :: A(:,:)
    type(GraphData) :: gd
    integer :: i, j, maxcolors, k, i2
    integer, allocatable, dimension(:) :: cntarr
    
    gd%ndim = size(A, 1)
    ! check input
    ! first check diagonal
    do i = 1, gd%ndim
        if(A(i, i) /= 0.D0) then
            write (*, *) "the main-diagonal must be zero!"
            stop
        endif
    enddo
    ! check symmetry of input matrix
    do i = 1, gd%ndim
        do j = 1, gd%ndim
            if(A(i, j) /= conjg(A(j, i))) then
                write (*, *) "Non-hermitian matrix encountered!"
                stop
            endif
        enddo
    enddo
    allocate(gd%verts(gd%ndim), cntarr(gd%ndim))
    associate(ndim => gd%ndim, verts => gd%verts)
        ! calculate Vertex degree of each vertex
        cntarr = 0
        do j = 1, ndim-1
            do i = j+1, ndim
                if(A(i, j) /= 0.D0) then
                    cntarr(i) = cntarr(i) + 1
                    cntarr(j) = cntarr(j) + 1
                endif
            enddo
        enddo
        gd%deltag = maxval(cntarr)
        write (*,*) "Delta(G) = ", gd%deltag
        maxcolors = gd%deltag + 1
        allocate(gd%elems(sum(cntarr)))
        i2 = 1
        do j = 1, ndim
            call colorvertex_init(gd%verts(j), cntarr(j), maxcolors)
!            call gd%verts(j)%init(cntarr(j), maxcolors)
            k = 1
            do i = 1, ndim
                if(A(i, j) /= 0.D0) then
                    verts(j)%nbrs(k) = i
                    k = k + 1
                endif
            enddo
            call quicksort(verts(j)%nbrs, 1, verts(j)%degree) !< it is probably sorted already...
            do i = 1, verts(j)%degree ! set up array of elements
                gd%elems(i2) = A(verts(j)%nbrs(i), j) 
                i2 = i2 + 1
            enddo
        enddo
        deallocate(cntarr)
    end associate
end function
end module colorvertex_mod
