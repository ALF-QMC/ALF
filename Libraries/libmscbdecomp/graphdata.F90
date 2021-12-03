! MIT License
! 
! Copyright (c) 2021 Florian Goth
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

module graphdata_mod
    use colorvertex_mod, only: ColorVertex, colorvertex_init
    use vertex_mod
    use Exponentials_mod, only: EulerExp, FullExp
    implicit none

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

subroutine dealloc_graphdata(gd)
    implicit none
    type(GraphData) :: gd
    
    deallocate(gd%elems, gd%verts)
end subroutine

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> An implementation of quicksort to sort an array of integers.
!
!> @param[in] a the array in which we search
!> @param[in] first where to start sorting
!> @param[in] last where to stop sorting
!--------------------------------------------------------------------
recursive subroutine quicksort(a, first, last)
  implicit none
  integer, dimension(:), intent(inout) :: a
  integer, intent(in) :: first, last
  integer :: x, t
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
!> A function that transforms a matrix into our internal data structure.
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

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> A function to fill the color information of a graph.
!
!> @param[inout] gd A graphdata object.
!--------------------------------------------------------------------
subroutine determine_used_colors_of_graph(gd)
    implicit none
    type(GraphData), intent(inout) :: gd
    integer :: i, k

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
end subroutine

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> A function that transforms the graphdata object to an array of nodes.
!
!> @param[in] gd A fully prepared graphdata object.
!> @result An allocated array of nodes.
!--------------------------------------------------------------------
function gd_to_nodes(gd) result(nodes)
    use node_mod
    implicit none
    type(GraphData), intent(in) :: gd
    type(node), allocatable, dimension(:) :: nodes
    
    integer :: k, elempos, i, l, mynbr, nbr1
    logical, allocatable, dimension(:) :: usedcols

    allocate(nodes(gd%nredges), usedcols(gd%usedcolors))
! set up data in an edges based layout
    k = 0
    elempos = 0
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
            if (mynbr > 0) then ! only do sth. if the color is associated with an edge
                nbr1 = gd%verts(i)%nbrs(mynbr)
                if (nbr1 > i) then ! nbr1 could be zero if there is no such edge
                    k = k+1
                    nodes(k)%x = i
                    nodes(k)%y = nbr1
                    nodes(k)%axy = gd%elems(elempos + mynbr)
                    nodes(k)%col = l
                endif
            endif
        enddo
        elempos = elempos + gd%verts(i)%degree
    enddo
    deallocate(usedcols)
end function

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> A function that findsfrom an array of colors 
!> The color that has the most edges associated to it,
!> according to the data in edges_per_color.
!
!> @param[in] edges_per_color An array containing for each color the number of edges
!> @param[in] cols An array of colors from which we rty to find the biggest
!> @result The color in col that has the most edges associated with it.
!--------------------------------------------------------------------
function find_biggest_color(edges_per_color, cols) result(ret)
    integer, allocatable, dimension(:), intent(in) :: edges_per_color, cols
    integer :: ret
    integer :: i, start
    
    start = 1
    do while (cols(start) == 0)
        start = start + 1
    enddo
    
    ret = cols(start)
    do i = start+1, size(cols)
        if (cols(i) > 0) then
           if (edges_per_color(cols(i)) > edges_per_color(ret)) ret = cols(i)
        endif
    enddo
end function

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function encodes the strategy that is ued to distribute the
!> main-diagonal ov the various colors.
!
!> Currently we try to put the diagonal into the color which hosts the most
!> edges. If the edges of touch the entire diagonal then only a single color
!> with a GeneralExp, or HomogeneousExp is generated. The rest will be the fast
!> ZeroDiag type exponentials.
!
!> @param gd  The matrix in the graphdata representation.
!> @param nodes The matrix in an array of nodes represntation.
!> @param diag The diagonal of the matrix
!> @result A nrofcolors x ndim array where for each color a diagonal is present according to the strategy.
!--------------------------------------------------------------------
function distribute_diagonal_over_colors(gd, nodes, diag) result(dcol)
    use node_mod
    implicit none
    type(GraphData), intent(in) :: gd
    type(node), allocatable, dimension(:), intent(in) :: nodes
    real(kind=kind(0.D0)), intent(in), allocatable, dimension(:) :: diag
    real(kind=kind(0.D0)), allocatable, dimension(:, :) :: dcol ! color separated diagonal
    
    integer :: i, tmpcol
    integer, allocatable, dimension(:) :: edgespercol


    allocate(edgespercol(gd%usedcolors), dcol(gd%usedcolors, gd%ndim) )
    edgespercol = 0
    do i = 1, size(nodes)
        edgespercol(nodes(i)%col) = edgespercol(nodes(i)%col) + 1
    enddo
    
    
    dcol = 0
    do i = 1, gd%ndim ! for every chemical potential on the diagonal, do
        ! Data is actually present.
        ! Scale-wise decisions are taken when setting up the *expbase classes
        if(abs(diag(i)) > 0.D0 ) then
            tmpcol = find_biggest_color(edgespercol, gd%verts(i)%cols)
            dcol(tmpcol, i) = diag(i) ! move diagonal to color
        endif
    enddo
    deallocate(edgespercol)
end function

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
!> @result A EulerExp object that can be utilized for matrix multiplications
!--------------------------------------------------------------------
function createEulerExponentialfromGraphData(gd, diags) result(ee)
    use node_mod
    implicit none
    type(GraphData) :: gd
    real(kind=kind(0.D0)), intent(in), allocatable, dimension(:) :: diags
    type(EulerExp) :: ee
    real(kind=kind(0.D0)) :: weight
    real(kind=kind(0.D0)), allocatable, dimension(:, :) :: dcol ! color separated diagonal
    type(node), allocatable, dimension(:) :: nodes
    
    if ((gd%usedcolors == 0) .or. (gd%nredges == 0)) then ! check that those are available
        call determine_used_colors_of_graph(gd)
    endif

    ! set up data in an edges based layout
    nodes = gd_to_nodes(gd)
    
    ! distribute the diagonal over the colors
    dcol = distribute_diagonal_over_colors(gd, nodes, diags)

    weight = 1.0
    call ee%init(nodes, gd%usedcolors, dcol, weight)
    deallocate(nodes, dcol)
end function

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function takes a graphdata object as e.g. determined with the
!> help of the MvG_decomp function and creates a FullExp(=product 
!> of checkerboard exponentials) object from it.
!
!> @param gd
!> @param method the method that determines the approximation type.
!> @result A FullExp object that can be utilized for matrix multiplications.
!--------------------------------------------------------------------
function createFullExponentialfromGraphData(gd, diags, method) result(fe)
    use node_mod
    implicit none
    type(GraphData) :: gd
    real(kind=kind(0.D0)), intent(in), allocatable, dimension(:) :: diags
    real(kind=kind(0.D0)), allocatable, dimension(:,:) :: dcol
    integer, intent(in) :: method
    type(FullExp) :: fe
    real(kind=kind(0.D0)) :: weight
    type(node), allocatable, dimension(:) :: nodes

    if ((gd%usedcolors == 0) .or. (gd%nredges == 0)) then ! check that those are available
        call determine_used_colors_of_graph(gd)
    endif

    ! set up data in an edges based layout
    nodes = gd_to_nodes(gd)
    
    ! distribute the diagonal over the colors
    dcol = distribute_diagonal_over_colors(gd, nodes, diags)
    
    weight = 1.0
    call fe%init(nodes, gd%usedcolors, dcol, method, weight)
    deallocate(nodes, dcol)
end function
end module graphdata_mod
