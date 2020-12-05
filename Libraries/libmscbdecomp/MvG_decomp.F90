! MIT License
! 
! Copyright (c) 2018 Florian Goth
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

module MvG_mod
    Use colorvertex_mod
    implicit none
    !Another helper type for the cd_X path construction
    type :: Path
        integer :: avamem
        integer :: tail ! Last valid Fortran index
        integer, allocatable, dimension(:) :: vertices
        integer, allocatable, dimension(:) :: nbrindex
    contains
        procedure :: init => Path_init
        procedure :: dealloc => Path_dealloc
        procedure :: pushback => Path_pushback
        procedure :: at => Path_at
        procedure :: back => Path_back
        procedure :: length => Path_length
    end type Path
    
contains

subroutine Path_init(this)
    class(Path) :: this
    integer :: temp
    this%tail = 0 ! when the vector is empty, this is invalid memory
    this%avamem = 4096/SIZEOF(temp) ! allocate a page of memory
    allocate(this%vertices(this%avamem), this%nbrindex(this%avamem))
end subroutine Path_init

subroutine Path_dealloc(this)
    class(Path) :: this
    deallocate(this%vertices, this%nbrindex)
end subroutine

subroutine Path_pushback(this, vert, nbridx)
    class(Path) :: this
    integer, intent(in) :: vert, nbridx
    integer, allocatable, dimension(:) :: temp1, temp2
    integer :: i
    if (this%tail == this%avamem) then
        ! reallocate the memory
        write (*,*) "not enough space!"
        call MOVE_ALLOC(this%vertices, temp1)
        call MOVE_ALLOC(this%nbrindex, temp2)
        allocate(this%vertices(2*this%avamem), this%nbrindex(2*this%avamem))
        do i = 1, this%avamem
            this%vertices(i) = temp1(i)
            this%nbrindex(i) = temp2(i)
        enddo
        deallocate(temp1, temp2)
        this%avamem = 2*this%avamem
    endif
    this%tail = this%tail + 1
    this%vertices(this%tail) = vert
    this%nbrindex(this%tail) = nbridx
end subroutine

subroutine Path_at(this, pos, vert, nbridx)
    class(Path) :: this
    integer, intent(in) :: pos
    integer, intent(out) :: vert, nbridx
    vert = this%vertices(pos)
    nbridx = this%nbrindex(pos)
end subroutine

subroutine Path_back(this, vert, nbridx)
    class(Path) :: this
    integer, intent(out) :: vert, nbridx
    vert = this%vertices(this%tail)
    nbridx = this%nbrindex(this%tail)
end subroutine

function Path_length(this) result(l)
    class(Path) :: this
    integer :: l
    l = this%tail
end function 

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function returns the smallest available color
!> that is free at both vertices
!
!> @param[in] v one vertex that we consider
!> @param[in] w the other vertex that we consider
!> @param[in] maxcols The available maximum number colors
!> @return the smallest common available color
!--------------------------------------------------------------------
function find_common_free_color(v, w, maxcols) result(col)
    implicit none
    class(ColorVertex), intent(in) :: v, w
    integer, intent(in) :: maxcols
    integer :: col, i

    col = 0
    i = 1
    do while (.not.((v%nbrbycol(i) == 0) .and. (w%nbrbycol(i) == 0)) .and. i <= maxcols)
        i = i + 1
    enddo
    if (i <= maxcols) col = i
!     do i = maxcols, 1, -1
!         if ((v%nbrbycol(i) == 0) .and. (w%nbrbycol(i) == 0)) col = i
!     enddo
end function find_common_free_color

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function returns the smallest available color
!> that is free at both vertices and tries to set that color if it
!> is free.
!>
!> @param[in] v one vertex that we consider
!> @param[in] w the other vertex that we consider
!> @param[in] maxcols The available maximum number colors
!> @return the smallest common available color
!--------------------------------------------------------------------
function find_and_try_to_set_common_free_color(i, j, verts, maxcols) result(col)
    implicit none
    type(ColorVertex), dimension(:), allocatable, intent(inout) :: verts
    integer, intent(in) :: maxcols, i, j
    integer :: col, k, l, pos, kf, ko

    col = 0
    k = 1
!     if (Modulo(maxcols, 2) == 0) then
!         pos = -1
!         do while ((pos < 0) .and. (k <= maxcols))
! !            write (*,*) "k = ", k, pos, verts(i)%nbrbycol(k), verts(j)%nbrbycol(k)
!             if ((verts(i)%nbrbycol(k) == 0) .and. (verts(j)%nbrbycol(k) == 0)) then
!             pos = 0
!             else
!             k = k + 1
!             endif
!         !    if ((verts(i)%nbrbycol(k+1) == 0) .and. (verts(j)%nbrbycol(k+1) == 0)) pos = 1
! !        write (*,*) "k = ", k, pos
!         enddo
!         k = k + pos
!         kf = k
! !        write (*,*) "kf = ", k, pos, verts(i)%nbrbycol(k), verts(j)%nbrbycol(k)
!         
!        k = 1 
!     do while (.not.((verts(i)%nbrbycol(k) == 0) .and. (verts(j)%nbrbycol(k) == 0)) .and. k <= maxcols)
! !    write (*,*) "k = ", k, verts(i)%nbrbycol(k), verts(j)%nbrbycol(k)
!         k = k + 1
!     enddo
!     ko = k
! !    write (*,*) "ko = ", k, pos, verts(i)%nbrbycol(k), verts(j)%nbrbycol(k)
!     if (kf .ne. ko) then
!     write (*,*) kf, ko, maxcols
!     endif
! !STOP
!     else
    do while (.not.((verts(i)%nbrbycol(k) == 0) .and. (verts(j)%nbrbycol(k) == 0)) .and. k <= maxcols)
        k = k + 1
    enddo
!     endif
    if (k <= maxcols) then
        col = k
        l = binarysearch(verts(i)%nbrs, j)
        verts(i)%cols(l) = col
        verts(i)%nbrbycol(col) = l
        l = binarysearch(verts(j)%nbrs, i)
        verts(j)%cols(l) = col
        verts(j)%nbrbycol(col) = l
    endif
end function find_and_try_to_set_common_free_color

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function constructs an alternating cd path at start
!> And then inverts the colors
!>
!> @param[in] col the two colors that make up the path
!> @param[inout] verts The array of vertices
!> @param[in] start the starting vertex
!--------------------------------------------------------------------
subroutine construct_and_invert_path(col, verts, start)
    implicit none
    type(ColorVertex), allocatable, dimension(:) :: verts
    integer, intent(in) :: col(2), start
    integer :: k, nbr1, ver, colctr, curver, nbr, tmpcol
    logical :: stoppath
    Type(Path) :: p

#ifndef NDEBUG
    write (*,*) "construct path of colors", col(1), col(2), "at", start
#endif
    call p%init()
    do k = 1, verts(start)%degree
        if (verts(start)%cols(k) == col(1)) nbr1 = k
    enddo
    ! We start the path at verts(i)
    ! (i, nbr1) is the first piece of the path. nbr1 is the position within the arrays of verts(i)
    ! The edge i -> nbr1 has color col(1)
#ifndef NDEBUG
    write (*,*) "starting path to ", nbr1, "which is", verts(start)%nbrs(nbr1)
#endif
    call p%pushback(start, nbr1)
    stoppath = .false.
    colctr = 2
    do while(stoppath .eqv. .false.)
        ! current vertex
        call p%back(ver, nbr)
        curver = verts(ver)%nbrs(nbr)
        ! find nbr of current color
        nbr1 = 0
        do k = 1, verts(curver)%degree
            if (verts(curver)%cols(k) == col(colctr)) nbr1 = k
        enddo
        if (nbr1 == 0) then
            ! presumably we have reached the end of the path
            stoppath = .true.
        else
            call p%pushback(curver, nbr1)
            colctr = colctr + 1
            if (colctr > 2) colctr = 1
        endif
    enddo
#ifndef NDEBUG
    write (*, *) "Length of path", p%length()
#endif
    ! Now we switch the colors on the path
    k = 1
    colctr = 2
    ! erase colors first
    do while(k<=p%length())
        call p%at(k, ver, nbr)
        tmpcol = verts(ver)%cols(nbr)
        verts(ver)%nbrbycol(tmpcol) = 0 ! to be on the safe side and create consistent data
                        
        call verts(verts(ver)%nbrs(nbr))%erase_edge_color(ver)
        k = k + 1
    enddo
    ! Now we perform the actual setting of the new path colors
    k = 1                   
    do while(k<=p%length())
        call p%at(k, ver, nbr)
!       tmpcol = verts(ver)%cols(nbr)
        verts(ver)%cols(nbr) = col(colctr)
        verts(ver)%nbrbycol(col(colctr)) = nbr
!       verts(ver)%nbrbycol(tmpcol) = 0 ! to be on the safe side and create consistent data
!       call verts(verts(ver)%nbrs(nbr))%erase_edge_color(ver)
        call verts(verts(ver)%nbrs(nbr))%set_edge_color(ver, col(colctr))
        k = k + 1
        colctr = colctr + 1
        if (colctr > 2) colctr = 1
    enddo
    call p%dealloc()
end subroutine

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> This function performs a downshift of a fan that was previously
!> constructed around fanpos. After that it sets the new available edge 
!> to the color col.
!
!> @param[in] fanpos the position where the vertex was constructed
!> @param[in] fan the actual length of the Vizing fan
!> @param[in] fannbr the positions of the neighbours in the array associated with verts(fanpos)
!> @param[in] fanedgecol the colors of the respective edges
!> @param[in] tail the length/end of the fan
!> @param[in] col the new color
!> @param[inout] the array of vertices.
!--------------------------------------------------------------------
subroutine downshift_fan_and_set_color(fanpos, fan, fannbr, fanedgecol, tail, col, verts)
    integer, allocatable, dimension(:), intent(in) :: fan, fannbr, fanedgecol
    type(ColorVertex), allocatable, dimension(:) :: verts
    integer, intent(in) :: tail, fanpos, col
    integer :: k, oldpos

    ! erase colors first that are about to become invalidated    
    do k = 2, tail
        ! erase current color at fanpos
        verts(fanpos)%cols(fannbr(k)) = 0
        verts(fanpos)%nbrbycol(fanedgecol(k)) = 0
    enddo
    
    do k = 2, tail-1
        ! set to new color
        verts(fanpos)%cols(fannbr(k)) = fanedgecol(k+1)
        verts(fanpos)%nbrbycol(fanedgecol(k+1)) = fannbr(k)
        
        ! set current color from fan(k) to fanpos. erase old entres first
        oldpos = verts(fan(k))%nbrbycol(fanedgecol(k))
        verts(fan(k))%nbrbycol(fanedgecol(k)) = 0

        verts(fan(k))%cols(oldpos) = fanedgecol(k+1)
!         
         verts(fan(k))%nbrbycol(fanedgecol(k+1)) = oldpos
    enddo
    
    verts(fanpos)%cols(fannbr(1)) = fanedgecol(2)
    verts(fanpos)%nbrbycol(fanedgecol(2)) = fannbr(1)
    call verts(fan(1))%set_edge_color(fanpos, fanedgecol(2))
    
    ! erase current color at fan(tail)
    verts(fan(tail))%cols( verts(fan(tail))%nbrbycol(fanedgecol(tail)) ) = 0
    verts(fan(tail))%nbrbycol(fanedgecol(tail)) = 0
    ! set to new color
    verts(fanpos)%cols(fannbr(tail)) = col
    verts(fanpos)%nbrbycol(col) = fannbr(tail)
    call verts(fan(tail))%set_edge_color(fanpos, col)
end subroutine

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> Binary search within the array
!
!> @param[in] a the array in which we search
!> @param[in] value the value that we are looking for
!> @result the position in the array
!--------------------------------------------------------------------
function binarySearch (a, value)
    integer                  :: binarySearch
    integer, intent(in), target :: a(:)
    integer, intent(in)         :: value
    integer, pointer            :: p(:)
    integer                  :: mid, offset
 
    p => a
    binarySearch = 0
    offset = 0
    do while (size(p) > 0)
        mid = size(p)/2 + 1
        if (p(mid) > value) then
            p => p(:mid-1)
        else if (p(mid) < value) then
            offset = offset + mid
            p => p(mid+1:)
        else
            binarySearch = offset + mid    ! SUCCESS!!
            return
        end if
    end do
end function binarySearch

!--------------------------------------------------------------------
!> @author
!> Florian Goth
!
!> @brief 
!> the Misra van-Gries graph decomposition
!> If everything works out, the col array that each vertex has is completely set afterwards.
!
!> @param[in] verts our internal data structure composed of vertices.
!--------------------------------------------------------------------
subroutine MvG_decomp(verts)
    implicit none
    type(ColorVertex), allocatable, dimension(:) :: verts
    integer :: maxcolors
    integer :: i, i2, j, k, availablecolor, fanlen, oldcol, ver
    integer, allocatable, dimension(:) :: fan, fannbr, fanedgecol
    integer :: col(2)
    maxcolors = 0
    ! determine the maximum degree and therefore the maximum allowed colors beforehand
    do i = 1, size(verts)
        maxcolors = max(maxcolors, verts(i)%degree)
    enddo
    maxcolors = maxcolors + 1
    ! Starting Vizings algorith as outlined in https://thorehusfeldt.files.wordpress.com/2010/08/gca.pdf
    ! and the help of the wikipedia page.
    do i = 1, size(verts)-1
        do i2 = 1, verts(i)%degree
        if ( verts(i)%cols(i2) == 0) then
            j = verts(i)%nbrs(i2)
        ! Edge found between vertex i and j
! ! ! ! A debugging check that the data is consistent:
! ! !  check = .false.
! ! ! do k = 1, verts(i)%degree
! ! !  if (verts(i)%nbrs(k) == j) check = .true.
! ! ! enddo
! ! ! if(check .eqv. .false.) then
! ! ! write(*,*) "inconsistent data!"
! ! ! STOP
! ! ! endif
        ! Let's check wether we have free edges at every vertex
            availablecolor = find_and_try_to_set_common_free_color(i, j, verts, maxcolors)
            if(availablecolor == 0) then
                ! Our starting vertex is verts(i), our target vertex is verts(j)
                ! Now we need to construct a Vizing fan around verts(i)
#ifndef NDEBUG
                write (*,*) "Out of colors. construct fan around", i
#endif
                allocate(fan(verts(i)%degree), fannbr(verts(i)%degree), fanedgecol(verts(i)%degree))
                
                fannbr = 0
                fanedgecol = 0
                do k = 1, verts(i)%degree
                    if(verts(i)%nbrs(k) == j) fannbr(1) = k
                enddo
                oldcol = verts(i)%find_maximal_fan(verts, j, maxcolors, fan, fannbr, fanedgecol, fanlen)
#ifndef NDEBUG
                write (*,*) "oldcol = ", oldcol, "fanlen = ", fanlen, fan
#endif
                if (oldcol .ne. 0) then
#ifndef NDEBUG
                write (*,*) "Free color available in fan -> downshift"
#endif
                    ! the end of the fan has a free color -> down-shifting sufficient
                    call downshift_fan_and_set_color(i, fan, fannbr, fanedgecol, fanlen, oldcol, verts)
                else
                    ! We would need to inverse a path
                    ! set up the start of the path at the fan
                    ! determine the two colors for the c-d path
                    col(1) = verts(fan(fanlen))%findfreecolor()
                    col(2) = verts(i)%findfreecolor()
                    call construct_and_invert_path(col, verts, i)
                    ! try again to obtain a color
                    availablecolor = find_and_try_to_set_common_free_color(i, j, verts, maxcolors)
                    if (availablecolor == 0) then
                        ! We have to do a proper downshifting
                        ! find w in fan()
                        k = 1
                        ver = 0
                        do while (k <= fanlen)
                            if(verts(fan(k))%iscolorfree(col(1))) then
#ifndef NDEBUG
                                write(*,*) "color found at ",k, "of ", fanlen
#endif
                                ver = k
                                k = fanlen
                            endif
                            k = k+1
                        enddo
#ifndef NDEBUG
                        write (*,*) "determined ", ver
#endif
                        ! the path inversion has modified the information contained in the auxiliary arrays
                        ! fannbr and fanedgecol -> We need to update that
                        ! The edge that has been modified had color col(1) and now has color col(2)
                        ! search the edge that has changed
                        k = 1
                        do while ((k <= fanlen) .and. (fanedgecol(k) .ne. col(1)))
                            k = k+1
                        enddo
#ifndef NDEBUG
                        if (k > fanlen) then
                            write (*,*) "changed color ", col(1), "not found"
                            STOP
                        endif
#endif
                        fanedgecol(k) = col(2)
                        ! downshift the fan until the position of w
                        call downshift_fan_and_set_color(i, fan, fannbr, fanedgecol, ver, col(1), verts)
                    endif
                endif
                deallocate(fan, fannbr, fanedgecol)
            endif
        endif
        enddo
    enddo
end subroutine

end module MvG_mod
