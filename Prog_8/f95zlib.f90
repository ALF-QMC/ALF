!
! Fortran interfaces for zlib and for POSIX pipes (popen and friends)
! Author:  David L Duffy 2009-2014
! 
! Copyright (c) 2014 David L Duffy
! 
! This software is provided 'as-is', without any express or implied
! warranty. In no event will the authors be held liable for any damages
! arising from the use of this software.
! 
! Permission is granted to anyone to use this software for any purpose,
! including commercial applications, and to alter it and redistribute it
! freely, subject to the following restrictions:
! 
!    1. The origin of this software must not be misrepresented; you must not
!    claim that you wrote the original software. If you use this software
!    in a product, an acknowledgment in the product documentation would be
!    appreciated but is not required.
! 
!    2. Altered source versions must be plainly marked as such, and must not be
!    misrepresented as being the original software.
! 
!    3. This notice may not be removed or altered from any source
!    distribution.
!
! AS182 was published in the journal Applied Statistics, under the copyright of the Royal
! Statistical Society, and I understand was thus made available subject 
! to the restriction that no fee is charged for redistribution
!
! C preprocessor directives support various Fortran compilers. Where
! iso_c_binding is not available, gzipped files can be decompressed
! using the gzip executable (providing the Fortran compiler supports the
! SYSTEM command)
!
! The readline subroutine may fail reading very large files
!
! 20140106: initial release
! 20140110: added in missing iocodes module, fixed typos 
!           pointed out by John David and yeonpil
!
! iostat codes needed for wrinline etc
!
module iocodes
#if OPEN64
  integer, parameter :: eofcode = -4001
  integer, parameter :: eolcode = -4006
  character (len=10), parameter :: stream_access = 'sequential'
  character (len=6), parameter :: stream_form = 'binary'
#else
  integer, parameter :: eofcode = -1
  integer, parameter :: eolcode = -2
  character (len=6), parameter :: stream_access = 'stream'
  character (len=11), parameter :: stream_form = 'unformatted'
#endif
end module iocodes
!
! Definition of a port
!   slots: associated file name
!          1=uncompressed 2=gzipped 3=unzipped copy 4=pipe
!          Fortran style logical unit number
!          gzip C-style file handle
!
module ioports
#if !OPEN64
  use, intrinsic :: iso_c_binding
#endif
  integer, parameter :: PORT_STANDARD = 1, PORT_GZIPPED = 2,  &
                        PORT_COPY = 3, PORT_PIPE = 4
  type, public :: ioport
    character (len=256) :: filnam
    character (len=1) :: stat = ' '
    integer :: filtyp
    integer :: fstream
#if !OPEN64
    type (c_ptr) :: handle = c_null_ptr
#endif
  end type ioport
end module ioports
!
! Fortran interface to zlib 
!   based on looking at fgzlib, fgsl and Janus Weil's example
!   on comp.lang.fortran May 2009
!   currently enough functionality to read gzipped text files
!
#if ZLIB
module f95zlib
  use, intrinsic :: iso_c_binding
  use ioports
! buffer for gzread
  integer, parameter :: ZBUFLEN = 65536
  character (len=ZBUFLEN), target :: zbuffer
! current character and end of zbuffer
  integer :: zbufpos=0, zbufend=ZBUFLEN
! gzopen  
  interface
    function gzopen(path, mode) bind(C, name='gzopen')
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: path, mode
      type (c_ptr) :: gzopen
    end function
  end interface
! gzread  
  interface
    function gzread(filehandle, buf, len) bind(C, name='gzread')
      use, intrinsic :: iso_c_binding
      integer (c_int) :: gzread
      type (c_ptr), value :: filehandle
      type (c_ptr), value :: buf
      integer(c_int), value :: len
    end function
  end interface
! gzwrite
  interface
    function gzwrite (filehandle, buf, len) bind(C, name='gzwrite')
      use, intrinsic :: iso_c_binding
      integer (c_int) :: gzwrite
      type (c_ptr), value :: filehandle
      type (c_ptr), value :: buf
      integer(c_int), value :: len
    end function
  end interface
! gzgetc  
  interface
    function gzgetc(filehandle) bind(C, name='gzgetc')
      use, intrinsic :: iso_c_binding
      integer (c_int) :: gzgetc
      type (c_ptr), value :: filehandle
    end function
  end interface
! gzrewind 
  interface
    function gzrewind(filehandle) bind(C, name='gzrewind')
      use, intrinsic :: iso_c_binding
      integer(c_int) :: gzrewind
      type (c_ptr), value :: filehandle
    end function
  end interface
! gzclose  
  interface
    function gzclose(filehandle) bind(C, name='gzclose')
    use, intrinsic :: iso_c_binding
    integer(c_int) :: gzclose
    type (c_ptr), value :: filehandle
    end function
  end interface
contains
!
! Wrapper for gzopen
!   also reinitializes gzread's buffer
!
  subroutine fgz_open(path, mode, fd, ios)
    use, intrinsic :: iso_c_binding
    character(kind=c_char, len=*), intent(in) :: path, mode
    type (ioport) :: fd
#if SUN
    character(kind=c_char, len=len_trim(path)+1) :: cpath
    character(kind=c_char, len=len_trim(mode)+1) :: cmode
    integer :: eos
#endif
    integer :: ios

    ios=0
    fd%filnam=path
    fd%filtyp=PORT_GZIPPED 
    fd%fstream=-1
#if SUN
    eos=len_trim(path)
    cpath=path
    cpath((eos+1):(eos+1))=c_null_char
    eos=len_trim(mode)
    cmode=mode
    cmode((eos+1):(eos+1))=c_null_char
    fd%handle = gzopen(cpath, cmode)
#else
    fd%handle = gzopen(trim(path) // c_null_char, trim(mode) // c_null_char)
#endif
    if (.not.c_associated(fd%handle)) ios=-1
    zbufpos=0
  end subroutine fgz_open
!
! Wrapper for gzrewind
!
  subroutine fgz_rewind(fd, ios)
    use, intrinsic :: iso_c_binding
    type(ioport) :: fd
    integer :: ios
    integer(c_int) :: ir
    ios = 0
    ir = gzrewind(fd%handle)
    if (ir /= 0) ios=ir
    zbufpos=0
  end subroutine fgz_rewind
!
! Wrapper for gzread
!   read one line of text from buffer
!
  subroutine fgz_read(fd, lin, advance, ios)
    use, intrinsic :: iso_c_binding
    use iocodes
    type(ioport) :: fd
    character(len=*) :: lin
    character(len=*), intent(in), optional :: advance
    integer, intent(out) :: ios

    integer :: i, j, linlen, nchar, newzpos, pos
    integer(c_int) :: blen, rlen
!
!  eol morez more
!   F    T    T    read buffer, copy to output
!   F    T    F    read buffer, output full
!   T    F    F    found <NL>
!  advancing
!   no             after output full, exit with buffer pos at end of text
!   yes            after output full, exit with buffer pos at next <NL>
!
    logical :: advancing, eol, more, morez

    type (c_ptr) :: buf = c_null_ptr

    advancing=.true.
    if (present(advance)) advancing=(advance == 'yes') 
    linlen=len(lin)
    ios=0
    lin=' '
    sta=1
    nchar=-1
    pos=0
    j=0
    eol=.false.
    more=.true.
    morez=.true.
    do while (morez)
      j=j+1
! refill buffer if necessary
      if (zbufpos == 0) then
        blen=ZBUFLEN
        buf=c_loc(zbuffer(1:1))
        rlen=gzread(fd%handle, buf, blen)
        if (rlen <= 0) then
          ios=-1
          return
        end if
        zbufpos=1
        zbufend=rlen
      end if
! place buffer index at <NL> or buffer end
! if <NL> will exit after updating output
      newzpos=zbufend+1
      nchar=zbufend-zbufpos+1
      do i=zbufpos, zbufend
        if (zbuffer(i:i) == achar(10)) then
          eol=.true.
          morez=.false.
          newzpos=i+1
          nchar=i-zbufpos
          exit
        end if
      end do
! read in min(buffer, remaining output)
! if not advancing move buffer idx back to last character read and exit
      if (more) then
        if (linlen < pos+nchar) then
          more=.false.
          nchar=linlen-pos
          if (.not.advancing) then
            newzpos=zbufpos+nchar
            morez=.false.
            eol=.false.
          end if
        end if
        lin((pos+1):(pos+nchar))=zbuffer(zbufpos:(zbufpos+nchar-1))
        pos=pos+nchar
      end if
      zbufpos=newzpos
      if (zbufpos > zbufend) then
        zbufpos=0
      end if
    end do
    if (.not.advancing .and. eol) ios=eolcode
  end subroutine fgz_read
!
! write one line of text to a gzipped textfile
!
  subroutine fgz_write(fd, lin, advance, ios)
    use, intrinsic :: iso_c_binding
    use iocodes
    type(ioport) :: fd
    character(len=*) :: lin
    character(len=*), intent(in), optional :: advance
    integer, intent(out) :: ios

    logical :: advancing
    integer :: ioerr, linlen, lsta, lpos
    integer(c_int) :: blen, wlen
    type (c_ptr) :: buf = c_null_ptr

    advancing=.true.
    if (present(advance)) then
      advancing=(advance == 'yes')
    end if
    ios=0
    lpos=0
    lenlin=len_trim(lin)
    do 
      lsta=lpos+1
      lpos=min(lsta+ZBUFLEN-1, lenlin)
      zbuffer=lin(lsta:lpos)
      buf=c_loc(zbuffer(1:1))
      blen=lpos-lsta+1
      wlen=gzwrite(fd%handle, buf, blen)
      ioerr=wlen
      if (ioerr == 0) exit
      if (lpos == lenlin) then
        if (advancing) then
          zbuffer=char(10)
          buf=c_loc(zbuffer(1:1))
          blen=1
          wlen=gzwrite(fd%handle, buf, blen)
          ioerr=wlen
        end if
        exit
      end if
    end do
    if (ioerr == 0) ios=-1
  end subroutine fgz_write
!
! Wrapper for gzclose
!
  subroutine fgz_close(fd, ios)
    use, intrinsic :: iso_c_binding
    type(ioport) :: fd
    integer :: ios
    integer(c_int) :: ic 
    ios = 0
    ic = gzclose(fd%handle)
    if (ic /= 0) ios = ic
  end subroutine fgz_close
end module f95zlib
#endif
