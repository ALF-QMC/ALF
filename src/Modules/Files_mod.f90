Module Files_mod
   contains

     Character (len=64) function File_i( file, I)
        character (len=64) :: file
        integer            :: i
        write(File_i,'(A,"_",I0)') trim(file),i
      end function File_i

     Character (len=64) function File_add( file, file1)
        character (len=64) :: file, file1
        write(File_add,'(A,A)') trim(file),Trim(file1)
      end function File_add

end Module Files_mod
