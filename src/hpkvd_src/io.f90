module io

  use iso_c_binding

  implicit none

  interface

     subroutine myfwrite(filename,buffer,size,offset) &
         bind(c,name="myfwrite") 

       use iso_c_binding

       character(C_CHAR), intent(in) :: filename
       type(C_PTR), value            :: buffer
       integer(C_LONG), value        :: size, offset

     end subroutine myfwrite

   end interface

contains

subroutine get_filename(ibox,filename)

  use mpivars
  use textlib
  use input_parameters

  character *128 filename,basename
  integer ibox

  basename = filein
  filename = basename
  if    (ibox<   10) then
     write(filename(indlnb(basename)+1:indlnb(basename)+2),40) ibox
  elseif(ibox<  100) then
     write(filename(indlnb(basename)+1:indlnb(basename)+3),41) ibox
  elseif(ibox< 1000) then
     write(filename(indlnb(basename)+1:indlnb(basename)+4),42) ibox
  elseif(ibox<10000) then
     write(filename(indlnb(basename)+1:indlnb(basename)+5),43) ibox
  else
     write(*,*) 'too many boxes...'
     call mpi_finalize(ierr)
     stop
  endif

40 format('.',i1)
41 format('.',i2)
42 format('.',i3)
43 format('.',i4)

  return

  end subroutine get_filename

end module io
