program pks2asci

  use pksc

  character *128 filein,fileout

  call getarg(1,filein)
  call getarg(2,fileout)

  call loadpksc(filein)

  open(unit=1,file=fileout)
  do i=1,nhalo
     write(1,*) (pos(j,i),j=1,3)
  enddo
  close(1)

end program pks2asci
